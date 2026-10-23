"""Batch-validate MotifCraft MPNN sequences with ESMFold2 on Modal (monomer, no MSA).

Image and model calls follow KrajTechne/Optimus (modal_common.py, RunESMFold2.py), but nothing is imported from it.

    modal run --detach modal_esmfold2_batch.py            # all sequences
    modal run modal_esmfold2_batch.py --limit 3           # quick test

Inputs (defaults are for the 1PRW run):
  --stats      MotifCraft mpnn_design_stats.csv (Design, Sequence, AF2 metrics)
  --log        MotifCraft run log; gives each trajectory's motif residues (input vs. in scaffold)
  --reference  starting PDB the motif came from
Outputs: experiments/<name>/data/esmfold2_results.csv locally; predicted .cif files in Volume `motifcraft-esmfold2`.
"""
import io
import re
from pathlib import Path

import modal

app = modal.App("motifcraft-esmfold2")
out_vol = modal.Volume.from_name("motifcraft-esmfold2", create_if_missing=True)
hf_vol = modal.Volume.from_name("esmfold2-hf-cache", create_if_missing=True)
OUT_DIR, HF_DIR = "/out", "/hf-cache"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "torch==2.11.0", "torchvision==0.26.0", "torchaudio==2.11.0",
        index_url="https://download.pytorch.org/whl/cu130",
    )
    .apt_install("git")
    .pip_install("esm@git+https://github.com/Biohub/esm.git@main")
    .pip_install("hf-transfer", "biotite", "numpy", "pandas")
    .env({"HF_HOME": HF_DIR})
)


def motif_rmsd(pred_atoms, pred_idx, ref_atoms, ref_res_ids):
    """Cα RMSD after superimposing the predicted motif residues (0-based positions in chain A) onto the reference motif."""
    import biotite.structure as struc
    import numpy as np

    pred_ca = pred_atoms[(pred_atoms.chain_id == "A") & (pred_atoms.atom_name == "CA")]
    ref_ca = ref_atoms[(ref_atoms.chain_id == "A") & (ref_atoms.atom_name == "CA")]
    ref_motif = ref_ca[np.isin(ref_ca.res_id, ref_res_ids)]
    pred_motif = pred_ca[pred_idx]
    fitted, _ = struc.superimpose(ref_motif, pred_motif)
    return float(struc.rmsd(ref_motif, fitted))


def expand(spec):
    """'A1-20,A21-40' -> [1..20, 21..40]"""
    out = []
    for seg in spec.split(","):
        a, b = re.fullmatch(r"[A-Za-z](\d+)-(\d+)", seg.strip()).groups()
        out += list(range(int(a), int(b) + 1))
    return out


@app.function(
    image=image, gpu="A100", timeout=2 * 3600,
    volumes={OUT_DIR: out_vol, HF_DIR: hf_vol},
)
def predict_all(jobs: list, reference_pdb: str, run_name: str, num_loops: int, num_sampling_steps: int, seed: int):
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io.pdbx as pdbx
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder, EsmFold2Model, ProteinInput, StructurePredictionInput,
    )

    model = EsmFold2Model.from_pretrained("biohub/ESMFold2", device="cuda").eval()
    hf_vol.commit()  # keep the downloaded weights for later runs
    ref = pdb.PDBFile.read(io.StringIO(reference_pdb)).get_structure(model=1)
    out = Path(OUT_DIR) / run_name
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for i, job in enumerate(jobs):
        try:
            spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=job["sequence"], msa=None)])
            pred = ESMFold2InputBuilder().fold(
                model, spi, num_loops=num_loops, num_sampling_steps=num_sampling_steps,
                num_diffusion_samples=1, seed=seed,
            )
            if isinstance(pred, list):
                pred = pred[0]
            cif = pred.complex.to_mmcif()
            (out / f"{job['design']}.cif").write_text(cif)

            row = {"Design": job["design"], "esm_pLDDT": float(pred.plddt.mean().item()), "esm_pTM": float(pred.ptm)}
            structure = pdbx.get_structure(pdbx.CIFFile.read(io.StringIO(cif)), model=1)
            sc_idx = [n - 1 for n in job["scaffold_res"]]  # scaffold residue numbers are 1-based and sequential
            row["esm_motif_ca_rmsd"] = motif_rmsd(structure, sc_idx, ref, job["input_res"])
            row["esm_error"] = ""
        except Exception as e:  # keep going; record the failure
            row = {"Design": job["design"], "esm_error": f"{type(e).__name__}: {e}"}
        rows.append(row)
        if (i + 1) % 25 == 0:
            print(f"{i + 1}/{len(jobs)} done", flush=True)
    out_vol.commit()
    return rows


def motif_map(log_path):
    """trajectory name -> (input residues, scaffold residues), from the MotifCraft log."""
    mapping, inp, scaf = {}, None, None
    for line in Path(log_path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("backbone_constrained_residues (input):"):
            inp = expand(line.split(":", 1)[1])
        elif line.startswith("backbone_constrained_residues (in scaffold):"):
            scaf = expand(line.split(":", 1)[1])
        elif line.startswith("Starting trajectory:") and inp and scaf:
            mapping[line.split(":", 1)[1].strip()] = (inp, scaf)
    return mapping


@app.local_entrypoint()
def main(
    stats: str = "experiments/1PRW/data/mpnn_design_stats.csv",
    log: str = "experiments/1PRW/logs/run.log",
    reference: str = "targets/1PRW.pdb",
    limit: int = 0,
    out: str = "",
    num_loops: int = 20,
    num_sampling_steps: int = 100,
    seed: int = 0,
):
    import pandas as pd

    mapping = motif_map(log)
    df = pd.read_csv(stats)
    jobs, skipped = [], []
    for _, r in df.iterrows():
        traj = r["Design"].rsplit("_mpnn", 1)[0]
        if traj not in mapping or len(mapping[traj][0]) != len(mapping[traj][1]):
            skipped.append(r["Design"])
            continue
        jobs.append({"design": r["Design"], "sequence": r["Sequence"],
                     "input_res": mapping[traj][0], "scaffold_res": mapping[traj][1]})
    if skipped:
        print(f"Skipping {len(skipped)} designs without a usable motif mapping in the log: {skipped[:3]}...")
    if limit:
        jobs = jobs[:limit]
    print(f"Predicting {len(jobs)} sequences")

    run_name = Path(stats).parent.parent.name  # experiments/<name>/data/... -> <name>
    rows = predict_all.remote(jobs, Path(reference).read_text(), run_name, num_loops, num_sampling_steps, seed)

    keep = [c for c in ["Design", "Sequence", "MPNN_score", "MPNN_seq_recovery", "Average_pLDDT", "1_pLDDT", "1_pTM"] if c in df.columns]
    merged = df[keep].merge(pd.DataFrame(rows), on="Design", how="right")
    out_csv = Path(out) if out else Path(stats).parent / "esmfold2_results.csv"
    merged.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")
    ok = merged[merged["esm_error"] == ""]
    print(f"{len(ok)}/{len(merged)} predicted OK")
    if len(ok):
        print(ok[["esm_pLDDT", "esm_pTM", "esm_motif_ca_rmsd"]].describe().round(2))
