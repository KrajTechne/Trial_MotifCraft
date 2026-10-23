"""Run MotifCraft on Modal.

One-time setup (downloads AF2 weights into a Volume):
    modal run modal_app.py::download_weights

Run a design job (local settings JSON + local PDB are uploaded; results stay in a Volume).
Use --detach so the job survives this terminal closing; run details go to outputs/<PDB ID>/run_info.json:
    modal run --detach modal_app.py --settings configs/a_1PRW_motif_modal.json --pdb targets/1PRW.pdb \
        --filters MotifCraft/settings_filters/default_filters_partial.json \
        --advanced configs/default_4stage_multimer_partial_af2.json [--gpu A100]
Save the log (printed as "Logs:" after spawning):
    modal app logs <app_id> -f > experiments/<PDB ID>/logs/run.log   # -f streams until the app stops

Fetch results:
    modal volume get motifcraft-results <binder_name> ./results
"""
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

import modal

APP_DIR = "/root/MotifCraft"
PARAMS_DIR = "/params"
RESULTS_DIR = "/results"

app = modal.App("motifcraft")
params_vol = modal.Volume.from_name("motifcraft-params", create_if_missing=True)
results_vol = modal.Volume.from_name("motifcraft-results", create_if_missing=True)

CONDA_PKGS = [
    "pandas", "matplotlib", "numpy<2.0.0", "biopython", "scipy", "pdbfixer", "seaborn",
    "libgfortran5", "tqdm", "fsspec", "py3dmol", "pyrosetta",
    "chex", "dm-haiku", "flax<0.10.0", "dm-tree", "joblib", "ml-collections",
    "immutabledict", "optax", "jax>=0.4,<=0.6.0", "jaxlib>=0.4,<=0.6.0=*cuda*",
    "cuda-nvcc", "cudnn",
]

image = (
    modal.Image.micromamba(python_version="3.10")
    .apt_install("git", "wget")
    .env({"CONDA_OVERRIDE_CUDA": "12.4"})
    .micromamba_install(
        *CONDA_PKGS,
        channels=["conda-forge", "nvidia", "https://conda.graylab.jhu.edu"],
    )
    # matplotlib needs ffmpeg to render trajectory animations (the install script gets it from conda)
    .apt_install("ffmpeg")
    # custom ColabDesign (motif scaffolding code) is a nested git submodule of MotifCraft, installed without deps.
    # The pinned commit has an over-indented line in colabdesign/af/inputs.py (IndentationError on import). It is
    # dedented with sed (a patch file breaks on the CRLF line endings of a Windows checkout), and py_compile makes
    # the build fail loudly if the file is still broken, e.g. after bumping the submodule.
    .add_local_dir("MotifCraft/ColabDesign", f"{APP_DIR}/ColabDesign", copy=True)
    .run_commands(
        f"sed -i 's/^        seq = soft_seq(/    seq = soft_seq(/' {APP_DIR}/ColabDesign/colabdesign/af/inputs.py",
        f"python -m py_compile {APP_DIR}/ColabDesign/colabdesign/af/inputs.py",
        f"pip install {APP_DIR}/ColabDesign --no-deps",
    )
    .add_local_dir(
        "MotifCraft", APP_DIR,
        ignore=["**/.git", "**/.git/**", "ColabDesign/**", "params/**", "**/__pycache__/**", "**/*.pyc"],
        copy=True,
    )
    # dssp / DAlphaBall binaries are gitignored in MotifCraft; BindCraft ships them
    .run_commands(
        *[
            f"wget -q -O {APP_DIR}/functions/{b} https://raw.githubusercontent.com/martinpacesa/BindCraft/main/functions/{b}"
            for b in ("dssp", "DAlphaBall.gcc")
        ],
        f"chmod +x {APP_DIR}/functions/dssp {APP_DIR}/functions/DAlphaBall.gcc",
    )
)


@app.function(image=image, volumes={PARAMS_DIR: params_vol}, timeout=3600)
def download_weights():
    """Download AlphaFold2 weights (~5.3 GB) into the params Volume."""
    if os.path.exists(f"{PARAMS_DIR}/params_model_5_ptm.npz"):
        print("AF2 weights already present.")
        return
    tar = f"{PARAMS_DIR}/alphafold_params.tar"
    subprocess.run(
        ["wget", "-O", tar, "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"],
        check=True,
    )
    subprocess.run(["tar", "-xf", tar, "-C", PARAMS_DIR], check=True)
    os.remove(tar)
    params_vol.commit()
    print("Done:", sorted(os.listdir(PARAMS_DIR)))


@app.function(
    image=image,
    gpu="A100",
    volumes={PARAMS_DIR: params_vol, RESULTS_DIR: results_vol},
    timeout=24 * 3600,
)
def run_motifcraft(settings: dict, pdb_bytes: bytes, pdb_name: str, filters: dict, advanced: dict):
    work = Path("/tmp/job")
    work.mkdir(parents=True, exist_ok=True)

    design_path = Path(RESULTS_DIR) / settings["binder_name"]
    design_path.mkdir(parents=True, exist_ok=True)
    pdb_path = work / pdb_name
    pdb_path.write_bytes(pdb_bytes)

    # the example settings point at cluster paths; point them at this container instead
    settings = {**settings, "design_path": str(design_path), "starting_pdb": str(pdb_path)}
    advanced = {**advanced, "af_params_dir": PARAMS_DIR}

    paths = {}
    for name, cfg in [("settings", settings), ("filters", filters), ("advanced", advanced)]:
        paths[name] = work / f"{name}.json"
        paths[name].write_text(json.dumps(cfg, indent=2))

    try:
        subprocess.run(
            ["python", "-u", "motifcraft.py",
             "--settings", str(paths["settings"]),
             "--filters", str(paths["filters"]),
             "--advanced", str(paths["advanced"])],
            cwd=APP_DIR, check=True,
        )
    finally:
        results_vol.commit()  # keep partial results even if the run fails/times out


@app.local_entrypoint()
def main(
    settings: str,
    pdb: str,
    filters: str = "settings_filters/default_filters.json",
    advanced: str = "settings_advanced/default_4stage_multimer.json",
    gpu: str = "",
):
    load = lambda p: json.loads(Path(p).read_text())
    fn = run_motifcraft.with_options(gpu=gpu) if gpu else run_motifcraft
    # spawn returns immediately; with `modal run --detach` the job keeps running after this process exits
    call = fn.spawn(load(settings), Path(pdb).read_bytes(), Path(pdb).name, load(filters), load(advanced))

    out_dir = Path("outputs") / Path(pdb).stem
    out_dir.mkdir(parents=True, exist_ok=True)
    info = {
        "app_id": app.app_id,
        "function_call_id": call.object_id,
        "started": datetime.now().isoformat(timespec="seconds"),
        "settings": settings, "pdb": pdb, "filters": filters, "advanced": advanced,
        "gpu": gpu or "default",
    }
    (out_dir / "run_info.json").write_text(json.dumps(info, indent=2))
    print(f"Spawned. app_id={app.app_id} call_id={call.object_id}")
    print(f"Run info: {out_dir / 'run_info.json'}")
    print(f"Logs:     modal app logs {app.app_id} -f > experiments/{Path(pdb).stem}/logs/run.log  (create the folder first)")
    print(f"Stop:     modal app stop {app.app_id}")
