# MotifCraft on Modal: notes

Running log of setup, gotchas and runs. Newest run at the bottom of the run history.

## Layout

```
Trial_MotifCraft/
├── MotifCraft/              git submodule of oceanefollonier/MotifCraft (pinned commit; keep pristine)
├── modal_app.py             MotifCraft on Modal
├── modal_esmfold2_batch.py  ESMFold2 re-validation of MPNN sequences
├── configs/                 our settings and preset copies
├── targets/                 input PDBs (1PRW, 6VW1)
├── outputs/<PDB ID>/        scratch (gitignored): run_info.json written at each launch
├── experiments/<PDB ID>/    one folder per experiment: write-up .md, data/ (CSVs, run_info), logs/ (gitignored)
└── MODAL_NOTES.md           this file
```

Fresh clone: `git clone --recurse-submodules <repo>`; MotifCraft's own `.gitmodules` points ColabDesign at SSH, so if that fails run
`git -C MotifCraft -c url."https://github.com/".insteadOf=git@github.com: submodule update --init`.

## Setup

- `modal_app.py` builds the image (micromamba, Python 3.10, JAX+CUDA, PyRosetta, custom ColabDesign, ffmpeg) and defines:
  - `download_weights`: one-time AF2 weights (~5.3 GB) into Volume `motifcraft-params` (done).
  - `run_motifcraft`: GPU job (A100, 24 h timeout). Overwrites `design_path`/`starting_pdb` in the settings with container paths. Results go to Volume `motifcraft-results/<binder_name>`.
- Windows: set `PYTHONUTF8=1 PYTHONIOENCODING=utf-8` before running `modal`, or the CLI crashes on a console encoding error.
- The image dedents the broken line in the ColabDesign copy with `sed` and checks it with `py_compile` before installing it (see issue 1).

## Launch a motif-only run (e.g. 1PRW)

Run from the repo root:

```
modal run --detach modal_app.py \
  --settings configs/a_1PRW_motif_modal.json \
  --pdb targets/1PRW.pdb \
  --filters MotifCraft/settings_filters/default_filters_partial.json \
  --advanced configs/default_4stage_multimer_partial_af2.json
```

- `--detach` is a CLI flag (cannot be set in the file). The entrypoint uses `.spawn()` and writes `outputs/<PDB>/run_info.json` (scratch; copy it into `experiments/<PDB>/data/` for a run you want to keep).
- Capture the log: `modal app logs <app_id> -f > experiments/<PDB>/logs/run.log` (without `-f` it returns only the last 100 lines).
- Stop: `modal app stop <app_id>`. Fetch results: `modal volume get motifcraft-results a_1PRW ./results`.

## Config rules (these caused failures)

| Run type | Advanced preset | Filters |
|---|---|---|
| Motif only (no target) | `default_4stage_multimer_partial*.json` | `default_filters_partial.json` |
| Binder with target | `default_4stage_multimer.json` | `default_filters.json` |

- `default_4stage_multimer_partial.json` sets `use_esmfold_prediction: true` with a cluster `esm-fold` path. ESMFold is not in the Modal image, so `configs/default_4stage_multimer_partial_af2.json` is a copy with it set to `false` (AF2 validates instead).
- The paths in the upstream example settings (`/scicore/...`) are cluster paths. They are ignored on Modal; `configs/a_1PRW_motif_modal.json` is a copy with local-looking paths and `number_of_final_designs: 10`.

## Issues hit and fixes

| # | Symptom | Cause | Fix |
|---|---|---|---|
| 1 | `IndentationError` in `colabdesign/af/inputs.py:19` | Bug in the pinned ColabDesign commit (also in latest upstream I could see) | One-line dedent applied with `sed` during the image build, followed by `py_compile` so the build fails if it is still broken (the submodule itself stays pristine). A patch file did not work: the Windows checkout has CRLF line endings and `patch` rejected the hunk |
| 2 | `Incompatible shapes ... (86,3),(40,3)` in `pseudo_beta_fn` | Used binder preset with `predict_initial_guess: true`; no target, so initial coordinates were motif-length only | Use the `_partial` advanced preset |
| 3 | `MovieWriter (ffmpeg) not available` | `ffmpeg` missing from image | `.apt_install("ffmpeg")` |
| 4 | `KeyError: 'i_pae'` after MPNN | Binder filters check interface metrics that don't exist without a target | Use `default_filters_partial.json` |

Missing from the repo: `functions/dssp` and `functions/DAlphaBall.gcc` are gitignored. The image fetches them from BindCraft.

## Run history (a_1PRW)

| Run | Outcome |
|---|---|
| 1 | Failed: issue 1 |
| 2 | Failed: issue 2 |
| 3 | Design stage ran (motif RMSD ~6 to ~1.2 Å); failed: issue 3 |
| 4 | Ran fine; stopped by user from dashboard |
| 5 (`ap-p1k4cIsPHKjaqw76hq5mRx`) | Detached. Trajectory passed (pLDDT 0.87), MPNN ran; failed: issue 4 |
| 6 (`ap-xDtCaxFEQPinJSQr6N4IGD`) | Detached, partial filters. Ran 14 h 37 min, 201 trajectories, then stopped itself (rejection check). **0 accepted designs.** Pipeline ran end to end with no crash. 113 trajectories rejected at the start (pLDDT too low), 26 finished all stages and reached MPNN (520 sequences). 374 of them passed the AF2 pLDDT/pTM cutoffs but all failed the motif RMSD cutoffs (`Hotspot_RMSD`/`Target_RMSD` <= 1 Å; best seen 1.25 / 0.74 Å); the other ~146 failed pLDDT/pTM. (An earlier version of this row blamed pLDDT/pTM for everything; that was wrong.) Log: `experiments/1PRW/logs/run.log`. Follow-up: `experiments/1PRW/EXP-001_1PRW_esmfold2_revalidation.md` |

## Open questions / ideas

- Does the pipeline count designs already in `a_1PRW/Accepted` toward `number_of_final_designs`? Not confirmed.
- ESMFold validation is not set up in the image.
- The `inputs.py` indentation bug could be reported to the MotifCraft authors.
- Restructure done 2026-10-23 (MotifCraft as a pinned submodule, our files at the repo root). No commits were made as part of it; `.gitmodules` and the submodule pointer are staged, everything else is untracked.
