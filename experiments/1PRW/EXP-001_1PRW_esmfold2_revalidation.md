# EXP-001: ESMFold2 re-validation of MotifCraft MPNN sequences (1PRW motif scaffolding)

- **Status:** ESMFold2 runs complete (520 of 521 sequences scored); analysis below. Follow-up (rerun MotifCraft with a looser motif cutoff) not started.
- **Started:** 2026-10-23
- **Related:** `MODAL_NOTES.md` (setup, issues, run history)

## 1. Background

A MotifCraft run on Modal scaffolding the 1PRW motif (chain A residues 1-40, two segments of 20) produced **0 accepted designs** after 14 h 37 min. The pipeline ran end to end, so the setup works; the question is why nothing passed the filters and whether the filters or the validation predictor are the cause.

MotifCraft run (app `ap-xDtCaxFEQPinJSQr6N4IGD`, log `experiments/1PRW/logs/run.log`):

| Stage | Count |
|---|---|
| Trajectories started | 201 |
| Rejected at start (initial pLDDT too low) | 113 |
| Completed all four design stages and reached MPNN | 26 |
| MPNN FASTAs in the results folder (26 x 20 = 520 from this run, plus 1 leftover from the earlier crashed run, which wrote to the same folder) | 521 |
| Sequences in `mpnn_design_stats.csv` (passed AF2 pLDDT/pTM at model 1) | 374 |
| Sequences not in the stats file | 147 |
| Designs accepted | 0 |

Why the 374 were rejected (partial filter set, model-1 values, `default_filters_partial.json`):

| Filter | Threshold | Fails (of 374) | Observed |
|---|---|---|---|
| `Hotspot_RMSD` (motif Cα RMSD) | <= 1 Å | 374 | best 1.25 Å, median 2.67 |
| `Target_RMSD` (motif Cα RMSD) | <= 1 Å | 359 | best 0.74 Å, median 1.71 |
| `Binder_RMSD` (AF2 prediction vs trajectory) | <= 2 Å | 255 | best 0.60 Å, median 2.98 |
| `Binder_pLDDT` | >= 0.7 | 101 | |

An earlier version of this analysis (and a message to the user) attributed the failures to AF2 pLDDT and pTM. That was incomplete: those filters explain only the 147 sequences outside the stats file. For the 374, confidence passes (median pLDDT 0.84, pTM 0.64) and the motif RMSD cutoffs are what reject them.

## 2. Goals and questions

1. Score **every** sequence that reached MPNN (521) with an independent predictor, ESMFold2, instead of stopping at AF2's verdict.
2. Does ESMFold2 agree with AF2 about **motif geometry**? Compare ESMFold2 motif Cα RMSD with AF2's `Target_RMSD` on the 374 sequences that have both.
3. Is the **1 Å motif RMSD cutoff too strict** (the user's working hypothesis), and where would a reasonable cutoff fall?
4. Does ESMFold2 rescue any of the 147 sequences that AF2 rejected on confidence?

**Hypotheses (stated before seeing the results):**
- H1: A meaningful number of sequences have ESMFold2 motif RMSD below 2 Å, so a looser cutoff would yield candidates.
- H2: AF2 validation was not the main bottleneck, since the AF2 confidence checks passed for 374 sequences. ESMFold2 will mostly confirm AF2's motif RMSD ranking rather than reverse it.

Not a goal: producing final designs. Sequences that look good here would still need to go through the interface-free filters in MotifCraft (or a relaxed filter set) before being called designs.

## 3. Inputs

| Input | Path |
|---|---|
| Reference structure | `targets/1PRW.pdb` (chain A, residues 1-40 = motif) |
| MPNN sequences with AF2 metrics (374) | `experiments/1PRW/data/mpnn_design_stats.csv` |
| All MPNN sequences as FASTA (521) | `experiments/1PRW/data/Sequences/` (copied from Volume `motifcraft-results`, `a_1PRW/MPNN/Sequences`) |
| Sequences not in stats (147) | `experiments/1PRW/data/mpnn_missing_seqs.csv` (derived: FASTA set minus stats file) |
| Motif position in each scaffold | parsed from `run.log` lines `backbone_constrained_residues (input / in scaffold)`; all 26 trajectories mapped, segment lengths equal |

Check: all 374 stats sequences are identical to their FASTA sequences.

## 4. Methods

- **Script:** `modal_esmfold2_batch.py` (Modal app `motifcraft-esmfold2`, A100).
- **Model:** `biohub/ESMFold2`, loaded once per job, monomer, **no MSA** (single sequence), `num_loops=20`, `num_sampling_steps=100`, 1 diffusion sample, `seed=0`. Image and calls follow `KrajTechne/Optimus` (`modal_common.py`, `RunESMFold2.py`), without importing it.
- **Motif RMSD:** Cα atoms of the motif residues in the ESMFold2 model are superimposed on the 1PRW motif Cα (Biotite `superimpose`, then `rmsd`), using the input-to-scaffold residue mapping from the log. The helper was checked locally: a rigidly moved copy gives 0.0 Å; 1 Å Gaussian noise per axis gives 1.6 Å.
- **Jobs:**
  - Test: 3 sequences (older gemmi-based RMSD, kept as `esmfold2_results_test3.csv`; not used in the analysis).
  - Job 1: the 374 stats-file sequences.
  - Job 2: the 147 remaining sequences (146 scored; one has no motif mapping).
- **Comparison with AF2:** joined on `Design`. AF2 columns used: `1_pLDDT`, `1_pTM`, `Average_pLDDT`, and (from the stats file) `1_Target_RMSD`, `1_Hotspot_RMSD`, `1_Binder_RMSD`.

## 5. Outputs

| Output | Location |
|---|---|
| ESMFold2 metrics, 374 sequences | `experiments/1PRW/data/esmfold2_results.csv` |
| ESMFold2 metrics, 147 sequences | `experiments/1PRW/data/esmfold2_results_missing147.csv` |
| Predicted structures (.cif) | Modal Volume `motifcraft-esmfold2`, folder `1PRW/` |
| Job logs | `experiments/1PRW/logs/esmfold2_full.log`, `experiments/1PRW/logs/esmfold2_missing.log` |

Columns added: `esm_pLDDT`, `esm_pTM`, `esm_motif_ca_rmsd`, `esm_error`.

## 6. Results

520 of 521 sequences were scored by ESMFold2 with no prediction errors. One sequence was skipped (`a_1PRW_l80_s578946_mpnn1`, a leftover from the earlier crashed run, with no motif mapping in the final run's log). Totals: 374 in the stats set, 146 in the missing set.

### 6.1 ESMFold2 scores

| | Stats set (n=374) | Missing set (n=146) |
|---|---|---|
| Motif Cα RMSD: min / median | 0.65 / 2.08 Å | 0.81 / 4.01 Å |
| Motif RMSD under 1 / 1.5 / 2 / 3 Å | 13 / 52 / 175 / 262 | 4 / 10 / 24 / 38 |
| ESMFold2 pLDDT: median (max) | 0.64 (0.78) | 0.58 (0.77) |
| ESMFold2 pTM: median | 0.48 | 0.36 |
| pLDDT >= 0.7 / pTM >= 0.5 | 84 / 168 | 3 / 7 |

### 6.2 ESMFold2 vs AF2 (stats set, n=374)

| Metric | AF2 | ESMFold2 | Correlation |
|---|---|---|---|
| pLDDT, median | 0.84 | 0.64 | Pearson 0.48 |
| pTM, median | 0.64 | 0.48 | Pearson 0.57 |
| Motif RMSD vs AF2 `Target_RMSD` | median 1.71 Å | median 2.08 Å | Spearman 0.49 (Pearson 0.33) |
| Motif RMSD vs AF2 `Hotspot_RMSD` | | | Spearman 0.50 (Pearson 0.55) |

Sequences under each motif RMSD cutoff:

| Cutoff | AF2 `Target_RMSD` | AF2 `Hotspot_RMSD` | ESMFold2 | ESMFold2 and AF2 `Target_RMSD` |
|---|---|---|---|---|
| < 1.0 Å | 14 | 0 | 13 | 8 |
| < 1.5 Å | 72 | 13 | 52 | 32 |
| < 2.0 Å | 280 | 135 | 175 | 166 |

18 of the 26 trajectories have at least one sequence with ESMFold2 motif RMSD under 2 Å.

### 6.3 What a looser motif cutoff would have accepted

Applying every other partial-filter threshold as in the run (AF2 pLDDT >= 0.7, pTM >= 0.5, `Binder_pLDDT` >= 0.7, loop% <= 90, `Binder_RMSD` <= 2 Å), 118 of the 374 pass. Then varying the cutoff on both `Target_RMSD` and `Hotspot_RMSD`:

| Motif cutoff | Designs passing | From trajectories | Also confirmed by ESMFold2 (pLDDT >= 0.7, pTM >= 0.5, motif <= 2 Å) |
|---|---|---|---|
| 1.0 Å (as run) | 0 | 0 | 0 |
| 1.25 Å | 1 | 1 | 1 |
| 1.5 Å | 8 | 3 | 4 |
| 1.75 Å | 69 | 8 | 43 |
| 2.0 Å | 104 | 12 | 47 |

The 8 designs at 1.5 Å:

| Design | ESM pLDDT | ESM pTM | ESM motif RMSD | AF2 Target | AF2 Hotspot | AF2 Binder_RMSD |
|---|---|---|---|---|---|---|
| a_1PRW_l101_s822755_mpnn12 | 0.74 | 0.71 | 0.65 | 0.96 | 1.25 | 1.99 |
| a_1PRW_l83_s780823_mpnn5 | 0.73 | 0.62 | 1.44 | 1.41 | 1.47 | 1.66 |
| a_1PRW_l83_s780823_mpnn13 | 0.72 | 0.61 | 1.48 | 1.41 | 1.41 | 1.92 |
| a_1PRW_l83_s780823_mpnn20 | 0.73 | 0.63 | 1.57 | 1.44 | 1.44 | 1.55 |
| a_1PRW_l100_s496913_mpnn5 | 0.65 | 0.57 | 2.53 | 1.42 | 1.43 | 1.29 |
| a_1PRW_l100_s496913_mpnn7 | 0.66 | 0.58 | 2.02 | 1.44 | 1.49 | 1.36 |
| a_1PRW_l100_s496913_mpnn10 | 0.62 | 0.51 | 1.49 | 1.36 | 1.43 | 1.53 |
| a_1PRW_l100_s496913_mpnn18 | 0.68 | 0.59 | 1.65 | 1.43 | 1.48 | 1.29 |

The strongest single sequence is `a_1PRW_l101_s822755_mpnn12`: both predictors put its motif under 1 Å (ESMFold2 0.65, AF2 `Target_RMSD` 0.96), with AF2 pLDDT 0.88 and ESMFold2 pLDDT 0.74 / pTM 0.71. It misses the original filter set only on `Hotspot_RMSD` (1.25 vs the 1 Å cutoff).

### 6.4 The 146 sequences AF2 rejected on confidence

ESMFold2 does not rescue them: median pTM 0.36, only 3 with pLDDT >= 0.7 and 7 with pTM >= 0.5. 24 have motif RMSD under 2 Å, but with low confidence.

## 7. Discussion

- **Goal 3 (is 1 Å too strict): supported.** No sequence passes at 1 Å on both measures, the best `Hotspot_RMSD` is 1.25 Å, and loosening to 1.5 Å yields 8 designs. The cutoff sits just below where the distribution starts, and the run's rejections were geometry, not confidence.
- **Count designs by trajectory, not by sequence.** Sequences from one trajectory share a scaffold and are not independent. The 8 designs at 1.5 Å come from 3 trajectories and the 69 at 1.75 Å from 8, so the effective number of distinct scaffolds is small.
- **H1 (some sequences under 2 Å by ESMFold2): supported.** 175 of 374.
- **H2 (ESMFold2 mostly confirms AF2's ranking): only weakly supported.** The rank correlation on motif RMSD is about 0.5. ESMFold2 is also much less confident (median pLDDT 0.64 vs 0.84). Agreement is best at the top: 8 sequences are under 1 Å on both ESMFold2 and AF2 `Target_RMSD`, and 32 under 1.5 Å. I did not test why the two predictors differ; lower single-sequence confidence for designed proteins is plausible but unverified.
- **Goal 4 (rescue the 146): no.**
- **A second bottleneck exists.** `Binder_RMSD` (AF2 prediction vs the design trajectory) <= 2 Å fails for 255 of 374. Some sequences with excellent motifs have a very different overall fold from the design, e.g. `a_1PRW_l101_s822755_mpnn17` (ESMFold2 motif 0.67 Å, AF2 `Binder_RMSD` 8.31 Å). Whether that cutoff is also too strict was not examined.
- **`Target_RMSD` and `Hotspot_RMSD` are two measurements, not one.** In `motifcraft.py` both are aligned motif RMSDs for a motif-only run, but computed on different AF2 predictions (the design model and the binder-only model). Both are held to the same 1 Å cutoff.

## 8. Caveats

- The Modal image falls back to pure-PyTorch kernels (no Transformer Engine, xformers or flash-attn). The model warns outputs "will differ numerically" from the reference implementation. Treat small differences as noise.
- Motif RMSD here is Cα only, computed by my own script, not by MotifCraft's `aligned_motif_rmsd`. The two should be similar but are not guaranteed identical; I have not compared them directly.
- ESMFold2 numbers are from one seed and one diffusion sample per sequence.
- The motif mapping assumes scaffold segments appear in the same order as the input segments.
- The MotifCraft run used AF2 (not ESMFold) for its own validation because ESMFold is not installed in the MotifCraft image; the authors' example uses ESMFold validation.

## 9. Next steps

1. Look at the 8 candidates at 1.5 Å (start with `a_1PRW_l101_s822755_mpnn12`): inspect the AF2 and ESMFold2 structures, and run the relaxation and Rosetta scoring that MotifCraft skips for rejected designs. The structures are in Volumes `motifcraft-results` (AF2) and `motifcraft-esmfold2` (ESMFold2).
2. Rerun MotifCraft with a copy of `default_filters_partial.json` in `configs/` with motif cutoffs of about 1.5 Å (and consider `Binder_RMSD`), rather than editing the upstream file. Decide the cutoff with the trajectory-level count in mind.
3. Optionally add ESMFold-based validation to the MotifCraft image to match the authors' setup.
4. Record the conclusion and any new run IDs in `MODAL_NOTES.md`.
