# Validation record

## Compact-core correction (2026-10-29)

The cluster log revealed that the first residue-136 model yielded 272 pre-pruning pairs. The corrected callback removes excluded sequence placeholders before sequence alignment while retaining original residue objects. Four regression tests cover all insertion sites, joint-chain API wiring, rejection of the original 272-pair failure, and missing post-insertion residues. These use test doubles for ChimeraX objects; they do not claim a local ChimeraX integration run. The corrected numerical fit must be verified on Longleaf. Rerun into a new directory with the prior manifest's reference ID; keep the original results untouched.

Validated locally on 2026-10-28:

- Five unit tests pass: all insertion boundaries; exact sequence/cassette validation; rigid-transform invariance and reflection rejection; distortion sensitivity; a complete 100-model manifest and missing-input/overwrite rejection.
- Python files compile successfully.
- Both Bash scripts pass `bash -n`.
- The plotting smoke test uses **synthetic** 100-model tables, exercises all nine figure types (18 PNG/PDF files), and verifies 24 nonreference controls and 25 models per insertion in summaries. The multipanel figure was visually inspected. Synthetic results are excluded from the repository.
- Plotting environment used Python 3.12, NumPy 2.5.4, pandas 2.3.3 and Matplotlib 3.11.2. Requirements allow compatible versions for the cluster.

The provided AlphaFold archive's file listing, request and model-0 confidence arrays were inspected. The archive became unavailable at its original local path before the standalone archive-inspection utility could be rerun; no full CIF parser or b-factor cross-check against that archive is claimed in this validation record.

Subsequent cluster confirmation (2026-10-29): the user supplied successful ChimeraX 1.12 smoke-test output, followed by an end-to-end run validating all 100 models and producing both `alignment_complete.json` and `plotting_complete.json` in `results_corefix_v3`. This supersedes the earlier pending cluster validation. The production runner checks every expected residue pair and verifies MatchMaker's full RMSD against direct coordinate calculations. The optional 3D plot has local test coverage but no confirmed cluster run. Repository reorganization changes file locations and documentation, not analysis code.
