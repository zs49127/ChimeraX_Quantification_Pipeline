# ChimeraX Quantification Pipeline

Compare TFAP4 homodimer AlphaFold Server models with DrONv2 insertions after TFAP4 residues 76, 80, and 136 against one fixed control model. Four conditions × five actual random seeds × five ranked poses = **100 models**. Single-TFAP4 runs are excluded.

## Status and next step

The user confirmed successful processing of all 100 models on Longleaf with ChimeraX 1.12 in `results_corefix_v3`, including both alignment and plotting completion markers. Local tests separately cover numerical metrics, input validation, and synthetic plotting. The optional 3D plot has been tested locally; its cluster execution has not been confirmed. Inputs and scientific results are not committed here.

This study now lives in the repository's `20261028_TFAP4-DrON/` subdirectory. Run commands from that subdirectory, or copy its contents into your existing cluster project directory. Keep the relative `pipeline`, `scripts`, and `tests` folders intact and preserve your cluster-specific `config.json`. Existing cluster files and results do not need to move. Upload extracted AlphaFold job folders into:

```text
/path/to/project/Inputs_from_AF
```

The existing cluster project directory has a trailing space in its name; quote its exact path and use `Path.cwd()` when setting `project_root`. The repository directory has no trailing space. Inputs and generated results are not tracked by Git.

First, run this in a **Jupyter terminal** (Bash, not a Python notebook cell), from the directory containing this README:

```bash
bash scripts/inspect_cluster.sh /path/to/ChimeraX-1.12 /path/to/project
```

Send back the resulting `handoff/chimerax_environment_*.txt` report. It inventories the supplied installation, launcher candidates, libraries, loaded modules, and input directories. It does not change the installation or launch a computational job. The actual executable and job commands can then be verified against that report.

## Files

| File | Purpose |
|---|---|
| `config.json` | Project paths, allowed input folders, core and cassette length, reference override |
| `pipeline/prepare.py` | Validate four sequences, 20 distinct seeded jobs, all 100 poses; freeze input hashes |
| `pipeline/align_chimerax.py` | Joint-chain MatchMaker fits, atom deviations, confidence summaries, aligned structures |
| `pipeline/metrics.py` | NumPy coordinate metrics and symmetry diagnostic |
| `pipeline/plot_results.py` | CSV summaries and PNG/PDF figures |
| `scripts/inspect_cluster.sh` | Installation handoff report |
| `scripts/smoke_chimerax.py` | Small real-installation API and export check |
| `scripts/run_pipeline.sbatch` | Sequential CPU Slurm workflow; no GPU required |
| `tests/` | Numbering, metric, manifest, archive-inspection, and plotting checks |
| `METHODS.md` | Metric definitions, plotting rationale, interpretation and references |

## Input selection

`config.json` selects these folder patterns, each with seed labels 1–5:

```text
*_2xtfap4_2xre_seed[1-5]
*_2xtfap4_76_dronv2_2xre_seed[1-5]
*_2xtfap4_80_dronv2_2xre_seed[1-5]
*_2xtfap4_136_dronv2_2xre_seed[1-5]
```

Folders named `*_tfap4_2xre_seed*` are deliberately excluded. Folder selection is explicit, but **condition assignment is validated by sequence**, not trusted from the folder name. Each mutant sequence must equal the control after removing exactly 257 residues at exactly one of the three requested sites. The cassette must be identical across all mutant runs. All other TFAP4 residues must match the control exactly.

Each extracted job must contain its original `*_job_request.json` and all five corresponding `*_model_0.cif` … `*_model_4.cif`, `*_summary_confidences_N.json`, and `*_full_data_N.json`. MSAs and template CIFs are ignored. The job-request seed is used as the actual seed identifier; the filename's `seed1` is only a label. Duplicate actual seeds within a condition cause an error. DNA sequences are checked against the reference during atom correspondence validation.

## Alignment policy

The default reference is rank 0 of the numerically smallest actual seed in the control condition. This is deterministic and does not optimize the result. Set `reference_model_id` to another control ID from `inventory.csv` to change it, and use a new output directory.

All controls are processed before mutant models. Each is aligned independently to the same unchanged reference. The reference self-comparison is retained in raw tables for auditing and excluded from plotted summaries: control n=24, each insertion n=25.

| Condition | Construct residues used, in both A and B |
|---|---|
| TFAP4 Control | 1–179 |
| Insertion at Residue 76 | 1–76, 334–436 |
| Insertion at Residue 80 | 1–80, 338–436 |
| Insertion at Residue 136 | 1–136, 394–436 |

The original command-style selection was:

```text
mmaker #2/A,B:1-76,334-436 to #1/A,B:1-179 pairing ss alg nw matrix BLOSUM-62 ssFraction 0 cutoffDistance 2.0 showAlignment false
```

**The production pipeline now uses the MatchMaker Python API with compact core sequences instead of this command.** ChimeraX's command restriction retains excluded sequence positions as placeholders; in a residue-136 insertion this produced only 272 pairs. The new alignment callback constructs two 179-letter core sequences referencing the original residues, then uses MatchMaker's Needleman–Wunsch alignment and iterative fit. Coordinates, numbering and full-model sequences are unchanged. All 358 pairs are still required before pruning. `alignment_method=matchmaker_compact_core_v2` and `mm_full_ca` identify the corrected results.

`pairing ss` fits A→A and B→B in one rigid-body transformation. Default MatchMaker `pairing bb` would select a single best chain pair. Sequence-only scoring (`ssFraction 0`) is intentional because the homologous core sequence is identical: predicted secondary-structure changes must not change residue correspondence. The cassette, tail after native residue 179, and DNA do not influence the fit. The complete complex moves with that fit.

## Runtime setup (after handoff)

Set `project_root` in `config.json` to your absolute project directory before running. Public files intentionally contain placeholders for private cluster paths. Your local `cluster_config.json`, if supplied separately, can replace this example configuration.

ChimeraX runs its own Python and NumPy. Do **not** pip-install into ChimeraX. A separate Python ≥3.10 environment runs preparation and plotting. An environment setup example is:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -v
```

If Python is supplied by a Longleaf module, load the verified module before creating the environment. The Slurm script requires `CHIMERAX_BIN` and `ANALYSIS_PYTHON` to point to verified executables. `CODE_DIR` defaults to the submission directory; `CONFIG` defaults to its `config.json`. We will supply the exact smoke-test and submission commands after the handoff. Resource defaults are one CPU, 12 GB RAM, and four hours, with the site's default CPU partition; adjust if observed usage requires it.

The workflow creates a **new** output directory and refuses to overwrite an existing one. A failed run preserves partial outputs and logs; fix the cause and change `output_directory` (e.g. `results_run02`) before resubmitting. Completion requires both `alignment_complete.json` and `plotting_complete.json`. A partial CSV is not a complete result. Plotting alone can be rerun with `python pipeline/plot_results.py /absolute/path/to/results` after successful alignment.

## Outputs

All outputs go under the configured project root and output directory:

- `manifest.json`, `inventory.csv`: exact inputs, seed IDs, sequences, SHA-256 hashes, excluded requests, reference, parameters.
- `model_metrics.csv`: all 100 model rows, with core RMSDs, pruning coverage, per-chain metrics, threshold fractions, CA distance preservation, DNA displacement, AlphaFold confidence and symmetry diagnostics.
- `per_residue.csv`, `residue_summary.csv`: native TFAP4 numbering, construct numbering, Cα deviation, Cα pLDDT, local CA distance preservation, fit retention.
- `per_atom/<model_id>.csv`: homologous core and DNA heavy-atom deviations with atom names and pLDDT; no nonexistent control counterpart is assigned to DrON atoms.
- `condition_summary.csv`: mean, sample SD (`ddof=1`), median and counts across poses; `seed_means.csv` / `seed_summary.csv` separately describe variation between seed means.
- `confidence_regions.csv`, `confidence_region_summary.csv`: directional PAE and contact probabilities between chain-specific core, cassette, tail and DNA regions.
- `aligned/*.cif`: all complexes in the reference coordinate frame; `transforms.json` stores 3×4 scene transformations.
- `figures/`: nine figure types exported as PNG and PDF; plot labels use “Insertion at Residue 76”, etc., with line breaks on compact axes.
- `chimerax_stdout.log`, `chimerax_log.html`, `matchmaker_commands.txt`: runtime diagnostics and API-call descriptions with runtime model IDs; not a standalone reload script.

Read [METHODS.md](METHODS.md) before interpreting the plots. SD describes predicted-model spread, not experimental uncertainty or evidence of functional effects.

## Combined 3D residue plot

After alignment completes, run `python pipeline/plot_residue_3d.py /path/to/results`.
This standalone script reads the existing `per_residue.csv`, manifest and alignment
completion marker; no alignment rerun is needed. It creates a new `figures_3d/`
directory containing PNG, PDF and an underlying summary CSV. Use `--output-dir`
to choose another new destination on subsequent runs. Existing results are not changed.

Axes are native TFAP4 residue (X), mean Cα pLDDT (Y), and mean Cα deviation in Å (Z).
Each condition has a separate trajectory; A and B have separate panels with common
axis ranges. Points represent residue-level means across poses, with the reference
self-comparison excluded. Lines connect consecutive residues, not spatial contacts.
These are the same measurements as the original 2D profiles, not a new alignment.

For a rotatable, offline HTML version, install the optional Plotly dependency with
`python -m pip install 'plotly>=5.18,<7'` in your plotting environment, then append
`--html` to the command. Alternatively install `requirements-3d.txt` in a dedicated
environment. The HTML includes hover statistics, legend toggles and a mean ± SD
button. SD bars are marginal pose variation on each axis, not confidence intervals
or joint uncertainty regions, and are clipped to physically possible bounds.
Mean trajectories alone cannot establish an individual-model correlation between
confidence and displacement. The two original 2D profiles remain useful alongside
the 3D view because perspective can obscure overlap.
