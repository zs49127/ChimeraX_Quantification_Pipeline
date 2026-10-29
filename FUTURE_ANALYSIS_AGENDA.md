# Agenda for future AlphaFold → ChimeraX analyses

This is the reusable planning document learned from the TFAP4–DrON study. It is
intended to be read by future chats working in this repository. It is a guide for
adapting the workflow, not a claim that the TFAP4 scripts already support arbitrary
proteins without modification.

## 1. Define the biological comparison

Record the protein, construct sequences, oligomeric state, chain roles, ligands or
nucleic acids, conditions, actual seeds and poses per seed. Clarify whether the
question concerns overall structure, a conserved domain, an interface, or a local
insertion. Distinguish the number of conditions from the number of predicted models.

Ask only for information that cannot be obtained from the supplied files or this
chat. Use job requests to recover sequences and actual random seed values; a folder
label such as `seed1` is not the random seed itself. Exclude unrelated monomer runs
from a homodimer study unless a separate comparison is requested.

For a new protein, explicitly decide:

- The reference construct/model and whether it is predicted or experimental.
- Which native residues define the alignment core and why other regions are excluded.
- Exact insertion/deletion mapping, linker lengths and any sequence substitutions.
- Chain correspondence, including symmetry or interchangeable subunits.
- Primary structural endpoints and confidence measurements.

## 2. Create an isolated project directory

Use a dated folder, such as `YYYYMMDD_PROTEIN-CONSTRUCT/`, inside this repository.
Keep code, methods, tests and example configuration together. Maintain a separate
project root on the cluster containing inputs and uniquely named result directories.

The TFAP4 code has assumptions in both configuration and Python: A/B protein chains,
C/D DNA, exact 257-residue insertions, core end 179, three sites, four labels, five
seeds, five poses and AF Server filename conventions. Audit `common.py`, `prepare.py`,
`align_chimerax.py`, `plot_results.py` and `plot_residue_3d.py` when adapting them.
Do not change only `config.json` and assume everything else is generic.

Keep real cluster paths in a local configuration. Do not replace an existing working
cluster config with the repository's placeholder `/path/to/project` example.

## 3. Inspect the environment with a lightweight handoff

The proven installation was ChimeraX **1.12**, with a Linux executable at this
relative path within the extracted installation:

```text
<CHIMERAX_INSTALL>/usr/libexec/UCSF-ChimeraX/bin/ChimeraX
```

It has bundled Python 3.11. Run analysis code through the ChimeraX executable; do
not pip-install into that bundled interpreter. Separate plotting/preparation used
the Jupyter environment's Python 3.11.7 (`anaconda/2024.02.ood`). Recheck availability
in a new session rather than assuming the same module remains loaded.

`CLUSTER_HANDOFF.local.md` in the local workspace holds exact previously verified
locations and the trailing-space detail. If that file is absent, ask for the current
installation/project paths or a fresh environment report. Do not guess them.

Use `scripts/inspect_cluster.sh INSTALL PROJECT` to collect the launcher, library
resolution, Python, loaded modules, partitions and input inventory. First check a
suspect path with `ls`; quote it. A directory-name trailing space caused the original
`cd` failure. The repository directory intentionally has **no trailing space**.

## 4. Validate input data before alignment

Recognize AlphaFold Server model CIFs, summary-confidence JSON, full-data JSON and
job-request JSON. MSAs and template CIFs are not additional predicted poses. Record
the file inventory and SHA-256 hashes. Check counts, duplicate seeds, chain roles,
sequences, residue numbering and required companion files.

For exact cassette comparisons, deleting the specified inserted segment from each
construct must recover the reference sequence. If it does not, investigate linkers,
mutations, deletions or alternate isoforms rather than shifting residue numbers by
guesswork. Define an explicit map for more general changes.

## 5. Align the intended biological unit

Use a single fixed reference for the primary analysis and store its model ID.
TFAP4 used a deterministic control seed/rank reference; it did not select the
reference to make insertion results look favorable. Preserve that identity on reruns.

For a homodimer, fit both chains together using specific pairing, A→A and B→B.
Default MatchMaker best-chain pairing does not enforce a joint dimer fit. Transform
the whole complex, including DNA and excluded regions, as a rigid body.

**Lesson from the residue-136 failure:** restricted sequence specifications can
leave excluded residues as `?` placeholders. A long insertion then disrupted the
post-insertion alignment: only 272 of 358 expected Cα pairs were paired. The successful
fix creates compact core `StructureSeq` objects referencing original residues and
passes them to MatchMaker's sequence-alignment callback. The whole model keeps its
original numbering, sequences and coordinates. Require the exact expected atom
correspondence before pruning; a small RMSD alone does not validate the mapping.

Use `ss_fraction=0.0`, not `0`: the latter failed in the native alignment API even
though mocked tests passed. Exercise the real installed API before a full job.

## 6. Separate fitting from measurement

TFAP4 used a 2 Å iterative cutoff. MatchMaker removes worst-fitting pairs gradually
and refits until retained pairs meet the threshold. Pruned atoms are not deleted.

Always distinguish:

- **Pruned RMSD:** agreement of the retained fitting subset.
- **Retained fraction/count:** how much of the intended core supports that fit.
- **All-core RMSD:** all intended pairs evaluated after the same pruned transform.
- **Unpruned Kabsch RMSD:** best overall fit using every core pair (`kabsch_identity_rmsd`).
- **Per-atom/per-residue deviation:** spatial displacement after the chosen fit.

A low pruned RMSD does not establish preservation of the entire core. The same
per-residue measurements can change when a different fit defines the coordinate
frame. If unpruned per-residue profiles are requested, explicitly calculate them:
the existing TFAP4 profiles use the pruned transform.

Use native protein numbering for comparable residue plots. Cα/backbone profiles
are less sensitive to side-chain conformation and symmetric atom naming than all
heavy-atom RMSD. Label DNA displacement as measured **after the protein fit**.
Homodimer/DNA label swapping can affect these values; retain symmetry diagnostics.

## 7. Keep confidence separate from structural agreement

- Cα pLDDT is AlphaFold local confidence (0–100), not observed displacement.
- pTM/ipTM describe predicted structural/interface confidence, not affinity.
- PAE reports directional predicted uncertainty between positions or regions, not
  deviation from the chosen control.
- The implemented CA-only distance-preservation score compares internal Cα distances
  with reference neighbors within 15 Å at 0.5/1/2/4 Å thresholds. It is not pLDDT and
  not the full all-atom lDDT implementation. It is invariant to rigid-body movement.

Do not assign cassette-to-control RMSD when the cassette lacks a homologous control
counterpart. Use its confidence and relative placement uncertainty instead.

## 8. Plot variation honestly

Show individual poses and seed means alongside condition mean ± sample SD. Poses
within a seed are nested computational samples, not independent biological
replicates. Do not treat atoms or chains as independent replicates. Avoid automatic
significance tests. Exclude the reference self-row from summaries but retain it in
raw tables; this left 24 controls versus 25 models per insertion in TFAP4.

Useful views: global metrics with fit coverage, residue deviation/confidence profiles,
model-by-residue heatmaps, heavy-atom ECDFs, region PAE, and confidence versus deviation.
The combined 3D view uses native residue, mean pLDDT and mean deviation; lines only,
tight shared A/B ranges, optional SD and offline HTML. Perspective can obscure
differences, so keep the 2D plots. Mean trajectories do not establish individual-model
correlation. Schematic plots must be labeled illustrative, not measured results.

## 9. Test, submit, and diagnose without losing results

Run mapping/metric tests and a real ChimeraX smoke test including the problematic
construct. Use one CPU initially; the tested workflow does not require GPU alignment.
The existing Slurm template requests 12 GB and four hours, but tune resources to the
new project's size and observed usage. Export `CHIMERAX_BIN`, `ANALYSIS_PYTHON`,
`CODE_DIR` and `CONFIG`; submit from the project's code directory.

A disappeared queue entry is not proof of success. Check the Slurm output plus
`alignment_complete.json` and `plotting_complete.json`. The alignment traceback is
in `chimerax_stdout.log`, not necessarily the Slurm log. Inspect its last lines first.

For a failed run, preserve partial outputs, choose a new output directory and keep
the original reference. Replot completed CSVs without rerunning alignment. Prefer
commit-pinned downloads for code updates; current main-branch URLs now include the
project subdirectory. Do not overwrite a private working config during downloads.

## 10. Leave a useful handoff

Update each project's README/VALIDATION and local cluster notes with: biological
setup, sequence mapping, reference, method/commit, executable/environment, commands,
output directory, job ID, completion evidence, unresolved limitations and next steps.
Document user-reported cluster success separately from locally executed tests.
For a new chat, read these files first; they are the durable record of this workflow.
