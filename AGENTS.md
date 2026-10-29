# Guidance for future coding agents

Read `FUTURE_ANALYSIS_AGENDA.md` before planning another protein analysis. Read the
selected project's README, METHODS and VALIDATION files before editing it. If the
local-only `CLUSTER_HANDOFF.local.md` exists, consult it for previously verified
cluster paths and run history; do not include its private contents in GitHub commits.

- Keep each protein study in its own dated project directory. Preserve existing
  study code and results unless the user explicitly requests changes.
- Treat TFAP4-specific residues, chain IDs, insertion lengths, labels, file patterns
  and model counts as project settings, not universal biological assumptions.
- Use the latest corrected compact-core MatchMaker approach when known exact
  correspondence is appropriate. Never restore placeholder-based restricted-chain
  sequence alignment without checking the actual matched atom pairs.
- ChimeraX's native alignment API requires `ss_fraction=0.0`, not integer `0`.
  Python scripts execute in a ChimeraX sandbox; a `__main__` guard alone skips work.
- Keep pre-pruning correspondence validation, independent full-RMSD verification,
  input hashes, fixed reference identity and completion markers.
- Before another full cluster run after changing alignment behavior, provide a
  small actual-ChimeraX integration check that exercises both a control and an
  insertion/deletion model. Test doubles are not a native API integration test.
- Generate new result directories and separate rerun configurations. Never delete
  partial results simply to satisfy the pipeline's overwrite protection.
- Quote every shell path, preserve trailing spaces, and guard failed `cd` commands.
  Prefer a direct `ls` or short diagnostic over broad recursive searches.
- Keep user-facing cluster commands short, copy-pastable and explicit about their
  working directory. Read logs before suggesting another submission.
- The user prefers clean line-only 3D profiles with tight, shared scales and
  optional interactive HTML. Avoid adding point markers without a reason.
- Report what was locally tested, what was actually run on the cluster, and what
  remains unverified. Do not claim that assistant access to the cluster exists.
- Update the selected project's methods/handoff notes after a meaningful change
  so future chats can recover decisions from files rather than inferred memory.
