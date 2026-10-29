# ChimeraX Quantification Pipeline

Project-specific AlphaFold Server → ChimeraX structural comparison workflows.

| Directory / document | Purpose |
|---|---|
| [20261028_TFAP4-DrON](20261028_TFAP4-DrON/README.md) | TFAP4 homodimer with DrON insertions; code, configuration, methods and tests |
| [FUTURE_ANALYSIS_AGENDA.md](FUTURE_ANALYSIS_AGENDA.md) | Reusable planning and execution guide for new proteins |
| [AGENTS.md](AGENTS.md) | Instructions for coding agents continuing this repository |

## Repository layout

```text
ChimeraX_Quantification_Pipeline/
├── README.md
├── AGENTS.md
├── FUTURE_ANALYSIS_AGENDA.md
└── 20261028_TFAP4-DrON/
    ├── README.md
    ├── METHODS.md
    ├── VALIDATION.md
    ├── config.json
    ├── requirements.txt
    ├── requirements-3d.txt
    ├── pipeline/
    ├── scripts/
    └── tests/
```

Run project commands **from inside the project directory**. For example:

```bash
cd 20261028_TFAP4-DrON
python -m unittest discover -s tests -v
```

This repository reorganization does not move anything on Longleaf or invalidate
existing results. If updating an existing cluster project, copy the **contents**
of the matching repository project directory into its existing code location.
Do not overwrite its working cluster configuration with the placeholder example.
Raw GitHub URLs on the new main branch now include `20261028_TFAP4-DrON/` before
`pipeline/`, `scripts/`, etc.; older commit-pinned URLs retain their original paths.

New protein analyses should receive their own dated project directory rather than
overwriting the TFAP4 study. Large inputs, generated results and environment reports
remain on the cluster. A git-ignored `CLUSTER_HANDOFF.local.md`, when present in the
local workspace, records private paths for continuity without publishing them.
