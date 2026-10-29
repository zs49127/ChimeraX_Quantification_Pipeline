#!/usr/bin/env bash
#SBATCH --job-name=tfap4-quant
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=tfap4-%j.out
# Submit from the repository root, after environment handoff and setup.
set -euo pipefail
: "${CHIMERAX_BIN:?Set CHIMERAX_BIN to the verified executable path}"
: "${ANALYSIS_PYTHON:?Set ANALYSIS_PYTHON to the plotting environment python executable}"
CODE_DIR="${CODE_DIR:-$SLURM_SUBMIT_DIR}"
CONFIG="${CONFIG:-$CODE_DIR/config.json}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -x "$CHIMERAX_BIN"
test -x "$ANALYSIS_PYTHON"
cd "$CODE_DIR"
"$ANALYSIS_PYTHON" "$CODE_DIR/pipeline/prepare.py" --config "$CONFIG"
OUT="$("$ANALYSIS_PYTHON" -c 'import json,pathlib,sys; c=json.load(open(sys.argv[1])); print(pathlib.Path(c["project_root"])/c["output_directory"])' "$CONFIG")"
export MPLCONFIGDIR="$OUT/.matplotlib"
"$CHIMERAX_BIN" --nogui --exit --script "\"$CODE_DIR/pipeline/align_chimerax.py\" \"$OUT/manifest.json\"" > "$OUT/chimerax_stdout.log" 2>&1
test -f "$OUT/alignment_complete.json"
"$ANALYSIS_PYTHON" "$CODE_DIR/pipeline/plot_results.py" "$OUT"
test -f "$OUT/plotting_complete.json"
printf 'Pipeline complete: %s\n' "$OUT"
