#!/usr/bin/env bash
# Read-only inspection of the supplied installation; writes only a handoff report.
set -uo pipefail
CX_DIR="${1:?Usage: bash scripts/inspect_cluster.sh /path/to/ChimeraX /path/to/project}"
PROJECT="${2:?Supply the absolute project directory as the second argument}"
mkdir -p "$PROJECT/handoff"
REPORT="$PROJECT/handoff/chimerax_environment_$(date +%Y%m%d_%H%M%S).txt"
{
  printf 'UTC: '; date -u
  printf '\nOS / architecture\n'; uname -a
  printf '\nInstallation root\n'; ls -ld "$CX_DIR"
  printf '\nTop-level listing\n'; ls -la "$CX_DIR"
  printf '\nDirectory tree (depth 3)\n'; find "$CX_DIR" -maxdepth 3 -type d -print
  printf '\nExecutable candidates and symlinks (depth 5)\n'
  find "$CX_DIR" -maxdepth 5 \( -iname '*chimerax*' -o -name 'python*' \) -print
  printf '\nLauncher candidates\n'
  while IFS= read -r -d '' candidate; do
    ls -l "$candidate"
    file "$candidate"
    resolved="$(readlink -f "$candidate")"
    if file "$resolved" | grep -q ELF; then ldd "$resolved" 2>&1; fi
  done < <(find "$CX_DIR" -maxdepth 5 -type f -iname '*chimerax*' -executable -print0)
  printf '\nAvailable commands\n'; command -v python3 python sbatch srun module || true
  printf '\nPython version\n'; python3 --version 2>&1 || true
  printf '\nLoaded modules\n'; module list 2>&1 || true
  printf '\nSlurm partitions\n'; sinfo -o '%P %a %l %D' 2>&1 || true
  printf '\nProject input directories\n'
  find "$PROJECT/Inputs_from_AF" -maxdepth 2 -type d -print 2>&1
  printf '\nInput model count (excludes template CIFs)\n'
  find "$PROJECT/Inputs_from_AF" -type f -name '*_model_*.cif' | wc -l
} > "$REPORT" 2>&1
printf 'Handoff report: %s\n' "$REPORT"
