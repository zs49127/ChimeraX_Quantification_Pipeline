"""Cluster preflight: run inside ChimeraX with one real CIF as argument.

Makes two copies and fits their identical A/B 1:76 coordinates, checking the
MatchMaker return API and mmCIF export. Writes beside the supplied output prefix.
Usage: chimerax --nogui --exit --script 'scripts/smoke_chimerax.py input.cif output_prefix'
"""
import sys
from pathlib import Path
from chimerax.core.commands import run, StringArg

path, prefix = sys.argv[1:3]
prefix = Path(prefix)
prefix.parent.mkdir(parents=True, exist_ok=True)
ref = run(session, "open " + StringArg.unparse(path))[0]
mov = run(session, "open " + StringArg.unparse(path))[0]
result = run(session, f"mmaker #{mov.id_string}/A,B:1-76 to #{ref.id_string}/A,B:1-76 "
                     "pairing ss ssFraction 0 cutoffDistance 2.0 showAlignment false")
assert len(result) == 1
assert len(result[0]["full match atoms"]) == 152
assert result[0]["full RMSD"] < 1e-4
run(session, f"save {StringArg.unparse(str(prefix) + '.cif')} models #{mov.id_string} relModel #{ref.id_string}")
Path(str(prefix) + '.passed.txt').write_text("ChimeraX joint-chain fit and mmCIF export passed.\n")
print("CHIMERAX_SMOKE_TEST_PASSED")
