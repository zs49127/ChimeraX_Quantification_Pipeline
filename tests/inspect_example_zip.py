"""Read an example archive without extracting it or executing archive contents.

Usage: python tests/inspect_example_zip.py /path/to/example.zip
Uses the standard library only. Outputs a compact factual inventory to stdout.
"""
import collections
import json
import shlex
import sys
import zipfile


def inspect(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        requests = [n for n in names if n.endswith("_job_request.json")]
        request = json.loads(z.read(requests[0]))[0]
        summary_name = next(n for n in names if n.endswith("_summary_confidences_0.json"))
        summary = json.loads(z.read(summary_name))
        full = json.loads(z.read(next(n for n in names if n.endswith("_full_data_0.json"))))
        cif = z.read(next(n for n in names if n.endswith("_model_0.cif"))).decode().splitlines()
        fields, rows, in_atoms = [], [], False
        for line in cif:
            if line.startswith("_atom_site."):
                fields.append(line.strip())
                in_atoms = True
            elif in_atoms and line.startswith("#"):
                break
            elif in_atoms and line.strip():
                rows.append(dict(zip(fields, shlex.split(line))))
        chains = collections.Counter(r["_atom_site.label_asym_id"] for r in rows)
        residue_counts = {c: len({r["_atom_site.label_seq_id"] for r in rows if r["_atom_site.label_asym_id"] == c}) for c in chains}
        assert len(rows) == len(full["atom_plddts"])
        for atom, chain, plddt in zip(rows, full["atom_chain_ids"], full["atom_plddts"]):
            assert atom["_atom_site.label_asym_id"] == chain
            assert abs(float(atom["_atom_site.B_iso_or_equiv"]) - plddt) <= 0.011
        result = dict(job=request["name"], actual_seed=request["modelSeeds"],
                      archive_files=len(names), model_cifs=sum("_model_" in n and n.endswith(".cif") for n in names),
                      atom_counts=dict(chains), residue_counts=residue_counts,
                      full_data_shapes={k: [len(v), len(v[0])] if isinstance(v[0], list) else [len(v)] for k, v in full.items()},
                      model0_confidence={k: v for k, v in summary.items() if isinstance(v, (int, float))},
                      cif_bfactor_matches_full_data_plddt=True)
        print(json.dumps(result, indent=2))
        return result


if __name__ == "__main__":
    inspect(sys.argv[1])
