"""Run INSIDE ChimeraX: chimerax --nogui --exit --script '.../align_chimerax.py manifest.json'."""
import csv
import json
import sys
import traceback
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_json, write_json, write_csv, sha256, construct_number, core_spec
from metrics import summarize, distances, rmsd, best_fit_rmsd


def atom_key(a):
    return a.residue.chain_id, a.residue.number, a.name


def index_atoms(model):
    atoms = {}
    for atom in model.atoms:
        key = atom_key(atom)
        if key in atoms:
            raise ValueError(f"Duplicate atom key / alternate location: {key}")
        if atom.residue.insertion_code.strip():
            raise ValueError("Insertion codes are unsupported; expected AF Server numbering")
        atoms[key] = atom
    return atoms


def validate_model(model, entry):
    chains = {c.chain_id: c for c in model.chains}
    if set(chains) != set("ABCD"):
        raise ValueError(f"Expected exactly chains A,B,C,D: {list(chains)}")
    for chain_id, sequence in zip("ABCD", entry["sequences"]):
        chain = chains[chain_id]
        if chain.characters != sequence:
            raise ValueError(f"CIF/request sequence mismatch in {entry['model_id']} /{chain_id}")
        residues = [r for r in chain.residues if r is not None]
        if [r.number for r in residues] != list(range(1, len(sequence) + 1)):
            raise ValueError(f"Missing or nonsequential residues in /{chain_id}")


def confidence(entry, cfg):
    summary = read_json(entry["summary"])
    full = read_json(entry["full_data"])
    scalars = {f"af_{k}": float(v) for k, v in summary.items() if isinstance(v, (int, float, bool))}
    chains = np.asarray(full["token_chain_ids"])
    numbers = np.asarray(full["token_res_ids"], dtype=int)
    pae = np.asarray(full["pae"], dtype=float)
    cp = np.asarray(full["contact_probs"], dtype=float)
    if pae.shape != (len(chains), len(chains)) or cp.shape != pae.shape:
        raise ValueError("PAE/contact matrix shape mismatch")
    if not np.isfinite(pae).all() or not np.isfinite(cp).all():
        raise ValueError("Nonfinite confidence arrays")
    groups = {}
    site, length, end = entry["site"], cfg["insertion_length"], cfg["core_end"]
    for chain in "AB":
        selected = chains == chain
        core = np.isin(numbers, [construct_number(r, site, length) for r in range(1, end + 1)])
        insert = (numbers > site) & (numbers <= site + length) if site else np.zeros(len(chains), bool)
        groups[f"{chain}_core"] = selected & core
        groups[f"{chain}_tail"] = selected & ~core & ~insert
        if site:
            groups[f"{chain}_DrON"] = selected & insert
    for chain in "CD":
        groups[f"{chain}_DNA"] = chains == chain
    regions = []
    for g1, mask1 in groups.items():
        for g2, mask2 in groups.items():
            if not mask1.any() or not mask2.any():
                raise ValueError(f"Empty confidence region: {g1}/{g2}")
            values = pae[np.ix_(mask1, mask2)]
            regions.append(dict(model_id=entry["model_id"], condition=entry["condition"],
                                seed=entry["seed"], frame_region=g1, target_region=g2,
                                mean_pae=float(values.mean()), median_pae=float(np.median(values)),
                                mean_contact_probability=float(cp[np.ix_(mask1, mask2)].mean())))
    # AF Server can store token-level chain_ids in summary; preserve order, do not assume length=4.
    order = list(dict.fromkeys(full["token_chain_ids"]))
    if order != list("ABCD"):
        raise ValueError(f"Unexpected confidence chain order: {order}")
    for key in ("chain_pair_iptm", "chain_pair_pae_min"):
        matrix = np.asarray(summary[key], dtype=float)
        if matrix.shape != (4, 4):
            raise ValueError(f"Unexpected {key} shape")
        for i, a in enumerate(order):
            for j, b in enumerate(order):
                scalars[f"af_{key}_{a}_{b}"] = float(matrix[i, j])
    for key in ("chain_ptm", "chain_iptm"):
        for chain, value in zip(order, summary[key]):
            scalars[f"af_{key}_{chain}"] = float(value)
    return scalars, regions


def execute(session, manifest_path):
    from chimerax.core.commands import run, StringArg
    from chimerax.atomic import AtomicStructure
    from chimerax.core import __version__ as cx_version
    manifest = read_json(manifest_path)
    cfg, out = manifest["config"], Path(manifest["output"])
    if (out / "alignment_started.json").exists():
        raise FileExistsError("Alignment already started here; prepare a fresh output_directory")
    write_json(out / "alignment_started.json", {"chimerax_version": cx_version, "python": sys.version})
    (out / "aligned").mkdir()
    (out / "per_atom").mkdir()
    for entry in manifest["models"]:
        for name, digest in entry["hashes"].items():
            if sha256(entry[name]) != digest:
                raise ValueError(f"Input changed since preparation: {entry[name]}")

    def open_model(entry):
        opened = run(session, "open " + StringArg.unparse(entry["cif"]))
        structures = [m for m in opened if isinstance(m, AtomicStructure)]
        if len(structures) != 1:
            raise ValueError("Expected one atomic structure per model CIF")
        validate_model(structures[0], entry)
        return structures[0]

    ref_entry = next(m for m in manifest["models"] if m["model_id"] == manifest["reference_model_id"])
    reference = open_model(ref_entry)
    ref_atoms = index_atoms(reference)
    end, length = cfg["core_end"], cfg["insertion_length"]
    ca_keys = [(c, r, "CA") for c in "AB" for r in range(1, end + 1)]
    ref_ca = np.array([ref_atoms[k].scene_coord for k in ca_keys])
    rows, residue_rows, confidence_rows, transforms, commands = [], [], [], [], []
    # Manifest order puts all controls before insertion models.
    for entry in manifest["models"]:
        model = open_model(entry)  # independent copy, including reference self comparison
        atoms = index_atoms(model)
        site = entry["site"]
        mapped_keys = [(c, construct_number(r, site, length), n) for c, r, n in ca_keys]
        for key in mapped_keys:
            if key not in atoms:
                raise ValueError(f"Missing core CA: {key}")
        command = (f"mmaker #{model.id_string}/A,B:{core_spec(site, end, length)} "
                   f"to #{reference.id_string}/A,B:1-{end} pairing ss "
                   f"alg nw matrix BLOSUM-62 ssFraction 0 "
                   f"cutoffDistance {cfg['matchmaker_cutoff_angstrom']} showAlignment false")
        commands.append(f"# {entry['model_id']}\n{command}")
        result = run(session, command)
        if not isinstance(result, list) or len(result) != 1:
            raise RuntimeError("MatchMaker did not produce exactly one joint dimer fit")
        match = result[0]
        actual = {(atom_key(a), atom_key(b)) for a, b in zip(match["full match atoms"], match["full ref atoms"])}
        expected = set(zip(mapped_keys, ca_keys))
        if actual != expected:
            raise RuntimeError("MatchMaker correspondence differs from exact 358-site TFAP4 map; refusing metrics")
        used = {atom_key(a) for a in match["final match atoms"]}
        coords = np.array([atoms[k].scene_coord for k in mapped_keys])
        measures, deviations, local = summarize(coords, ref_ca)
        row = dict(model_id=entry["model_id"], condition=entry["condition"], site=site,
                   seed=entry["seed"], rank=entry["rank"],
                   is_reference=entry["model_id"] == ref_entry["model_id"],
                   mm_pruned_rmsd=float(match["final RMSD"]),
                   mm_full_rmsd=float(match["full RMSD"]), mm_retained_ca=len(used),
                   mm_retained_fraction=len(used) / len(ca_keys), **measures)
        if not np.isclose(row["mm_full_rmsd"], row["core_ca_rmsd"], atol=1e-4):
            raise RuntimeError("MatchMaker RMSD disagrees with scene-coordinate calculation")
        row["kabsch_identity_rmsd"] = best_fit_rmsd(coords, ref_ca)
        swapped = np.concatenate([ref_ca[end:], ref_ca[:end]])
        row["kabsch_swapped_AB_rmsd"] = best_fit_rmsd(coords, swapped)
        row["swap_improvement_angstrom"] = row["kabsch_identity_rmsd"] - row["kabsch_swapped_AB_rmsd"]
        for chain, sl in (("A", slice(0, end)), ("B", slice(end, 2 * end))):
            row[f"core_ca_rmsd_{chain}"] = rmsd(deviations[sl])
            row[f"core_ca_plddt_{chain}"] = float(np.mean([atoms[k].bfactor for k in mapped_keys[sl]]))
            row[f"mm_retained_ca_{chain}"] = sum(k[0] == chain for k in used)
        for k, mk, dev, lddt in zip(ca_keys, mapped_keys, deviations, local):
            residue_rows.append(dict(model_id=entry["model_id"], condition=entry["condition"], seed=entry["seed"],
                                     chain=k[0], tfap4_residue=k[1], construct_residue=mk[1],
                                     ca_deviation=float(dev), ca_lddt=float(lddt),
                                     ca_plddt=float(atoms[mk].bfactor), retained_in_fit=mk in used))
        atom_rows, heavy, backbone, dna = [], [], [], []
        for key, ra in ref_atoms.items():
            c, r, name = key
            if ra.element.number == 1 or (c in "AB" and r > end):
                continue
            if c not in "ABCD":
                continue
            mk = (c, construct_number(r, site, length) if c in "AB" else r, name)
            if mk not in atoms:
                raise ValueError(f"Missing homologous heavy atom {mk}")
            a = atoms[mk]
            if ra.residue.name != a.residue.name:
                raise ValueError(f"Residue identity mismatch: {key} -> {mk}")
            dev = float(np.linalg.norm(a.scene_coord - ra.scene_coord))
            region = "TFAP4_core" if c in "AB" else "DNA"
            atom_rows.append(dict(chain=c, reference_residue=r, construct_residue=mk[1],
                                  residue_name=a.residue.name, atom=name, region=region,
                                  deviation=dev, plddt=float(a.bfactor)))
            if c in "AB":
                heavy.append(dev)
                if name in ("N", "CA", "C", "O"):
                    backbone.append(dev)
            else:
                dna.append(dev)
        row.update(core_heavy_atom_rmsd=rmsd(heavy), core_backbone_rmsd=rmsd(backbone),
                   dna_heavy_atom_rmsd_after_protein_fit=rmsd(dna),
                   n_core_heavy_atoms=len(heavy), n_dna_heavy_atoms=len(dna))
        for chain in "AB":
            if site:
                cassette_atoms = [a for k, a in atoms.items() if k[0] == chain and site < k[1] <= site + length]
                row[f"dron_mean_atom_plddt_{chain}"] = float(np.mean([a.bfactor for a in cassette_atoms]))
        scalar_conf, region_conf = confidence(entry, cfg)
        row.update(scalar_conf)
        rows.append(row)
        confidence_rows.extend(region_conf)
        transforms.append(dict(model_id=entry["model_id"], matrix_3x4=model.scene_position.matrix.tolist()))
        write_csv(out / "per_atom" / (entry["model_id"] + ".csv"), atom_rows)
        run(session, f"save {StringArg.unparse(str(out / 'aligned' / (entry['model_id'] + '.cif')))} "
                     f"models #{model.id_string} relModel #{reference.id_string}")
        session.models.close([model])
        write_csv(out / "model_metrics.partial.csv", rows)
        print(f"Completed {len(rows)}/{len(manifest['models'])}: {entry['model_id']}", flush=True)
    write_csv(out / "model_metrics.csv", rows)
    write_csv(out / "per_residue.csv", residue_rows)
    write_csv(out / "confidence_regions.csv", confidence_rows)
    write_json(out / "transforms.json", transforms)
    (out / "matchmaker_commands.txt").write_text("\n".join(commands) + "\n", encoding="utf-8")
    run(session, "log save " + StringArg.unparse(str(out / "chimerax_log.html")))
    write_json(out / "alignment_complete.json", {"models": len(rows), "reference": ref_entry["model_id"]})


if "session" in globals():  # ChimeraX executes in a named sandbox, not __main__.
    try:
        execute(session, sys.argv[1])
    except Exception:
        traceback.print_exc()
        # Nonzero exit is necessary for Slurm; completion marker also checked by wrapper.
        raise SystemExit(1)
