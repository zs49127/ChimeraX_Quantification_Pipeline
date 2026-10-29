"""Validate all extracted AlphaFold Server jobs and freeze a reproducible manifest."""
import argparse
import fnmatch
from collections import Counter, defaultdict
from pathlib import Path
from common import read_json, write_json, write_csv, sha256, classify, label, request_sequences


def prepare(config_path):
    cfg = read_json(config_path)
    root = Path(cfg["project_root"]).resolve()
    inputs = root / cfg["input_directory"]
    all_requests = sorted(inputs.rglob("*_job_request.json"))
    patterns = cfg.get("include_folder_patterns", ["*"])
    requests = [p for p in all_requests if any(fnmatch.fnmatchcase(p.relative_to(inputs).parts[0], pattern)
                                             for pattern in patterns)]
    if not requests:
        raise ValueError(f"No extracted *_job_request.json files under {inputs}")
    jobs = [(p, *request_sequences(p)) for p in requests]
    if len({tuple(seqs[2:]) for _, _, seqs in jobs}) != 1:
        raise ValueError("DNA sequences/order differ across jobs; a fixed C/D correspondence is required")
    sequences = sorted({s[0] for _, _, s in jobs}, key=len)
    control = sequences[0]
    if len(sequences) != len(cfg["insertion_sites"]) + 1:
        raise ValueError("Expected exactly four unique protein sequences including a control")
    models, seen, cassettes = [], set(), set()
    counts = defaultdict(list)
    for request_path, request, seqs in jobs:
        site = classify(seqs[0], control, cfg["insertion_sites"], cfg["insertion_length"])
        if site:
            cassettes.add(seqs[0][site:site + cfg["insertion_length"]])
        seed = str(request["modelSeeds"][0])
        if (site, seed) in seen:
            raise ValueError(f"Repeated actual seed {seed} within condition {site}")
        seen.add((site, seed))
        counts[site].append(seed)
        stem = request_path.name.removesuffix("_job_request.json")
        expected = cfg["expected_models_per_seed"]
        files = sorted(request_path.parent.glob(stem + "_model_*.cif"))
        if len(files) != expected:
            raise ValueError(f"Expected {expected} model CIFs for {request_path}, got {len(files)}")
        for rank in range(expected):
            paths = {"cif": request_path.parent / f"{stem}_model_{rank}.cif",
                     "summary": request_path.parent / f"{stem}_summary_confidences_{rank}.json",
                     "full_data": request_path.parent / f"{stem}_full_data_{rank}.json",
                     "request": request_path}
            for path in paths.values():
                if not path.is_file():
                    raise ValueError(f"Missing input: {path}")
            models.append(dict(model_id=f"site{site}_seed{seed}_rank{rank}", condition=label(site),
                               site=site, seed=seed, rank=rank, job_name=request["name"], sequences=seqs,
                               **{k: str(v.resolve()) for k, v in paths.items()},
                               hashes={k: sha256(v) for k, v in paths.items()}))
    if len(cassettes) != 1:
        raise ValueError("Inserted cassette sequence differs across conditions")
    for site in [0] + cfg["insertion_sites"]:
        if len(counts[site]) != cfg["expected_seeds_per_condition"]:
            raise ValueError(f"Condition {site}: found {len(counts[site])} seeds, expected 5")
    models.sort(key=lambda r: (r["site"], int(r["seed"]), r["rank"]))
    reference_id = cfg.get("reference_model_id") or models[0]["model_id"]
    refs = [m for m in models if m["model_id"] == reference_id and m["site"] == 0]
    if len(refs) != 1:
        raise ValueError("reference_model_id must name one control model")
    out = root / cfg["output_directory"]
    if out.exists():
        raise FileExistsError(f"Refusing to mix results; choose a NEW output_directory: {out}")
    out.mkdir(parents=True)
    manifest = dict(config=cfg, reference_model_id=reference_id, models=models,
                    excluded_requests=[str(p) for p in all_requests if p not in requests],
                    control_sequence=control, output=str(out),
                    counts=dict(Counter(m["condition"] for m in models)))
    write_json(out / "manifest.json", manifest)
    write_csv(out / "inventory.csv", [{k: v for k, v in m.items() if k not in ("hashes", "sequences")} for m in models])
    print(f"Validated {len(models)} models. Reference: {reference_id}\n{out / 'manifest.json'}")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.json")
    prepare(parser.parse_args().config)
