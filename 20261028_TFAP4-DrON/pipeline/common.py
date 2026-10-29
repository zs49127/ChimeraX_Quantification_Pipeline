"""Shared I/O and explicit TFAP4 numbering; only the standard library required."""
import csv
import hashlib
import json
from pathlib import Path


def read_json(path):
    with open(path, encoding="utf-8-sig") as stream:
        return json.load(stream)


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows for {path}")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def construct_number(number, site, length=257):
    return number + length if site and number > site else number


def core_spec(site, end=179, length=257):
    return f"1-{site},{site + length + 1}-{end + length}" if site else f"1-{end}"


def classify(sequence, control, sites, length):
    if sequence == control:
        return 0
    candidates = [s for s in sites if len(sequence) == len(control) + length
                  and sequence[:s] + sequence[s + length:] == control]
    if len(candidates) != 1:
        raise ValueError(f"Sequence is not a unique exact {length}-aa insertion at {sites}: {candidates}")
    return candidates[0]


def label(site):
    return f"Insertion at Residue {site}" if site else "TFAP4 Control"


def request_sequences(path):
    request = read_json(path)
    if isinstance(request, list):
        if len(request) != 1:
            raise ValueError(f"Expected one job: {path}")
        request = request[0]
    chains = []
    for entry in request["sequences"]:
        if len(entry) != 1:
            raise ValueError(f"Unexpected sequence entry: {path}")
        kind, info = next(iter(entry.items()))
        if kind not in ("proteinChain", "dnaSequence"):
            raise ValueError(f"Unsupported chain type {kind}: {path}")
        chains.extend([(kind, info["sequence"].upper())] * info.get("count", 1))
    if len(chains) != 4 or [c[0] for c in chains] != ["proteinChain"] * 2 + ["dnaSequence"] * 2:
        raise ValueError(f"Expected two proteins then two DNA chains: {path}")
    if chains[0][1] != chains[1][1]:
        raise ValueError(f"A/B proteins differ: {path}")
    if len(request["modelSeeds"]) != 1:
        raise ValueError(f"Expected a single seed: {path}")
    return request, [c[1] for c in chains]
