"""Exercise every plot with explicitly synthetic data, never scientific results.

Usage: python tests/smoke_plot.py .local/plot_smoke
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
from common import write_json
from plot_results import ORDER, plot


def smoke(folder):
    out = Path(folder)
    out.mkdir(parents=True, exist_ok=False)
    (out / "per_atom").mkdir()
    rng = np.random.default_rng(91)
    models, residues, regions = [], [], []
    for site, condition in zip((0, 76, 80, 136), ORDER):
        for seed in range(5):
            for rank in range(5):
                model_id = f"site{site}_seed{seed}_rank{rank}"
                value = float(rng.uniform(0.5, 4))
                models.append(dict(model_id=model_id, condition=condition, site=site, seed=str(seed), rank=rank,
                                   is_reference=site == seed == rank == 0,
                                   core_ca_rmsd=value, mm_pruned_rmsd=value / 3, mm_retained_fraction=0.6,
                                   core_backbone_rmsd=value, core_heavy_atom_rmsd=value * 1.3,
                                   core_ca_lddt=0.8, dna_heavy_atom_rmsd_after_protein_fit=value * 2,
                                   af_iptm=float(rng.uniform(0.3, 0.8)), af_ptm=0.6))
                for chain in "AB":
                    for res in range(1, 180):
                        residues.append(dict(model_id=model_id, condition=condition, seed=str(seed), chain=chain,
                                             tfap4_residue=res, ca_deviation=float(rng.uniform(0, value)),
                                             ca_plddt=80.0, ca_lddt=0.8))
                for a in ("A_core", "B_core", "C_DNA", "D_DNA"):
                    for b in ("A_core", "B_core", "C_DNA", "D_DNA"):
                        regions.append(dict(model_id=model_id, condition=condition, frame_region=a,
                                            target_region=b, mean_pae=value, mean_contact_probability=0.1))
                pd.DataFrame(dict(region=["TFAP4_core"] * 100, deviation=rng.uniform(0, 10, 100))).to_csv(out / "per_atom" / f"{model_id}.csv", index=False)
    write_json(out / "manifest.json", dict(models=models, reference_model_id=models[0]["model_id"], config={"core_end": 179}))
    write_json(out / "alignment_complete.json", {"models": 100})
    pd.DataFrame(models).to_csv(out / "model_metrics.csv", index=False)
    pd.DataFrame(residues).to_csv(out / "per_residue.csv", index=False)
    pd.DataFrame(regions).to_csv(out / "confidence_regions.csv", index=False)
    plot(out)
    summary = pd.read_csv(out / "condition_summary.csv", header=[0, 1], index_col=0)
    assert summary.loc["TFAP4 Control", ("core_ca_rmsd", "count")] == 24
    assert summary.loc["Insertion at Residue 76", ("core_ca_rmsd", "count")] == 25
    assert len(list((out / "figures").glob("*.png"))) == 9
    print("Synthetic plotting smoke test passed; these are NOT biological results.")


if __name__ == "__main__":
    smoke(sys.argv[1])
