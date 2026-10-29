"""Descriptive plots: model-level SD and separate seed-level summaries; no p-values."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_json, write_json

ORDER = ["TFAP4 Control"] + [f"Insertion at Residue {s}" for s in (76, 80, 136)]
COLORS = ["#444444", "#0072B2", "#D55E00", "#009E73"]


def save(fig, folder, name):
    fig.tight_layout()
    fig.savefig(folder / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot(out):
    out = Path(out)
    manifest = read_json(out / "manifest.json")
    complete = read_json(out / "alignment_complete.json")
    df = pd.read_csv(out / "model_metrics.csv", dtype={"seed": str})
    if complete["models"] != len(manifest["models"]) or len(df) != complete["models"]:
        raise ValueError("Incomplete alignment results")
    if set(df.model_id) != {m["model_id"] for m in manifest["models"]} or not df.model_id.is_unique:
        raise ValueError("Model identity mismatch")
    # Keep reference self-row in raw metrics but exclude it from descriptive comparisons.
    data = df.loc[~df.is_reference].copy()
    figs = out / "figures"
    figs.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    metrics = [c for c in data.select_dtypes(include="number").columns if c not in ("site", "rank")]
    data.groupby("condition")[metrics].agg(["count", "mean", "std", "median"]).to_csv(out / "condition_summary.csv")
    seeds = data.groupby(["condition", "seed"])[metrics].mean()
    seeds.to_csv(out / "seed_means.csv")
    seeds.groupby("condition").agg(["count", "mean", "std"]).to_csv(out / "seed_summary.csv")
    panels = [("core_ca_rmsd", "All-core Cα RMSD (Å)"),
              ("mm_pruned_rmsd", "Pruned MatchMaker RMSD (Å)"),
              ("mm_retained_fraction", "Fraction of Cα retained in fit"),
              ("core_backbone_rmsd", "Core backbone RMSD (Å)"),
              ("core_heavy_atom_rmsd", "Core heavy-atom RMSD (Å)"),
              ("core_ca_lddt", "Cα distance preservation (0–1)"),
              ("dna_heavy_atom_rmsd_after_protein_fit", "DNA RMSD after protein fit (Å)"),
              ("af_iptm", "AlphaFold ipTM"), ("af_ptm", "AlphaFold pTM")]
    fig, axes = plt.subplots(3, 3, figsize=(17, 13))
    for ax, (metric, title) in zip(axes.flat, panels):
        for i, condition in enumerate(ORDER):
            sub = data[data.condition == condition]
            for j, (_, group) in enumerate(sub.groupby("seed", sort=True)):
                offset = (j - 2) * 0.11
                ax.scatter(i + offset + np.linspace(-0.025, 0.025, len(group)), group[metric],
                           s=15, alpha=0.5, color=COLORS[i])
                ax.scatter(i + offset, group[metric].mean(), s=40, marker="D", color=COLORS[i],
                           edgecolor="white", linewidth=0.5)
            ax.errorbar(i + 0.37, sub[metric].mean(), yerr=sub[metric].std(ddof=1),
                        fmt="o", color="black", capsize=3, markersize=4)
        ax.set_title(title)
        ax.set_xticks(range(4), ["TFAP4\nControl", "Insertion at\nResidue 76", "Insertion at\nResidue 80", "Insertion at\nResidue 136"], rotation=20)
        ax.grid(axis="y", alpha=0.15)
    fig.suptitle("TFAP4 core comparison | dots: poses; diamonds: seed means; black: pose mean ± SD\n"
                 "Control n=24 (reference excluded); each insertion n=25; five seeds per condition", y=1.01)
    save(fig, figs, "model_metric_superplots")

    residues = pd.read_csv(out / "per_residue.csv", dtype={"seed": str})
    residues = residues[residues.model_id != manifest["reference_model_id"]]
    aggregates = residues.groupby(["condition", "chain", "tfap4_residue"])[["ca_deviation", "ca_plddt", "ca_lddt"]].agg(["mean", "std", "count"])
    aggregates.to_csv(out / "residue_summary.csv")
    for metric, ylabel in [("ca_deviation", "Cα deviation (Å)"), ("ca_plddt", "Cα pLDDT"),
                           ("ca_lddt", "Cα distance preservation")]:
        fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
        for ax, chain in zip(axes, "AB"):
            for condition, color in zip(ORDER, COLORS):
                a = aggregates.loc[(condition, chain)][metric]
                x, mu, sd = a.index.to_numpy(), a["mean"].to_numpy(), a["std"].to_numpy()
                ax.plot(x, mu, label=condition, color=color)
                ax.fill_between(x, np.maximum(0, mu - sd), mu + sd, alpha=0.12, color=color)
            for site in (76, 80, 136):
                ax.axvline(site + 0.5, ls=":", color="gray", lw=0.8)
            ax.set_ylabel(ylabel)
            ax.set_title(f"Chain {chain}: pose mean ± SD; original TFAP4 numbering")
        axes[0].legend(ncol=2, fontsize=8)
        axes[-1].set_xlabel("TFAP4 residue (cassette excluded)")
        save(fig, figs, f"residue_{metric}")

    vmax = residues.ca_deviation.max()
    for chain in "AB":
        fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, layout="constrained")
        for ax, condition in zip(axes, ORDER):
            sub = residues[(residues.chain == chain) & (residues.condition == condition)]
            matrix = sub.pivot(index="model_id", columns="tfap4_residue", values="ca_deviation")
            im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=vmax, cmap="magma",
                           extent=[0.5, manifest['config']['core_end'] + 0.5, len(matrix) - 0.5, -0.5])
            ax.set_title(condition)
            ax.set_ylabel("Model index")
        axes[-1].set_xlabel("TFAP4 residue")
        fig.colorbar(im, ax=axes, label="Cα deviation (Å)", shrink=0.8)
        fig.savefig(figs / f"deviation_heatmap_{chain}.png", dpi=200)
        fig.savefig(figs / f"deviation_heatmap_{chain}.pdf")
        plt.close(fig)

    # Atom distributions are averaged within model first: no pooling atoms as replicates.
    thresholds = np.linspace(0, max(10, float(vmax)), 160)
    ecdfs = []
    for m in data.itertuples():
        atoms = pd.read_csv(out / "per_atom" / f"{m.model_id}.csv")
        dev = np.sort(atoms.loc[atoms.region == "TFAP4_core", "deviation"].to_numpy())
        values = np.searchsorted(dev, thresholds, side="right") / len(dev)
        ecdfs.extend(dict(model_id=m.model_id, condition=m.condition, threshold=t, fraction=f)
                     for t, f in zip(thresholds, values))
    ecdf = pd.DataFrame(ecdfs)
    ecdf.to_csv(out / "heavy_atom_ecdf.csv", index=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    for condition, color in zip(ORDER, COLORS):
        a = ecdf[ecdf.condition == condition].groupby("threshold").fraction.agg(["mean", "std"])
        ax.plot(a.index, a["mean"], label=condition, color=color)
        ax.fill_between(a.index, np.clip(a["mean"] - a["std"], 0, 1), np.clip(a["mean"] + a["std"], 0, 1), color=color, alpha=0.15)
    ax.set(xlabel="Core heavy-atom deviation threshold (Å)", ylabel="Fraction of atoms within threshold",
           title="Within-model ECDFs: mean ± SD across models")
    ax.legend()
    save(fig, figs, "heavy_atom_ecdf")

    regions = pd.read_csv(out / "confidence_regions.csv")
    regions = regions[regions.model_id != manifest["reference_model_id"]]
    regions.groupby(["condition", "frame_region", "target_region"])[["mean_pae", "mean_contact_probability"]].agg(["mean", "std"]).to_csv(out / "confidence_region_summary.csv")
    fig, axes = plt.subplots(1, 4, figsize=(19, 5), layout="constrained")
    limit = regions.mean_pae.max()
    for ax, condition in zip(axes, ORDER):
        sub = regions[regions.condition == condition]
        matrix = sub.pivot_table(index="frame_region", columns="target_region", values="mean_pae", aggfunc="mean")
        im = ax.imshow(matrix, vmin=0, vmax=limit, cmap="viridis")
        ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=90)
        ax.set_yticks(range(len(matrix)), matrix.index)
        ax.set_title(condition)
        ax.set_xlabel("Target region")
        ax.set_ylabel("Alignment-frame region")
    fig.colorbar(im, ax=axes, label="Mean predicted aligned error (Å)", shrink=0.8)
    fig.savefig(figs / "region_pae.png", dpi=200)
    fig.savefig(figs / "region_pae.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for condition, color in zip(ORDER, COLORS):
        sub = data[data.condition == condition]
        axes[0].scatter(sub.af_iptm, sub.core_ca_rmsd, c=color, label=condition, alpha=0.7)
        axes[1].scatter(sub.mm_retained_fraction, sub.core_ca_rmsd, c=color, alpha=0.7)
    axes[0].set(xlabel="AlphaFold ipTM", ylabel="All-core Cα RMSD (Å)")
    axes[1].set(xlabel="Fraction retained by MatchMaker", ylabel="All-core Cα RMSD (Å)")
    axes[0].legend(fontsize=7)
    save(fig, figs, "confidence_and_fit_coverage")
    write_json(out / "plotting_complete.json", {"models_plotted": len(data), "figure_files": len(list(figs.iterdir()))})
    print(f"Plots and summaries saved to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory")
    plot(parser.parse_args().output_directory)
