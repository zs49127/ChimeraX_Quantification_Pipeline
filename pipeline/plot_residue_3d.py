"""Combine existing C-alpha deviation and pLDDT profiles without rerunning alignment.

python pipeline/plot_residue_3d.py results_corefix_v3
Optional: --html requires plotly and writes a standalone, offline HTML figure.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ORDER = ["TFAP4 Control"] + [f"Insertion at Residue {s}" for s in (76, 80, 136)]
COLORS = ["#444444", "#0072B2", "#D55E00", "#009E73"]


def summarize_residues(folder):
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    completed = json.loads((folder / "alignment_complete.json").read_text())
    expected = {m["model_id"] for m in manifest["models"]}
    if completed["models"] != len(expected):
        raise ValueError("Alignment is incomplete")
    data = pd.read_csv(folder / "per_residue.csv")
    required = {"model_id", "condition", "chain", "tfap4_residue", "ca_deviation", "ca_plddt"}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing columns: {required - set(data.columns)}")
    if set(data.model_id) != expected:
        raise ValueError("Residue table does not contain exactly the manifest models")
    if data.duplicated(["model_id", "chain", "tfap4_residue"]).any():
        raise ValueError("Duplicate model/chain/residue rows")
    end = manifest["config"]["core_end"]
    for (model, chain), sub in data.groupby(["model_id", "chain"]):
        if set(sub.tfap4_residue) != set(range(1, end + 1)):
            raise ValueError(f"Incomplete residue coverage: {model}/{chain}")
    if set(data.chain) != {"A", "B"} or not (data.groupby("model_id").size() == 2 * end).all():
        raise ValueError("Each model must have both complete protein chains")
    values = data[["ca_deviation", "ca_plddt"]].to_numpy()
    if not np.isfinite(values).all() or (data.ca_deviation < 0).any() or not data.ca_plddt.between(0, 100).all():
        raise ValueError("Invalid deviation or pLDDT values")
    ref = manifest["reference_model_id"]
    data = data[data.model_id != ref]
    if set(data.condition) != set(ORDER):
        raise ValueError("Expected the control and three insertion conditions")
    summary = data.groupby(["condition", "chain", "tfap4_residue"], as_index=False).agg(
        n_models=("model_id", "nunique"),
        deviation_mean=("ca_deviation", "mean"), deviation_sd=("ca_deviation", "std"),
        plddt_mean=("ca_plddt", "mean"), plddt_sd=("ca_plddt", "std"))
    return summary, end, ref


def draw_static(summary, end, destination):
    zmax = max(1.0, float(summary.deviation_mean.max()) * 1.1)
    fig = plt.figure(figsize=(15, 7.5), facecolor="white")
    for panel, chain in enumerate("AB", 1):
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        for condition, color in zip(ORDER, COLORS):
            sub = summary[(summary.chain == chain) & (summary.condition == condition)].sort_values("tfap4_residue")
            ax.plot(sub.tfap4_residue, sub.plddt_mean, sub.deviation_mean, color=color,
                    linewidth=1.5, alpha=0.85)
            # Every point represents one native TFAP4 residue, not a new replicate.
            ax.scatter(sub.tfap4_residue, sub.plddt_mean, sub.deviation_mean, color=color, s=7, alpha=0.75)
        ax.set(xlim=(1, end), ylim=(0, 100), zlim=(0, zmax), title=f"Chain {chain}")
        ax.set_xlabel("Native TFAP4 residue", labelpad=10)
        ax.set_ylabel("Mean Cα pLDDT", labelpad=10)
        ax.set_zlabel("Mean Cα deviation (Å)", labelpad=12)
        ax.view_init(elev=24, azim=-62)
        ax.set_box_aspect((1.5, 1, 1))
        ax.tick_params(labelsize=8)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
    fig.suptitle("TFAP4 residue confidence and structural deviation", fontsize=17, y=0.98)
    handles = [Line2D([0], [0], color=c, lw=2, label=l) for l, c in zip(ORDER, COLORS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.92), ncol=4, frameon=False)
    fig.text(0.5, 0.04, "Each point: a residue’s pose-mean pLDDT and deviation. Lines follow sequence order.\n"
             "Reference self-comparison excluded; variation is available in the CSV and interactive HTML.",
             ha="center", fontsize=10, color="#444444")
    fig.subplots_adjust(left=0.02, right=0.94, top=0.85, bottom=0.14, wspace=0.08)
    fig.savefig(destination / "residue_confidence_deviation_3d.png", dpi=220)
    fig.savefig(destination / "residue_confidence_deviation_3d.pdf")
    plt.close(fig)


def draw_html(summary, end, destination):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    zmax = max(1.0, float((summary.deviation_mean + summary.deviation_sd.fillna(0)).max()) * 1.1)
    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=["Chain A", "Chain B"])
    for col, chain in enumerate("AB", 1):
        for condition, color in zip(ORDER, COLORS):
            sub = summary[(summary.chain == chain) & (summary.condition == condition)].sort_values("tfap4_residue")
            ps, ds = sub.plddt_sd.fillna(0), sub.deviation_sd.fillna(0)
            fig.add_trace(go.Scatter3d(
                x=sub.tfap4_residue, y=sub.plddt_mean, z=sub.deviation_mean,
                mode="lines+markers", name=condition, legendgroup=condition, showlegend=col == 1,
                line=dict(color=color, width=4), marker=dict(color=color, size=3),
                customdata=np.column_stack([ps, ds, sub.n_models]),
                error_y=dict(type="data", array=np.minimum(ps, 100 - sub.plddt_mean),
                             arrayminus=np.minimum(ps, sub.plddt_mean), symmetric=False, visible=False),
                error_z=dict(type="data", array=ds, arrayminus=np.minimum(ds, sub.deviation_mean),
                             symmetric=False, visible=False),
                hovertemplate=(f"{condition} · Chain {chain}<br>TFAP4 residue %{{x}}"
                               "<br>Mean pLDDT: %{y:.2f} ± %{customdata[0]:.2f} SD"
                               "<br>Mean deviation: %{z:.2f} ± %{customdata[1]:.2f} Å SD"
                               "<br>Models: %{customdata[2]:.0f}<extra></extra>")), row=1, col=col)
    scene = dict(xaxis=dict(title="Native TFAP4 residue", range=[1, end]),
                 yaxis=dict(title="Mean Cα pLDDT", range=[0, 100]),
                 zaxis=dict(title="Mean Cα deviation (Å)", range=[0, zmax]),
                 aspectmode="manual", aspectratio=dict(x=1.5, y=1, z=1))
    fig.update_layout(template="plotly_white", height=780, scene=scene, scene2=scene,
                      title="TFAP4 residue confidence and structural deviation",
                      legend=dict(orientation="h", y=1.12, x=0.5, xanchor="center", groupclick="togglegroup"),
                      margin=dict(t=150, b=110, l=15, r=15),
                      updatemenus=[dict(type="buttons", direction="right", x=0, y=1.03, buttons=[
                          dict(label="Mean only", method="restyle", args=[{"error_y.visible": False, "error_z.visible": False}]),
                          dict(label="Mean ± SD", method="restyle", args=[{"error_y.visible": True, "error_z.visible": True}])])])
    fig.add_annotation(text="Drag to rotate · Hover for residue statistics · Click legend to hide/show a condition<br>"
                       "SD describes pose variation, not uncertainty of the mean; bars clipped to physical bounds.<br>"
                       "Mean trajectories do not show individual-model correlations. Reference self-comparison excluded.",
                       x=0.5, y=-0.12, xref="paper", yref="paper", showarrow=False)
    fig.write_html(destination / "residue_confidence_deviation_3d.html", include_plotlyjs=True, full_html=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_directory", type=Path)
    parser.add_argument("--output-dir", type=Path, help="Default: RESULTS/figures_3d; must be new")
    parser.add_argument("--html", action="store_true", help="Also generate offline rotatable HTML (requires plotly)")
    args = parser.parse_args()
    if args.html:
        try:
            import plotly
        except ImportError:
            parser.error("--html requires plotly; install it in your plotting environment or omit --html")
    summary, end, ref = summarize_residues(args.results_directory)
    destination = args.output_dir or args.results_directory / "figures_3d"
    destination.mkdir(parents=True, exist_ok=False)
    summary.to_csv(destination / "residue_confidence_deviation_3d_summary.csv", index=False)
    draw_static(summary, end, destination)
    if args.html:
        draw_html(summary, end, destination)
    (destination / "README.txt").write_text(
        f"Reference excluded: {ref}\nX: native TFAP4 residue; Y: mean CA pLDDT; Z: mean CA deviation (angstrom).\n"
        "Both means use the same model rows at each residue; chains remain separate.\n"
        "Deviation comes from the existing alignment; no refitting is performed.\n"
        "SD uses ddof=1 across poses, not independent biological replicates.\n"
        "A relationship between condition means does not establish a within-model correlation.\n", encoding="utf-8")
    print(f"3D plots and underlying summary saved to {destination}")


if __name__ == "__main__":
    main()
