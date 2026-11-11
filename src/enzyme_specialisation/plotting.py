"""Publication panels from recorded states; variability is sample SD."""

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
CHIS = [0.0, 0.1, 0.325, 1.0]
R0 = [5.0, 7.5, 10.0, 12.5, 15.0]
RC = {
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.1,
    "svg.fonttype": "none",
    "svg.hashsalt": "enzyme-specialisation-v1",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
}


def save(fig, name, folder="figures/main"):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        metadata = (
            {"CreationDate": None, "ModDate": None}
            if ext == "pdf"
            else ({"Date": None} if ext == "svg" else {})
        )
        fig.savefig(folder / f"{name}.{ext}", dpi=600, bbox_inches="tight", metadata=metadata)
    preview = Path("figures/previews")
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=10)
    ax.text(-0.15, 1.06, letter, transform=ax.transAxes, fontweight="bold", fontsize=11)
    ax.tick_params(labelsize=7)


def response(ax, df, key):
    for chi, color in zip(CHIS, COLORS):
        for kernel, marker, style in [("loss_biased", "o", "-"), ("symmetric", "s", "--")]:
            sub = df[(df.chi == chi) & (df.mutation_kernel == kernel)]
            g = sub.groupby("redundancy_r0")[key].agg(["mean", "std"])
            if g.empty:
                continue
            ax.scatter(sub.redundancy_r0, sub[key], s=8, color=color, alpha=0.16, linewidth=0)
            ax.errorbar(
                g.index,
                g["mean"],
                yerr=g["std"],
                fmt=marker,
                ls=style,
                color=color,
                ms=3,
                capsize=2,
                elinewidth=0.6,
            )
    ax.set_xticks(R0)
    ax.set_xlabel(r"Initial reaction redundancy, $r_0$")


def legend(fig, ncol=6):
    handles = [Line2D([], [], color=c, label=rf"$\chi={x:g}$") for x, c in zip(CHIS, COLORS)]
    handles += [
        Line2D([], [], color="#333333", marker=m, ls=s, label=k)
        for k, m, s in [("Loss-biased", "o", "-"), ("Symmetric", "s", "--")]
    ]
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=7)


def grid(df, name, limits=None):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    fig.subplots_adjust(hspace=0.53, wspace=0.38, bottom=0.19)
    keys = [
        "delta_mean_breadth",
        "delta_strongest_share",
        "delta_weakest_P",
        "delta_mean_redundancy",
    ]
    titles = [
        "Functional breadth",
        "Catalytic concentration",
        "Whole-organism coverage",
        "Distributed redundancy",
    ]
    labels = [
        r"$\Delta B_{\mathrm{eff}}$",
        r"$\Delta q_{\max}$",
        r"$\Delta P_{\min}$",
        r"$\Delta r_{\mathrm{eff}}$",
    ]
    for i, (ax, key, title, label) in enumerate(zip(axes.flat, keys, titles, labels)):
        response(ax, df, key)
        ax.axhline(0, color=".65", lw=0.6, zorder=0)
        ax.set_ylabel(label)
        panel(ax, "ABCD"[i], title)
        if limits:
            ax.set_ylim(limits[key])
    legend(fig)
    save(fig, name)


def dynamics(df, h, name, exposure=False, end=250):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
    fig.subplots_adjust(wspace=0.48, bottom=0.29)
    response(axes[0], df, "final_mean_total_activity")
    axes[0].axhline(0.3, color=".6", lw=0.6)
    axes[0].set_ylabel("Mean total activity")
    axes[0].set_xlabel(r"Initial redundancy, $r_0$")
    panel(axes[0], "A", "Catalytic output")
    mask = np.asarray(h["fixed"]) <= end
    x = np.asarray(h["fixed"])[mask] / (100 if exposure else 1)
    for key, label, color in [
        ("mean_breadth", r"$B_{\mathrm{eff}}$", COLORS[0]),
        ("mean_strongest_share", r"$q_{\max}$", COLORS[1]),
        ("weakest_P", r"$P_{\min}$", COLORS[2]),
    ]:
        y = np.asarray(h[key])[mask]
        axes[1].step(x, y / y[0], where="post", color=color, label=label)
    axes[1].set_ylabel("Relative to ancestor")
    axes[1].legend(frameon=False, fontsize=7)
    panel(axes[1], "B", "Reference dynamics")
    axes[2].step(x, np.asarray(h["mean_breadth"])[mask], where="post", color=COLORS[0])
    axes[2].set_ylabel(r"$B_{\mathrm{eff}}$")
    panel(axes[2], "C", "Absolute breadth")
    for ax in axes[1:]:
        ax.set_xticks(np.linspace(0, end / 100 if exposure else end, 6))
        ax.set_xlabel("Substitutions per protein" if exposure else "Fixed substitutions")
    legend(fig)
    save(fig, name)


def exposure_curves(output, delta=False):
    groups = {}
    gridx = np.linspace(0, 5, 201)
    for p in (Path(output) / "trajectories").glob("*.npz"):
        with np.load(p, allow_pickle=False) as h:
            c = json.loads(str(h["metadata"]))["config"]
            n = c["n_proteins"]
            if h["fixed"][-1] < 5 * n:
                continue
            x = h["fixed"] / n
            idx = np.searchsorted(x, gridx, side="right") - 1
            y = h["mean_breadth"][idx]
            if delta:
                y = y - y[0]
            groups.setdefault((c["mutation_kernel"], c["chi"], n / 8), []).append(y)
    if not groups:
        return
    fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.1), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.15, hspace=0.42, bottom=0.19)
    palette = ["#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7"]
    for row, kernel in enumerate(["loss_biased", "symmetric"]):
        for col, chi in enumerate(CHIS):
            ax = axes[row, col]
            panel(ax, chr(65 + row * 4 + col), rf"$\chi={chi:g}$")
            for r, color in zip(R0, palette):
                ys = groups.get((kernel, chi, r))
                if ys is None:
                    continue
                a = np.array(ys)
                mean = a.mean(0)
                sd = a.std(0, ddof=1) if len(a) > 1 else np.zeros_like(mean)
                ax.step(gridx, mean, where="post", color=color, label=rf"$r_0={r:g}$")
                ax.fill_between(
                    gridx, mean - sd, mean + sd, step="post", color=color, alpha=0.10, lw=0
                )
            if col == 0:
                ax.set_ylabel(
                    ("Loss-biased" if row == 0 else "Symmetric")
                    + "\n"
                    + (r"$\Delta B_{\mathrm{eff}}$" if delta else r"$B_{\mathrm{eff}}$")
                )
            if row == 1:
                ax.set_xlabel(r"Exposure, $x=n_{\mathrm{fix}}/N$")
            ax.set_xlim(0, 5)
            ax.set_xticks(np.arange(6))
    handles = [Line2D([], [], color=c, label=rf"$r_0={r:g}$") for c, r in zip(palette, R0)]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False)
    save(fig, "Figure_5_breadth_vs_exposure" + ("_change" if delta else ""))


def alternatives(a, b=None):
    frames = [("250", a)] + ([("5N", b)] if b is not None else [])
    fig, axes = plt.subplots(len(frames), 2, figsize=(7, 2.7 * len(frames)), squeeze=False)
    fig.subplots_adjust(hspace=0.65, wspace=0.45, right=0.84, bottom=0.13, top=0.92)
    lo = min(df.delta_mean_breadth.min() for _, df in frames)
    for row, (label, df) in enumerate(frames):
        for col, kernel in enumerate(["loss_biased", "symmetric"]):
            matrix = (
                df[df.mutation_kernel == kernel]
                .pivot_table(index="chi", columns="redundancy_r0", values="delta_mean_breadth")
                .reindex(index=CHIS, columns=R0)
            )
            ax = axes[row, col]
            im = ax.pcolormesh(
                np.arange(6) - 0.5,
                np.arange(5) - 0.5,
                matrix.to_numpy(),
                cmap="cividis",
                vmin=lo,
                vmax=0,
                rasterized=False,
            )
            ax.invert_yaxis()
            ax.set_xticks(range(5), R0)
            ax.set_yticks(range(4), CHIS)
            ax.set_xlabel(r"$r_0$")
            ax.set_ylabel(r"$\chi$")
            panel(ax, chr(65 + row * 2 + col), f"{kernel.replace('_', ' ')} · {label}")
    cax = fig.add_axes([0.89, 0.23, 0.025, 0.54])
    colorbar = fig.colorbar(im, cax=cax, label=r"$\Delta B_{\mathrm{eff}}$")
    colorbar.solids.set_rasterized(False)
    save(fig, "Alternative_A_response_heatmaps", "figures/alternative")
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    fig.subplots_adjust(bottom=0.30)
    for chi, color in zip(CHIS, COLORS):
        for kernel, marker in [("loss_biased", "o"), ("symmetric", "s")]:
            sub = (
                a[(a.chi == chi) & (a.mutation_kernel == kernel)]
                .groupby("redundancy_r0")[["delta_mean_breadth", "delta_weakest_P"]]
                .mean()
            )
            ax.scatter(-sub.delta_mean_breadth, sub.delta_weakest_P, c=color, marker=marker, s=20)
    ax.set_xlabel(r"Specialisation, $-\Delta B_{\mathrm{eff}}$")
    ax.set_ylabel(r"$\Delta P_{\min}$")
    legend(fig, 3)
    save(fig, "Alternative_B_specialisation_coverage", "figures/alternative")
    if b is not None:
        keys = [
            "n_proteins",
            "mutation_kernel",
            "chi",
            "replicate",
            "ancestry_seed",
            "evolution_seed",
        ]
        paired = a.merge(b, on=keys, suffixes=("_250", "_5N"), validate="one_to_one")
        fig, ax = plt.subplots(figsize=(4.4, 3.8))
        fig.subplots_adjust(bottom=0.30)
        for chi, color in zip(CHIS, COLORS):
            for kernel, marker in [("loss_biased", "o"), ("symmetric", "s")]:
                sub = paired[(paired.chi == chi) & (paired.mutation_kernel == kernel)]
                ax.scatter(
                    sub.delta_mean_breadth_250,
                    sub.delta_mean_breadth_5N,
                    s=10,
                    color=color,
                    marker=marker,
                    alpha=0.5,
                )
        bounds = [min(ax.get_xlim()[0], ax.get_ylim()[0]), 0]
        ax.plot(bounds, bounds, color=".6", lw=0.6)
        ax.set_xlabel(r"$\Delta B_{\mathrm{eff}}^{(250)}$")
        ax.set_ylabel(r"$\Delta B_{\mathrm{eff}}^{(5N)}$")
        ax.set_xlim(bounds)
        ax.set_ylim(bounds)
        legend(fig, 3)
        save(fig, "Alternative_C_paired_endpoints", "figures/alternative")


def plot_all(output, experiment="all"):
    if experiment not in {"all", "equal_total_250", "equal_exposure_c5"}:
        raise ValueError(f"Unknown experiment: {experiment}")
    with plt.rc_context(RC):
        a = pd.read_csv("data/reference/replicate_results.csv")
        regenerated = Path(output) / "equal_total_250/replicate_results.csv"
        if regenerated.exists():
            a = pd.read_csv(regenerated)
        h = pd.read_csv("data/reference/reference_history.csv")
        path = Path(output) / "equal_exposure_c5/replicate_results.csv"
        b = pd.read_csv(path) if path.exists() else None
        limits = None
        if b is not None:
            counts = b.groupby(["mutation_kernel", "chi", "redundancy_r0"]).size()
            if len(counts) != 40 or counts.min() < 5:
                raise ValueError("New figure export requires the complete five-seed grid")
            limits = {}
            for key in [
                "delta_mean_breadth",
                "delta_strongest_share",
                "delta_weakest_P",
                "delta_mean_redundancy",
            ]:
                values = pd.concat([a[key], b[key]])
                lo = min(0.0, values.min())
                hi = max(0.0, values.max())
                pad = 0.07 * max(hi - lo, 0.01)
                limits[key] = (lo - pad, hi + pad)
        caption = Path("figures/captions.md")
        caption.parent.mkdir(exist_ok=True)
        sample_a = sorted(
            a.groupby(["mutation_kernel", "chi", "redundancy_r0"]).size().unique().tolist()
        )
        sample_b = (
            sorted(b.groupby(["mutation_kernel", "chi", "redundancy_r0"]).size().unique().tolist())
            if b is not None
            else []
        )
        if b is not None and sample_a != sample_b:
            raise ValueError(
                "Stopping regimes have unequal replicate counts; complete paired aggregation before plotting"
            )
        caption.write_text(
            f"# Export sample counts\n\n250-substitution analysis: {sample_a} replicates per cell.\n\n5N analysis and exposure bands: {sample_b} replicates per cell.\n\nError bars and bands are sample SD (ddof=1). Full panel descriptions: [captions](../docs/figure_captions.md).\n"
        )
        if experiment in {"all", "equal_total_250"}:
            grid(a, "Figure_1_specialisation_grid", limits)
            ref = Path(output) / "trajectories/N100_loss_biased_chi0.325_rep000.npz"
            if ref.exists():
                with np.load(ref) as z:
                    dynamics(a, z, "Figure_2_output_and_dynamics")
            else:
                dynamics(a, h, "Figure_2_output_and_dynamics")
        if experiment in {"all", "equal_exposure_c5"}:
            if b is None:
                raise ValueError(
                    "Run and aggregate the complete paper grid before exposure figures"
                )
            grid(b, "Figure_3_equal_exposure_grid", limits)
            ref = Path(output) / "trajectories/N100_loss_biased_chi0.325_rep000.npz"
            with np.load(ref) as z:
                dynamics(b, z, "Figure_4_equal_exposure_dynamics", True, 500)
            exposure_curves(output)
            exposure_curves(output, True)
        if experiment == "all":
            alternatives(a, b)
