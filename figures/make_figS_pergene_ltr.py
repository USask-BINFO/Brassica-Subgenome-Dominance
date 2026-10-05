#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
ST = os.path.join(HERE, "..", "..", "expression_rebuild", "08_structure")

SAMPLES = ["napus_pollen","napus_stigma","carinata_pollen","carinata_stigmaE",
           "carinata_stigmaL","juncea_pollen","juncea_stigma"]
PRETTY = {"napus_pollen":"napus pollen","napus_stigma":"napus stigma",
          "carinata_pollen":"carinata pollen","carinata_stigmaE":"carinata stigma (early)",
          "carinata_stigmaL":"carinata stigma (late)","juncea_pollen":"juncea pollen",
          "juncea_stigma":"juncea stigma"}

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches=None)
    print("wrote", name)

M = pd.read_csv(os.path.join(ST, "pergene_bias_matrix.csv"), index_col=0).loc[SAMPLES, SAMPLES]
C = pd.read_csv(os.path.join(ST, "pergene_bias_correlation.csv"))

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), layout="constrained",
                         width_ratios=[1.25, 1])
fig.get_layout_engine().set(w_pad=0.08, wspace=0.10)

ax = axes[0]
im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-0.75, vmax=0.75)
ax.set_xticks(range(len(SAMPLES)))
ax.set_xticklabels([style.ital(PRETTY[s]) for s in SAMPLES], rotation=45, ha="right",
                   fontsize=6.2)
ax.set_yticks(range(len(SAMPLES)))
ax.set_yticklabels([style.ital(PRETTY[s]) for s in SAMPLES], fontsize=6.2)
for i in range(len(SAMPLES)):
    for j in range(len(SAMPLES)):
        v = M.iat[i, j]
        ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.4,
                color="white" if abs(v) > 0.45 else "black")
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("Per-gene bias groups by species, not tissue", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("Spearman $\\rho$", fontsize=7); cb.ax.tick_params(labelsize=6)

ax = axes[1]
C["kind"] = np.where(C.same_species, "same species", "different species")
for k, col, off in [("same species", style.A, -0.16), ("different species", style.GREY, 0.16)]:
    v = C[C.kind == k].rho.to_numpy()
    jit = (np.random.default_rng(1).random(len(v)) - 0.5) * 0.18
    ax.plot(np.full(len(v), 0 if k == "same species" else 1) + off + jit, v, "o",
            ms=4.5, color=col, alpha=0.85, mec="none")
    ax.plot([(0 if k == "same species" else 1) + off - 0.14,
             (0 if k == "same species" else 1) + off + 0.14], [v.mean()] * 2,
            "-", lw=1.5, color="black", zorder=4)
ax.axhline(0, color="black", lw=0.7, ls="--")
ax.set_xticks([0 - 0.16, 1 + 0.16])
ax.set_xticklabels(["same\nspecies", "different\nspecies"])
ax.set_xlim(-0.6, 1.6)
ax.set_ylabel("Spearman $\\rho$ between samples")
ax.set_title("Shared direction, different genes", fontsize=8)
ax.legend(handles=[Line2D([], [], color="black", lw=1.5, label="mean")],
          loc="upper right", fontsize=6.5)
fig.canvas.draw()
for ax_, L in zip(axes, "ab"):
    p = ax_.get_position()
    fig.text(p.x0 - 0.055, p.y1 + 0.05, L, fontsize=10, fontweight="bold", va="top", ha="left")
save(fig, "FigS14_pergene_bias")

L = pd.read_csv(os.path.join(ST, "ltr_by_layer.csv"))
LAY = ["LF", "MF1", "MF2"]
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), layout="constrained", sharey=True)
fig.get_layout_engine().set(w_pad=0.06, wspace=0.07)
for ax, sp in zip(axes, ["napus", "juncea", "carinata"]):
    g = L[L.species == sp].set_index("layer").loc[LAY]
    x = np.arange(len(LAY)); w = 0.36
    ax.bar(x - w/2, g.mean_ltr1 * 1000, width=w, color=style.GENOME[g.sub1.iloc[0]],
           edgecolor="white", linewidth=0.4, label=g.sub1.iloc[0])
    ax.bar(x + w/2, g.mean_ltr2 * 1000, width=w, color=style.GENOME[g.sub2.iloc[0]],
           edgecolor="white", linewidth=0.4, label=g.sub2.iloc[0])
    for xi, p in zip(x, g.p):
        if p < 0.05:
            top = max(g.mean_ltr1.iloc[xi], g.mean_ltr2.iloc[xi]) * 1000
            ax.text(xi, top + 0.22, "*", ha="center", va="bottom", fontsize=9,
                    color=style.ANNOT)
    ax.set_xticks(x); ax.set_xticklabels(LAY)
    ax.set_xlabel("Br-subgenome")
    ax.set_title(style.ital(sp), fontsize=8)
    ax.legend(fontsize=6.5, loc="upper right", handlelength=1.1, ncol=2, columnspacing=0.8)
axes[0].set_ylabel("flanking LTR density\n(per kb, $\\times10^{-3}$)")
axes[0].set_ylim(0, 10.6)
fig.canvas.draw()
for ax_, Lt in zip(axes, "abc"):
    p = ax_.get_position()
    fig.text(p.x0 - 0.045, p.y1 + 0.055, Lt, fontsize=10, fontweight="bold",
             va="top", ha="left")
save(fig, "FigS15_ltr_layers")
