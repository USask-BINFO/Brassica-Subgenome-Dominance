#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
ST = os.path.join(HERE, "..", "..", "expression_rebuild", "08_structure")

D = pd.read_csv(os.path.join(ST, "he_dosage.csv"))
LAB = {"Bna_A":"An","Bna_C":"Cn","Bju_A":"Aj","Bju_B":"Bj","Bca_B":"Bc","Bca_C":"Cc"}
GEN = {"Bna_A":"A","Bna_C":"C","Bju_A":"A","Bju_B":"B","Bca_B":"B","Bca_C":"C"}
ORDER = ["Bna_A","Bna_C","Bju_A","Bju_B","Bca_B","Bca_C"]
D = D[D.retained.isin(ORDER)]

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")
fig.get_layout_engine().set(w_pad=0.07, wspace=0.09)

ax = axes[0]
for i, t in enumerate(ORDER):
    v = np.log2(D[D.retained == t].dose.clip(lower=0.05).values)
    if not len(v): continue
    p = ax.violinplot([v], positions=[i], widths=0.8, showextrema=False, showmedians=False)
    b = p["bodies"][0]
    b.set_facecolor(style.GENOME[GEN[t]]); b.set_alpha(0.55)
    b.set_edgecolor(style.GENOME[GEN[t]]); b.set_linewidth(0.6)
    ax.plot([i], [np.median(v)], "o", ms=4, mfc="white", mec="black", mew=0.9, zorder=5)
ax.axhline(0, ls="-",  lw=0.8, color=style.GREY, zorder=1)
ax.axhline(1, ls="--", lw=0.9, color="black", zorder=2)
ax.text(5.6, 1.0, " 2x, expected\n if the segment\n were replaced", fontsize=6,
        va="center", ha="left", color=style.ANNOT)
ax.text(5.6, 0.0, " 1x, expected\n if it were\n simply lost", fontsize=6,
        va="center", ha="left", color=style.ANNOT)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([LAB[t] for t in ORDER])
ax.set_ylabel(r"dose of the retained copy (log$_2$)")
ax.set_xlabel("retained Allo-subgenome")
ax.set_ylim(-3.2, 3.2); ax.set_xlim(-0.7, 7.6)
ax.set_title("No second copy where the partner is missing", fontsize=8)

ax = axes[1]
for t in ORDER:
    d = D[D.retained == t]
    ax.scatter(d.anchors, d.dose.clip(upper=8), s=9, alpha=0.5,
               color=style.GENOME[GEN[t]], edgecolor="none", zorder=3)
ax.axhline(1, lw=0.8, color=style.GREY, zorder=1)
ax.axhline(2, ls="--", lw=0.9, color="black", zorder=2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("run length (anchors missing on the partner)")
ax.set_ylabel("dose of the retained copy")
ax.set_title("Longer runs are no closer to a double dose", fontsize=8)
h = [Line2D([], [], marker="o", ls="", color=style.GENOME[g], ms=5, label=g + " genome")
     for g in ("A", "B", "C")]
ax.legend(handles=h, loc="upper right", fontsize=6.3, frameon=False)

fig.canvas.draw()
for ax_, L in zip(axes, "ab"):
    p = ax_.get_position()
    fig.text(p.x0 - 0.048, p.y1 + 0.05, L, fontsize=10, fontweight="bold", va="top", ha="left")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "FigS20_he_dosage.%s" % ext), bbox_inches="tight")
print("wrote FigS20_he_dosage")
