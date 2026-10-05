#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
RB = os.path.join(HERE, "..", "..", "expression_rebuild")
ST = os.path.join(RB, "08_structure")

ORDER  = ["An", "Aj", "Bj", "Bc", "Cn", "Cc"]
GEN    = {"An": "A", "Aj": "A", "Bj": "B", "Bc": "B", "Cn": "C", "Cc": "C"}
SPEC   = {"An": "napus", "Aj": "juncea", "Bj": "juncea", "Bc": "carinata",
          "Cn": "napus", "Cc": "carinata"}
COL    = {k: style.GENOME[GEN[k]] for k in ORDER}

def letters(fig, axes, dx=0.052, dy=0.052):
    fig.canvas.draw()
    for ax, L in zip(axes, "abcdef"):
        p = ax.get_position()
        fig.text(p.x0 - dx, p.y1 + dy, L, fontsize=10, fontweight="bold",
                 va="top", ha="left")

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches=None)
    print("wrote", name)

D = pd.read_csv(os.path.join(ST, "divergence_pairs.csv.gz"))
S = pd.read_csv(os.path.join(ST, "divergence_summary.csv")).set_index("subg")
rec = D[D.recent]

fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), layout="constrained")
fig.get_layout_engine().set(w_pad=0.06, wspace=0.08)

ax = axes[0]
data = []
for s in ORDER:
    v = rec[(rec.subg == s)].ks.dropna()
    v = v[(v > 0) & (v <= 1.0)]
    data.append(np.log10(v.to_numpy()))
parts = ax.violinplot(data, positions=range(len(ORDER)), widths=0.8,
                      showextrema=False, showmedians=False)
for b, s in zip(parts["bodies"], ORDER):
    b.set_facecolor(COL[s]); b.set_alpha(0.55); b.set_edgecolor(COL[s]); b.set_linewidth(0.6)
for i, s in enumerate(ORDER):
    ax.plot([i], [np.log10(S.loc[s, "median_ks"])], "o", ms=4, mfc="white",
            mec="black", mew=0.9, zorder=5)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_ylabel(r"$K_s$ (log$_{10}$)")
ax.set_xlabel("Allo-subgenome")
ax.set_title("Divergence from the progenitor", fontsize=8)

ax = axes[1]
data = [rec[rec.subg == s].pid.dropna().to_numpy() for s in ORDER]
parts = ax.violinplot(data, positions=range(len(ORDER)), widths=0.8,
                      showextrema=False, showmedians=False)
for b, s in zip(parts["bodies"], ORDER):
    b.set_facecolor(COL[s]); b.set_alpha(0.55); b.set_edgecolor(COL[s]); b.set_linewidth(0.6)
for i, s in enumerate(ORDER):
    ax.plot([i], [S.loc[s, "median_pid"]], "o", ms=4, mfc="white", mec="black", mew=0.9, zorder=5)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_ylim(90, 100.5)
ax.set_ylabel("identity to progenitor (%)")
ax.set_xlabel("Allo-subgenome")
ax.set_title("Sequence identity", fontsize=8)

ax = axes[2]
FD = os.path.join(ST, "4dtv")
med4 = {}
for s_ in ORDER:
    v = np.loadtxt(os.path.join(FD, s_ + ".txt"))
    v = np.sort(v)
    med4[s_] = float(np.median(v))
    y = np.arange(1, len(v) + 1) / len(v)
    ax.step(v, y, where="post", color=COL[s_], lw=1.2,
            ls="-" if s_ in ("An", "Bj", "Cn") else (0, (3.5, 1.5)),
            label="%s  %.3f" % (s_, med4[s_]), zorder=3)
ax.axhline(0.5, ls=":", lw=0.8, color=style.GREY, zorder=1)
ax.set_xlim(0, 0.08); ax.set_ylim(0, 1.0)
ax.set_xlabel("4DTv per gene pair")
ax.set_ylabel("cumulative share of pairs")
ax.set_title("4DTv is zero-inflated", fontsize=8)
leg = ax.legend(loc="lower right", fontsize=6, frameon=False,
                handlelength=1.6, labelspacing=0.25, borderpad=0.2,
                title="median", title_fontsize=6)
leg._legend_box.align = "left"
ax.text(0.03, 0.955, "intercept = share with no transversion",
        transform=ax.transAxes, fontsize=6.2, color=style.ANNOT, va="top")

letters(fig, axes)
save(fig, "FigS11_divergence")

R = pd.read_csv(os.path.join(ST, "retention_by_layer.csv"))
T = pd.read_csv(os.path.join(ST, "retention_total.csv"))
SY = pd.read_csv(os.path.join(ST, "retention_syntenic.csv"))
TR = ["Bra_A", "Bol_C", "Bni_B", "Bna_A", "Bna_C", "Bju_A", "Bju_B", "Bca_B", "Bca_C"]
LAYER = ["LF", "MF1", "MF2"]

fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.0), layout="constrained",
                         width_ratios=[1.6, 1, 1])
fig.get_layout_engine().set(w_pad=0.07, wspace=0.09)

ax = axes[0]
w = 0.26
x = np.arange(len(TR))
HATCH = {"LF": "", "MF1": "//", "MF2": "xx"}
for j, lay in enumerate(LAYER):
    v = [R[(R.track == t) & (R.layer == lay)].pct.iloc[0] for t in TR]
    ax.bar(x + (j - 1) * w, v, width=w,
           color=[style.GENOME[style.SUB[t]] for t in TR],
           alpha=1.0 - 0.28 * j, edgecolor="white", linewidth=0.4, hatch=HATCH[lay])
ax.set_xticks(x)
ax.set_xticklabels([style.ital(R[R.track == t].label.iloc[0]) for t in TR],
                   rotation=45, ha="right", fontsize=6.3)
ax.set_ylabel("Arabidopsis anchors retained (%)")
ax.set_title("LF > MF1 > MF2 in all nine tracks", fontsize=8)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=style.GREY, alpha=1.0 - 0.28 * j,
                                 hatch=HATCH[l], ec="white") for j, l in enumerate(LAYER)],
          labels=LAYER, loc="upper right", ncol=3, fontsize=6.3, handlelength=1.3)
ax.set_ylim(0, 62)

ax = axes[1]
SY = SY.set_index("sub").loc[["An", "Cn", "Aj", "Bj", "Bc", "Cc"]].reset_index()
y = np.arange(len(SY))
ax.barh(y, SY.med_bin_any, height=0.6,
        xerr=[SY.med_bin_any - SY.q1, SY.q3 - SY.med_bin_any],
        error_kw=dict(lw=0.7, ecolor=style.GREY, capsize=1.8),
        color=[style.GENOME[g] for g in SY.genome])
ax.set_yticks(y)
ax.set_yticklabels([style.ital("%s %s" % (sp, g)) for sp, g in zip(SY.species, SY.genome)],
                   fontsize=6.3)
ax.invert_yaxis()
ax.set_xlim(0, 1.0)
ax.set_xlabel("progenitor genes retained")
ax.set_title("Progenitor gene set", fontsize=8)

ax = axes[2]
allo = T[T.progenitor.notna()].copy()
allo["subg"] = allo.track.map({"Bna_A": "An", "Bna_C": "Cn", "Bju_A": "Aj",
                               "Bju_B": "Bj", "Bca_B": "Bc", "Bca_C": "Cc"})
allo = allo.set_index("subg").loc[["An", "Cn", "Aj", "Bj", "Bc", "Cc"]].reset_index()
y = np.arange(len(allo))
ax.barh(y, allo.delta_pct, height=0.6,
        color=[style.GENOME[s[0]] for s in allo.subg])
for yi, (_, r) in zip(y, allo.iterrows()):
    ax.text(r.delta_pct - 0.1, yi, "%+.2f%s" % (r.delta_pct, "" if r.p < 0.05 else " n.s."),
            va="center", ha="right", fontsize=6.0, color=style.ANNOT)
ax.set_yticks(y); ax.set_yticklabels([])
ax.invert_yaxis()
ax.axvline(0, color="black", lw=0.7)
ax.set_xlim(-6.6, 0.9)
ax.set_xlabel("retention vs progenitor (pp)")
ax.set_title("Conserved anchors", fontsize=8)

letters(fig, axes, dx=0.052)
save(fig, "FigS12_retention")

C = pd.read_csv(os.path.join(ST, "collinearity_by_chr.csv"))
fig, ax = plt.subplots(figsize=(4.2, 2.7), layout="constrained")
for i, s in enumerate(ORDER):
    v = C[C.subg == s].top_share.to_numpy() * 100
    jit = (np.random.default_rng(0).random(len(v)) - 0.5) * 0.34
    ax.plot(np.full(len(v), i) + jit, v, "o", ms=3.2, color=COL[s], alpha=0.75,
            mec="none", zorder=3)
    ax.plot([i - 0.28, i + 0.28], [v.mean()] * 2, "-", lw=1.4, color="black", zorder=4)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_xlabel("Allo-subgenome")
ax.set_ylabel("pairs on the modal\nprogenitor chromosome (%)")
ax.set_ylim(94, 100.4)
ax.set_title("Collinearity with the progenitor is near complete", fontsize=8)
ax.legend(handles=[Line2D([], [], color="black", lw=1.4, label="subgenome mean"),
                   Line2D([], [], marker="o", ls="none", ms=3.2, color=style.GREY,
                          label="one chromosome")],
          loc="lower left", fontsize=6.5)
save(fig, "FigS13_collinearity")
