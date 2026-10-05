#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
REB = os.path.join(HERE, "..", "..", "expression_rebuild")
NTOP = 2000

D = pd.read_csv(os.path.join(REB, "01_counts", "gref_counts.tsv"), sep="\t")
cnt = D.drop(columns=["AT_gene", "br_sub"])
cnt.index = D.AT_gene
tracks = [c.split("|")[0] for c in cnt.columns]
libs   = [c.split("|")[1] for c in cnt.columns]

SPEC = {"Bra_A":"rapa","Bni_B":"nigra","Bol_C":"oleracea","Bna_A":"napus","Bna_C":"napus",
        "Bju_A":"juncea","Bju_B":"juncea","Bca_B":"carinata","Bca_C":"carinata"}
GEN  = {"Bra_A":"A","Bni_B":"B","Bol_C":"C","Bna_A":"A","Bna_C":"C",
        "Bju_A":"A","Bju_B":"B","Bca_B":"B","Bca_C":"C"}
LAB  = {"Bra_A":"rapa A","Bni_B":"nigra B","Bol_C":"oleracea C","Bna_A":"An","Bna_C":"Cn",
        "Bju_A":"Aj","Bju_B":"Bj","Bca_B":"Bc","Bca_C":"Cc"}
def tissue(l): return "pollen" if "Pol" in l else "stigma"

M = cnt.groupby(libs, axis=1).sum()
lib_spec = {}
for t, l in zip(tracks, libs): lib_spec[l] = SPEC[t]
cpm = M / M.sum(axis=0) * 1e6
lg  = np.log2(cpm + 1.0)
keep = lg.loc[lg.var(axis=1).sort_values(ascending=False).index[:NTOP]]
X = keep.T.values; X = X - X.mean(axis=0)
pca = PCA(n_components=5).fit(X); Y = pca.transform(X)
ev = pca.explained_variance_ratio_ * 100
print("libraries %d, genes %d, PC1 %.1f%%, PC2 %.1f%%" % (M.shape[1], M.shape[0], ev[0], ev[1]))

fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.5), layout="constrained")
fig.get_layout_engine().set(w_pad=0.08, wspace=0.10)

ax = axes[0]
MK = {"pollen": "o", "stigma": "^"}
seen = []
for i, l in enumerate(M.columns):
    sp = lib_spec[l]; ti = tissue(l)
    ax.scatter(Y[i, 0], Y[i, 1], s=24, marker=MK[ti],
               facecolor=style.SPECIES[sp], edgecolor="black", linewidth=0.4, zorder=3)
    if sp not in seen: seen.append(sp)
ax.set_xlabel("PC1  (%.1f%%)" % ev[0]); ax.set_ylabel("PC2  (%.1f%%)" % ev[1])
ax.set_title("All 39 libraries, Gref gene set", fontsize=8)
ht = [Line2D([], [], marker=MK[t], ls="", mfc="white", mec="black", ms=5, label=t)
      for t in ("pollen", "stigma")]
hs = [Line2D([], [], marker="o", ls="", mfc=style.SPECIES[sp], mec="black", mew=0.4, ms=5,
             label=style.ital(sp)) for sp in style.SPNAMES if sp in seen]
l1 = ax.legend(handles=hs, loc="center", fontsize=6.2, frameon=False,
               handletextpad=0.4, labelspacing=0.3)
ax.add_artist(l1)
ax.legend(handles=ht, loc="lower right", fontsize=6.5, frameon=False,
          handletextpad=0.4, labelspacing=0.3)

prof = {}
for t, l, c in zip(tracks, libs, cnt.columns):
    prof.setdefault((t, tissue(l)), []).append(c)
P = pd.DataFrame({ "%s %s" % (LAB[t], ti): np.log2(
        (cnt[cs].sum(axis=1) / cnt[cs].sum(axis=1).sum() * 1e6) + 1.0)
    for (t, ti), cs in prof.items() })
R = P.corr(method="spearman")
Z = linkage(squareform((1 - R).values, checks=False), method="average")
order = dendrogram(Z, no_plot=True)["leaves"]
R = R.iloc[order, order]
ax = axes[1]
im = ax.imshow(R.values, cmap="viridis", vmin=float(R.values.min()), vmax=1.0)
ax.set_xticks(range(len(R))); ax.set_xticklabels(R.columns, rotation=90, fontsize=5.2)
ax.set_yticks(range(len(R))); ax.set_yticklabels(R.index, fontsize=5.2)
ax.set_title("Subgenome profiles, Spearman", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=6)

fig.canvas.draw()
for ax_, L in zip(axes, "ab"):
    p = ax_.get_position()
    fig.text(p.x0 - 0.045, p.y1 + 0.045, L, fontsize=10, fontweight="bold", va="top", ha="left")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "FigS18_ordination.%s" % ext), bbox_inches="tight")
print("wrote FigS18_ordination")
