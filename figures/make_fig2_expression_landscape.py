#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style; style.use()

REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild")
OUT=os.path.dirname(os.path.abspath(__file__))
NTOP=2000

SPECIES_OF={"Bra":"rapa","Bni":"nigra","Bol":"oleracea","Bna":"napus","Bju":"juncea","Bca":"carinata"}
def tissue_of(lib):
    if "PolA" in lib or "PolB" in lib or "PolC" in lib: return "pollen"
    if "_E_Stig" in lib: return "stigmaE"
    if "_L_Stig" in lib: return "stigmaL"
    return "stigma"

d=pd.read_csv(os.path.join(REB,"01_counts/gref_counts.tsv"), sep="\t")
meta=[c for c in d.columns if "|" in c]
libs=sorted({c.split("|")[1] for c in meta})
M=pd.DataFrame(index=d.index)
for L in libs:
    cols=[c for c in meta if c.split("|")[1]==L]
    M[L]=d[cols].sum(axis=1)
sp=[SPECIES_OF[L.split("_")[0]] for L in M.columns]
ti=[tissue_of(L) for L in M.columns]
print("libraries: %d  genes: %d" % (M.shape[1], M.shape[0]))
print("per species:", {s:sp.count(s) for s in dict.fromkeys(sp)})

cpm = M / M.sum(axis=0) * 1e6
lg  = np.log2(cpm + 1.0)
keep = lg.loc[lg.var(axis=1).sort_values(ascending=False).index[:NTOP]]
X = keep.T.values
X = X - X.mean(axis=0)
pca = PCA(n_components=5).fit(X); Y = pca.transform(X)
ev = pca.explained_variance_ratio_*100
print("variance explained: " + ", ".join("PC%d %.1f%%"%(i+1,v) for i,v in enumerate(ev)))

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.0], width_ratios=[1.3, 1.0],
                      hspace=0.46, wspace=0.40)
ax = fig.add_subplot(gs[0, :])
xmin,xmax = Y[:,0].min(), Y[:,0].max(); pad=(xmax-xmin)*0.10
ymin,ymax = Y[:,1].min(), Y[:,1].max(); ypad=(ymax-ymin)*0.30
ax.set_xlim(xmin-pad, xmax+pad*1.6); ax.set_ylim(ymin-ypad*0.45, ymax+ypad)
for lab, xs in (("stigma", Y[[i for i,t in enumerate(ti) if t!="pollen"],0]),
                ("pollen", Y[[i for i,t in enumerate(ti) if t=="pollen"],0])):
    ax.text(xs.mean(), ax.get_ylim()[1]*0.97, lab, ha="center", va="top",
            fontsize=13, color="#E3E3E3", fontweight="bold", zorder=0)
ax.axhline(0,color=style.LIGHT,lw=0.5,zorder=0); ax.axvline(0,color=style.LIGHT,lw=0.5,zorder=0)
for s_ in dict.fromkeys(sp):
    for t in dict.fromkeys(ti):
        idx=[i for i in range(len(sp)) if sp[i]==s_ and ti[i]==t]
        if not idx: continue
        ax.scatter(Y[idx,0], Y[idx,1], s=42, marker=style.TISSUE_MARKER[t],
                   facecolor=style.SPECIES[s_], edgecolor="white", linewidth=0.7,
                   alpha=0.95, zorder=3)
ax.set_xlabel("PC1 (%.1f%% of variance)"%ev[0]); ax.set_ylabel("PC2 (%.1f%%)"%ev[1])
order_sp=["rapa","nigra","oleracea","napus","juncea","carinata"]
x0 = 0.455
for j,s_ in enumerate(order_sp):
    ax.text(x0, 0.93-j*0.085, s_, transform=ax.transAxes, fontsize=8, fontstyle="italic",
            fontweight="bold", color=style.SPECIES[s_], ha="left", va="top")
ax.text(x0, 0.985, "species", transform=ax.transAxes, fontsize=7,
        color=style.ANNOT, ha="left", va="top")
h=[Line2D([],[],marker=style.TISSUE_MARKER[t],color="none",markerfacecolor=style.GREY,
          markeredgecolor="white",markersize=6,label=t) for t in ["pollen","stigma","stigmaE","stigmaL"]]
ax.legend(handles=h, loc="lower center", ncol=4, columnspacing=0.9, handletextpad=0.2,
          bbox_to_anchor=(0.5,-0.04), fontsize=6.8)
ax.set_title("PC1 separates tissue, not species "
             "(%d most variable Gref genes, log$_2$(CPM+1), %d libraries)"%(NTOP,M.shape[1]),
             fontsize=7.5, color=style.ANNOT, pad=4)
style.panel(ax,"a",dx=-0.072,dy=1.06)

ax2 = fig.add_subplot(gs[1,0])
ax2.bar(np.arange(1,6), ev, color=style.GREY, width=0.62)
for i,v in enumerate(ev): ax2.text(i+1, v+1.2, "%.1f"%v, ha="center", fontsize=6.5)
ax2.set_xticks(np.arange(1,6)); ax2.set_xlabel("principal component")
ax2.set_ylabel("variance explained (%)"); ax2.set_ylim(0, max(ev)*1.26)
style.panel(ax2,"b",dx=-0.16)

ax3 = fig.add_subplot(gs[1,1])
R = np.corrcoef(keep.T.rank(axis=1).values)
Z = linkage(1-R, method="average")
order = dendrogram(Z, no_plot=True)["leaves"]
im = ax3.imshow(R[np.ix_(order,order)], cmap="viridis",
                vmin=np.percentile(R,2), vmax=1, extent=[0,len(order),len(order),0])
ax3.set_xticks([]); ax3.set_yticks([])
ax3.set_xlabel("libraries, hierarchically clustered")
cb=fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.03)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0)
cb.set_label("Spearman $r$", fontsize=7)
sp_o=[sp[i] for i in order]
for i,s_ in enumerate(sp_o):
    ax3.add_patch(plt.Rectangle((i,-1.9),1,1.4,color=style.SPECIES[s_],lw=0,clip_on=False))
ax3.set_ylim(len(order),-2.1)
style.panel(ax3,"c",dx=-0.1,dy=1.10)

for f in ("Fig2_expression_landscape.png","Fig2_expression_landscape.pdf"):
    fig.savefig(os.path.join(OUT,f))
print("wrote Fig2_expression_landscape.png / .pdf")
pd.DataFrame({"library":M.columns,"species":sp,"tissue":ti,
              "PC1":Y[:,0],"PC2":Y[:,1],"PC3":Y[:,2]}).to_csv(
    os.path.join(OUT,"Fig2_pca_coordinates.csv"), index=False)
print("wrote Fig2_pca_coordinates.csv")
