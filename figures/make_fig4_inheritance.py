#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=os.path.dirname(os.path.abspath(__file__))

res=pd.read_csv(os.path.join(REB,"05_validation/inherited_vs_novel.csv"))
res["lab"]=(res.species+" "+res.tissue
            .replace({"stigmaE":"stigma (early)","stigmaL":"stigma (late)"}))
res=res.sort_values(["species","tissue"]).reset_index(drop=True)

W,H = 7.2, 5.25
fig=plt.figure(figsize=(W,H))
gs1=fig.add_gridspec(1,3, left=0.085, right=0.985, top=0.945, bottom=0.565, wspace=0.16)
gs2=fig.add_gridspec(1,2, left=0.085, right=0.985, top=0.405, bottom=0.095,
                     width_ratios=[1.0,1.42], wspace=0.62)

PANELS=[("napus","stigma","napus stigma"),("juncea","stigma","juncea stigma"),
        ("carinata","stigmaL","carinata stigma (late)")]
lim=6
for j,(sp,ti,ttl) in enumerate(PANELS):
    ax=fig.add_subplot(gs1[0,j])
    d=pd.read_csv(os.path.join(REB,"02_heb/theta1_a05/pairs_%s_%s.csv"%(sp,ti)))
    d=d[np.isfinite(d.allo_lfc)&np.isfinite(d.par_lfc)]
    ax.hexbin(d.par_lfc.clip(-lim,lim), d.allo_lfc.clip(-lim,lim), gridsize=32,
              cmap="Greys", bins="log", mincnt=1, linewidths=0)
    ax.plot([-lim,lim],[-lim,lim],color="#9A9A9A",lw=0.8,ls=(0,(4,3)),zorder=4)
    r=res[(res.species==sp)&(res.tissue==ti)].iloc[0]
    xs=np.array([-lim,lim]); ax.plot(xs,r.slope*xs,color=style.SPECIES[sp],lw=1.8,zorder=5)
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect("equal","box")
    ax.set_xticks([-5,0,5]); ax.set_yticks([-5,0,5])
    ax.set_title(style.ital(ttl),fontsize=7.6,pad=3)
    ax.text(0.95,0.05,"slope %.2f\n$R^2$ %.2f"%(r.slope,r.R2),transform=ax.transAxes,
            fontsize=6.8,ha="right",va="bottom",color=style.SPECIES[sp],fontweight="bold")
    if j>0: ax.set_yticklabels([])
    else:
        ax.set_ylabel("bias inside the\nallotetraploid (log$_2$)",fontsize=7.3)
        style.panel(ax,"a",dx=-0.32,dy=1.15)
    ax.set_xlabel("")
_axes_a=fig.axes[:3]
_pos=[a.get_position() for a in _axes_a]
_y0=min(q.y0 for q in _pos); _xmid=(min(q.x0 for q in _pos)+max(q.x1 for q in _pos))/2
fig.text(_xmid, _y0-0.062, "bias between the diploid progenitors (log$_2$)",
         ha="center", va="top", fontsize=7.3)

axb=fig.add_subplot(gs2[0,0])
for _,r in res.iterrows():
    m="o" if "pollen" in r.tissue else "s"
    axb.scatter(r.slope,r.R2,s=52,marker=m,color=style.SPECIES[r.species],
                edgecolor="white",linewidth=0.8,zorder=3)
axb.set_xlabel("slope   (1 = fully inherited)",fontsize=7.3)
axb.set_ylabel("$R^2$",fontsize=7.3)
axb.set_xlim(0.17,0.49); axb.set_ylim(0.055,0.235)
axb.set_xticks([0.2,0.3,0.4]); axb.set_yticks([0.10,0.15,0.20])
axb.legend(handles=[Line2D([],[],marker="o",color="none",markerfacecolor=style.GREY,
                           markeredgecolor="white",markersize=6,label="pollen"),
                    Line2D([],[],marker="s",color="none",markerfacecolor=style.GREY,
                           markeredgecolor="white",markersize=6,label="stigma")],
           fontsize=6.6,loc="upper left",bbox_to_anchor=(-0.02,1.02),
           handletextpad=0.3,borderpad=0.2)
style.panel(axb,"b",dx=-0.30,dy=1.13)

axc=fig.add_subplot(gs2[0,1])
cats=["pct_inherited","pct_novel","pct_lost","pct_reversed"]
names=["inherited","novel","lost","reversed"]
cols=["#3C3C3C","#0072B2","#D55E00","#B8B8B8"]
y=np.arange(len(res))[::-1]; left=np.zeros(len(res))
for c,nm,col in zip(cats,names,cols):
    axc.barh(y,res[c],left=left,color=col,height=0.66,label=nm,edgecolor="white",linewidth=0.5)
    left=left+res[c].values
axc.set_yticks(y); axc.set_yticklabels([style.ital(v) for v in res.lab],fontsize=6.8)
axc.set_xlabel("% of tested homoeolog pairs",fontsize=7.3)
axc.set_xlim(0,100); axc.set_xticks([0,25,50,75,100])
axc.legend(ncol=4,fontsize=6.6,loc="lower center",bbox_to_anchor=(0.5,1.01),
           columnspacing=1.5,handlelength=1.1,handletextpad=0.4,frameon=False)
style.panel(axc,"c",dx=-0.40,dy=1.13)

for f in ("Fig4_inheritance.png","Fig4_inheritance.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig4_inheritance  figsize %.2f x %.2f"%(W,H))
