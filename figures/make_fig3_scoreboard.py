#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0,os.path.join(HERE,"..","tables"))
import style; style.use()
from canonical_values import H, SENS, ORDER, SUBS, heb

LABEL={"juncea_pollen":"juncea pollen","juncea_stigma":"juncea stigma",
       "napus_pollen":"napus pollen","napus_stigma":"napus stigma",
       "carinata_pollen":"carinata pollen","carinata_stigmaE":"carinata stigma (early)",
       "carinata_stigmaL":"carinata stigma (late)"}
H=H.set_index("sample").loc[ORDER].reset_index()

fig=plt.figure(figsize=(7.2,4.5))
gs=fig.add_gridspec(2,3, height_ratios=[1.0,0.80], width_ratios=[2.25,0.95,1.25],
                    left=0.155, right=0.985, top=0.90, bottom=0.10,
                    hspace=0.62, wspace=0.10)

ax=fig.add_subplot(gs[0,0:2])
y=np.arange(len(H))[::-1]
for i,r in H.iterrows():
    col=style.GENOME[r.sub1] if r.ratio>1 else style.GENOME[r.sub2]
    a=1.0 if r.resolved else 0.40
    ax.plot([r.ci_lo,r.ci_hi],[y[i],y[i]],color=col,lw=2.4,alpha=a,solid_capstyle="round",zorder=2)
    ax.plot([r.ratio],[y[i]],"o",ms=6.5,color=col,mec="white",mew=0.9,alpha=a,zorder=3)
ax.axvline(1.0,color="black",lw=0.9,zorder=1)
ax.set_yticks(y); ax.set_yticklabels([style.ital(LABEL[s]) for s in H["sample"]],fontsize=7.4)
ax.set_xlim(0.74,1.26); ax.set_xticks([0.8,0.9,1.0,1.1,1.2])
ax.set_xlabel("ratio of biased pairs, subgenome 1 : subgenome 2   (95% CI)",fontsize=7.4)
ax.set_ylim(-0.7,len(H)-0.3)
style.panel(ax,"a",dx=-0.235,dy=1.10)

axt=fig.add_subplot(gs[0,2]); axt.axis("off")
axt.set_ylim(ax.get_ylim()); axt.set_xlim(0,1)
for i,r in H.iterrows():
    if r.resolved:
        axt.text(0.0,y[i],"%s favoured"%r.favoured,fontsize=7.0,va="center",
                 color=style.GENOME[r.favoured])
    else:
        axt.text(0.0,y[i],"not resolved",fontsize=7.0,va="center",color=style.ANNOT,style="italic")
    p=r.p
    axt.text(0.62,y[i],("p = %.2g"%p if p>=1e-4 else "p = %.0e"%p),fontsize=6.6,
             va="center",color=style.ANNOT)
fig.text(0.5,0.955,"two-fold change, padj $\\leq$ 0.05, CDS counting",
         ha="center",fontsize=7.2,color=style.ANNOT)

axb=fig.add_subplot(gs[1,0])
mk=["o","s","^","D","v"]
grey=["#1A1A1A","#4D4D4D","#7A7A7A","#A5A5A5","#C8C8C8"]
for k,(tag,lab) in enumerate(SENS):
    try: d=heb(tag)
    except Exception: continue
    xs=[];ys=[]
    for i,s in enumerate(H["sample"]):
        if s in d.index:
            xs.append(d.loc[s].toward1/d.loc[s].toward2); ys.append(y[i])
    axb.scatter(xs,ys,s=23,marker=mk[k],label=lab,color=grey[k],alpha=0.95,
            edgecolor="white",linewidth=0.35)
axb.axvline(1.0,color="black",lw=0.9)
axb.set_yticks(y); axb.set_yticklabels([style.ital(LABEL[s]) for s in H["sample"]],fontsize=6.3)
axb.set_xlabel("ratio across thresholds",fontsize=7.2)
axb.set_xlim(0.74,1.26); axb.set_xticks([0.8,1.0,1.2])
axb.set_ylim(-0.7,len(H)-0.3)
axb.legend(fontsize=5.9,ncol=1,loc="center left",bbox_to_anchor=(1.02,0.5),
           handletextpad=0.25,borderpad=0.2,labelspacing=0.35)
style.panel(axb,"b",dx=-0.335,dy=1.14)

axc=fig.add_subplot(gs[1,2]); axc.axis("off")
def chain(yy,items,rel):
    x=0.02
    for j,it in enumerate(items):
        axc.text(x,yy,it,fontsize=13,fontweight="bold",color=style.GENOME[it],
                 transform=axc.transAxes,va="center",ha="left")
        if j<len(items)-1:
            axc.text(x+0.105,yy,rel[j],fontsize=10,color=style.ANNOT,
                     transform=axc.transAxes,va="center",ha="left")
            x+=0.215
axc.text(0.02,0.97,"resulting order",fontsize=7.2,color=style.ANNOT,transform=axc.transAxes,va="top")
axc.text(0.02,0.78,"pollen",fontsize=7.4,transform=axc.transAxes,va="center")
chain(0.60,["B","A","C"],[">","="])
axc.text(0.02,0.38,"stigma",fontsize=7.4,transform=axc.transAxes,va="center")
chain(0.20,["B","A","C"],[">",">"])
style.panel(axc,"c",dx=-0.06,dy=1.14)

for f in ("Fig3_scoreboard.png","Fig3_scoreboard.pdf"): fig.savefig(os.path.join(HERE,f))
print("wrote Fig3_scoreboard (canonical primary values)")
print(H[["sample","ratio","ci_lo","ci_hi","resolved","favoured"]].to_string(index=False))
