#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=os.path.dirname(os.path.abspath(__file__))

d=pd.read_csv(os.path.join(REB,"05_validation/cytonuclear.csv"))
d=d[d.species.isin(["juncea","carinata"])]
SETS=["plastid","chloroplast","thylakoid","photosynthesis","mitochondrion"]
d=d[d.geneset.isin(SETS)]
MAT={"juncea":"A (maternal)","carinata":"B (maternal)"}

fig,axes=plt.subplots(1,2,figsize=(7.2,2.9),sharey=True,gridspec_kw={"wspace":0.10})
for ax,sp in zip(axes,["juncea","carinata"]):
    sub=d[d.species==sp]
    tiss=sorted(sub.tissue.unique())
    for i,gs_ in enumerate(SETS):
        for j,t in enumerate(tiss):
            r=sub[(sub.geneset==gs_)&(sub.tissue==t)]
            if r.empty: continue
            r=r.iloc[0]
            x=r.delta; y=len(SETS)-1-i+(j-(len(tiss)-1)/2)*0.22
            sig = r.p<0.05 and abs(r.delta)>1e-6
            col = style.GENOME["A"] if sp=="juncea" else style.GENOME["B"]
            ax.scatter(x,y,s=40 if sig else 22,
                       color=col if x>0 else style.LIGHT,
                       edgecolor="black" if sig else "white",linewidth=0.7 if sig else 0.4,zorder=3)
            if sig: ax.text(x+0.006,y,"p=%.3g"%r.p,fontsize=5.6,va="center")
    ax.axvline(0,color="black",lw=0.8)
    ax.set_yticks(range(len(SETS))); ax.set_yticklabels(SETS[::-1],fontsize=7)
    ax.set_xlim(-0.16,0.16)
    ax.set_xlabel("median shift toward the maternal subgenome\n(log$_2$)")
    ax.set_title(style.ital(sp)+", maternal = %s"%MAT[sp],fontsize=8,pad=4)
axes[0].text(0.5,-0.42,"a positive shift in BOTH would support maternal inheritance; it appears only in juncea",
             transform=axes[0].transAxes,ha="left",fontsize=6.6,color=style.ANNOT)
style.panel(axes[0],"a",dx=-0.22,dy=1.14); style.panel(axes[1],"b",dx=-0.06,dy=1.14)
for f in ("Fig7_cytonuclear.png","Fig7_cytonuclear.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig7_cytonuclear")
