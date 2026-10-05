#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=HERE
TH=1.0
SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
ORDER=[("napus","pollen"),("napus","stigma"),("carinata","pollen"),
       ("carinata","stigmaE"),("carinata","stigmaL"),("juncea","pollen"),("juncea","stigma")]
NICE={"stigmaE":"stigma\n(early)","stigmaL":"stigma\n(late)","stigma":"stigma","pollen":"pollen"}
CLS=["maintained","novel","switched"]
CLSLAB={"maintained":"parental bias maintained","novel":"truly novel (parents unbiased)",
        "switched":"switched (anti-parental by definition)"}

rows=[]
for sp,ti in ORDER:
    d=pd.read_csv(os.path.join(REB,"02_heb/theta1_a05/pairs_%s_%s.csv"%(sp,ti)))
    d=d[np.isfinite(d.allo_lfc)&np.isfinite(d.par_lfc)]
    b=d[d.allo.isin(["toward1","toward2"])].copy()
    pB=b.par_lfc.abs()>=TH; same=np.sign(b.allo_lfc)==np.sign(b.par_lfc)
    b["cls"]=np.where(pB&same,"maintained",np.where(pB&~same,"switched","novel"))
    g1,g2=SUBS[sp]
    for c in CLS:
        s=b[b.cls==c]
        k1=int((s.allo=="toward1").sum()); k2=int((s.allo=="toward2").sum())
        rows.append(dict(species=sp,tissue=ti,cls=c,t1=k1,t2=k2,n=k1+k2,
                         ratio=(k1/k2 if k2 else np.nan),g1=g1,g2=g2))
R=pd.DataFrame(rows); R.to_csv(os.path.join(REB,"05_validation/novel_switched_split.csv"),index=False)

fig=plt.figure(figsize=(7.2,5.3))
gs=fig.add_gridspec(2,1,height_ratios=[1.0,0.92],left=0.105,right=0.985,
                    top=0.885,bottom=0.145,hspace=0.80)
x=np.arange(len(ORDER))
HATCH={"maintained":None,"novel":"..","switched":"//"}
ALPHA={"maintained":1.0,"novel":0.65,"switched":0.30}

ax=fig.add_subplot(gs[0])
for i,(sp,ti) in enumerate(ORDER):
    base=0; tot=R[(R.species==sp)&(R.tissue==ti)].n.sum()
    for c in CLS:
        r=R[(R.species==sp)&(R.tissue==ti)&(R.cls==c)].iloc[0]
        ax.bar(x[i],r.n,bottom=base,width=0.6,color="#4D4D4D",alpha=ALPHA[c],
               hatch=HATCH[c],edgecolor="white",linewidth=0.6)
        base+=r.n
    mn=R[(R.species==sp)&(R.tissue==ti)&(R.cls=="maintained")].iloc[0].n
    ax.text(x[i],base+50,"%.0f%% maintained"%(100*mn/tot),ha="center",fontsize=6.4)
ax.set_xticks(x); ax.set_xticklabels([style.ital(sp)+"\n"+NICE[ti] for sp,ti in ORDER],fontsize=6.8)
ax.set_ylabel("biased homoeolog pairs",fontsize=7.4)
ax.set_title("every biased pair, split by what its diploid progenitors did",
             fontsize=7.4,color=style.ANNOT,pad=20)
ax.legend(handles=[mpatches.Patch(facecolor="#4D4D4D",alpha=ALPHA[c],hatch=HATCH[c],
                                  edgecolor="white",label=CLSLAB[c]) for c in CLS],
          fontsize=6.3,ncol=3,loc="lower center",bbox_to_anchor=(0.5,1.005),
          frameon=False,columnspacing=1.2,handlelength=1.5)
style.panel(ax,"a",dx=-0.082,dy=1.20)

axb=fig.add_subplot(gs[1])
MK={"maintained":"o","novel":"s","switched":"^"}
for i,(sp,ti) in enumerate(ORDER):
    g1,g2=SUBS[sp]
    for j,c in enumerate(CLS):
        r=R[(R.species==sp)&(R.tissue==ti)&(R.cls==c)].iloc[0]
        if not np.isfinite(r.ratio): continue
        xx=x[i]+(j-1)*0.21
        if c=="switched":
            axb.scatter(xx,r.ratio,s=34,marker=MK[c],facecolor="none",
                        edgecolor="#AAAAAA",linewidth=1.1,zorder=3)
        else:
            col=style.GENOME[g1] if r.ratio>1 else style.GENOME[g2]
            axb.scatter(xx,r.ratio,s=52,marker=MK[c],color=col,
                        edgecolor="white",linewidth=0.8,zorder=4)
axb.axhline(1.0,color="black",lw=0.9)
axb.axhspan(0.95,1.05,color="#F2F2F2",zorder=0)
axb.set_xticks(x)
axb.set_xticklabels([style.ital(sp)+"\n"+NICE[ti] for sp,ti in ORDER],fontsize=6.8)
axb.set_ylabel("ratio, first : second genome",fontsize=7.4)
axb.set_yscale("log"); axb.set_yticks([0.6,0.8,1.0,1.2,1.4])
axb.get_yaxis().set_major_formatter(plt.ScalarFormatter())
axb.set_ylim(0.54,1.55); axb.set_xlim(-0.6,len(ORDER)-0.4)
axb.set_title("the lean is carried by MAINTAINED bias; truly novel bias is near balance",
              fontsize=7.4,color=style.ANNOT,pad=20)
axb.legend(handles=[Line2D([],[],marker="o",color="none",markerfacecolor=style.GREY,
                           markeredgecolor="white",markersize=6,label="maintained"),
                    Line2D([],[],marker="s",color="none",markerfacecolor=style.GREY,
                           markeredgecolor="white",markersize=6,label="truly novel"),
                    Line2D([],[],marker="^",color="none",markerfacecolor="none",
                           markeredgecolor="#AAAAAA",markersize=6,label="switched (not interpretable)")],
            fontsize=6.3,ncol=3,loc="lower center",bbox_to_anchor=(0.5,1.005),
            frameon=False,columnspacing=1.2)
fig.text(0.985,0.018,"shaded band: within 5% of balance",ha="right",va="bottom",
         fontsize=6.3,color=style.ANNOT)
style.panel(axb,"b",dx=-0.082,dy=1.22)

for f in ("Fig6b_novel_bias.png","Fig6b_novel_bias.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig6b_novel_bias")
