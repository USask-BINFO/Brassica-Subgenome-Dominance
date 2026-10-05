#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binomtest
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=HERE
LAYER={1:"LF",2:"MF1",3:"MF2"}
SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
SAMPLES=[("juncea","pollen"),("juncea","stigma"),("napus","pollen"),("napus","stigma"),
         ("carinata","pollen"),("carinata","stigmaE"),("carinata","stigmaL")]
LAB={"stigmaE":"stigma (early)","stigmaL":"stigma (late)"}

rows=[]
for sp,ti in SAMPLES:
    f=os.path.join(REB,"02_heb/theta1_a05/pairs_%s_%s.csv"%(sp,ti))
    if not os.path.exists(f): continue
    d=pd.read_csv(f)
    d["layer"]=d["key"].str.split("|").str[1].astype(int)
    for L in (1,2,3):
        s=d[(d.layer==L)&(d.allo.isin(["toward1","toward2"]))]
        k1=int((s.allo=="toward1").sum()); k2=int((s.allo=="toward2").sum())
        if k1+k2<30: continue
        p=binomtest(k1,k1+k2,0.5).pvalue
        rows.append(dict(species=sp,tissue=ti,layer=LAYER[L],n=k1+k2,
                         toward1=k1,toward2=k2,ratio=k1/k2 if k2 else np.nan,p=p,
                         sub1=SUBS[sp][0],sub2=SUBS[sp][1]))
R=pd.DataFrame(rows)
R.to_csv(os.path.join(REB,"05_validation/heb_by_br_subgenome.csv"),index=False)
print(R.to_string(index=False))

fig,ax=plt.subplots(figsize=(7.2,3.5))
order=[(sp,ti) for sp,ti in SAMPLES if ((R.species==sp)&(R.tissue==ti)).any()]
ygap=0.28
yticks=[];ylabs=[]
for i,(sp,ti) in enumerate(order[::-1]):
    base=i*1.0
    if i%2==0: ax.axhspan(base-0.46,base+0.46,color='#F4F4F4',zorder=0)
    for k,L in enumerate(["LF","MF1","MF2"]):
        r=R[(R.species==sp)&(R.tissue==ti)&(R.layer==L)]
        if r.empty: continue
        r=r.iloc[0]
        col=style.GENOME[r.sub1] if r.ratio>1 else style.GENOME[r.sub2]
        yy=base+(1-k)*ygap
        sig = r.p<0.05
        ax.scatter(r.ratio,yy,s=52 if sig else 30,marker=["o","s","^"][k],color=col,
                   edgecolor="black" if sig else "white",linewidth=0.7 if sig else 0.5,
                   alpha=1.0 if sig else 0.55,zorder=3)
    yticks.append(base); ylabs.append("%s %s"%(sp,LAB.get(ti,ti)))
ax.axvline(1.0,color="black",lw=0.9,zorder=1)
ax.set_yticks(yticks); ax.set_yticklabels([style.ital(v) for v in ylabs],fontsize=7.4)
ax.set_xlabel("homoeolog-bias ratio within each Br-subgenome layer, subgenome 1 : subgenome 2",fontsize=7.4)
ax.set_xlim(0.55,1.85)
from matplotlib.lines import Line2D
h=[Line2D([],[],marker=m,color="none",markerfacecolor=style.GREY,markeredgecolor="white",
          markersize=6,label=l) for m,l in zip(["o","s","^"],["LF","MF1","MF2"])]
h.append(Line2D([],[],marker="o",color="none",markerfacecolor="none",markeredgecolor="black",
                markersize=6,label="p < 0.05"))
ax.legend(handles=h,fontsize=6.8,ncol=4,loc="lower center",bbox_to_anchor=(0.5,1.01),
          frameon=False,columnspacing=1.4,handletextpad=0.35)
fig.text(0.985,0.018,"a black outline marks a layer whose lean differs from 1:1",
         ha="right",va="bottom",fontsize=6.4,color=style.ANNOT)
fig.subplots_adjust(left=0.20,right=0.985,top=0.86,bottom=0.235)
for f in ("Fig7_br_subgenome.png","Fig7_br_subgenome.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig7_br_subgenome")
