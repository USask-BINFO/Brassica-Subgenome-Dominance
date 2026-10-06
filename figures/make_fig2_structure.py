#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); import style; style.use()
OUT=HERE

LTR={"A":[("napus",0.005541),("juncea",0.005320)],
     "B":[("juncea",0.004500),("carinata",0.004164)],
     "C":[("napus",0.008167),("carinata",0.008923)]}
DIP={"rapa":0.007733,"oleracea":0.011238}
FEAT=pd.DataFrame({
 "sub":["An","Aj","Bj","Bc","Cn","Cc"], "genome":list("AABBCC"),
 "fractionation":[0.382,0.328,0.375,0.355,0.434,0.448],
 "DCJ":[36,49,43,26,42,38],
 "Ka/Ks":[0.1825,0.1843,0.1967,0.1921,0.2783,0.2893],
 "LTR":[0.005541,0.005320,0.004500,0.004164,0.008167,0.008923]})

fig=plt.figure(figsize=(7.2,4.9))
gs=fig.add_gridspec(2,2, height_ratios=[1.0,1.0], width_ratios=[1.0,1.0],
                    left=0.100, right=0.985, top=0.905, bottom=0.095,
                    hspace=0.46, wspace=0.38)

ax=fig.add_subplot(gs[0,0])
HOSTM={0:"o",1:"^"}
for i,(g,hosts) in enumerate(LTR.items()):
    vals=[v for _,v in hosts]
    for j,(h,v) in enumerate(hosts):
        ax.scatter(i+(j-0.5)*0.30, v*100, s=58, color=style.GENOME[g], marker=HOSTM[j],
                   edgecolor="white", linewidth=0.8, zorder=3)
    ax.plot([i-0.30,i+0.30],[np.mean(vals)*100]*2,color=style.GENOME[g],lw=1.0,alpha=0.35,zorder=1)
    ax.annotate("%.1f%% apart"%(100*abs(vals[0]-vals[1])/np.mean(vals)),
                (i, min(vals)*100), xytext=(0,-13), textcoords="offset points",
                ha="center", fontsize=6.3, color=style.ANNOT)
ax.set_xticks(range(3)); ax.set_xticklabels(["A","B","C"],fontweight="bold",fontsize=9)
ax.set_ylabel("gene-proximal LTR density\n(% of the 5 kb flank)",fontsize=7.3)
ax.set_xlim(-0.55,2.55); ax.set_ylim(0.28,1.03)
ax.set_title("each genome keeps its own density in both hybrids it enters",
             fontsize=7.3,color=style.ANNOT,pad=5)
ax.legend(handles=[Line2D([],[],marker="o",color="none",markerfacecolor=style.GREY,
                          markeredgecolor="white",markersize=6,label="first host"),
                   Line2D([],[],marker="^",color="none",markerfacecolor=style.GREY,
                          markeredgecolor="white",markersize=6,label="second host")],
          fontsize=6.3,loc="upper left",bbox_to_anchor=(-0.02,1.0),
          handletextpad=0.3,borderpad=0.25,labelspacing=0.3)
style.panel(ax,"a",dx=-0.26,dy=1.18)

axb=fig.add_subplot(gs[0,1])
NIGRA=0.006876
bars=[(style.ital("rapa")+"\n(A progenitor)",   DIP["rapa"]*100,     style.GENOME["A"], False),
      (style.ital("nigra")+"\n(B progenitor)",  NIGRA*100,           style.GENOME["B"], True),
      (style.ital("oleracea")+"\n(C progenitor)",DIP["oleracea"]*100, style.GENOME["C"], False)]
for i,(nm,v,c,sep) in enumerate(bars):
    axb.bar(i,v,color=("white" if sep else c),width=0.52,
            edgecolor=c,linewidth=(1.2 if sep else 0),hatch=("///" if sep else None))
    axb.text(i,v+0.03,"%.2f%%"%v,ha="center",fontsize=7.2)
axb.set_xticks(range(3)); axb.set_xticklabels([b[0] for b in bars],fontsize=6.4)
axb.set_ylabel("gene-proximal LTR density\n(% of the 5 kb flank)",fontsize=7.3)
axb.set_ylim(0,1.95); axb.set_xlim(-0.6,2.6)
r_dip=DIP["oleracea"]/DIP["rapa"]
axb.annotate("",xy=(2,1.62),xytext=(0,1.62),arrowprops=dict(arrowstyle="<->",lw=0.9,color="black"))
axb.text(1.0,1.67,style.ital("oleracea")+" / "+style.ital("rapa")+" = %.2f"%r_dip,
         ha="center",fontsize=7.4,fontweight="bold")
axb.text(1.0,1.37,"C$_n$ / A$_n$ inside "+style.ital("napus")+" = 1.47",
         ha="center",fontsize=7.0,color="black")
axb.set_title("the two comparable progenitors already differ by that factor",
              fontsize=7.0,color=style.ANNOT,pad=5)
style.panel(axb,"b",dx=-0.26,dy=1.18)

axc=fig.add_subplot(gs[1,:])
mets=["fractionation","DCJ","Ka/Ks","LTR"]; w=0.2
alphas=[1.0,0.78,0.56,0.34]
for k,m in enumerate(mets):
    v=FEAT[m].values.astype(float)
    norm=(v-v.min())/(v.max()-v.min())*0.92+0.08
    for i,(s_,g) in enumerate(zip(FEAT["sub"],FEAT["genome"])):
        axc.bar(i+(k-1.5)*w, norm[i], width=w*0.88, color=style.GENOME[g],
                alpha=alphas[k], edgecolor="white", linewidth=0.35)
axc.set_xticks(range(6)); axc.set_xticklabels(FEAT["sub"],fontweight="bold",fontsize=8.5)
axc.set_ylabel("scaled to the range\nof each measure",fontsize=7.3)
axc.set_ylim(0,1.12); axc.set_yticks([0,0.5,1.0])
axc.legend(handles=[mpatches.Patch(facecolor=style.GREY,alpha=a,label=m)
                    for m,a in zip(mets,alphas)],
           ncol=4,fontsize=6.6,loc="lower center",bbox_to_anchor=(0.5,1.01),
           columnspacing=1.6,handlelength=1.1,frameon=False)
style.panel(axc,"c",dx=-0.075,dy=1.16)

for f in ("Fig2_structure.png","Fig2_structure.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig2_structure")
