#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=os.path.dirname(os.path.abspath(__file__))

ltr=pd.read_csv(os.path.join(REB,"05_validation/ltr_vs_bias_all_species.csv"))
cen=pd.read_csv(os.path.join(REB,"05_validation/centromere_vs_bias.csv"))
kk =pd.read_csv(os.path.join(REB,"05_validation/kaks_per_gene.csv"))
er =pd.read_csv(os.path.join(REB,"05_validation/er_control.csv"))
def key(d):
    t=d.tissue.replace("stigmaE","stigma (early)").replace("stigmaL","stigma (late)")
    return d.species+" "+t
for d in (ltr,cen,kk,er): d["lab"]=key(d)
order=ltr.sort_values(["species","tissue"]).lab.tolist()

SERIES=[("flanking LTR density", ltr.set_index("lab").rho, "#0072B2"),
        ("centromere distance",  cen.set_index("lab").rho, "#009E73"),
        ("Ka/Ks, raw",           kk.set_index("lab").rho,  "#D55E00"),
        ("Ka/Ks, expression-controlled", er.set_index("lab").residual_rho, "#555555")]

fig,(ax,ax2)=plt.subplots(1,2,figsize=(7.2,3.1),gridspec_kw={"width_ratios":[1.55,1.0],"wspace":0.42})
y=np.arange(len(order))[::-1]
off=[0.24,0.08,-0.08,-0.24]
for (nm,s,col),o in zip(SERIES,off):
    v=[s.get(l,np.nan) for l in order]
    ax.scatter(v,y+o,s=26,color=col,label=nm,edgecolor="white",linewidth=0.5,zorder=3)
ax.axvline(0,color="black",lw=0.8,zorder=1)
ax.axvspan(-0.05,0.05,color="#F0F0F0",zorder=0)
ax.set_yticks(y); ax.set_yticklabels([style.ital(v) for v in order],fontsize=6.6)
ax.set_xlabel("Spearman $\\rho$ with homoeolog expression bias")
ax.set_xlim(-0.16,0.10)
ax.legend(fontsize=6.3,loc="lower left",bbox_to_anchor=(0.0,1.045),ncol=2,
          columnspacing=1.0,handletextpad=0.3)
ax.text(0.5,1.005,"shaded band: |$\\rho$| < 0.05",transform=ax.transAxes,
        fontsize=6.2,color=style.ANNOT,ha="center",va="bottom")
style.panel(ax,"a",dx=-0.30,dy=1.20)

sub=pd.DataFrame({"metric":[style.ital("LTR density\nC vs A (napus)"),style.ital("LTR density\nC vs B (carinata)"),
                            style.ital("Ka/Ks\nC vs A (napus)"),style.ital("Ka/Ks\nC vs B (carinata)")],
                  "ratio":[0.008063/0.005524, 0.0085/0.0041, 0.2783/0.1825, 0.2893/0.1921]})
ax2.barh(np.arange(len(sub))[::-1], sub.ratio, color=["#0072B2","#0072B2","#D55E00","#D55E00"],
         height=0.6, edgecolor="white", linewidth=0.5)
ax2.axvline(1.0,color="black",lw=0.8)
for i,(m,r) in enumerate(zip(sub.metric,sub.ratio)):
    ax2.text(r+0.04, len(sub)-1-i, "%.2f$\\times$"%r, va="center", fontsize=6.6)
ax2.set_yticks(np.arange(len(sub))[::-1]); ax2.set_yticklabels(sub.metric,fontsize=6.4)
ax2.set_xlabel("subgenome-level ratio"); ax2.set_xlim(0,2.5)
ax2.set_title("the same measures BETWEEN subgenomes differ 1.5-2.1$\times$",
              fontsize=7,color=style.ANNOT,pad=4)
style.panel(ax2,"b",dx=-0.44,dy=1.20)

for f in ("Fig5_genelevel.png","Fig5_genelevel.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig5_genelevel")
