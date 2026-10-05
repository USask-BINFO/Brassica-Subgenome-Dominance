#!/usr/bin/env python3
import os, sys, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binomtest
from scipy.cluster.hierarchy import linkage, dendrogram
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=os.path.dirname(os.path.abspath(__file__))
SPECIES_OF={"Bra":"rapa","Bni":"nigra","Bol":"oleracea","Bna":"napus","Bju":"juncea","Bca":"carinata"}
def tissue_of(l):
    if "Pol" in l: return "pollen"
    if "_E_Stig" in l: return "stigmaE"
    if "_L_Stig" in l: return "stigmaL"
    return "stigma"
def save(fig,name):
    for e in ("png","pdf"): fig.savefig(os.path.join(OUT,"%s.%s"%(name,e)))
    plt.close(fig); print("wrote",name)

cnt=pd.read_csv(os.path.join(REB,"01_counts/gref_counts.tsv"),sep="\t")
meta=[c for c in cnt.columns if "|" in c]
libs=sorted({c.split("|")[1] for c in meta})
M=pd.DataFrame({L:cnt[[c for c in meta if c.split("|")[1]==L]].sum(axis=1) for L in libs})
sp=[SPECIES_OF[L.split("_")[0]] for L in M.columns]; ti=[tissue_of(L) for L in M.columns]

import matplotlib.transforms as mtransforms
fig,ax=plt.subplots(figsize=(7.6,3.4))
expressed=(M>=3).sum(axis=0)
x=np.arange(len(M.columns))
ax.bar(x,expressed.values,color=[style.SPECIES[s] for s in sp],width=0.72,
       edgecolor="white",linewidth=0.3)
ax.set_xticks(x); ax.set_xticklabels(M.columns,rotation=90,fontsize=5.0)
ax.set_ylabel("Gref gene copies with $\geq$ 3 reads",fontsize=7.6)
ax.set_title("expressed gene copies per library; stigma exceeds pollen in every species",
             fontsize=7.6,color=style.ANNOT,pad=5)
ax.set_ylim(0,expressed.max()*1.06); ax.set_xlim(-1.0,len(x)-0.0)
tr=mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
for s_ in dict.fromkeys(sp):
    idx=[i for i in range(len(sp)) if sp[i]==s_]
    ax.plot([min(idx)-0.38,max(idx)+0.38],[-0.255,-0.255],transform=tr,
            color=style.SPECIES[s_],lw=3.2,solid_capstyle="butt",clip_on=False,zorder=5)
    ax.text(np.mean(idx),-0.30,s_,transform=tr,ha="center",va="top",fontsize=7.0,
            color=style.SPECIES[s_],fontstyle="italic",fontweight="bold",clip_on=False)
fig.subplots_adjust(left=0.095,right=0.99,top=0.90,bottom=0.30)
save(fig,"FigS1_expressed_copies")

cpm=M/M.sum(axis=0)*1e6; lg=np.log2(cpm+1)
keep=lg.loc[lg.var(axis=1).sort_values(ascending=False).index[:2000]]
R=np.corrcoef(keep.T.rank(axis=1).values)
Z=linkage(1-R,method="average"); order=dendrogram(Z,no_plot=True)["leaves"]
fig,ax=plt.subplots(figsize=(6.8,6.3))
im=ax.imshow(R[np.ix_(order,order)],cmap="viridis",vmin=np.percentile(R,2),vmax=1)
ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
ax.set_xticklabels([M.columns[i] for i in order],rotation=90,fontsize=4.6)
ax.set_yticklabels([M.columns[i] for i in order],fontsize=4.6)
for i,o in enumerate(order):
    ax.add_patch(plt.Rectangle((i-0.5,-2.4),1,1.4,color=style.SPECIES[sp[o]],lw=0,clip_on=False))
cb=fig.colorbar(im,ax=ax,fraction=0.044,pad=0.02); cb.set_label("Spearman $r$",fontsize=7,color=style.ANNOT)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0)
fig.suptitle("library correlation, hierarchically clustered; the strip marks species",
             fontsize=8.5,y=0.985)
save(fig,"FigS2_replicate_correlation")

THRESH=[("theta0585_a05","1.5-fold"),("theta1_a05","2-fold"),("span_theta1_a05","2-fold span"),
        ("theta1585_a05","3-fold"),("Q10_theta1_a05","2-fold MAPQ10")]
ORDER=["juncea_pollen","juncea_stigma","napus_pollen","napus_stigma",
       "carinata_pollen","carinata_stigmaE","carinata_stigmaL"]
SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
D={t:pd.read_csv(os.path.join(REB,"02_heb",t,"heb_summary.csv")).set_index("sample")
   for t,_ in THRESH if os.path.exists(os.path.join(REB,"02_heb",t,"heb_summary.csv"))}
fig,ax=plt.subplots(figsize=(7.2,3.5))
mk=["o","s","^","D","v"]
grey=["#1A1A1A","#4D4D4D","#7A7A7A","#A5A5A5","#C8C8C8"]
NICE={"stigmaE":"stigma (early)","stigmaL":"stigma (late)"}
def pretty(s_):
    sp,ti=s_.split("_",1); return "%s %s"%(sp,NICE.get(ti,ti))
for k,(t,lab) in enumerate(THRESH):
    if t not in D: continue
    xs=[];ys=[]
    for i,s_ in enumerate(ORDER):
        if s_ not in D[t].index: continue
        r=D[t].loc[s_]; xs.append(r.toward1/r.toward2); ys.append(len(ORDER)-1-i)
    ax.scatter(xs,ys,s=30,marker=mk[k],label=lab,color=grey[k],alpha=0.95,
               edgecolor="white",linewidth=0.4)
ax.axvline(1,color="black",lw=0.9)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([style.ital(pretty(s_)) for s_ in ORDER][::-1],fontsize=7)
ax.set_xlabel("ratio of biased pairs, subgenome 1 : subgenome 2",fontsize=7.4)
ax.legend(fontsize=6.4,ncol=5,loc="lower center",bbox_to_anchor=(0.5,1.015),
          frameon=False,columnspacing=1.3,handletextpad=0.35)
fig.text(0.5,0.965,"the direction of every contest is unchanged by all five settings",
         ha="center",va="bottom",fontsize=7.2,color=style.ANNOT)
fig.subplots_adjust(left=0.20,right=0.985,top=0.80,bottom=0.155)
save(fig,"FigS3_threshold_sensitivity")

MB=pd.read_csv(os.path.join(REB,"09_mapbias","mapping_bias.csv"))
MB=MB.rename(columns={"fwd_pct":"fwd","rev_pct":"rev"})
MB["contrast"]=MB["contrast"].str.replace("<->","$\\leftrightarrow$",regex=False)
MB["asym"]=MB.fwd-MB.rev
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.2,2.7),gridspec_kw={"wspace":0.34})
w=0.34; x=np.arange(3)
a1.bar(x-w/2,MB.fwd,width=w,label="copy1 $\\to$ copy2",color=style.GREY)
a1.bar(x+w/2,MB.rev,width=w,label="copy2 $\\to$ copy1",color=style.LIGHT)
a1.set_xticks(x); a1.set_xticklabels([style.ital(a)+"\n"+b for a,b in zip(MB.species,MB.contrast)],fontsize=6.6)
a1.set_ylabel("cross-mapping (% of reads)"); a1.legend(fontsize=6.4)
style.panel(a1,"a",dx=-0.20)
a2.bar(x,MB.asym,color=[style.SPECIES[s] for s in MB.species],width=0.5)
a2.axhline(0,color="black",lw=0.8)
a2.set_xticks(x); a2.set_xticklabels([style.ital(v) for v in MB.species],fontsize=6.8)
a2.set_ylabel("asymmetry (percentage points)"); a2.set_ylim(0,0.55)
a2.text(1,0.50,"largest imbalance 0.47 pp",ha="center",fontsize=6.6,color=style.ANNOT)
style.panel(a2,"b",dx=-0.22)
save(fig,"FigS4_mapping_bias")

res=pd.read_csv(os.path.join(REB,"05_validation/inherited_vs_novel.csv"))
fig,axes=plt.subplots(2,4,figsize=(7.2,3.9))
for ax,(_,r) in zip(axes.ravel(),res.iterrows()):
    f=os.path.join(REB,"02_heb/theta2_a05/pairs_%s_%s.csv"%(r.species,r.tissue))
    d=pd.read_csv(f); d=d[np.isfinite(d.allo_lfc)&np.isfinite(d.par_lfc)]
    lim=6
    ax.hexbin(d.par_lfc.clip(-lim,lim),d.allo_lfc.clip(-lim,lim),gridsize=28,cmap="Greys",
              bins="log",mincnt=1,linewidths=0)
    ax.plot([-lim,lim],[-lim,lim],color=style.LIGHT,lw=0.7,ls="--")
    xs=np.array([-lim,lim]); ax.plot(xs,r.slope*xs,color=style.SPECIES[r.species],lw=1.3)
    ax.set_title(style.ital("%s %s"%(r.species,r.tissue)),fontsize=6.4,pad=2)
    ax.text(0.04,0.96,"%.2f / %.2f"%(r.slope,r.R2),transform=ax.transAxes,fontsize=5.8,va="top")
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.tick_params(labelsize=5.5)
axes.ravel()[-1].axis("off")
axes.ravel()[-1].text(0.0,0.5,"slope / $R^2$ on\neach panel\n\ndashed line 1:1",fontsize=6.4,va="center")
fig.supxlabel("bias between the diploid progenitors (log$_2$)",fontsize=7.5)
fig.supylabel("bias inside the allotetraploid (log$_2$)",fontsize=7.5)
fig.tight_layout()
save(fig,"FigS7_inheritance_all")

a=pd.read_csv(os.path.join(REB,"05_validation/ltr_vs_bias_all_species.csv"))
b=pd.read_csv(os.path.join(REB,"05_validation/ltr_vs_bias_coge_sensitivity.csv"))
a["lab"]=a.species+" "+a.tissue; b["lab"]=b.species+" "+b.tissue
m=a.merge(b,on="lab",suffixes=("_bed","_coge"))
fig,ax=plt.subplots(figsize=(4.6,3.0))
ax.scatter(m.rho_bed,m.rho_coge,s=44,color=[style.SPECIES[s] for s in m.species_bed],
           edgecolor="white",linewidth=0.6,zorder=3)
lo,hi=-0.045,0.025
ax.plot([lo,hi],[lo,hi],color=style.LIGHT,ls="--",lw=0.8)
ax.axhline(0,color=style.LIGHT,lw=0.5); ax.axvline(0,color=style.LIGHT,lw=0.5)
ax.set_xlabel("$\\rho$, Br-subgenome BED coordinates"); ax.set_ylabel("$\\rho$, complete CoGe coordinates")
ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
ax.set_title("raising the join rate to 92-98% changes nothing",fontsize=7,color=style.ANNOT,pad=4)
save(fig,"FigS8_ltr_join_sensitivity")

er=pd.read_csv(os.path.join(REB,"05_validation/er_control.csv"))
er["lab"]=er.species+" "+er.tissue
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.2,3.1),gridspec_kw={"wspace":0.34})
y=np.arange(len(er))[::-1]
a1.barh(y,er.global_ER_rho,color=style.GREY,height=0.6)
a1.set_yticks(y); a1.set_yticklabels([style.ital(v) for v in er.lab],fontsize=6.2)
a1.set_xlabel("global expression-rate $\\rho$\n(all gene copies)")
a1.axvline(0,color="black",lw=0.8)
style.panel(a1,"a",dx=-0.52)
a2.scatter(er.raw_rho,y+0.12,s=30,color="#D55E00",label="raw",zorder=3)
a2.scatter(er.residual_rho,y-0.12,s=30,color=style.GREY,label="after control",zorder=3)
for i,r in enumerate(er.itertuples()):
    a2.plot([r.raw_rho,r.residual_rho],[y[i]+0.12,y[i]-0.12],color=style.LIGHT,lw=0.7,zorder=1)
a2.axvline(0,color="black",lw=0.8); a2.axvspan(-0.05,0.05,color="#F0F0F0",zorder=0)
a2.set_yticks(y); a2.set_yticklabels([]); a2.set_xlabel("within-pair $\\rho$ (Ka/Ks vs bias)")
a2.set_xticks([-0.12,-0.08,-0.04,0.0,0.04])
a2.legend(fontsize=6.4,loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2,
          frameon=False,columnspacing=1.4,handletextpad=0.35)
style.panel(a2,"b",dx=-0.10)
save(fig,"FigS9_er_control")
print("done")
