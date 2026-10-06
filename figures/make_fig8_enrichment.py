#!/usr/bin/env python3
import os, sys, pickle, collections
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); import style; style.use()
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=HERE
TH=1.0
MINFG=8
BG_MIN,BG_MAX=10,400

GO=pickle.load(open(os.path.join(REB,"00_inputs/go_propagated.pkl"),"rb"))
G2GO, NAME, NS = GO["g2go"], GO["name"], GO["ns"]
SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
COLS=[("napus","A","An"),("juncea","A","Aj"),("juncea","B","Bj"),
      ("carinata","B","Bc"),("napus","C","Cn"),("carinata","C","Cc")]
TISSUE={"stigma":[("napus","stigma"),("juncea","stigma"),("carinata","stigmaL")],
        "pollen":[("napus","pollen"),("juncea","pollen"),("carinata","pollen")]}

def enrich(fg, bg):
    cf=collections.Counter(); cb=collections.Counter()
    for g in fg:
        for t in G2GO.get(g,()): cf[t]+=1
    for g in bg:
        for t in G2GO.get(g,()): cb[t]+=1
    NF,NB=len(fg),len(bg); out=[]
    for t,a in cf.items():
        if a<MINFG: continue
        b=cb.get(t,0)
        if b<a: continue
        odds,p=fisher_exact([[a,NF-a],[b-a,NB-NF-(b-a)]],alternative="greater")
        out.append((t,a,b,p))
    out.sort(key=lambda x:x[3])
    m=len(out); res=[]
    for i,(t,a,b,p) in enumerate(out,1):
        res.append(dict(term=t,term_name=NAME.get(t,t),ns=NS.get(t,""),fg=a,bg=b,p=p,fdr=min(1,p*m/i)))
    return pd.DataFrame(res)

records=[]
for tis,samples in TISSUE.items():
    for sp,ti in samples:
        f=os.path.join(REB,"02_heb/theta1_a05/pairs_%s_%s.csv"%(sp,ti))
        if not os.path.exists(f): continue
        d=pd.read_csv(f)
        d=d[np.isfinite(d.allo_lfc)&np.isfinite(d.par_lfc)]
        d["at_gene"]=d["key"].str.split("|").str[0]
        aB=d.allo_lfc.abs()>=TH; pB=d.par_lfc.abs()>=TH
        same=np.sign(d.allo_lfc)==np.sign(d.par_lfc)
        novel_or_switched = (aB & ~pB) | (aB & pB & ~same)
        bg=set(d["at_gene"])
        g1,g2=SUBS[sp]
        for g,sign in ((g1,1),(g2,-1)):
            fg=set(d.loc[novel_or_switched & (np.sign(d.allo_lfc)==sign),"at_gene"])
            if len(fg)<25: continue
            e=enrich(fg,bg)
            if e.empty: continue
            e["tissue"]=tis; e["species"]=sp; e["genome"]=g
            e["col"]=[c[2] for c in COLS if c[0]==sp and c[1]==g][0]
            records.append(e)
E=pd.concat(records,ignore_index=True)
E.to_csv(os.path.join(REB,"05_validation/enrichment_novel_switched.csv"),index=False)
print("tested terms:",len(E),"| significant (FDR<0.05):",int((E.fdr<0.05).sum()))

sig=E[(E.fdr<0.05)&(E.tissue=="stigma")&(E.bg>=BG_MIN)&(E.bg<=BG_MAX)]
top=(sig.groupby("term_name").fdr.min().sort_values().head(14).index.tolist())
colnames=[c[2] for c in COLS]
M=pd.DataFrame(np.nan,index=top,columns=colnames)
for _,r in E[(E.tissue=="stigma")&(E.term_name.isin(top))&(E.bg>=BG_MIN)&(E.bg<=BG_MAX)].iterrows():
    M.loc[r["term_name"],r["col"]]=-np.log10(max(r.fdr,1e-16))
fig,ax=plt.subplots(figsize=(7.2,max(3.1,0.30*len(top)+1.9)))
cmap=plt.get_cmap("Reds").copy(); cmap.set_bad("#F2F2F2")
im=ax.imshow(M.values,cmap=cmap,aspect="auto",vmin=0,
             vmax=np.nanmax(M.values) if np.isfinite(np.nanmax(M.values)) else 1)
ax.set_xticks(range(len(colnames)))
ax.set_xticklabels(["A$_n$","A$_j$","B$_j$","B$_c$","C$_n$","C$_c$"],fontsize=8)
for i,c in enumerate(colnames):
    ax.get_xticklabels()[i].set_color(style.GENOME[c[0]])
ax.set_yticks(range(len(top))); ax.set_yticklabels(top,fontsize=6.8)
ax.set_title("stigma: enriched terms among genes with novel or switched bias\n"
             "propagated GO, Fisher's exact test, Benjamini-Hochberg, terms of 10 to 400 genes",
             fontsize=7.6,pad=6)
cb=fig.colorbar(im,ax=ax,fraction=0.030,pad=0.015)
cb.set_label("$-$log$_{10}$ FDR",fontsize=7); cb.ax.tick_params(labelsize=6)
cb.outline.set_linewidth(0)
ax.set_xticks(np.arange(-0.5,len(colnames),1),minor=True)
ax.set_yticks(np.arange(-0.5,len(top),1),minor=True)
ax.grid(which="minor",color="white",linewidth=1.1); ax.tick_params(which="minor",length=0)
fig.subplots_adjust(left=0.42,right=0.93,top=0.86,bottom=0.07)
for f in ("Fig8_enrichment.png","Fig8_enrichment.pdf"): fig.savefig(os.path.join(OUT,f))
print("wrote Fig8_enrichment with %d terms"%len(top))
print(sig.groupby("col").size().to_string())
