#!/usr/bin/env python3
import os, bisect, collections, csv, re, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon, binomtest

LTRDIR="00_inputs/ltr"; GFF="00_inputs/coge_gff"; HEB="02_heb/theta2_a05"; FLANK=5000
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")
SP={
 "juncea": dict(ltr="Bjuncea_genome_v1.fa.out.LTR.gff",
   g1="A", g2="B",
   gff1="Brassica_juncea_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68103.gff",
   gff2="Brassica_juncea_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68108.gff",
   col1="Bju_A_subgenome_%d", col2="Bju_B_subgenome_%d", tissues=["pollen","stigma"]),
 "carinata": dict(ltr="Bcarinata.v1.genome.fasta.out.LTR.gff",
   g1="B", g2="C",
   gff1="Brassica_carinata_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68110.gff",
   gff2="Brassica_carinata_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68111.gff",
   col1="Bca_B_subgenome_%d", col2="Bca_C_subgenome_%d", tissues=["pollen","stigmaE","stigmaL"]),
}
def load_ltr(p):
    iv=collections.defaultdict(list)
    for line in open(p):
        if line.startswith("#"): continue
        f=line.split("\t")
        if len(f)<5: continue
        try: iv[f[0]].append((int(f[3]),int(f[4])))
        except ValueError: continue
    return {c:([x[0] for x in sorted(l)], sorted(l)) for c,l in iv.items()}
def ltr_bp(idx,c,s,e):
    if c not in idx or e<s: return 0
    starts,l=idx[c]; i=bisect.bisect_left(starts,s)-1; tot=0
    while i<len(l):
        if i>=0:
            a,b=l[i]
            if a>e: break
            tot+=max(0,min(b,e)-max(a,s)+1)
        i+=1
    return tot
def genes_gff(p):
    g={}
    for line in open(p):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="gene": continue
        c=f[0]
        if c.startswith(("Scaffold","utg","contig")): continue
        gid=None
        for kv in f[8].split(";"):
            if kv.startswith("ID="): gid=kv[3:]; break
        if gid is None: continue
        g[STRIP.sub("",gid)]=(c,int(f[3]),int(f[4]))
    return g
rows=list(csv.reader(open("00_inputs/Subgenomes_Brassica.txt"),delimiter="\t"))
H={x.strip():i for i,x in enumerate(rows[0])}
print("%-10s %-9s %7s %8s | %8s %10s | %-26s" % ("species","tissue","n","join%","rho","p","biased: suppressed copy"))
print("-"*94)
out=[]
for sp,c in SP.items():
    idx=load_ltr(os.path.join(LTRDIR,c["ltr"]))
    gene={}; gene.update(genes_gff(os.path.join(GFF,c["gff1"]))); gene.update(genes_gff(os.path.join(GFF,c["gff2"])))
    rec={}
    for r in rows[1:]:
        for s in (1,2,3):
            a=r[H[c["col1"]%s]].strip(); b=r[H[c["col2"]%s]].strip()
            if a in ABSENT or b in ABSENT: continue
            rec["%s|%d"%(r[H["AT_geneid"]],s)]=(STRIP.sub("",a),STRIP.sub("",b))
    cache={}
    def dens(g):
        if g in cache: return cache[g]
        if g not in gene: cache[g]=np.nan; return np.nan
        ch,s,e=gene[g]
        v=(ltr_bp(idx,ch,max(0,s-FLANK),s-1)+ltr_bp(idx,ch,e+1,e+FLANK))/(2.0*FLANK)
        cache[g]=v; return v
    k=[];d1=[];d2=[]
    for key,(a,b) in rec.items():
        x,y=dens(a),dens(b)
        if np.isnan(x) or np.isnan(y): continue
        k.append(key); d1.append(x); d2.append(y)
    pr=pd.DataFrame({"key":k,"ltr1":d1,"ltr2":d2}); pr["ltr_diff"]=pr["ltr1"]-pr["ltr2"]
    jr=100.0*len(pr)/max(1,len(rec))
    print("%-10s %-9s %7s %7.1f%% | %-8s %10s | mean %s %.4f  %s %.4f"
          % (sp,"(subgenome)",len(pr),jr,"","",c["g1"],np.mean(d1),c["g2"],np.mean(d2)))
    for t in c["tissues"]:
        f="%s/pairs_%s_%s.csv"%(HEB,sp,t)
        if not os.path.exists(f): continue
        h=pd.read_csv(f)
        m=pr.merge(h[["key","allo","allo_lfc"]],on="key",how="inner").dropna(subset=["ltr_diff","allo_lfc"])
        if len(m)<50: continue
        rho,p=spearmanr(m["ltr_diff"],m["allo_lfc"])
        b=m[m["allo"].isin(["toward1","toward2"])].copy()
        b["sup"]=np.where(b["allo"]=="toward1",b["ltr2"],b["ltr1"])
        b["dom"]=np.where(b["allo"]=="toward1",b["ltr1"],b["ltr2"])
        d=(b["sup"]-b["dom"]); nz=d[d!=0]
        txt=("more LTR in %.1f%% (n=%d, p=%.2g)"%(100*(nz>0).mean(),len(nz),binomtest(int((nz>0).sum()),len(nz)).pvalue)) if len(nz)>10 else "too few"
        print("%-10s %-9s %7d %7.1f%% | %8.4f %10.3g | %-26s"%(sp,t,len(m),jr,rho,p,txt))
        out.append(dict(species=sp,tissue=t,n=len(m),join_pct=jr,rho=rho,p=p,source="coge_gff"))
pd.DataFrame(out).to_csv("05_validation/ltr_vs_bias_coge_sensitivity.csv",index=False)
print("\nCompare with scripts/12 (BWA sub BEDs). Wrote 05_validation/ltr_vs_bias_coge_sensitivity.csv")
