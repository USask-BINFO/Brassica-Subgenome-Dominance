#!/usr/bin/env python3
import os, csv, re, glob, collections
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

BED="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new/Research/BWA"
CEN="00_inputs/Brassica_centromere"
HEB="02_heb/theta2_a05"
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")

SP={
 "napus":   dict(cen="napus_AACC_centromere.txt", bed="Bna", pre="Bna_sub", g=["A","C"],
                 col1="Bna_A_subgenome%d",   col2="Bna_C_subgenome%d",
                 remap=lambda c: ("A"+c[1:]) if c.startswith("N") and int(c[1:])<=10
                                 else (("C"+c[1:]) if c.startswith("N") else c),
                 tissues=["pollen","stigma"]),
 "juncea":  dict(cen="juncea_AABB_centromere.txt", bed="Bju", pre="Bju_sub", g=["A","B"],
                 col1="Bju_A_subgenome_%d", col2="Bju_B_subgenome_%d",
                 remap=lambda c: c, tissues=["pollen","stigma"]),
 "carinata":dict(cen="carinata_BBCC_centromere.txt", bed="Bca", pre="Bca_sub", g=["B","C"],
                 col1="Bca_B_subgenome_%d", col2="Bca_C_subgenome_%d",
                 remap=lambda c: c, tissues=["pollen","stigmaE","stigmaL"]),
}
def centromeres(p):
    d={}
    for line in open(p):
        f=line.rstrip("\n").split("\t")
        if len(f)<3: continue
        try: d[f[0]]=(int(f[1]),int(f[2]))
        except ValueError: continue
    return d
def genes(bdir,pre,groups,remap):
    g={}
    for s in (1,2,3):
        for grp in groups:
            p="%s/%s/%s%d%s.bed"%(BED,bdir,pre,s,grp)
            if not os.path.exists(p): continue
            for line in open(p):
                f=line.rstrip("\n").split("\t")
                if len(f)<6: continue
                g[f[3]]=(remap(f[0]), (int(f[1])+int(f[2]))//2)
    return g
rows=list(csv.reader(open("00_inputs/Subgenomes_Brassica.txt"),delimiter="\t"))
H={x.strip():i for i,x in enumerate(rows[0])}
print("Prediction: POSITIVE rho (copy further from the centromere is the higher expressed one)\n")
print("%-10s %-9s %7s | %8s %9s | %-34s" %
      ("species","tissue","n","rho","p","subgenome mean distance (Mb)"))
print("-"*92)
out=[]
for sp,c in SP.items():
    CM=centromeres(os.path.join(CEN,c["cen"]))
    G=genes(c["bed"],c["pre"],c["g"],c["remap"])
    def dist(gid):
        if gid not in G: return np.nan
        ch,mid=G[gid]
        if ch not in CM: return np.nan
        a,b=CM[ch]
        return 0.0 if a<=mid<=b else float(min(abs(mid-a),abs(mid-b)))
    rec={}
    for r in rows[1:]:
        for s in (1,2,3):
            a=r[H[c["col1"]%s]].strip(); b=r[H[c["col2"]%s]].strip()
            if a in ABSENT or b in ABSENT: continue
            rec["%s|%d"%(r[H["AT_geneid"]],s)]=(STRIP.sub("",a),STRIP.sub("",b))
    keys=[];d1=[];d2=[]
    for k,(a,b) in rec.items():
        x,y=dist(a),dist(b)
        if np.isnan(x) or np.isnan(y): continue
        keys.append(k); d1.append(x); d2.append(y)
    if not keys:
        print("%-10s  no joined genes (check chromosome naming)"%sp); continue
    df=pd.DataFrame({"key":keys,"d1":d1,"d2":d2}); df["dd"]=df["d1"]-df["d2"]
    m1,m2=np.mean(d1)/1e6,np.mean(d2)/1e6
    for t in c["tissues"]:
        f="%s/pairs_%s_%s.csv"%(HEB,sp,t)
        if not os.path.exists(f): continue
        h=pd.read_csv(f)
        m=df.merge(h[["key","allo","allo_lfc"]],on="key",how="inner").dropna(subset=["dd","allo_lfc"])
        if len(m)<50: continue
        rho,p=spearmanr(m["dd"],m["allo_lfc"])
        print("%-10s %-9s %7d | %8.4f %9.3g | %s %.2f vs %s %.2f" %
              (sp,t,len(m),rho,p,c["g"][0],m1,c["g"][1],m2))
        out.append(dict(species=sp,tissue=t,n=len(m),rho=rho,p=p,
                        mean_dist_sub1_Mb=m1,mean_dist_sub2_Mb=m2))
pd.DataFrame(out).to_csv("05_validation/centromere_vs_bias.csv",index=False)
print("\nWrote 05_validation/centromere_vs_bias.csv")
