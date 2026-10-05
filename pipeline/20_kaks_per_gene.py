#!/usr/bin/env python3
import glob, os, csv, re, statistics
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon

LAB={'68107_68114':('napus','A'),'68103_68114':('juncea','A'),'68108_68112':('juncea','B'),
     '68110_68112':('carinata','B'),'68109_68113':('napus','C'),'68111_68113':('carinata','C')}
SP={"napus":dict(col1="Bna_A_subgenome%d",col2="Bna_C_subgenome%d",g1="A",g2="C",
                 tissues=["pollen","stigma"]),
    "juncea":dict(col1="Bju_A_subgenome_%d",col2="Bju_B_subgenome_%d",g1="A",g2="B",
                 tissues=["pollen","stigma"]),
    "carinata":dict(col1="Bca_B_subgenome_%d",col2="Bca_C_subgenome_%d",g1="B",g2="C",
                 tissues=["pollen","stigmaE","stigmaL"])}
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")

def per_gene(path, lo=91, hi=100):
    out={}; blk=[]; sim=[]
    for num,line in enumerate(open(path)):
        if num in (0,1,2): continue
        if line[0]!='#':
            p=line.rstrip("\n").split("\t")
            ks,kn=p[0],p[1]
            sim.append(float(p[3].split('||')[8]))
            gid=p[3].split('||')[3]
            if ks not in ('NA','undef') and kn not in ('NA','undef'):
                k,n=float(ks),float(kn)
                if 0.001<k<=4.1 and n>=0:
                    r=n/k
                    if np.isfinite(r) and r<10: blk.append((gid,r))
        else:
            if blk and sim and lo<=statistics.mean(sim)<=hi:
                for g,r in blk: out.setdefault(STRIP.sub("",g),[]).append(r)
            blk=[]; sim=[]
    return {g:float(np.median(v)) for g,v in out.items()}

KK={}
for f in sorted(glob.glob("00_inputs/synmap/*.ks.txt")):
    k=os.path.basename(f).split('.')[0]
    if k in LAB: KK[LAB[k]]=per_gene(f)
for key,d in KK.items(): print("  %-9s %-2s  genes with Ka/Ks: %d" % (key[0],key[1],len(d)))
rows=list(csv.reader(open("00_inputs/Subgenomes_Brassica.txt"),delimiter="\t"))
H={x.strip():i for i,x in enumerate(rows[0])}
print()
print("Prediction: NEGATIVE rho (the less constrained copy is the lower expressed one)\n")
print("%-10s %-9s %7s | %8s %9s | %-30s" % ("species","tissue","n","rho","p","biased pairs: suppressed copy"))
print("-"*94)
out=[]
for sp,c in SP.items():
    d1=KK.get((sp,c["g1"]),{}); d2=KK.get((sp,c["g2"]),{})
    if not d1 or not d2: print("%-10s missing Ka/Ks"%sp); continue
    keys=[];v1=[];v2=[]
    for r in rows[1:]:
        for s in (1,2,3):
            a=r[H[c["col1"]%s]].strip(); b=r[H[c["col2"]%s]].strip()
            if a in ABSENT or b in ABSENT: continue
            a=STRIP.sub("",a); b=STRIP.sub("",b)
            if a in d1 and b in d2:
                keys.append("%s|%d"%(r[H["AT_geneid"]],s)); v1.append(d1[a]); v2.append(d2[b])
    df=pd.DataFrame({"key":keys,"k1":v1,"k2":v2}); df["dk"]=df["k1"]-df["k2"]
    for t in c["tissues"]:
        f="02_heb/theta2_a05/pairs_%s_%s.csv"%(sp,t)
        if not os.path.exists(f): continue
        h=pd.read_csv(f)
        m=df.merge(h[["key","allo","allo_lfc"]],on="key",how="inner").dropna(subset=["dk","allo_lfc"])
        if len(m)<50: continue
        rho,p=spearmanr(m["dk"],m["allo_lfc"])
        b=m[m["allo"].isin(["toward1","toward2"])].copy()
        b["sup"]=np.where(b["allo"]=="toward1",b["k2"],b["k1"])
        b["dom"]=np.where(b["allo"]=="toward1",b["k1"],b["k2"])
        dd=(b["sup"]-b["dom"]); nz=dd[dd!=0]
        txt=("higher Ka/Ks in %.1f%% (n=%d)"%(100*(nz>0).mean(),len(nz))) if len(nz)>10 else "too few"
        print("%-10s %-9s %7d | %8.4f %9.3g | %-30s" % (sp,t,len(m),rho,p,txt))
        out.append(dict(species=sp,tissue=t,n=len(m),rho=rho,p=p,
                        frac_sup_higher_kaks=float((nz>0).mean()) if len(nz)>10 else None))
pd.DataFrame(out).to_csv("05_validation/kaks_per_gene.csv",index=False)
print("\nWrote 05_validation/kaks_per_gene.csv")
