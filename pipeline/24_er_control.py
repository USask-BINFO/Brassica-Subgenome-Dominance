#!/usr/bin/env python3
import glob, os, csv, re, statistics
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

LAB={'68107_68114':('napus','A'),'68103_68114':('juncea','A'),'68108_68112':('juncea','B'),
     '68110_68112':('carinata','B'),'68109_68113':('napus','C'),'68111_68113':('carinata','C')}
SP={"napus":dict(c1="Bna_A",c2="Bna_C",col1="Bna_A_subgenome%d",col2="Bna_C_subgenome%d",
        tis={"pollen":["Bna_PolA","Bna_PolB","Bna_PolC"],"stigma":["Bna_StigA","Bna_StigB","Bna_StigC"]}),
    "juncea":dict(c1="Bju_A",c2="Bju_B",col1="Bju_A_subgenome_%d",col2="Bju_B_subgenome_%d",
        tis={"pollen":["Bju_PolA","Bju_PolB","Bju_PolC"],"stigma":["Bju_StigA","Bju_StigB","Bju_StigC"]}),
    "carinata":dict(c1="Bca_B",c2="Bca_C",col1="Bca_B_subgenome_%d",col2="Bca_C_subgenome_%d",
        tis={"pollen":["Bca_PolA","Bca_PolB","Bca_PolC"],
             "stigmaE":["Bca_E_StigA","Bca_E_StigB","Bca_E_StigC"],
             "stigmaL":["Bca_L_StigA","Bca_L_StigB","Bca_L_StigC"]})}
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")

def per_gene_kaks(path, lo=91, hi=100):
    out={}; blk=[]; sim=[]
    for num,line in enumerate(open(path)):
        if num in (0,1,2): continue
        if line[0]!='#':
            p=line.rstrip("\n").split("\t"); ks,kn=p[0],p[1]
            sim.append(float(p[3].split('||')[8])); gid=p[3].split('||')[3]
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
    if k in LAB: KK[LAB[k]]=per_gene_kaks(f)
CNT=pd.read_csv("01_counts/gref_counts.tsv",sep="\t")
CNT["key"]=CNT["AT_gene"].astype(str)+"|"+CNT["br_sub"].astype(str)
CNT=CNT.set_index("key")
rows=list(csv.reader(open("00_inputs/Subgenomes_Brassica.txt"),delimiter="\t"))
H={x.strip():i for i,x in enumerate(rows[0])}

print("%-10s %-9s %7s | %-24s | %8s %8s %9s" %
      ("species","tissue","n","global E-R rho (all copies)","raw rho","resid rho","change"))
print("-"*96)
out=[]
for sp,c in SP.items():
    d1=KK.get((sp,c["c1"].split('_')[1] if False else {"Bna_A":"A","Bna_C":"C","Bju_A":"A","Bju_B":"B","Bca_B":"B","Bca_C":"C"}[c["c1"]]),{})
    d2=KK.get((sp,{"Bna_A":"A","Bna_C":"C","Bju_A":"A","Bju_B":"B","Bca_B":"B","Bca_C":"C"}[c["c2"]]),{})
    if not d1 or not d2: continue
    pairs=[]
    for r in rows[1:]:
        for s in (1,2,3):
            a=r[H[c["col1"]%s]].strip(); b=r[H[c["col2"]%s]].strip()
            if a in ABSENT or b in ABSENT: continue
            a=STRIP.sub("",a); b=STRIP.sub("",b)
            if a in d1 and b in d2:
                pairs.append(("%s|%d"%(r[H["AT_geneid"]],s), d1[a], d2[b]))
    P=pd.DataFrame(pairs, columns=["key","k1","k2"])
    for t,libs in c["tis"].items():
        cc1=[c["c1"]+"|"+l for l in libs]; cc2=[c["c2"]+"|"+l for l in libs]
        if not all(x in CNT.columns for x in cc1+cc2): continue
        e1=CNT.reindex(P["key"])[cc1].mean(axis=1).values
        e2=CNT.reindex(P["key"])[cc2].mean(axis=1).values
        df=P.copy(); df["e1"]=e1; df["e2"]=e2
        df=df[(df.e1>0)&(df.e2>0)].dropna()
        if len(df)<200: continue
        ex=np.log10(np.concatenate([df.e1.values,df.e2.values]))
        kk=np.concatenate([df.k1.values,df.k2.values])
        er_rho,_=spearmanr(ex,kk)
        sm_fit=lowess(kk,ex,frac=0.3,return_sorted=True)
        pred=np.interp(ex,sm_fit[:,0],sm_fit[:,1])
        res=kk-pred
        n=len(df); r1=res[:n]; r2=res[n:]
        f="02_heb/theta2_a05/pairs_%s_%s.csv"%(sp,t)
        if not os.path.exists(f): continue
        h=pd.read_csv(f)
        df["dk"]=df.k1-df.k2; df["dres"]=r1-r2
        m=df.merge(h[["key","allo_lfc"]],on="key",how="inner").dropna(subset=["dk","dres","allo_lfc"])
        if len(m)<200: continue
        raw,_=spearmanr(m.dk,m.allo_lfc); rsd,p2=spearmanr(m.dres,m.allo_lfc)
        print("%-10s %-9s %7d | %-24.3f | %8.4f %8.4f %+9.1f%%" %
              (sp,t,len(m),er_rho,raw,rsd,100*(abs(rsd)-abs(raw))/abs(raw)))
        out.append(dict(species=sp,tissue=t,n=len(m),global_ER_rho=er_rho,
                        raw_rho=raw,residual_rho=rsd,residual_p=p2))
pd.DataFrame(out).to_csv("05_validation/er_control.csv",index=False)
print("\nglobal E-R rho: negative confirms the confounder exists (highly expressed copies evolve slower)")
print("if residual rho stays close to raw rho, section 15 is NOT merely the E-R anticorrelation")
