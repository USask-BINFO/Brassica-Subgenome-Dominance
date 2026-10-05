#!/usr/bin/env python3
import os, sys, glob
import numpy as np, pandas as pd
from scipy.stats import binomtest
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild"); OUT=os.path.dirname(os.path.abspath(__file__))

def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    lo,hi=c-h,c+h; return lo/(1-lo), hi/(1-hi)

def write(df,name,note=""):
    df.to_csv(os.path.join(OUT,name+".tsv"),sep="\t",index=False)
    with open(os.path.join(OUT,name+".md"),"w") as f:
        f.write("| "+" | ".join(df.columns)+" |\n")
        f.write("|"+"|".join([":-"]*len(df.columns))+"|\n")
        for _,r in df.iterrows():
            f.write("| "+" | ".join(str(x) for x in r.values)+" |\n")
        if note: f.write("\n"+note+"\n")
    print("wrote %s.tsv / .md  (%d rows)"%(name,len(df)))

SUB={"rapa":["A"],"nigra":["B"],"oleracea":["C"],
     "napus":["A","C"],"juncea":["A","B"],"carinata":["B","C"]}
PLOIDY={"rapa":"diploid (2n=20)","nigra":"diploid (2n=16)","oleracea":"diploid (2n=18)",
        "napus":"allotetraploid (2n=38)","juncea":"allotetraploid (2n=36)",
        "carinata":"allotetraploid (2n=34)"}
ASM={"rapa":"Brapa v3.0","nigra":"B. nigra NI100 v2","oleracea":"B. oleracea v2.1",
     "napus":"B. napus 3DH v3.1","juncea":"B. juncea 3DH v1","carinata":"B. carinata 3DH v1"}
TIS={"rapa":"pollen, stigma","nigra":"pollen, stigma","oleracea":"pollen, stigma",
     "napus":"pollen, stigma","juncea":"pollen, stigma",
     "carinata":"pollen, stigma (early), stigma (late)"}
cnt=pd.read_csv(os.path.join(REB,"01_counts/gref_counts.tsv"),sep="\t",nrows=1)
libs=[c.split("|")[1] for c in cnt.columns if "|" in c]
PRE={"rapa":"Bra","nigra":"Bni","oleracea":"Bol","napus":"Bna","juncea":"Bju","carinata":"Bca"}
rows=[]
for sp in ["rapa","nigra","oleracea","napus","juncea","carinata"]:
    n=len({l for l in libs if l.startswith(PRE[sp]+"_")})
    rows.append({"Species":"B. "+sp,"Ploidy":PLOIDY[sp],"Subgenomes":" + ".join(SUB[sp]),
                 "Reference assembly":ASM[sp],"Tissues":TIS[sp],"Libraries":n})
t1=pd.DataFrame(rows)
write(t1,"Table1_design",
 "Three replicates per tissue. Libraries total %d. All species are scored on the conserved\n"
 "homoeologous gene set (Gref), so cross-species comparisons use the same orthologous genes."%sum(t1.Libraries))

SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
ORDER=["juncea_pollen","juncea_stigma","napus_pollen","napus_stigma",
       "carinata_pollen","carinata_stigmaE","carinata_stigmaL"]
THRESH=[("theta0585_a05","1.5-fold"),("theta1_a05","2-fold (primary)"),
        ("theta1585_a05","3-fold"),("span_theta1_a05","2-fold, gene span"),
        ("Q10_theta1_a05","2-fold, MAPQ>=10")]
main=pd.read_csv(os.path.join(REB,"02_heb/theta1_a05/heb_summary.csv")).set_index("sample")
alt={t:pd.read_csv(os.path.join(REB,"02_heb",t,"heb_summary.csv")).set_index("sample")
     for t,_ in THRESH if os.path.exists(os.path.join(REB,"02_heb",t,"heb_summary.csv"))}
rows=[]
for s in ORDER:
    if s not in main.index: continue
    r=main.loc[s]; sp=s.split("_")[0]; g1,g2=SUBS[sp]
    k1,k2=int(r.toward1),int(r.toward2); n=k1+k2
    ratio=k1/k2; lo,hi=wilson(k1,n); p=binomtest(k1,n,0.5).pvalue
    span=[alt[t].loc[s].toward1/alt[t].loc[s].toward2 for t,_ in THRESH if t in alt and s in alt[t].index]
    rows.append({"Sample":s.replace("_"," "),"Contest":"%s vs %s"%(g1,g2),
                 "Pairs tested":int(r.tested),"Toward sub 1":k1,"Toward sub 2":k2,
                 "Ratio":"%.3f"%ratio,"95% CI":"%.3f-%.3f"%(lo,hi),
                 "p":("%.1g"%p if p>=1e-4 else "%.0e"%p),
                 "Resolved":"yes" if not (lo<1<hi) else "no",
                 "Range over the settings":"%.3f-%.3f"%(min(span),max(span)),
                 "Favoured":(g1 if ratio>1 else g2) if not (lo<1<hi) else "-"})
t2=pd.DataFrame(rows)
write(t2,"Table2_scoreboard",
 "Ratio is biased pairs favouring subgenome 1 over subgenome 2 at 2-fold and padj <= 0.05.\n"
 "95% CI is the Wilson interval on the proportion, re-expressed as a ratio; p is an exact\n"
 "binomial test against 1:1. A contest is resolved when the interval excludes 1.")

FRAC={"An":0.382,"Aj":0.328,"Bj":0.375,"Bc":0.355,"Cn":0.434,"Cc":0.448}
DCJ ={"An":36,"Aj":49,"Bj":43,"Bc":26,"Cn":42,"Cc":38}
KS  ={"An":0.0376,"Aj":0.0355,"Bj":0.0398,"Bc":0.0405,"Cn":0.0222,"Cc":0.0220}
KAKS={"An":0.1825,"Aj":0.1843,"Bj":0.1967,"Bc":0.1921,"Cn":0.2783,"Cc":0.2893}
LTRS={"An":0.005541,"Aj":0.005320,"Bj":0.004500,"Bc":0.004164,"Cn":0.008167,"Cc":0.008923}
HOST={"An":"napus","Cn":"napus","Aj":"juncea","Bj":"juncea","Bc":"carinata","Cc":"carinata"}
PROG={"A":"rapa","B":"nigra","C":"oleracea"}
rows=[]
for s in ["An","Aj","Bj","Bc","Cn","Cc"]:
    rows.append({"Allo-subgenome":s,"Host":"B. "+HOST[s],"Progenitor":"B. "+PROG[s[0]],
                 "Fractionation ratio":FRAC[s],"DCJ distance":DCJ[s],
                 "Median Ks":"%.4f"%KS[s],"Ka/Ks (median)":"%.4f"%KAKS[s],
                 "Gene-proximal LTR":"%.5f"%LTRS[s]})
t3=pd.DataFrame(rows)
write(t3,"Table3_genome_features",
 "Fractionation ratio and DCJ distance are measured on the synteny blocks between each\n"
 "Allo-subgenome and its diploid progenitor. Ks is the median of per-pair Ks over the\n"
 "allotetraploidization blocks with Ks <= 1. Gene-proximal LTR is the fraction of the 5 kb\n"
 "either side of a gene covered by an annotated LTR, on one common TE library.")

inh=pd.read_csv(os.path.join(REB,"05_validation/inherited_vs_novel.csv"))
inh=inh.sort_values(["species","tissue"])
t4=pd.DataFrame({"Species":"B. "+inh.species,"Tissue":inh.tissue,"Pairs":inh.n,
                 "Spearman rho":inh["rho"].round(3),"R2":inh.R2.round(3),
                 "Slope":inh.slope.round(3),
                 "Inherited %":inh.pct_inherited.round(1),"Novel %":inh.pct_novel.round(1),
                 "Lost %":inh.pct_lost.round(1),"Reversed %":inh.pct_reversed.round(1)})
write(t4,"Table4_inheritance",
 "Regression of homoeolog bias inside the allotetraploid on the bias between the same two\n"
 "genes in the diploid progenitors. Slope below 1 means parental bias is dampened in the\n"
 "hybrid. Categories use a 2-fold cut on both axes. par_lfc carries annotation noise, which\n"
 "biases the slope downward, so slopes are a lower bound on what is retained.")
