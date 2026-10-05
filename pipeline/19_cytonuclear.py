#!/usr/bin/env python3
import pickle, glob, os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

GO=pickle.load(open("00_inputs/go_propagated.pkl","rb")); G2GO=GO["g2go"]; NAME=GO["name"]
SETS={
 "chloroplast"      : "GO:0009507",
 "plastid"          : "GO:0009536",
 "mitochondrion"    : "GO:0005739",
 "photosynthesis"   : "GO:0015979",
 "thylakoid"        : "GO:0009579",
}
MAT={"juncea":("A","subgenome1 = Bju_A"),"carinata":("B","subgenome1 = Bca_B"),
     "napus":(None,"unclear - no prediction")}
TIS={"juncea":["pollen","stigma"],"carinata":["pollen","stigmaE","stigmaL"],
     "napus":["pollen","stigma"]}
rows=[]
print("Positive shift = bias toward subgenome 1. For juncea and carinata that is the MATERNAL copy.\n")
for sp in ("juncea","carinata","napus"):
    mat,desc=MAT[sp]
    print("=== %-9s maternal: %s" % (sp,desc))
    for t in TIS[sp]:
        f="02_heb/theta2_a05/pairs_%s_%s.csv"%(sp,t)
        if not os.path.exists(f): continue
        d=pd.read_csv(f)
        d["at"]=d["key"].str.split("|").str[0]
        d=d[np.isfinite(d["allo_lfc"])]
        for label,term in SETS.items():
            inset=d["at"].map(lambda g: term in G2GO.get(g,()))
            a=d.loc[inset,"allo_lfc"].values; b=d.loc[~inset,"allo_lfc"].values
            if len(a)<20: continue
            u,p=mannwhitneyu(a,b,alternative="two-sided")
            sgn = "toward sub1" if np.median(a)>np.median(b) else "toward sub2"
            star = "*" if p<0.05 else " "
            print("    %-9s %-15s n=%-5d median %+.3f vs %+.3f  d=%+.3f  p=%-9.3g %s %s"
                  % (t,label,len(a),np.median(a),np.median(b),
                     np.median(a)-np.median(b),p,star,sgn))
            rows.append(dict(species=sp,tissue=t,geneset=label,n=len(a),
                             median_set=float(np.median(a)),median_bg=float(np.median(b)),
                             delta=float(np.median(a)-np.median(b)),p=float(p),
                             maternal_is_sub1=(mat is not None)))
    print()
pd.DataFrame(rows).to_csv("05_validation/cytonuclear.csv",index=False)
print("Wrote 05_validation/cytonuclear.csv")
print("\nRead the juncea vs carinata contrast: a maternal effect requires a POSITIVE delta in")
print("BOTH, since maternal is subgenome 1 in both. A positive delta only in juncea means the")
print("signal follows the A genome, matching the paper's enrichment-based conclusion.")
