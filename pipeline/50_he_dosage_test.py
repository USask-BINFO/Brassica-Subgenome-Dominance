#!/usr/bin/env python3
import os, csv
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = {"Bna_A":("Bra_A",["Bna_A_subgenome1","Bna_A_subgenome2","Bna_A_subgenome3"],
                        ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"]),
       "Bna_C":("Bol_C",["Bna_C_subgenome1","Bna_C_subgenome2","Bna_C_subgenome3"],
                        ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"]),
       "Bju_A":("Bra_A",["Bju_A_subgenome_1","Bju_A_subgenome_2","Bju_A_subgenome_3"],
                        ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"]),
       "Bju_B":("Bni_B",["Bju_B_subgenome_1","Bju_B_subgenome_2","Bju_B_subgenome_3"],
                        ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"]),
       "Bca_B":("Bni_B",["Bca_B_subgenome_1","Bca_B_subgenome_2","Bca_B_subgenome_3"],
                        ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"]),
       "Bca_C":("Bol_C",["Bca_C_subgenome_1","Bca_C_subgenome_2","Bca_C_subgenome_3"],
                        ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"])}

CNT = pd.read_csv(os.path.join(HERE,"01_counts","gref_counts.tsv"), sep="\t")
SUB = pd.read_csv(os.path.join(HERE,"00_inputs","Subgenomes_Brassica.txt"), sep="\t",
                  dtype=str).fillna("x")
RUNS = pd.read_csv(os.path.join(HERE,"08_structure","he_presence_runs.csv"))
key = CNT.set_index(["AT_gene","br_sub"])

def cpm(track):
    libs=[c for c in CNT.columns if c.startswith(track+"|")]
    m=key[libs]
    return (m / m.sum(axis=0) * 1e6).mean(axis=1)

CP = {t: cpm(t) for t in set(list(REF)+[v[0] for v in REF.values()])}
IDX = {}
for t,(ref,cols,rcols) in REF.items():
    d={}
    for layer,c in enumerate(cols, start=1):
        for at,g in zip(SUB.AT_geneid, SUB[c]):
            if g!="x": d[g]=(at,layer)
    IDX[t]=d

out=[]
print("dosage of the retained subgenome where its partner's copies are missing\n")
print("%-9s %-7s %6s %10s %12s %10s" %
      ("species","retained","runs","median dose","runs >= 1.5","runs >= 1.8"))
print("-"*62)
for (sp,keep), d in RUNS.groupby(["species","retained"]):
    ref,cols,rcols = REF[keep]
    a, b = CP[keep], CP[ref]
    common = a.index.intersection(b.index)
    base = np.median([a[k]/b[k] for k in common if b[k] > 1 and a[k] > 0]) or 1.0
    doses=[]
    for _,r in d.iterrows():
        ks=[(t.split(":")[0], int(t.split(":")[1])) for t in str(r.anchor_list).split(";") if ":" in t]
        sel=[k for k in ks if k in a.index and k in b.index]
        vals=[(a[k]/b[k])/base for k in sel if b[k] > 1 and a[k] > 0]
        if vals: doses.append(float(np.median(vals))); out.append(
            dict(species=sp, retained=keep, chrom=r.chrom, anchors=r.anchors,
                 first_gene=r.first_gene, dose=round(float(np.median(vals)),2)))
    if doses:
        doses=np.array(doses)
        print("%-9s %-7s %6d %10.2f %12d %10d" % (sp, keep, len(doses), np.median(doses),
              int((doses>=1.5).sum()), int((doses>=1.8).sum())))

with open(os.path.join(HERE,"08_structure","he_dosage.csv"),"w",newline="") as fh:
    if out:
        w=csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print("\nwrote 08_structure/he_dosage.csv (%d runs scored)"%len(out))
