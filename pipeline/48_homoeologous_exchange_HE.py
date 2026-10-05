#!/usr/bin/env python3
import os, sys, collections, csv
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAIRS = {
 "napus":    ("Bna_A", ["Bna_A_subgenome1","Bna_A_subgenome2","Bna_A_subgenome3"],
              "Bna_C", ["Bna_C_subgenome1","Bna_C_subgenome2","Bna_C_subgenome3"], "Bna"),
 "juncea":   ("Bju_A", ["Bju_A_subgenome_1","Bju_A_subgenome_2","Bju_A_subgenome_3"],
              "Bju_B", ["Bju_B_subgenome_1","Bju_B_subgenome_2","Bju_B_subgenome_3"], "Bju"),
 "carinata": ("Bca_B", ["Bca_B_subgenome_1","Bca_B_subgenome_2","Bca_B_subgenome_3"],
              "Bca_C", ["Bca_C_subgenome_1","Bca_C_subgenome_2","Bca_C_subgenome_3"], "Bca"),
}
MINRUN  = 5
MAXGAP  = 2
OFFMAX  = 1
ONMIN   = 10

SUB = pd.read_csv(os.path.join(HERE,"00_inputs","Subgenomes_Brassica.txt"), sep="\t",
                  dtype=str).fillna("x")
CNT = pd.read_csv(os.path.join(HERE,"01_counts","gref_counts.tsv"), sep="\t")

def coords(code):
    d = pd.read_csv(os.path.join(HERE,"00_inputs",code+".saf"), sep="\t")
    g = d.groupby("GeneID").agg(chrom=("Chr","first"), start=("Start","min")).reset_index()
    return dict(zip(g.GeneID, zip(g.chrom, g.start)))

rows = []
print("homoeologous exchange between the two subgenomes of one nucleus\n")
print("%-9s %-7s %7s %8s %9s %9s %9s" % ("species","dir","pairs","silent","runs","genes","median x"))
print("-"*62)
for sp,(t1,c1,t2,c2,code) in PAIRS.items():
    CO = coords(code)
    libs1 = [c for c in CNT.columns if c.startswith(t1+"|")]
    libs2 = [c for c in CNT.columns if c.startswith(t2+"|")]
    key = CNT.set_index(["AT_gene","br_sub"])
    M1 = key[libs1]; M2 = key[libs2]
    for direction,(ta,ca,tb,cb,Ma,Mb) in (("%s lost"%t1,(t1,c1,t2,c2,M1,M2)),
                                          ("%s lost"%t2,(t2,c2,t1,c1,M2,M1))):
        recs=[]
        for layer,(colA,colB) in enumerate(zip(ca,cb), start=1):
            for at,ga,gb in zip(SUB.AT_geneid, SUB[colA], SUB[colB]):
                if ga=="x" or gb=="x" or ga not in CO: continue
                k=(at,layer)
                if k not in Ma.index or k not in Mb.index: continue
                a=Ma.loc[k].values; b=Mb.loc[k].values
                if np.ndim(a)>1: continue
                recs.append((CO[ga][0], CO[ga][1], at, layer, ga, gb, a.max(), b.sum(), b))
        recs.sort(key=lambda r:(r[0],r[1]))
        silent = sum(1 for r in recs if r[6]<=OFFMAX and r[7]>=ONMIN)
        runs=[]; cur=[]; gap=0
        last_chr=None
        for r in recs:
            if r[0]!=last_chr: 
                if len(cur)>=MINRUN: runs.append(cur)
                cur=[]; gap=0; last_chr=r[0]
            if r[6]<=OFFMAX and r[7]>=ONMIN:
                cur.append(r); gap=0
            else:
                gap+=1
                if gap>MAXGAP:
                    if len(cur)>=MINRUN: runs.append(cur)
                    cur=[]; gap=0
        if len(cur)>=MINRUN: runs.append(cur)
        base = {}
        for r in recs:
            base.setdefault(r[0],[]).append(r[7])
        ratios=[]
        for run in runs:
            med_out = np.median([v for v in base[run[0][0]] if v>0]) or 1
            ratios.append(np.median([r[7] for r in run]) / med_out)
            rows.append(dict(species=sp, direction=direction, chrom=run[0][0],
                             pairs=len(run), first_gene=run[0][4], last_gene=run[-1][4],
                             partner_fold=round(float(np.median([r[7] for r in run])/med_out),2)))
        print("%-9s %-7s %7d %8d %9d %9d %9s" % (
            sp, direction.split()[0], len(recs), silent, len(runs),
            sum(len(r) for r in runs),
            "%.2f"%np.median(ratios) if ratios else "-"))

os.makedirs("08_structure", exist_ok=True)
with open("08_structure/he_candidates.csv","w",newline="") as fh:
    if rows:
        w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote 08_structure/he_candidates.csv (%d runs)"%len(rows))
