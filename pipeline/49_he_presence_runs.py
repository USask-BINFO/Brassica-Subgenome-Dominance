#!/usr/bin/env python3
import os, sys, csv
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
MINRUN, MAXGAP = 5, 1

SUB = pd.read_csv(os.path.join(HERE,"00_inputs","Subgenomes_Brassica.txt"), sep="\t",
                  dtype=str).fillna("x")
CNT = pd.read_csv(os.path.join(HERE,"01_counts","gref_counts.tsv"), sep="\t")

def coords(code):
    d = pd.read_csv(os.path.join(HERE,"00_inputs",code+".saf"), sep="\t")
    g = d.groupby("GeneID").agg(chrom=("Chr","first"), start=("Start","min")).reset_index()
    return dict(zip(g.GeneID, zip(g.chrom, g.start)))

rows=[]
print("missing homoeolog runs along the retained subgenome\n")
print("%-9s %-14s %8s %9s %8s %7s %8s %10s" %
      ("species","retained","anchors","missing","missing%","runs","genes","median fold"))
print("-"*78)
for sp,(t1,c1,t2,c2,code) in PAIRS.items():
    CO = coords(code)
    key = CNT.set_index(["AT_gene","br_sub"])
    for keep,(ck,co),other in ((t1,(c1,c2),t2),(t2,(c2,c1),t1)):
        libs=[c for c in CNT.columns if c.startswith(keep+"|")]
        M=key[libs]
        recs=[]
        for layer,(colK,colO) in enumerate(zip(ck,co), start=1):
            for at,gk,go in zip(SUB.AT_geneid, SUB[colK], SUB[colO]):
                if gk=="x" or gk not in CO: continue
                k=(at,layer)
                if k not in M.index: continue
                v=M.loc[k].values
                if np.ndim(v)>1: continue
                recs.append((CO[gk][0], CO[gk][1], at, layer, gk, go=="x", float(np.sum(v))))
        recs.sort(key=lambda r:(r[0],r[1]))
        miss=sum(1 for r in recs if r[5])
        runs=[]; cur=[]; gap=0; last=None
        for r in recs:
            if r[0]!=last:
                if len(cur)>=MINRUN: runs.append(cur)
                cur=[]; gap=0; last=r[0]
            if r[5]: cur.append(r); gap=0
            else:
                gap+=1
                if gap>MAXGAP:
                    if len(cur)>=MINRUN: runs.append(cur)
                    cur=[]; gap=0
        if len(cur)>=MINRUN: runs.append(cur)
        base={}
        for r in recs: base.setdefault(r[0],[]).append(r[6])
        folds=[]
        for run in runs:
            ref=np.median([v for v in base[run[0][0]] if v>0]) or 1.0
            f=float(np.median([r[6] for r in run])/ref)
            folds.append(f)
            rows.append(dict(species=sp, retained=keep, missing_from=other, chrom=run[0][0],
                             anchors=len(run), first_gene=run[0][4], last_gene=run[-1][4],
                             retained_fold=round(f,2),
                             anchor_list=";".join("%s:%d"%(r[2],r[3]) for r in run)))
        print("%-9s %-14s %8d %9d %7.1f%% %7d %8d %10s" %
              (sp, keep, len(recs), miss, 100*miss/max(len(recs),1), len(runs),
               sum(len(r) for r in runs), "%.2f"%np.median(folds) if folds else "-"))

with open("08_structure/he_presence_runs.csv","w",newline="") as fh:
    if rows:
        w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote 08_structure/he_presence_runs.csv (%d runs)"%len(rows))
