#!/usr/bin/env python3
import pandas as pd, numpy as np, sys
g = pd.read_csv("01_counts/gref_counts.tsv", sep="\t")
cols = lambda p, t: [c for c in g.columns if c.startswith(p + "|") and t in c]
S = [("napus","Bna_A","Bna_C","Pol"),("napus","Bna_A","Bna_C","Stig"),
     ("carinata","Bca_B","Bca_C","Pol"),("carinata","Bca_B","Bca_C","E_Stig"),
     ("carinata","Bca_B","Bca_C","L_Stig"),
     ("juncea","Bju_A","Bju_B","Pol"),("juncea","Bju_A","Bju_B","Stig")]
rows=[]
print("%-10s %-7s %-6s %9s %9s %9s %11s %12s" % (
      "species","tissue","pair","summed","median","genes 1>2","genes 2>1","n to erase"))
for sp,g1,g2,tis in S:
    A=g[cols(g1,tis)]; C=g[cols(g2,tis)]
    ok=A.notna().all(axis=1)&C.notna().all(axis=1)
    a=A[ok].sum(axis=1).values; c=C[ok].sum(axis=1).values
    m=(a+c)>=20
    lfc=np.log2((a[m]+1)/(c[m]+1))
    d=a-c; tot=d.sum()
    pos=np.sort(d[d>0])[::-1]; cum=np.cumsum(pos)
    k=int(np.searchsorted(cum, abs(tot)))+1 if tot>0 else np.nan
    if tot<0:
        neg=np.sort(d[d<0]); cum=np.cumsum(neg); k=int(np.searchsorted(-cum, -tot))+1
    lab=g1.split('_')[1]+":"+g2.split('_')[1]
    print("%-10s %-7s %-6s %9.3f %+9.3f %9d %11d %7d (%.2f%%)" % (
          sp,tis,lab,a.sum()/c.sum(),np.median(lfc),(lfc>0).sum(),(lfc<0).sum(),k,100*k/len(d)))
    rows.append(dict(species=sp,tissue=tis,pair=lab,summed_ratio=round(a.sum()/c.sum(),4),
        median_log2=round(float(np.median(lfc)),4), genes_first=int((lfc>0).sum()),
        genes_second=int((lfc<0).sum()), pairs_to_erase=int(k),
        pct_to_erase=round(100*k/len(d),3), n_pairs=int(len(d))))
pd.DataFrame(rows).to_csv("05_validation/summed_expression.csv", index=False)
print("\nwritten to 05_validation/summed_expression.csv")
