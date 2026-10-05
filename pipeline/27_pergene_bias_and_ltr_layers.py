#!/usr/bin/env python3
import os, itertools
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon

HEB = "02_heb/theta1_a05"
VAL = "05_validation"
OUT = "08_structure"
os.makedirs(OUT, exist_ok=True)

SAMPLES = ["napus_pollen","napus_stigma","carinata_pollen","carinata_stigmaE",
           "carinata_stigmaL","juncea_pollen","juncea_stigma"]
SUBS = {"napus":("A","C"), "juncea":("A","B"), "carinata":("B","C")}
LAYNAME = {1:"LF", 2:"MF1", 3:"MF2"}

lfc = {}
for s in SAMPLES:
    d = pd.read_csv(f"{HEB}/pairs_{s}.csv")
    lfc[s] = d.set_index("key").allo_lfc
M = pd.DataFrame(lfc)
print("pairs keyed across samples: %d rows, %d samples" % M.shape)

rows = []
for a, b in itertools.combinations(SAMPLES, 2):
    sub = M[[a, b]].dropna()
    rho, p = spearmanr(sub[a], sub[b])
    rows.append(dict(a=a, b=b, n=len(sub), rho=rho, p=p,
                     same_species=a.split("_")[0] == b.split("_")[0]))
C = pd.DataFrame(rows)
C.to_csv(f"{OUT}/pergene_bias_correlation.csv", index=False)

print("\n=== per-gene bias correlation, within vs between species ===")
print("  within  species: rho %.3f to %.3f (n=%d pairs of samples)"
      % (C[C.same_species].rho.min(), C[C.same_species].rho.max(), C.same_species.sum()))
print("  between species: rho %.3f to %.3f (n=%d pairs of samples)"
      % (C[~C.same_species].rho.min(), C[~C.same_species].rho.max(), (~C.same_species).sum()))

print("\n=== the two C-carrying allotetraploids, napus against carinata ===")
for a in ["napus_pollen", "napus_stigma"]:
    for b in ["carinata_pollen", "carinata_stigmaE", "carinata_stigmaL"]:
        r = C[(C.a == a) & (C.b == b)].iloc[0]
        print("  %-16s vs %-18s n=%5d  rho=%+.3f  p=%.3g" % (a, b, r.n, r.rho, r.p))

Mat = pd.DataFrame(np.eye(len(SAMPLES)), index=SAMPLES, columns=SAMPLES)
for _, r in C.iterrows():
    Mat.loc[r.a, r.b] = Mat.loc[r.b, r.a] = r.rho
Mat.to_csv(f"{OUT}/pergene_bias_matrix.csv")

print("\n=== gene-proximal LTR density per Br-subgenome layer ===")
out = []
for sp in ["napus", "juncea", "carinata"]:
    d = pd.read_csv(f"{VAL}/ltr_per_gene_{sp}.csv")
    d["layer"] = d.key.str.split("|").str[1].astype(int).map(LAYNAME)
    s1, s2 = SUBS[sp]
    for lay in ["LF", "MF1", "MF2"]:
        g = d[d.layer == lay].dropna(subset=["ltr1", "ltr2"])
        if len(g) < 30:
            continue
        try:
            stat, p = wilcoxon(g.ltr1, g.ltr2)
        except ValueError:
            p = np.nan
        out.append(dict(species=sp, layer=lay, n=len(g), sub1=s1, sub2=s2,
                        mean_ltr1=g.ltr1.mean(), mean_ltr2=g.ltr2.mean(),
                        median_diff=float(g.ltr_diff.median()), p=p))
        print("  %-9s %-4s n=%5d   %s %.5f vs %s %.5f   median diff %+.5f  p=%.3g"
              % (sp, lay, len(g), s1, g.ltr1.mean(), s2, g.ltr2.mean(),
                 g.ltr_diff.median(), p))
L = pd.DataFrame(out)
L.to_csv(f"{OUT}/ltr_by_layer.csv", index=False)

print("\n=== does LTR density follow the LF > MF1 > MF2 retention hierarchy? ===")
for sp in ["napus", "juncea", "carinata"]:
    g = L[L.species == sp].set_index("layer")
    if len(g) < 3:
        continue
    m = (g.mean_ltr1 + g.mean_ltr2) / 2
    print("  %-9s LF %.5f  MF1 %.5f  MF2 %.5f   %s"
          % (sp, m.LF, m.MF1, m.MF2,
             "monotone increasing" if m.LF < m.MF1 < m.MF2 else
             "monotone decreasing" if m.LF > m.MF1 > m.MF2 else "not monotone"))
print("\nwrote %s/pergene_bias_correlation.csv, pergene_bias_matrix.csv, ltr_by_layer.csv" % OUT)
