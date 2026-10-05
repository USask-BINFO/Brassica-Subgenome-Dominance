#!/usr/bin/env python3
import pandas as pd, numpy as np
from scipy.stats import fisher_exact

IN  = "00_inputs/Subgenomes_Brassica.txt"
OUT = "08_structure"

TRACK = {
  "Bra_A": ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"],
  "Bol_C": ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"],
  "Bni_B": ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"],
  "Bna_A": ["Bna_A_subgenome1","Bna_A_subgenome2","Bna_A_subgenome3"],
  "Bna_C": ["Bna_C_subgenome1","Bna_C_subgenome2","Bna_C_subgenome3"],
  "Bju_A": ["Bju_A_subgenome_1","Bju_A_subgenome_2","Bju_A_subgenome_3"],
  "Bju_B": ["Bju_B_subgenome_1","Bju_B_subgenome_2","Bju_B_subgenome_3"],
  "Bca_B": ["Bca_B_subgenome_1","Bca_B_subgenome_2","Bca_B_subgenome_3"],
  "Bca_C": ["Bca_C_subgenome_1","Bca_C_subgenome_2","Bca_C_subgenome_3"],
}
LAYER    = ["LF","MF1","MF2"]
PROGEN   = {"Bna_A":"Bra_A","Bju_A":"Bra_A","Bna_C":"Bol_C","Bca_C":"Bol_C",
            "Bju_B":"Bni_B","Bca_B":"Bni_B"}
PRETTY   = {"Bra_A":"rapa A","Bol_C":"oleracea C","Bni_B":"nigra B","Bna_A":"napus A",
            "Bna_C":"napus C","Bju_A":"juncea A","Bju_B":"juncea B",
            "Bca_B":"carinata B","Bca_C":"carinata C"}

df = pd.read_csv(IN, sep="\t", dtype=str).fillna("x")
N  = len(df)
print("Arabidopsis genes in the table: %d" % N)

rows = []
for tr, cols in TRACK.items():
    for lay, c in zip(LAYER, cols):
        present = int((df[c].str.strip().str.lower() != "x").sum())
        rows.append(dict(track=tr, label=PRETTY[tr], layer=lay,
                         retained=present, total=N, pct=100.0*present/N))
ret = pd.DataFrame(rows)
ret.to_csv(f"{OUT}/retention_by_layer.csv", index=False)

print("\n=== retention per layer (%% of %d Arabidopsis genes) ===" % N)
piv = ret.pivot(index="label", columns="layer", values="pct")[LAYER]
cnt = ret.pivot(index="label", columns="layer", values="retained")[LAYER]
order = [PRETTY[t] for t in TRACK]
print(piv.loc[order].round(1).to_string())

tot = {}
for tr, cols in TRACK.items():
    any_copy = (df[cols].apply(lambda s: s.str.strip().str.lower() != "x")).any(axis=1)
    tot[tr] = int(any_copy.sum())

print("\n=== total retention, and loss relative to the diploid progenitor ===")
out = []
for tr in TRACK:
    r = dict(track=tr, label=PRETTY[tr], retained_any=tot[tr], total=N,
             pct_any=100.0*tot[tr]/N)
    if tr in PROGEN:
        p = PROGEN[tr]
        odds, pv = fisher_exact([[tot[tr], N-tot[tr]], [tot[p], N-tot[p]]])
        r.update(progenitor=PRETTY[p], pct_progenitor=100.0*tot[p]/N,
                 delta_pct=100.0*(tot[tr]-tot[p])/N, p=pv)
    out.append(r)
tt = pd.DataFrame(out)
tt.to_csv(f"{OUT}/retention_total.csv", index=False)
for _, r in tt.iterrows():
    if "progenitor" in tt.columns and isinstance(r.get("progenitor"), str):
        print("  %-12s %5.1f%%   vs %-12s %5.1f%%   delta %+5.2f pp   p=%.3g"
              % (r.label, r.pct_any, r.progenitor, r.pct_progenitor, r.delta_pct, r.p))
    else:
        print("  %-12s %5.1f%%   (diploid)" % (r.label, r.pct_any))

print("\n=== is the LF > MF1 > MF2 ordering intact in every track? ===")
ok = True
for tr in TRACK:
    v = ret[ret.track == tr].set_index("layer").loc[LAYER, "retained"]
    mono = v.LF > v.MF1 > v.MF2
    ok &= bool(mono)
    _, p_lf_mf1 = fisher_exact([[v.LF, N-v.LF], [v.MF1, N-v.MF1]])
    print("  %-12s LF %5d  MF1 %5d  MF2 %5d   %s  (LF vs MF1 p=%.3g)"
          % (PRETTY[tr], v.LF, v.MF1, v.MF2, "monotone" if mono else "NOT monotone", p_lf_mf1))
print("\nordering holds in every track: %s" % ok)
print("\nwrote %s/retention_by_layer.csv and %s/retention_total.csv" % (OUT, OUT))
