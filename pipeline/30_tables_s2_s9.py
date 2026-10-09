#!/usr/bin/env python3
import os
import pandas as pd

ABSENT_CELLS = {"", "x", "-", "na", "#n/a", "0", "#ref!", "#value!", "null"}
def present_mask(s):
    return ~s.astype(str).str.strip().str.lower().isin(ABSENT_CELLS)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(HERE, "08_structure")

GREF_ANCHORS = 27203
LAYERS       = 3
POSSIBLE     = GREF_ANCHORS * LAYERS

SPECIES = {"Bra": "rapa", "Bni": "nigra", "Bol": "oleracea",
           "Bna": "napus", "Bju": "juncea", "Bca": "carinata"}
ORDER = ["rapa", "nigra", "oleracea", "napus", "juncea", "carinata"]

df = pd.read_csv(os.path.join(HERE, "01_counts", "gref_counts.tsv"), sep="\t")
libs = [c for c in df.columns if "|" in c]

def tissue_of(lib):
    s = lib.split("|")[1]
    if "Pol" in s: return "pollen"
    if "E_Stig" in s: return "stigma (early)"
    if "L_Stig" in s: return "stigma (late)"
    return "stigma"

rows = []
tracks = sorted({c.split("|")[0] for c in libs})
for tr in tracks:
    sp, gen = tr.split("_")
    cols = [c for c in libs if c.startswith(tr + "|")]
    for tis in sorted({tissue_of(c) for c in cols}):
        tc = [c for c in cols if tissue_of(c) == tis]
        present = df[tc].notna().all(axis=1)
        tot = df.loc[present, tc].sum(axis=1)
        rows.append(dict(species=SPECIES[sp], subgenome=gen, tissue=tis,
                         replicates=len(tc),
                         possible=POSSIBLE,
                         copies_present=int(present.sum()),
                         expressed=int((tot >= 3).sum()),
                         pct_of_possible=100.0 * (tot >= 3).sum() / POSSIBLE,
                         pct_of_present=100.0 * (tot >= 3).sum() / present.sum()))
S2 = pd.DataFrame(rows)
S2["_s"] = pd.Categorical(S2.species, ORDER)
S2 = S2.sort_values(["_s", "subgenome", "tissue"]).drop(columns="_s")
S2.to_csv(os.path.join(OUT, "table_s2_expressed.csv"), index=False)
print("=== Table S2: expressed gene copies (>= 3 reads over the three replicates) ===")
print(S2.to_string(index=False, float_format=lambda v: "%.1f" % v))
print("\nstigma expresses more copies than pollen in every species:")
for sp in ORDER:
    g = S2[S2.species == sp]
    pol = g[g.tissue == "pollen"].expressed.sum()
    sti = g[g.tissue != "pollen"].expressed.sum() / max(g[g.tissue != "pollen"].tissue.nunique(), 1)
    print("  %-9s pollen %6d   stigma %8.0f   %s" % (sp, pol, sti, "yes" if sti > pol else "NO"))

SUB = pd.read_csv(os.path.join(HERE, "00_inputs", "Subgenomes_Brassica.txt"),
                  sep="\t", dtype=str).fillna("x")
TRACK = {
  ("rapa","A"):     ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"],
  ("oleracea","C"): ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"],
  ("nigra","B"):    ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"],
  ("napus","A"):    ["Bna_A_subgenome1","Bna_A_subgenome2","Bna_A_subgenome3"],
  ("napus","C"):    ["Bna_C_subgenome1","Bna_C_subgenome2","Bna_C_subgenome3"],
  ("juncea","A"):   ["Bju_A_subgenome_1","Bju_A_subgenome_2","Bju_A_subgenome_3"],
  ("juncea","B"):   ["Bju_B_subgenome_1","Bju_B_subgenome_2","Bju_B_subgenome_3"],
  ("carinata","B"): ["Bca_B_subgenome_1","Bca_B_subgenome_2","Bca_B_subgenome_3"],
  ("carinata","C"): ["Bca_C_subgenome_1","Bca_C_subgenome_2","Bca_C_subgenome_3"],
}
N = len(SUB)
rows = []
for (sp, gen), cols in TRACK.items():
    assigned = {}
    for lay, c in zip(["LF", "MF1", "MF2"], cols):
        assigned[lay] = int(present_mask(SUB[c]).sum())
    anyc = SUB[cols].apply(present_mask).any(axis=1).sum()
    rows.append(dict(species=sp, subgenome=gen, anchors=N,
                     LF=assigned["LF"], MF1=assigned["MF1"], MF2=assigned["MF2"],
                     total_assigned=int(anyc),
                     pct_assigned=100.0 * anyc / N))
S9 = pd.DataFrame(rows)
S9["_s"] = pd.Categorical(S9.species, ORDER)
S9 = S9.sort_values(["_s", "subgenome"]).drop(columns="_s")
S9.to_csv(os.path.join(OUT, "table_s9_synteny.csv"), index=False)
print("\n=== Table S9: Br-subgenome assignment, over %d Arabidopsis anchors ===" % N)
print(S9.to_string(index=False, float_format=lambda v: "%.1f" % v))
print("\nLF > MF1 > MF2 in every track: %s"
      % bool((S9.LF > S9.MF1).all() and (S9.MF1 > S9.MF2).all()))
print("\nwrote table_s2_expressed.csv and table_s9_synteny.csv")

HEB = pd.read_csv(os.path.join(HERE, "02_heb", "theta1_a05", "heb_summary.csv"))
PAIR = {"napus": ("Bna_A", "Bna_C"), "juncea": ("Bju_A", "Bju_B"),
        "carinata": ("Bca_B", "Bca_C")}
TIS = {"pollen": "Pol", "stigma": "Stig", "stigmaE": "E_Stig", "stigmaL": "L_Stig"}
slots_with_counts = len(df)
rows = []
for _, r in HEB.iterrows():
    sp, tis = r["sample"].split("_", 1)
    t1, t2 = PAIR[sp]
    tag = TIS[tis]
    c1 = [c for c in libs if c.startswith(t1 + "|") and tag in c.split("|")[1]]
    c2 = [c for c in libs if c.startswith(t2 + "|") and tag in c.split("|")[1]]
    if tis == "stigma":
        c1 = [c for c in c1 if "E_Stig" not in c and "L_Stig" not in c]
        c2 = [c for c in c2 if "E_Stig" not in c and "L_Stig" not in c]
    both = (df[c1].notna().all(axis=1) & df[c2].notna().all(axis=1)).sum()
    biased = int(r.toward1 + r.toward2)
    rows.append(dict(sample=r["sample"].replace("_", " "), possible_pairs=POSSIBLE,
                     both_copies_present=int(both), tested=int(r.tested),
                     biased=biased, balanced=int(r.balanced),
                     indeterminate=int(r.indeterminate),
                     pct_biased_of_tested=100.0 * biased / r.tested,
                     pct_tested_of_possible=100.0 * r.tested / POSSIBLE))
S10 = pd.DataFrame(rows)
S10.to_csv(os.path.join(OUT, "table_s10_flow.csv"), index=False)
print("\n=== Table S10: gene-set flow ===")
print(S10.to_string(index=False, float_format=lambda v: "%.1f" % v))
print("\ncheck, biased + balanced + indeterminate == tested:",
      bool(((S10.biased + S10.balanced + S10.indeterminate) == S10.tested).all()))
print("wrote table_s10_flow.csv")
