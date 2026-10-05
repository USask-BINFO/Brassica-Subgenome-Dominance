#!/usr/bin/env python3
import sys, bisect, collections
import numpy as np, pandas as pd
LTR="00_inputs/ltr/Bna_genome_v3.1.fa.19K_repeat.LTR.gff"
BED="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new/Research/BWA/Bna"
FLANK=5000

iv=collections.defaultdict(list)
for line in open(LTR):
    if line.startswith("#"): continue
    f=line.split("\t")
    if len(f)<5: continue
    iv[f[0]].append((int(f[3]),int(f[4])))
idx={}
for c,l in iv.items():
    l.sort(); idx[c]=([x[0] for x in l], l)
sys.stderr.write("LTR features: %d on %d sequences\n" % (sum(len(v) for v in iv.values()), len(iv)))

def ltr_bp(c, s, e):
    if c not in idx: return 0
    starts, l = idx[c]
    i = bisect.bisect_left(starts, s) - 1
    tot = 0
    while i < len(l):
        if i >= 0:
            a,b = l[i]
            if a > e: break
            tot += max(0, min(b,e) - max(a,s) + 1)
        i += 1
    return tot

gene={}
for s in (1,2,3):
    for g in ("A","C"):
        for line in open("%s/Bna_sub%d%s.bed" % (BED,s,g)):
            f=line.rstrip("\n").split("\t")
            if len(f)<6: continue
            gene[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])
sys.stderr.write("genes with coordinates: %d\n" % len(gene))

def flank_density(gid):
    if gid not in gene: return np.nan
    c,s,e,st = gene[gid]
    up = ltr_bp(c, max(0,s-FLANK), s-1)
    dn = ltr_bp(c, e+1, e+FLANK)
    return (up+dn) / (2.0*FLANK)

import csv, re
SYN = "00_inputs/Subgenomes_Brassica.txt"
rows=list(csv.reader(open(SYN),delimiter="\t")); h={x.strip():i for i,x in enumerate(rows[0])}
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")
rec={}
for r in rows[1:]:
    for s in (1,2,3):
        a=r[h["Bna_A_subgenome%d"%s]].strip(); c=r[h["Bna_C_subgenome%d"%s]].strip()
        if a in ABSENT or c in ABSENT: continue
        rec["%s|%d"%(r[h["AT_geneid"]],s)] = (STRIP.sub("",a), STRIP.sub("",c))
sys.stderr.write("napus pairs with both ids: %d\n" % len(rec))

cache={}
def dens(g):
    if g not in cache: cache[g]=flank_density(g)
    return cache[g]

out=[]
for k,(ga,gc) in rec.items():
    da,dc = dens(ga), dens(gc)
    if np.isnan(da) or np.isnan(dc): continue
    out.append((k, da, dc, da-dc))
df=pd.DataFrame(out, columns=["key","ltr_A","ltr_C","ltr_diff"])
df.to_csv("05_validation/napus_ltr_per_gene.csv", index=False)
sys.stderr.write("pairs with LTR density for both copies: %d\n" % len(df))
print(df.head().to_string())
