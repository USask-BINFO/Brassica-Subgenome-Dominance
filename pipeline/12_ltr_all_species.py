#!/usr/bin/env python3
import sys, os, bisect, collections, csv, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon, binomtest

LTRDIR = "00_inputs/ltr"
BED    = "/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new/Research/BWA"
HEB    = "02_heb/theta2_a05"
FLANK  = 5000
OUT    = "05_validation"

SPECIES = {
 "napus": dict(ltr="Bna_genome_v3.1.fa.19K_repeat.LTR.gff", bed="Bna", pre="Bna_sub",
               g1="A", g2="C", col1="Bna_A_subgenome%d",   col2="Bna_C_subgenome%d",
               tissues=["pollen","stigma"], lib="19K"),
 "juncea": dict(ltr="Bjuncea_genome_v1.fa.out.LTR.gff", bed="Bju", pre="Bju_sub",
               g1="A", g2="B", col1="Bju_A_subgenome_%d", col2="Bju_B_subgenome_%d",
               tissues=["pollen","stigma"], lib="19K"),
 "carinata": dict(ltr="Bcarinata.v1.genome.fasta.out.LTR.gff", bed="Bca", pre="Bca_sub",
               g1="B", g2="C", col1="Bca_B_subgenome_%d",  col2="Bca_C_subgenome_%d",
               tissues=["pollen","stigmaE","stigmaL"], lib="19K"),
}
GFF="00_inputs/coge_gff"
DIPLOID = {
 "rapa":     dict(ltr="Brapa_sequence_v3.0.fasta.19K_repeat.LTR.gff",
                  gff="Brassica_rapa_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68114.gff", lib="19K"),
 "oleracea": dict(ltr="Boleracea.v2.1.genome.fasta.19K_repeat.LTR.gff",
                  gff="Brassica_oleracea_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68113.gff", lib="19K"),
 "nigra":    dict(ltr="Bnigra_NI100.v2.genome.fasta.out.LTR.gff",
                  gff="Brassica_nigra_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68112.gff", lib="OTHER"),
}
STRIP=re.compile(r"\.\d+$"); ABSENT=("","x","-","NA","#N/A","0")

def load_ltr(path):
    iv=collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"): continue
        f=line.split("\t")
        if len(f)<5: continue
        try: iv[f[0]].append((int(f[3]),int(f[4])))
        except ValueError: continue
    idx={}
    for c,l in iv.items():
        l.sort(); idx[c]=([x[0] for x in l], l)
    n=sum(len(v) for v in iv.values())
    return idx, n

def ltr_bp(idx, c, s, e):
    if c not in idx or e<s: return 0
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

def load_genes(bdir, pre, groups):
    gene={}
    for s in (1,2,3):
        for g in groups:
            p="%s/%s/%s%d%s.bed" % (BED,bdir,pre,s,g)
            if not os.path.exists(p):
                sys.stderr.write("  missing %s\n"%p); continue
            for line in open(p):
                f=line.rstrip("\n").split("\t")
                if len(f)<6: continue
                gene[f[3]]=(f[0], int(f[1]), int(f[2]))
    return gene

def density(idx, gene, gid):
    if gid not in gene: return np.nan
    c,s,e = gene[gid]
    return (ltr_bp(idx,c,max(0,s-FLANK),s-1) + ltr_bp(idx,c,e+1,e+FLANK)) / (2.0*FLANK)

rows=list(csv.reader(open("00_inputs/Subgenomes_Brassica.txt"),delimiter="\t"))
H={x.strip():i for i,x in enumerate(rows[0])}

print("="*96)
print("SUBGENOME-LEVEL MEAN FLANKING LTR DENSITY  (fraction of the 5 kb flanks covered)")
print("="*96)
sub_summary={}; PAIRS={}
for sp,cfg in SPECIES.items():
    idx,nf = load_ltr(os.path.join(LTRDIR,cfg["ltr"]))
    gene   = load_genes(cfg["bed"], cfg["pre"], [cfg["g1"],cfg["g2"]])
    rec={}
    for r in rows[1:]:
        for s in (1,2,3):
            a=r[H[cfg["col1"]%s]].strip(); b=r[H[cfg["col2"]%s]].strip()
            if a in ABSENT or b in ABSENT: continue
            rec["%s|%d"%(r[H["AT_geneid"]],s)]=(STRIP.sub("",a), STRIP.sub("",b))
    hit=sum(1 for a,b in rec.values() if a in gene and b in gene)
    rate=hit/max(1,len(rec))
    d1=[];d2=[];keys=[]
    for k,(a,b) in rec.items():
        x,y = density(idx,gene,a), density(idx,gene,b)
        if np.isnan(x) or np.isnan(y): continue
        keys.append(k); d1.append(x); d2.append(y)
    PAIRS[sp]=pd.DataFrame({"key":keys,"ltr1":d1,"ltr2":d2})
    PAIRS[sp]["ltr_diff"]=PAIRS[sp]["ltr1"]-PAIRS[sp]["ltr2"]
    m1,m2=float(np.mean(d1)),float(np.mean(d2))
    sub_summary[sp]=(cfg["g1"],m1,cfg["g2"],m2)
    print("%-9s LTR feats %-8d genes %-7d pairs %-6d joined %.1f%%  | %s %.4f  %s %.4f  ratio %s/%s = %.2f"
          % (sp,nf,len(gene),len(rec),100*rate,cfg["g1"],m1,cfg["g2"],m2,cfg["g2"],cfg["g1"],
             (m2/m1 if m1 else float('nan'))))
    if rate < 0.80:
        print("   WARNING: only %.1f%% of pairs joined to coordinates" % (100*rate))
print()
def load_genes_gff(path):
    gene={}
    for line in open(path):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="gene": continue
        c=f[0]
        if c.startswith("Scaffold") or c.startswith("utg") or c.startswith("contig"): continue
        gid=None
        for kv in f[8].split(";"):
            if kv.startswith("ID="): gid=kv[3:]; break
        if gid is None: gid="%s:%s-%s"%(c,f[3],f[4])
        gene[gid]=(c,int(f[3]),int(f[4]))
    return gene

for sp,cfg in DIPLOID.items():
    idx,nf = load_ltr(os.path.join(LTRDIR,cfg["ltr"]))
    gene = load_genes_gff(os.path.join(GFF,cfg["gff"]))
    vals=[density(idx,gene,g) for g in gene]
    vals=[v for v in vals if not np.isnan(v)]
    flag = "" if cfg["lib"]=="19K" else "   <-- DIFFERENT TE LIBRARY, absolute value not comparable"
    print("%-9s LTR feats %-8d genes %-7d  mean density %.4f  lib=%-5s%s"
          % (sp,nf,len(gene),float(np.mean(vals)),cfg["lib"],flag))
    sub_summary["diploid_"+sp]=float(np.mean(vals))

print()
print("="*96)
print("PER-PAIR TEST: does the copy with more flanking LTR have lower expression?")
print("Freeling/Woodhouse predicts rho < 0.  copy1 minus copy2 on both axes.")
print("="*96)
print("%-10s %-9s %7s | %8s %10s | %-28s" % ("species","tissue","n","rho","p","biased pairs: suppressed copy"))
print("-"*96)
res=[]
for sp,cfg in SPECIES.items():
    pr=PAIRS[sp]
    for t in cfg["tissues"]:
        f="%s/pairs_%s_%s.csv" % (HEB,sp,t)
        if not os.path.exists(f):
            print("%-10s %-9s   MISSING %s"%(sp,t,f)); continue
        h=pd.read_csv(f)
        m=pr.merge(h[["key","allo","allo_lfc"]], on="key", how="inner").dropna(subset=["ltr_diff","allo_lfc"])
        if len(m)<50: continue
        rho,p = spearmanr(m["ltr_diff"], m["allo_lfc"])
        b=m[m["allo"].isin(["toward1","toward2"])].copy()
        b["ltr_sup"]=np.where(b["allo"]=="toward1", b["ltr2"], b["ltr1"])
        b["ltr_dom"]=np.where(b["allo"]=="toward1", b["ltr1"], b["ltr2"])
        d=b["ltr_sup"]-b["ltr_dom"]
        nz=d[d!=0]
        if len(nz)>10:
            w=wilcoxon(nz)[1]; frac=float((nz>0).mean()); bt=binomtest(int((nz>0).sum()),len(nz)).pvalue
            txt="more LTR in %.1f%% (n=%d, p=%.2g)" % (100*frac,len(nz),bt)
        else:
            w=float('nan'); frac=float('nan'); txt="too few"
        print("%-10s %-9s %7d | %8.4f %10.3g | %-28s" % (sp,t,len(m),rho,p,txt))
        res.append(dict(species=sp,tissue=t,n=len(m),rho=rho,p=p,
                        n_biased=len(nz),frac_sup_more_ltr=frac,
                        wilcoxon_p=w,
                        mean_sup=float(b["ltr_sup"].mean()),mean_dom=float(b["ltr_dom"].mean())))
pd.DataFrame(res).to_csv(os.path.join(OUT,"ltr_vs_bias_all_species.csv"), index=False)
for sp in PAIRS: PAIRS[sp].to_csv(os.path.join(OUT,"ltr_per_gene_%s.csv"%sp), index=False)
print()
print("Expected under the mechanism: rho clearly negative, and suppressed copy with more LTR")
print("in well over 50%% of biased pairs. Wrote %s/ltr_vs_bias_all_species.csv" % OUT)
