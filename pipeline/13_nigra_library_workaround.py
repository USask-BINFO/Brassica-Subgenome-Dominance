#!/usr/bin/env python3
import os, bisect, collections, csv, re, sys
import numpy as np

LTR="00_inputs/ltr"; GFF="00_inputs/coge_gff"; FLANK=5000
BED="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new/Research/BWA"
TA="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new"

def load_ltr(p):
    iv=collections.defaultdict(list)
    for line in open(p):
        if line.startswith("#"): continue
        f=line.split("\t")
        if len(f)<5: continue
        try: iv[f[0]].append((int(f[3]),int(f[4])))
        except ValueError: continue
    out={}
    for c,l in iv.items():
        l.sort()
        m=[]
        for a,b in l:
            if m and a<=m[-1][1]+1: m[-1]=(m[-1][0],max(m[-1][1],b))
            else: m.append((a,b))
        out[c]=([x[0] for x in m], m)
    return out
def bp(idx,c,s,e):
    if c not in idx or e<s: return 0
    st,l=idx[c]; i=bisect.bisect_left(st,s)-1; t=0
    while i<len(l):
        if i>=0:
            a,b=l[i]
            if a>e: break
            t+=max(0,min(b,e)-max(a,s)+1)
        i+=1
    return t
def dens(idx,c,s,e):
    return (bp(idx,c,max(0,s-FLANK),s-1)+bp(idx,c,e+1,e+FLANK))/(2.0*FLANK)

print("="*92)
print("ROUTE A: each subgenome measured in TWO independent allotetraploids, common 19K library")
print("="*92)
ALLO={"napus":    ("Bna_genome_v3.1.fa.19K_repeat.LTR.gff","Bna","Bna_sub",["A","C"]),
      "juncea":   ("Bjuncea_genome_v1.fa.out.LTR.gff","Bju","Bju_sub",["A","B"]),
      "carinata": ("Bcarinata.v1.genome.fasta.out.LTR.gff","Bca","Bca_sub",["B","C"])}
per={}
for sp,(lt,bd,pre,groups) in ALLO.items():
    idx=load_ltr(os.path.join(LTR,lt))
    for g in groups:
        vals=[]
        for s in (1,2,3):
            p="%s/%s/%s%d%s.bed"%(BED,bd,pre,s,g)
            if not os.path.exists(p): continue
            for line in open(p):
                f=line.rstrip("\n").split("\t")
                if len(f)<6: continue
                vals.append(dens(idx,f[0],int(f[1]),int(f[2])))
        per[(g,sp)]=(len(vals),float(np.mean(vals)))
        print("  %s in %-9s n=%-6d mean %.6f" % (g,sp,len(vals),float(np.mean(vals))))
print()
print("%-3s %-26s %-26s %9s" % ("sub","host 1","host 2","agreement"))
for g,(h1,h2) in (("A",("napus","juncea")),("B",("juncea","carinata")),("C",("napus","carinata"))):
    n1,v1=per[(g,h1)]; n2,v2=per[(g,h2)]
    print("%-3s %-26s %-26s %8.1f%%" % (g,"%s %.6f"%(h1,v1),"%s %.6f"%(h2,v2),
          100*abs(v1-v2)/((v1+v2)/2)))
print()
mean={g:np.mean([per[(g,h)][1] for h in ("napus","juncea","carinata") if (g,h) in per]) for g in "ABC"}
print("pooled per subgenome:  A %.6f   B %.6f   C %.6f" % (mean["A"],mean["B"],mean["C"]))
print("ratios:  C/A = %.2f   C/B = %.2f   A/B = %.2f" % (mean["C"]/mean["A"],mean["C"]/mean["B"],mean["A"]/mean["B"]))

print()
print("="*92)
print("ROUTE B: background-normalised enrichment in the three diploids")
print("="*92)
DIP={"rapa":    ("Brapa_sequence_v3.0.fasta.19K_repeat.LTR.gff",
                 "Brassica_rapa_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68114.gff",
                 TA+"/brapa/Brapa_sequence_v3.0.fasta.fai","19K"),
     "oleracea":("Boleracea.v2.1.genome.fasta.19K_repeat.LTR.gff",
                 "Brassica_oleracea_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68113.gff",
                 TA+"/boleracea/to1000/Boleracea.v2.1.genome.fasta.fai","19K"),
     "nigra":   ("Bnigra_NI100.v2.genome.fasta.out.LTR.gff",
                 "Brassica_nigra_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68112.gff",
                 TA+"/bnigra/Bnigra_NI100.v2.genome.fasta.fai","OTHER")}
print("%-9s %-6s %7s %12s %12s %12s %10s" %
      ("genome","lib","genes","gene-prox","genome-wide","enrichment","chr Mb"))
print("-"*78)
E={}
for sp,(lt,gf,fai,lib) in DIP.items():
    idx=load_ltr(os.path.join(LTR,lt))
    L={}
    for line in open(fai):
        f=line.split("\t"); L[f[0]]=int(f[1])
    chroms=[c for c in idx if c in L and not c.startswith(("Scaffold","utg","contig"))]
    tot_len=sum(L[c] for c in chroms)
    tot_ltr=sum(sum(b-a+1 for a,b in idx[c][1]) for c in chroms)
    gw=tot_ltr/tot_len
    vals=[]
    for line in open(os.path.join(GFF,gf)):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="gene": continue
        if f[0] not in idx or f[0] not in L: continue
        if f[0].startswith(("Scaffold","utg","contig")): continue
        vals.append(dens(idx,f[0],int(f[3]),int(f[4])))
    gp=float(np.mean(vals)); E[sp]=(gp,gw,gp/gw)
    print("%-9s %-6s %7d %12.6f %12.6f %12.3f %10.1f" %
          (sp,lib,len(vals),gp,gw,gp/gw,tot_len/1e6))
print()
print("ABSOLUTE gene-proximal density (library-dependent, nigra NOT comparable):")
print("   oleracea/rapa = %.3f    oleracea/nigra = %.3f  <-- unusable" %
      (E["oleracea"][0]/E["rapa"][0], E["oleracea"][0]/E["nigra"][0]))
print("ENRICHMENT near genes relative to genome background (library largely cancels):")
print("   oleracea/rapa = %.3f    oleracea/nigra = %.3f" %
      (E["oleracea"][2]/E["rapa"][2], E["oleracea"][2]/E["nigra"][2]))
print()
print("Compare with the allotetraploid partner ratios measured on one library:")
print("   napus    C/A = %.2f" % (mean["C"]/mean["A"]))
print("   carinata C/B = %.2f" % (mean["C"]/mean["B"]))
