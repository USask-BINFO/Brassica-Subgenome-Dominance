#!/usr/bin/env python3
import os, bisect, collections
import numpy as np

RM="06_repeatmasker/diploids"; GFF="00_inputs/coge_gff"; FLANK=5000
SP={"rapa":     ("rapa/rapa.19K.LTR.gff",         "Brassica_rapa_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68114.gff"),
    "nigra":    ("nigra/nigra.19K.LTR.gff",       "Brassica_nigra_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68112.gff"),
    "oleracea": ("oleracea/oleracea.19K.LTR.gff", "Brassica_oleracea_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68113.gff")}

def feats(path):
    out=collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<5: continue
        try: out[f[0]].append((int(f[3]),int(f[4])))
        except ValueError: continue
    return out

def merged(d):
    o={}
    for c,l in d.items():
        l.sort(); m=[]
        for a,b in l:
            if m and a<=m[-1][1]+1: m[-1]=(m[-1][0],max(m[-1][1],b))
            else: m.append((a,b))
        o[c]=([x[0] for x in m], m)
    return o

def bp_in(idx,c,s,e):
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
    return (bp_in(idx,c,max(0,s-FLANK),s-1)+bp_in(idx,c,e+1,e+FLANK))/(2.0*FLANK)

D={}
for sp,(gffltr,gene_gff) in SP.items():
    idx=merged(feats(os.path.join(RM,gffltr)))
    genes=[]
    for line in open(os.path.join(GFF,gene_gff)):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="gene": continue
        if f[0].startswith(("Scaffold","utg","contig")): continue
        genes.append((f[0],int(f[3]),int(f[4])))
    v=[dens(idx,c,s,e) for c,s,e in genes if c in idx]
    D[sp]=float(np.mean(v))
    print("  %-9s %6d genes on assembled chromosomes   density %.6f" % (sp,len(v),D[sp]))

print("\nAll three from one run, so the pipeline offset cancels in every ratio below.\n")
ALLO={"napus (C/A)":(1.460,1.474), "carinata (C/B)":(2.073,2.143), "juncea (A/B)":(1.108,1.182)}
PR={"napus (C/A)":("oleracea / rapa",  D["oleracea"]/D["rapa"]),
    "carinata (C/B)":("oleracea / nigra",D["oleracea"]/D["nigra"]),
    "juncea (A/B)":("rapa / nigra",     D["rapa"]/D["nigra"])}
print("%-16s %-18s %10s %18s   %s" % ("pairing","progenitor ratio","progenitors","inside the hybrid","verdict"))
print("-"*92)
for k in ("napus (C/A)","carinata (C/B)","juncea (A/B)"):
    lab,pr=PR[k]; lo,hi=ALLO[k]
    inside = lo*0.95 <= pr <= hi*1.05
    print("%-16s %-18s %10.2f %8.2f-%-8.2f   %s"
          % (k,lab,pr,lo,hi,"inherited" if inside else "WIDENED in the hybrid"))
pr=D["oleracea"]/D["nigra"]; lo,hi=ALLO["carinata (C/B)"]
print("\ncarinata: progenitors differ by %.2f, subgenomes by %.2f-%.2f, a gap %.0f-%.0f%% wider."
      % (pr,lo,hi,100*(lo/pr-1),100*(hi/pr-1)))
