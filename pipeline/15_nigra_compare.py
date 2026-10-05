#!/usr/bin/env python3
import sys, os, bisect, collections
import numpy as np

LTR="00_inputs/ltr"; GFF="00_inputs/coge_gff"; FLANK=5000
RM="06_repeatmasker"
TARGET=dict(n=970, bp=332379, chrlen=29591584)

def feats(path, chrom=None):
    out=collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<5: continue
        if chrom and f[0]!=chrom: continue
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

stage=sys.argv[1] if len(sys.argv)>1 else "validate"

if stage=="validate":
    mine_p=os.path.join(RM,"validate_rapa","rapa_A01.mine.LTR.gff")
    if not os.path.exists(mine_p): sys.exit("not found: %s - run scripts/14 validate first"%mine_p)
    mine=feats(mine_p)
    n=sum(len(v) for v in mine.values()); b=sum(y-x+1 for v in mine.values() for x,y in v)
    print("VALIDATION on rapa chromosome A01, same library, my RepeatMasker run vs Sampath's")
    print("%-26s %10s %12s %10s" % ("","features","LTR bp","% of chr"))
    print("%-26s %10d %12d %9.4f%%" % ("Sampath",TARGET["n"],TARGET["bp"],100*TARGET["bp"]/TARGET["chrlen"]))
    print("%-26s %10d %12d %9.4f%%" % ("mine",n,b,100*b/TARGET["chrlen"]))
    dn=100*(n-TARGET["n"])/TARGET["n"]; db=100*(b-TARGET["bp"])/TARGET["bp"]
    print("%-26s %9.1f%% %11.1f%%" % ("difference",dn,db))
    print()
    if abs(dn)<=10 and abs(db)<=10:
        print("PASS: settings reproduce Sampath's calls. nigra run will be comparable.")
    else:
        print("FAIL: settings differ from Sampath's. Do NOT treat a nigra run as comparable")
        print("      to the other five until this is reconciled (try -s, or -q, or ask Sampath")
        print("      for the exact command).")
    sys.exit(0)

mine_p=os.path.join(RM,"nigra","nigra.19K.LTR.gff")
if not os.path.exists(mine_p): sys.exit("not found: %s - run scripts/14 nigra first"%mine_p)
idx=merged(feats(mine_p))
old=merged(feats(os.path.join(LTR,"Bnigra_NI100.v2.genome.fasta.out.LTR.gff")))
gff=os.path.join(GFF,"Brassica_nigra_annos0-cds0-id_typename-nu1-upa1-add_chr0.gid68112.gff")
genes=[]
for line in open(gff):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="gene": continue
    if f[0].startswith(("Scaffold","utg","contig")): continue
    genes.append((f[0],int(f[3]),int(f[4])))
new=[dens(idx,c,s,e) for c,s,e in genes if c in idx]
oldv=[dens(old,c,s,e) for c,s,e in genes if c in old]
NG=float(np.mean(new))
print("nigra gene-proximal LTR density, %d genes on chromosomes B1-B8" % len(new))
print("  old calls (different library) : %.6f" % float(np.mean(oldv)))
print("  NEW, common 19K library       : %.6f" % NG)
print()
RAPA, OLER = 0.007733, 0.011238
print("Diploid progenitors, all now on the common library:")
print("  rapa %.6f   nigra %.6f   oleracea %.6f" % (RAPA,NG,OLER))
print()
print("%-34s %8s %8s" % ("comparison","diploids","allotetraploid partner ratio"))
print("%-34s %8.2f   napus    C/A = 1.46-1.47" % ("oleracea / rapa", OLER/RAPA))
print("%-34s %8.2f   carinata C/B = 2.07-2.14" % ("oleracea / nigra", OLER/NG))
print("%-34s %8.2f   juncea   A/B = 1.11-1.18" % ("rapa / nigra", RAPA/NG))
print()
print("The inheritance claim holds for carinata if oleracea/nigra lands near 2.1, and for")
print("juncea if rapa/nigra lands near 1.1.")
