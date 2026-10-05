#!/usr/bin/env python3
import csv, re, sys, random, numpy as np
from itertools import product
D="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new"
A="/scratch2/data/Brassica/SyntenyLink/Compare"
CDS={"Bra":D+"/brapa/Brapa_genome_v3.0_cds.fasta",
     "Bol":A+"/Bol/Boleracea.v2.1.cds.fasta",
     "Bni":A+"/Bni/Bnigra_NI100.v2.cds.fasta"}
SYN = "00_inputs/Subgenomes_Brassica.txt"
ABSENT=("","x","-","NA","#N/A","0"); random.seed(1)

def load(path):
    seqs={}; name=None; buf=[]
    for line in open(path):
        if line.startswith(">"):
            if name: seqs[name]="".join(buf)
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip().upper())
    if name: seqs[name]="".join(buf)
    return seqs
S={k:load(v) for k,v in CDS.items()}
for k in S: sys.stderr.write("%s CDS: %d\n" % (k,len(S[k])))

CODON={}
bases="TCAG"
aas="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i,c in enumerate(product(bases,repeat=3)): CODON["".join(c)]=aas[i]
def syn_sites(cod):
    if cod not in CODON or CODON[cod]=="*": return None
    s=0.0
    for i in range(3):
        for b in "ACGT":
            if b==cod[i]: continue
            alt=cod[:i]+b+cod[i+1:]
            if alt in CODON and CODON[alt]==CODON[cod]: s+=1/3.0
    return s
PEP={"Bra":D+"/brapa/Brapa_genome_v3.0_pep.fasta","Bol":A+"/Bol/Boleracea.v2.1.pep.fasta","Bni":A+"/Bni/Bnigra_NI100.v2.pep.fasta"}
P={k:load(v) for k,v in PEP.items()}
from Bio import Align
_al=Align.PairwiseAligner(); _al.open_gap_score=-11; _al.extend_gap_score=-1; _al.substitution_matrix=None
_al.match_score=2; _al.mismatch_score=-1; _al.mode="global"
def codon_align(na, nb, pa, pb):
    try: aln = _al.align(pa, pb)[0]
    except Exception: return None
    A1,A2 = aln[0], aln[1]
    oa=ob=0; ca=[]; cb=[]
    for x,y in zip(A1,A2):
        if x!="-" and y!="-":
            if (oa+1)*3<=len(na) and (ob+1)*3<=len(nb):
                ca.append(na[oa*3:oa*3+3]); cb.append(nb[ob*3:ob*3+3])
        if x!="-": oa+=1
        if y!="-": ob+=1
    return "".join(ca), "".join(cb)

def ng86(a,b):
    N=Ssites=0.0; Nd=Sd=0.0; n=0
    for i in range(0,min(len(a),len(b))-2,3):
        ca,cb=a[i:i+3],b[i:i+3]
        if len(ca)<3 or len(cb)<3: break
        if any(x not in "ACGT" for x in ca+cb): continue
        if ca not in CODON or cb not in CODON: continue
        if CODON[ca]=="*" or CODON[cb]=="*": continue
        sa=syn_sites(ca); sb=syn_sites(cb)
        if sa is None or sb is None: continue
        s=(sa+sb)/2.0; Ssites+=s; N+=3-s; n+=1
        if ca==cb: continue
        diff=[j for j in range(3) if ca[j]!=cb[j]]
        if len(diff)==1:
            j=diff[0]
            if CODON[ca]==CODON[cb]: Sd+=1
            else: Nd+=1
        else:
            tot=sd=nd=0
            for order in ([diff,list(reversed(diff))] if len(diff)==2 else [diff]):
                cur=ca; ok=True; s_,n_=0,0
                for j in order:
                    nxt=cur[:j]+cb[j]+cur[j+1:]
                    if nxt not in CODON or CODON[nxt]=="*": ok=False; break
                    if CODON[cur]==CODON[nxt]: s_+=1
                    else: n_+=1
                    cur=nxt
                if ok: tot+=1; sd+=s_; nd+=n_
            if tot: Sd+=sd/tot; Nd+=nd/tot
    if n<30 or Ssites<=0 or N<=0: return None
    pS,pN=Sd/Ssites, Nd/N
    if pS>=0.75 or pN>=0.75: return None
    try:
        dS=-0.75*np.log(1-4*pS/3); dN=-0.75*np.log(1-4*pN/3)
    except Exception: return None
    if not np.isfinite(dS) or dS<=0.001 or not np.isfinite(dN): return None
    return dN/dS

rows=list(csv.reader(open(SYN),delimiter="\t")); h={x.strip():i for i,x in enumerate(rows[0])}
pairs={("Bra","Bol"):[], ("Bni","Bol"):[], ("Bra","Bni"):[]}
cand=[]
for r in rows[1:]:
    for s in (1,2,3):
        g={k:r[h["%s_subgenome%d"%(k,s)]].strip() for k in ("Bra","Bol","Bni")}
        if any(v in ABSENT for v in g.values()): continue
        cand.append(g)
random.shuffle(cand); cand=cand[:4000]
sys.stderr.write("ortholog triplets sampled: %d\n" % len(cand))
for g in cand:
    for a,b in pairs:
        sa,sb=S[a].get(g[a]), S[b].get(g[b])
        pa,pb=P[a].get(g[a]), P[b].get(g[b])
        if not sa or not sb or not pa or not pb: continue
        al=codon_align(sa,sb,pa.rstrip("*"),pb.rstrip("*"))
        if al is None: continue
        v=ng86(al[0],al[1])
        if v is not None and v<10: pairs[(a,b)].append(v)
print("\n%-16s %8s %10s %10s" % ("comparison","n","median","mean"))
for (a,b),v in pairs.items():
    v=np.array(v)
    if len(v): print("%-16s %8d %10.4f %10.4f" % ("%s vs %s"%(a,b),len(v),np.median(v),v.mean()))
print("\nIf oleracea inflates Ka/Ks, both oleracea comparisons should exceed Bra vs Bni.")
