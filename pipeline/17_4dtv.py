#!/usr/bin/env python3
import glob, os, sys, random, statistics, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio import Align

TA="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new"
CMP={
 'An':('68107_68114', TA+"/bnapus/Bnapus_3DH.cds_20211001.fasta",    TA+"/bnapus/Bnapus_3DH.pep_20211001.fasta",
        TA+"/brapa/Brapa_genome_v3.0_cds.fasta", TA+"/brapa/Brapa_genome_v3.0_pep.fasta"),
 'Aj':('68103_68114', TA+"/bjuncea/Bjuncea_3DH.cds_20211001.fasta",  TA+"/bjuncea/Bjuncea_3DH.pep_20211001.fasta",
        TA+"/brapa/Brapa_genome_v3.0_cds.fasta", TA+"/brapa/Brapa_genome_v3.0_pep.fasta"),
 'Bj':('68108_68112', TA+"/bjuncea/Bjuncea_3DH.cds_20211001.fasta",  TA+"/bjuncea/Bjuncea_3DH.pep_20211001.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bni/Bnigra_NI100.v2.cds.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bni/Bnigra_NI100.v2.pep.fasta"),
 'Bc':('68110_68112', TA+"/bcarinata/Bcarinata_3DH.cds_20211001.fasta", TA+"/bcarinata/Bcarinata_3DH.pep_20211001.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bni/Bnigra_NI100.v2.cds.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bni/Bnigra_NI100.v2.pep.fasta"),
 'Cn':('68109_68113', TA+"/bnapus/Bnapus_3DH.cds_20211001.fasta",    TA+"/bnapus/Bnapus_3DH.pep_20211001.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bol/Boleracea.v2.1.cds.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bol/Boleracea.v2.1.pep.fasta"),
 'Cc':('68111_68113', TA+"/bcarinata/Bcarinata_3DH.cds_20211001.fasta", TA+"/bcarinata/Bcarinata_3DH.pep_20211001.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bol/Boleracea.v2.1.cds.fasta",
        "/scratch2/data/Brassica/SyntenyLink/Compare/Bol/Boleracea.v2.1.pep.fasta"),
}
PUB={'An':0.018,'Aj':0.018,'Bj':0.020,'Bc':0.020,'Cn':0.013,'Cc':0.013}
WIN=(int(sys.argv[1]) if len(sys.argv)>1 else 91, 100)
NBLOCK=int(sys.argv[2]) if len(sys.argv)>2 else 60
random.seed(1)

CP12={'GC','CG','GG','CT','CC','TC','AC','GT'}
CP3TV={'AT','AC','TA','CA','GT','GC','TG','CG'}
def cal_4dtv(ca,cb):
    n4=n4t=0
    for i in range(0,min(len(ca),len(cb))-2,3):
        a=ca[i:i+3]; b=cb[i:i+3]
        if a[0:2]==b[0:2] and a[0:2] in CP12:
            n4+=1
            if a[2]+b[2] in CP3TV: n4t+=1
    return 0.0 if n4==0 else float("%.3f"%(n4t/n4))

def fasta(path):
    d={}; k=None; buf=[]
    for line in open(path):
        if line[0]=='>':
            if k: d[k]="".join(buf)
            k=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if k: d[k]="".join(buf)
    return d

_al=Align.PairwiseAligner(); _al.open_gap_score=-11; _al.extend_gap_score=-1
_al.substitution_matrix=None; _al.match_score=2; _al.mismatch_score=-1; _al.mode="global"
def codon_align(na,nb,pa,pb):
    try: aln=_al.align(pa,pb)[0]
    except Exception: return None
    oa=ob=0; ca=[]; cb=[]
    for x,y in zip(aln[0],aln[1]):
        if x!="-" and y!="-":
            if (oa+1)*3<=len(na) and (ob+1)*3<=len(nb):
                ca.append(na[oa*3:oa*3+3]); cb.append(nb[ob*3:ob*3+3])
        if x!="-": oa+=1
        if y!="-": ob+=1
    return "".join(ca),"".join(cb)

def blocks(path, lo, hi):
    out=[]; cur=[]; sim=[]
    hits=open(path).readlines(); flen=len(hits)-1
    for num,line in enumerate(hits):
        if num in (0,1,2) or line[0:3]=='#Ks': continue
        if line[0]!='#':
            f=line.split('\t')
            sim.append(float(f[3].split('||')[8]))
            cur.append((f[3].split('||')[3], f[7].split('||')[3]))
        if line[0]=='#' or num==flen:
            if cur and sim and lo<=statistics.mean(sim)<=hi: out.append(cur[:])
            cur=[]; sim=[]
    return out

print("4DTv, recent blocks identity %d-%d, %d blocks sampled per comparison" % (WIN[0],WIN[1],NBLOCK))
print("%-4s %7s %7s | %9s %9s %9s | %9s" %
      ("sub","blocks","pairs","pair mean","pair med","blockmean","published"))
print("-"*74)
for s in ['An','Aj','Bj','Bc','Cn','Cc']:
    key,acds,apep,dcds,dpep=CMP[s]
    f=[p for p in glob.glob("00_inputs/synmap/*.ks.txt") if os.path.basename(p).startswith(key)]
    if not f: print("%-4s  synmap missing"%s); continue
    B=blocks(f[0],*WIN)
    random.shuffle(B); B=B[:NBLOCK]
    AC=fasta(acds); AP=fasta(apep); DC=fasta(dcds); DP=fasta(dpep)
    per=[]; bmeans=[]
    for blk in B:
        buf=[]
        for ga,gb in blk:
            na=AC.get(ga); pa=AP.get(ga); nb=DC.get(gb); pb=DP.get(gb)
            if not(na and pa and nb and pb): continue
            r=codon_align(na,nb,pa,pb)
            if not r: continue
            buf.append(cal_4dtv(*r))
        if buf:
            bmeans.append(statistics.mean(buf)); per+=buf
    if not per: print("%-4s  no usable pairs"%s); continue
    print("%-4s %7d %7d | %9.4f %9.4f %9.4f | %9.3f" %
          (s,len(bmeans),len(per),statistics.mean(per),statistics.median(per),
           statistics.mean(bmeans),PUB[s]), flush=True)
    os.makedirs("08_structure/4dtv", exist_ok=True)
    with open("08_structure/4dtv/%s.txt" % s, "w") as fh:
        fh.write("\n".join("%.4f" % v for v in per) + "\n")
