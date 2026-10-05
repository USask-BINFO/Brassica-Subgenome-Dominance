#!/usr/bin/env python3
import glob, os, sys, random, statistics, subprocess, tempfile, shutil
from Bio import Align

HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ["PATH"]=HERE+"/envs/rm/bin:"+os.environ["PATH"]
PAL2NAL=HERE+"/envs/tools/pal2nal.v14/pal2nal.pl"
TA="/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new"
SL="/scratch2/data/Brassica/SyntenyLink/Compare"
CMP={
 'An':('68107_68114',TA+"/bnapus/Bnapus_3DH.cds_20211001.fasta",TA+"/bnapus/Bnapus_3DH.pep_20211001.fasta",
       TA+"/brapa/Brapa_genome_v3.0_cds.fasta",TA+"/brapa/Brapa_genome_v3.0_pep.fasta"),
 'Aj':('68103_68114',TA+"/bjuncea/Bjuncea_3DH.cds_20211001.fasta",TA+"/bjuncea/Bjuncea_3DH.pep_20211001.fasta",
       TA+"/brapa/Brapa_genome_v3.0_cds.fasta",TA+"/brapa/Brapa_genome_v3.0_pep.fasta"),
 'Bj':('68108_68112',TA+"/bjuncea/Bjuncea_3DH.cds_20211001.fasta",TA+"/bjuncea/Bjuncea_3DH.pep_20211001.fasta",
       SL+"/Bni/Bnigra_NI100.v2.cds.fasta",SL+"/Bni/Bnigra_NI100.v2.pep.fasta"),
 'Bc':('68110_68112',TA+"/bcarinata/Bcarinata_3DH.cds_20211001.fasta",TA+"/bcarinata/Bcarinata_3DH.pep_20211001.fasta",
       SL+"/Bni/Bnigra_NI100.v2.cds.fasta",SL+"/Bni/Bnigra_NI100.v2.pep.fasta"),
 'Cn':('68109_68113',TA+"/bnapus/Bnapus_3DH.cds_20211001.fasta",TA+"/bnapus/Bnapus_3DH.pep_20211001.fasta",
       SL+"/Bol/Boleracea.v2.1.cds.fasta",SL+"/Bol/Boleracea.v2.1.pep.fasta"),
 'Cc':('68111_68113',TA+"/bcarinata/Bcarinata_3DH.cds_20211001.fasta",TA+"/bcarinata/Bcarinata_3DH.pep_20211001.fasta",
       SL+"/Bol/Boleracea.v2.1.cds.fasta",SL+"/Bol/Boleracea.v2.1.pep.fasta"),
}
PUB={'An':0.018,'Aj':0.018,'Bj':0.020,'Bc':0.020,'Cn':0.013,'Cc':0.013}
NPAIR=int(sys.argv[1]) if len(sys.argv)>1 else 400
random.seed(1)
CP12={'GC','CG','GG','CT','CC','TC','AC','GT'}
CP3TV={'AT','AC','TA','CA','GT','GC','TG','CG'}
def cal_4dtv(ca,cb):
    n4=n4t=0
    for i in range(0,min(len(ca),len(cb))-2,3):
        a,b=ca[i:i+3],cb[i:i+3]
        if a[0:2]==b[0:2] and a[0:2] in CP12:
            n4+=1
            if a[2]+b[2] in CP3TV: n4t+=1
    return None if n4==0 else float("%.3f"%(n4t/n4))
def fasta(p):
    d={};k=None;buf=[]
    for line in open(p):
        if line[0]=='>':
            if k: d[k]="".join(buf)
            k=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if k: d[k]="".join(buf)
    return d
_al=Align.PairwiseAligner(); _al.open_gap_score=-11; _al.extend_gap_score=-1
_al.substitution_matrix=None; _al.match_score=2; _al.mismatch_score=-1; _al.mode="global"
def fast_align(na,nb,pa,pb):
    try: aln=_al.align(pa,pb)[0]
    except Exception: return None
    oa=ob=0;ca=[];cb=[]
    for x,y in zip(aln[0],aln[1]):
        if x!="-" and y!="-":
            if (oa+1)*3<=len(na) and (ob+1)*3<=len(nb):
                ca.append(na[oa*3:oa*3+3]); cb.append(nb[ob*3:ob*3+3])
        if x!="-": oa+=1
        if y!="-": ob+=1
    return "".join(ca),"".join(cb)
def pub_align(ga,gb,na,nb,pa,pb):
    cwd=os.getcwd(); td=tempfile.mkdtemp()
    try:
        os.chdir(td)
        with open("p.fa","w") as f: f.write(">%s\n%s\n>%s\n%s\n"%(ga,pa,gb,pb))
        r=subprocess.run(["clustalw","-infile=p.fa"],stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL,timeout=120)
        if not os.path.exists("p.aln"): return None
        with open("c.fa","w") as f: f.write(">%s\n%s\n>%s\n%s\n"%(ga,na,gb,nb))
        out=subprocess.run(["perl",PAL2NAL,"p.aln","c.fa","-nogap","-output","fasta"],
                           capture_output=True,text=True,timeout=120).stdout
        if not out or out.lstrip().startswith("#---  ERROR"): return None
        seqs={};k=None;buf=[]
        for line in out.splitlines():
            if line.startswith(">"):
                if k: seqs[k]="".join(buf)
                k=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if k: seqs[k]="".join(buf)
        if len(seqs)<2: return None
        v=list(seqs.values())
        return v[0],v[1]
    except Exception: return None
    finally:
        os.chdir(cwd); shutil.rmtree(td,ignore_errors=True)
def pairs_of(path,lo=91,hi=100):
    out=[];cur=[];sim=[]
    hits=open(path).readlines(); flen=len(hits)-1
    for num,line in enumerate(hits):
        if num in (0,1,2) or line[0:3]=='#Ks': continue
        if line[0]!='#':
            f=line.split('\t'); sim.append(float(f[3].split('||')[8]))
            cur.append((f[3].split('||')[3],f[7].split('||')[3]))
        if line[0]=='#' or num==flen:
            if cur and sim and lo<=statistics.mean(sim)<=hi: out+=cur
            cur=[];sim=[]
    return out
print("4DTv on IDENTICAL pairs, two aligners. %d pairs attempted per comparison." % NPAIR)
print("%-4s %7s | %9s %9s | %9s %9s | %9s" %
      ("sub","n ok","fast mean","fast med","pub mean","pub med","published"))
print("-"*76)
for s in ['An','Aj','Bj','Bc','Cn','Cc']:
    key,acds,apep,dcds,dpep=CMP[s]
    g=[p for p in glob.glob("00_inputs/synmap/*.ks.txt") if os.path.basename(p).startswith(key)]
    if not g: continue
    P=pairs_of(g[0]); random.shuffle(P)
    AC=fasta(acds);AP=fasta(apep);DC=fasta(dcds);DP=fasta(dpep)
    fv=[];pv=[];n=0
    for ga,gb in P:
        if n>=NPAIR: break
        na,pa,nb,pb=AC.get(ga),AP.get(ga),DC.get(gb),DP.get(gb)
        if not(na and pa and nb and pb): continue
        fa=fast_align(na,nb,pa,pb); pb_=pub_align(ga,gb,na,nb,pa,pb)
        if not fa or not pb_: continue
        a=cal_4dtv(*fa); b=cal_4dtv(*pb_)
        if a is None or b is None: continue
        fv.append(a); pv.append(b); n+=1
    if not fv: print("%-4s   none usable"%s); continue
    print("%-4s %7d | %9.4f %9.4f | %9.4f %9.4f | %9.3f" %
          (s,len(fv),statistics.mean(fv),statistics.median(fv),
           statistics.mean(pv),statistics.median(pv),PUB[s]))
