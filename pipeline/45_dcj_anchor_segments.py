#!/usr/bin/env python3
import os, sys, glob, itertools, statistics, subprocess, tempfile, shutil, collections

HERE   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYN    = os.path.join(HERE, "00_inputs", "synmap")
UNIMOG = "/binfo-nas4/data/AncestralGenome/scripts/shuffle-analyze/UniMoG-java11.jar"
LAB={'68107_68114':('An','rapa'),'68103_68114':('Aj','rapa'),'68108_68112':('Bj','nigra'),
     '68110_68112':('Bc','nigra'),'68109_68113':('Cn','oleracea'),'68111_68113':('Cc','oleracea')}
ORDER=['An','Aj','Bj','Bc','Cn','Cc']
PUB={'An':36,'Aj':49,'Bj':43,'Bc':26,'Cn':42,'Cc':38}
IDLO=91.0

def anchors(path):
    out=[]; cur=None; buf=[]
    def flush():
        if cur and buf:
            ident=statistics.mean(x[4] for x in buf)
            if IDLO<=ident<=100:
                out.extend([(cur[0],x[0],cur[1],x[1],x[2],x[3]) for x in buf])
    for line in open(path):
        if line.startswith('#'):
            f=line.rstrip("\n").split("\t")
            if len(f)>=5 and f[0][1:].isdigit():
                flush(); buf=[]
                ca=f[2].split('_',1)[1]; cb=f[3].split('_',1)[1]
                cur=(ca,cb) if not (ca.lower().startswith(("scaffold","utg","contig")) or
                                    cb.lower().startswith(("scaffold","utg","contig"))) else None
            continue
        if cur is None: continue
        f=line.rstrip("\n").split("\t")
        if len(f)<10: continue
        p1=f[3].split('||'); p2=f[7].split('||')
        try:
            oa=int(p1[7]); ob=int(p2[7]); pid=float(p1[8])
        except (IndexError,ValueError): continue
        buf.append((oa,ob,p1[3],p2[3],pid))
    flush()
    return out

def segments(anc, gap, mn):
    seen_a=set(); seen_b=set(); clean=[]
    for a in anc:
        if a[4] in seen_a or a[5] in seen_b: continue
        seen_a.add(a[4]); seen_b.add(a[5]); clean.append(a)
    clean.sort(key=lambda x:(x[0],x[1]))
    segs=[]; cur=None
    for a in clean:
        ca,oa,cb,ob = a[0],a[1],a[2],a[3]
        if cur and cur["ca"]==ca and cur["cb"]==cb:
            step = ob - cur["last"]
            if step!=0 and abs(step)<=gap and (cur["dir"]==0 or (step>0)==(cur["dir"]>0)):
                cur["n"]+=1; cur["last"]=ob
                cur["b0"]=min(cur["b0"],ob); cur["b1"]=max(cur["b1"],ob); cur["a1"]=oa
                if cur["dir"]==0: cur["dir"]=1 if step>0 else -1
                continue
        cur=dict(ca=ca,cb=cb,a0=oa,a1=oa,b0=ob,b1=ob,last=ob,dir=0,n=1); segs.append(cur)
    return [s for s in segs if s["n"]>=mn]

def dcj(segs):
    if not segs: return None
    for i,s in enumerate(sorted(segs,key=lambda x:(x["ca"],x["a0"])),1): s["bid"]=i
    L=[">allo"]
    for c,g in itertools.groupby(sorted(segs,key=lambda x:(x["ca"],x["a0"])),key=lambda x:x["ca"]):
        L.append(" ".join(str(x["bid"]) for x in g)+" |")
    L.append(">diploid")
    for c,g in itertools.groupby(sorted(segs,key=lambda x:(x["cb"],x["b0"])),key=lambda x:x["cb"]):
        L.append(" ".join(("%d" if x["dir"]>=0 else "-%d")%x["bid"] for x in g)+" |")
    td=tempfile.mkdtemp(); p=os.path.join(td,"u.unimog"); open(p,"w").write("\n".join(L)+"\n")
    r=subprocess.run(["java","-jar",UNIMOG,"-m=6","-d",p],capture_output=True,text=True,timeout=3600)
    shutil.rmtree(td)
    for ln in (r.stdout+r.stderr).splitlines():
        if "istance" in ln:
            t=[x for x in ln.replace(":"," ").split() if x.lstrip('-').isdigit()]
            if t: return int(t[-1])
    return None

if __name__=="__main__":
    ANC={}
    for path in sorted(glob.glob(SYN+"/*.ks.txt")):
        sub,prog=LAB[os.path.basename(path).split('.')[0]]
        ANC[sub]=anchors(path)
    for s in ORDER: print("%s: %d recent anchors"%(s,len(ANC[s])), flush=True)
    print("\n%-6s %-5s"%("GAP","MIN")+"".join("%13s"%s for s in ORDER)+"   An<Cn Bc<Cc Bj<Aj")
    print("%-12s"%"published"+"".join("%13d"%PUB[s] for s in ORDER)+"      OK    OK    OK")
    print("-"*(12+13*6+26))
    for gap in (5,10,20,50):
        for mn in (2,3,5,8):
            row="%-6d %-5d"%(gap,mn); got={}
            for sub in ORDER:
                sg=segments(ANC[sub],gap,mn); got[sub]=dcj(sg)
                row+="%13s"%("%s/%d"%(got[sub],len(sg)))
            pat=""
            for a,b in (("An","Cn"),("Bc","Cc"),("Bj","Aj")):
                pat+="%6s"%("OK" if (got[a] is not None and got[b] is not None and got[a]<got[b]) else "no")
            print(row+pat, flush=True)
