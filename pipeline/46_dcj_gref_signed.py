#!/usr/bin/env python3
import os, sys, itertools, subprocess
import pandas as pd

HERE   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNIMOG = "/binfo-nas4/data/AncestralGenome/scripts/shuffle-analyze/UniMoG-java11.jar"
OUT    = os.path.join(HERE, "11_dcj_gref"); os.makedirs(OUT, exist_ok=True)

SUB = pd.read_csv(os.path.join(HERE,"00_inputs","Subgenomes_Brassica.txt"), sep="\t",
                  dtype=str).fillna("x")
COLS = {
 "An": (["Bna_A_subgenome1","Bna_A_subgenome2","Bna_A_subgenome3"], "Bna",
        ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"], "Bra", "rapa"),
 "Aj": (["Bju_A_subgenome_1","Bju_A_subgenome_2","Bju_A_subgenome_3"], "Bju",
        ["Bra_subgenome1","Bra_subgenome2","Bra_subgenome3"], "Bra", "rapa"),
 "Bj": (["Bju_B_subgenome_1","Bju_B_subgenome_2","Bju_B_subgenome_3"], "Bju",
        ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"], "Bni", "nigra"),
 "Bc": (["Bca_B_subgenome_1","Bca_B_subgenome_2","Bca_B_subgenome_3"], "Bca",
        ["Bni_subgenome1","Bni_subgenome2","Bni_subgenome3"], "Bni", "nigra"),
 "Cn": (["Bna_C_subgenome1","Bna_C_subgenome2","Bna_C_subgenome3"], "Bna",
        ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"], "Bol", "oleracea"),
 "Cc": (["Bca_C_subgenome_1","Bca_C_subgenome_2","Bca_C_subgenome_3"], "Bca",
        ["Bol_subgenome1","Bol_subgenome2","Bol_subgenome3"], "Bol", "oleracea"),
}
ORDER=["An","Aj","Bj","Bc","Cn","Cc"]
PUB={'An':36,'Aj':49,'Bj':43,'Bc':26,'Cn':42,'Cc':38}

def coords(code):
    d=pd.read_csv(os.path.join(HERE,"00_inputs",code+".saf"), sep="\t")
    g=d.groupby("GeneID").agg(chrom=("Chr","first"), start=("Start","min")).reset_index()
    g=g[~g.chrom.str.lower().str.startswith(("scaffold","utg","contig"))]
    return dict(zip(g.GeneID, zip(g.chrom, g.start)))
CO={c:coords(c) for c in ("Bna","Bju","Bca","Bra","Bni","Bol")}

def runs_for(sub, minrun):
    acols,acode,dcols,dcode,prog = COLS[sub]
    A, D = CO[acode], CO[dcode]
    pairs=[]
    for ac,dc in zip(acols,dcols):
        for ag,dg in zip(SUB[ac], SUB[dc]):
            if ag=="x" or dg=="x": continue
            if ag in A and dg in D:
                pairs.append((A[ag][0],A[ag][1],D[dg][0],D[dg][1]))
    pairs.sort(key=lambda t:(t[0],t[1]))
    runs=[]; cur=None
    def close():
        if cur and cur["n"]>=minrun: runs.append(cur)
    for ca,sa,cb,sb in pairs:
        if cur and ca==cur["ca"] and cb==cur["cb"]:
            step = sb - cur["lb"]
            if step!=0 and (cur["dir"]==0 or (step>0)==(cur["dir"]>0)):
                if cur["dir"]==0: cur["dir"] = 1 if step>0 else -1
                cur["n"]+=1; cur["lb"]=sb
                cur["b0"]=min(cur["b0"],sb); cur["b1"]=max(cur["b1"],sb); cur["a1"]=sa
                continue
        close()
        cur=dict(ca=ca,cb=cb,a0=sa,a1=sa,b0=sb,b1=sb,lb=sb,dir=0,n=1)
    close()
    return pairs, runs, prog

def dcj(runs, prog, sub):
    if not runs: return None
    for i,r in enumerate(sorted(runs,key=lambda x:(x["ca"],x["a0"])),1): r["bid"]=i
    L=[">%s"%sub]
    for c,g in itertools.groupby(sorted(runs,key=lambda x:(x["ca"],x["a0"])),key=lambda x:x["ca"]):
        L.append(" ".join(str(x["bid"]) for x in g)+" |")
    L.append(">%s"%prog)
    for c,g in itertools.groupby(sorted(runs,key=lambda x:(x["cb"],x["b0"])),key=lambda x:x["cb"]):
        L.append(" ".join(("%d" if x["dir"]>=0 else "-%d")%x["bid"] for x in g)+" |")
    p=os.path.join(OUT,"%s.unimog"%sub); open(p,"w").write("\n".join(L)+"\n")
    r=subprocess.run(["java","-jar",UNIMOG,"-m=6","-d",p],capture_output=True,text=True,timeout=3600)
    for ln in (r.stdout+r.stderr).splitlines():
        if "istance" in ln:
            t=[x for x in ln.replace(":"," ").split() if x.lstrip('-').isdigit()]
            if t: return int(t[-1])
    return None

if __name__=="__main__":
    print("DCJ-indel on Arabidopsis-anchored markers, orientation kept\n")
    print("%-7s"%"MINRUN"+"".join("%13s"%s for s in ORDER)+"   An<Cn Bc<Cc Bj<Aj")
    print("%-7s"%"pub"+"".join("%13d"%PUB[s] for s in ORDER)+"      OK    OK    OK")
    print("-"*(7+13*6+26))
    for mr in (2,3,5,8,12,20,30):
        row="%-7d"%mr; got={}
        for s in ORDER:
            pairs,runs,prog = runs_for(s,mr)
            got[s]=dcj(runs,prog,s)
            row+="%13s"%("%s/%d"%(got[s],len(runs)))
        pat=""
        for a,b in (("An","Cn"),("Bc","Cc"),("Bj","Aj")):
            pat+="%6s"%("OK" if (got[a] is not None and got[b] is not None and got[a]<got[b]) else "no")
        print(row+pat, flush=True)
