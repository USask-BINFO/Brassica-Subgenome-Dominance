#!/usr/bin/env python3
import os, sys, glob, statistics
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYN  = os.path.join(HERE, "00_inputs", "synmap")
LAB  = {'68107_68114':('An','rapa'),'68103_68114':('Aj','rapa'),'68108_68112':('Bj','nigra'),
        '68110_68112':('Bc','nigra'),'68109_68113':('Cn','oleracea'),'68111_68113':('Cc','oleracea')}
ORDER = ['An','Aj','Bj','Bc','Cn','Cc']
PUB   = {'An':0.382,'Aj':0.328,'Bj':0.375,'Bc':0.355,'Cn':0.434,'Cc':0.448}
MAXGAP = 10

def blocks(path):
    out=[]; cur=None; pid=[]; anc=[]
    def flush():
        if cur and anc: out.append((statistics.mean(pid) if pid else 0.0, cur, list(anc)))
    for line in open(path):
        if line.startswith('#'):
            f=line.rstrip("\n").split("\t")
            if len(f)>=5 and f[0][1:].isdigit():
                flush(); pid=[]; anc=[]
                ca=f[2].split('_',1)[1]; cb=f[3].split('_',1)[1]
                cur=None if (ca.lower().startswith(("scaffold","utg","contig")) or
                             cb.lower().startswith(("scaffold","utg","contig"))) else (ca,cb)
            continue
        if cur is None: continue
        f=line.rstrip("\n").split("\t")
        if len(f)<10: continue
        p1=f[3].split('||'); p2=f[7].split('||')
        try:
            oa=int(p1[7]); ob=int(p2[7]); pid.append(float(p1[8]))
        except (IndexError,ValueError): continue
        anc.append((oa,ob))
    flush()
    return out

def ratio(anchors, maxgap=MAXGAP):
    a=sorted(set(anchors), key=lambda t:t[1])
    n=len(a)
    if n<2: return None
    g=0
    for (x0,y0),(x1,y1) in zip(a,a[1:]):
        d=abs(y1-y0)-1
        if 0 < d <= maxgap: g+=d
    den=(n-1)+g
    return g/den if den>0 else None

DATA={}
for path in sorted(glob.glob(SYN+"/*.ks.txt")):
    sub,prog=LAB[os.path.basename(path).split('.')[0]]
    DATA[sub]=blocks(path)

def ratio2(anchors, side, maxgap, strict):
    a=sorted(set(anchors), key=lambda t:t[1] if side=="prog" else t[0])
    if len(a)<2: return None,None
    g=0
    for p,q in zip(a,a[1:]):
        d=abs((q[1]-p[1]) if side=="prog" else (q[0]-p[0]))-1
        if 0<d and ((d<maxgap) if strict else (d<=maxgap)): g+=d
    return g, len(a)-1

def table(side, maxgap, strict, agg, lo):
    row={}
    for s in ORDER:
        gs=[]
        for ident,cur,a in DATA[s]:
            if not (lo<=ident<=100): continue
            g,d=ratio2(a,side,maxgap,strict)
            if g is not None: gs.append((g,d))
        row[s]=float(np.mean([g/(d+g) for g,d in gs if d+g>0])) if agg=="mean" \
               else sum(g for g,_ in gs)/(sum(d for _,d in gs)+sum(g for g,_ in gs))
    return row

if __name__ == "__main__":
    import csv, itertools
    CLAIMS = {
      "the two C subgenomes are the two highest":
          lambda r: sorted(r,key=r.get)[-2:]==sorted(["Cn","Cc"],key=r.get),
      "napus An < Cn":    lambda r: r["An"]<r["Cn"],
      "carinata Bc < Cc": lambda r: r["Bc"]<r["Cc"],
      "juncea Aj < Bj":   lambda r: r["Aj"]<r["Bj"],
      "Aj is the lowest of the six": lambda r: min(r,key=r.get)=="Aj",
    }
    grid=list(itertools.product((9,10,11),(False,True),("mean","pooled"),(90,91,92,93)))
    print("fractionation ratio = sum(gaps) / [(n-1) + sum(gaps)]\n")
    print("%-34s"%"reported"+"".join("%9.3f"%PUB[s] for s in ORDER))
    best=None; store={}
    for side in ("prog","allo"):
        rows=[table(side,*g) for g in grid]
        store[side]=rows
        for g,row in zip(grid,rows):
            d=float(np.mean([abs(row[s]-PUB[s]) for s in ORDER]))
            if side=="prog" and (best is None or d<best[0]): best=(d,g,row)
    d,g,row=best
    print("%-34s"%("closest: gaps on progenitor, maxgap %d, identity %d-100"%(g[0],g[3]))
          +"".join("%9.3f"%row[s] for s in ORDER)+"   mean |diff| %.4f"%d)
    print("\nhow often each claim holds, over %d settings per side"%len(grid))
    print("%-42s %10s %10s"%("claim","progenitor","Allo side"))
    print("-"*64)
    for name,fn in CLAIMS.items():
        print("%-42s %9.0f%% %9.0f%%"%(name,
            100*np.mean([fn(r) for r in store["prog"]]),
            100*np.mean([fn(r) for r in store["allo"]])))
    os.makedirs("08_structure", exist_ok=True)
    with open("08_structure/fractionation_reconstruction.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["setting"]+ORDER)
        w.writerow(["reported"]+[PUB[s] for s in ORDER])
        for gg,rr in zip(grid,store["prog"]):
            w.writerow(["prog maxgap=%d strict=%s agg=%s id=%d-100"%gg]+["%.4f"%rr[s] for s in ORDER])
    print("\nwrote 08_structure/fractionation_reconstruction.csv")
