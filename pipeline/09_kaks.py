#!/usr/bin/env python3
import glob, os, statistics, numpy as np
D=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "00_inputs", "synmap")
LAB={'68107_68114':('An','rapa'),'68103_68114':('Aj','rapa'),'68108_68112':('Bj','nigra'),
     '68110_68112':('Bc','nigra'),'68109_68113':('Cn','oleracea'),'68111_68113':('Cc','oleracea')}
def ratios(path, lo=91, hi=100):
    out=[]; blk=[]; sim=[]
    for num,line in enumerate(open(path)):
        if num in (0,1,2): continue
        if line[0]!='#':
            p=line.rstrip("\n").split("\t")
            ks,kn=p[0],p[1]
            sim.append(float(p[3].split('||')[8]))
            if ks not in ('NA','undef') and kn not in ('NA','undef'):
                k,n=float(ks),float(kn)
                if 0.001<k<=4.1 and n>=0: blk.append(n/k if k>0 else np.nan)
        else:
            if blk and sim and lo<=statistics.mean(sim)<=hi: out+=blk
            blk=[]; sim=[]
    return np.array([x for x in out if np.isfinite(x) and x<10])
res={}
print("%-6s %-10s %8s %10s %10s %10s" % ("sub","progenitor","n pairs","median","mean","Q1-Q3"))
for f in sorted(glob.glob(os.path.join(D,"*.ks.txt"))):
    k=os.path.basename(f).split('.')[0]
    if k not in LAB: continue
    sub,prog=LAB[k]; r=ratios(f); res[sub]=r
    q=np.percentile(r,[25,75])
    print("%-6s %-10s %8d %10.4f %10.4f  %.3f-%.3f" % (sub,prog,len(r),np.median(r),r.mean(),q[0],q[1]))
print()
from scipy.stats import mannwhitneyu
print("Partner comparison (prediction: the DOMINANT subgenome has the LOWER Ka/Ks)")
for lab,a,b,dom in (("napus    An vs Cn","An","Cn","A on most measures"),
                    ("carinata Bc vs Cc","Bc","Cc","B"),
                    ("juncea   Aj vs Bj","Aj","Bj","B")):
    x,y=res[a],res[b]
    u,p=mannwhitneyu(x,y)
    lower = a if np.median(x)<np.median(y) else b
    print("  %-20s median %.4f vs %.4f   lower = %-3s  p=%.2e   | expression favours %s"
          % (lab,np.median(x),np.median(y),lower,p,dom))
