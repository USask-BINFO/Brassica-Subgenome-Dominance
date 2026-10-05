#!/usr/bin/env python3
import pandas as pd, numpy as np, glob, os
from scipy.stats import binomtest

def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h,c+h

for th in sorted(glob.glob("02_heb/*/heb_summary.csv")):
    tag=os.path.basename(os.path.dirname(th))
    d=pd.read_csv(th)
    print("="*96); print("threshold set: %s"%tag); print("="*96)
    print("%-18s %7s %7s %7s | %7s %15s | %9s %s" %
          ("sample","sub1","sub2","n","ratio","ratio 95% CI","p","resolved?"))
    print("-"*96)
    for _,r in d.iterrows():
        k1,k2=int(r["toward1"]),int(r["toward2"]); n=k1+k2
        if n==0: continue
        ratio=k1/k2 if k2 else np.inf
        lo,hi=wilson(k1,n)
        rlo,rhi=lo/(1-lo), hi/(1-hi)
        p=binomtest(k1,n,0.5).pvalue
        res = "YES" if (rlo>1 or rhi<1) else "no - CI spans 1"
        print("%-18s %7d %7d %7d | %7.3f  [%5.3f, %5.3f] | %9.3g %s" %
              (r["sample"],k1,k2,n,ratio,rlo,rhi,p,res))
    print()
