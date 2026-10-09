#!/usr/bin/env python3
import os, json
import numpy as np, pandas as pd
from scipy.stats import binomtest
REB=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","expression_rebuild")
PRIMARY="theta1_a05"
SENS=[("theta0585_a05","1.5-fold"),("theta1_a05","2-fold (primary)"),
      ("theta1585_a05","3-fold"),("span_theta1_a05","2-fold, gene span"),
      ("Q10_theta1_a05","2-fold, MAPQ>=10")]
SUBS={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
ORDER=["juncea_pollen","juncea_stigma","napus_pollen","napus_stigma",
       "carinata_pollen","carinata_stigmaE","carinata_stigmaL"]

def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    lo,hi=c-h,c+h; return lo/(1-lo), hi/(1-hi)

def heb(tag): return pd.read_csv(os.path.join(REB,"02_heb",tag,"heb_summary.csv")).set_index("sample")

main=heb(PRIMARY)
alt={t:heb(t) for t,_ in SENS if os.path.exists(os.path.join(REB,"02_heb",t,"heb_summary.csv"))}
rows=[]
for s in ORDER:
    r=main.loc[s]; sp=s.split("_")[0]; g1,g2=SUBS[sp]
    k1,k2=int(r.toward1),int(r.toward2); n=k1+k2
    ratio=k1/k2; lo,hi=wilson(k1,n); p=binomtest(k1,n,0.5).pvalue
    span=[alt[t].loc[s].toward1/alt[t].loc[s].toward2 for t,_ in SENS if t in alt and s in alt[t].index]
    rows.append(dict(sample=s, species=sp, contest="%s vs %s"%(g1,g2), sub1=g1, sub2=g2,
                     tested=int(r.tested), toward1=k1, toward2=k2, ratio=round(ratio,3),
                     ci_lo=round(lo,3), ci_hi=round(hi,3), p=p,
                     resolved=bool(not (lo<1<hi)),
                     favoured=(g1 if ratio>1 else g2) if not (lo<1<hi) else None,
                     sens_min=round(min(span),3), sens_max=round(max(span),3)))
H=pd.DataFrame(rows)

eld=pd.read_csv(os.path.join(REB,"03_eld",PRIMARY,"eld_summary.csv")) \
    if os.path.exists(os.path.join(REB,"03_eld",PRIMARY,"eld_summary.csv")) \
    else pd.read_csv(os.path.join(REB,"03_eld","theta2_a05","eld_summary.csv"))
eld["margin"]=(eld.ELD_P1_pct-eld.ELD_P2_pct).round(1)

INH=pd.read_csv(os.path.join(REB,"05_validation/inherited_vs_novel.csv"))
LTR={"An":0.005541,"Aj":0.005320,"Bj":0.004500,"Bc":0.004164,"Cn":0.008167,"Cc":0.008923}
KAKS={"An":0.1825,"Aj":0.1843,"Bj":0.1967,"Bc":0.1921,"Cn":0.2783,"Cc":0.2893}
KS={"An":0.0376,"Aj":0.0355,"Bj":0.0398,"Bc":0.0405,"Cn":0.0222,"Cc":0.0220}
FRAC={"An":0.382,"Aj":0.328,"Bj":0.375,"Bc":0.355,"Cn":0.434,"Cc":0.448}
DCJ={"An":10,"Aj":34,"Bj":33,"Bc":38,"Cn":23,"Cc":49}
DIP_LTR={"rapa":0.007733,"oleracea":0.011238}

if __name__=="__main__":
    pd.set_option("display.width",200)
    print("PRIMARY THRESHOLD:",PRIMARY,"(two-fold, CDS)\n")
    print(H[["sample","contest","tested","toward1","toward2","ratio","ci_lo","ci_hi","p","resolved","favoured","sens_min","sens_max"]].to_string(index=False))
    print("\nELD margins (percentage points toward subgenome 1):")
    print(eld[["sample","ELD_P1_pct","ELD_P2_pct","margin"]].to_string(index=False))
    H.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"canonical_heb.csv"),index=False)
    print("\nwrote tables/canonical_heb.csv")
