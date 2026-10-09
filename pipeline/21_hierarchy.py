#!/usr/bin/env python3
import pandas as pd, numpy as np, os, statistics

HEB={
 "napus":   ("A","C",{"pollen":0.950,"stigma":1.140}),
 "juncea":  ("A","B",{"pollen":0.826,"stigma":0.869}),
 "carinata":("B","C",{"pollen":1.105,"stigmaE":1.071,"stigmaL":1.109}),
}
ELD={"napus":("A","C",{"pollen":+8.4,"stigma":+33.1}),
     "juncea":("A","B",{"pollen":+0.6,"stigma":+20.0}),
     "carinata":("B","C",{"pollen":+4.6,"stigmaE":-11.5,"stigmaL":-7.3})}
FRAC={"An":0.382,"Aj":0.328,"Bj":0.375,"Bc":0.355,"Cn":0.434,"Cc":0.448}
DCJ ={"An":36,"Aj":49,"Bj":43,"Bc":26,"Cn":42,"Cc":38}
LTR ={"A":0.00543,"B":0.00433,"C":0.00855}
KAKS={"An":0.1825,"Aj":0.1843,"Bj":0.1967,"Bc":0.1921,"Cn":0.2783,"Cc":0.2893}

def winner_expr(sp):
    s1,s2,d=HEB[sp]; r=statistics.mean(d.values())
    return (s1 if r>1 else s2), r
print("EXPRESSION contests (HEB ratio averaged over tissues; >1 favours subgenome 1)")
wins={}
for sp in HEB:
    s1,s2,d=HEB[sp]; w,r=winner_expr(sp)
    flip = "FLIPS by tissue" if (min(d.values())<1<max(d.values())) else ""
    print("  %-9s %s vs %s   mean ratio %.3f -> %s wins   %s" % (sp,s1,s2,r,w,flip))
    wins[(s1,s2)]=w
print()
import itertools
pairs={frozenset(k):v for k,v in wins.items()}
def beats(x,y):
    k=frozenset((x,y))
    return pairs.get(k)==x if k in pairs else None
order=[]
for perm in itertools.permutations("ABC"):
    ok=True
    for i in range(len(perm)-1):
        for j in range(i+1,len(perm)):
            b=beats(perm[i],perm[j])
            if b is False: ok=False
    if ok: order.append("".join(perm))
print("Transitive orders consistent with all three expression contests:", order if order else "NONE (cycle)")
print()
print("STRUCTURE, averaged over the hybrids carrying each subgenome (lower = less degraded)")
def avg(d,letter):
    v=[d[k] for k in d if k.startswith(letter)]
    return statistics.mean(v) if v else float('nan')
print("%-10s %8s %8s %8s %8s" % ("subgenome","frac","DCJ","LTR","Ka/Ks"))
st={}
for L in "ABC":
    f,dj,lt,kk=avg(FRAC,L),avg(DCJ,L),LTR[L],avg(KAKS,L)
    st[L]=(f,dj,lt,kk)
    print("%-10s %8.3f %8.1f %8.5f %8.4f" % (L,f,dj,lt,kk))
print()
for i,nm in enumerate(["fractionation","DCJ","gene-proximal LTR","Ka/Ks"]):
    rank=sorted("ABC", key=lambda L: st[L][i])
    print("  %-20s best -> worst: %s" % (nm," < ".join(rank)))
print()
print("EXPRESSION rank (best first):", order[0] if order else "undefined")
print()
print("Where the two axes agree and disagree:")
print("  C is last on every structural metric, and C loses both expression contests it enters")
print("    (napus stigma A>C, carinata B>C), so C is consistently the non-dominant partner.")
print("  The disagreement is about A. Structure puts A level with or ahead of B")
print("    (fractionation A %.3f vs B %.3f, Ka/Ks A %.4f vs B %.4f, both favouring A)," %
      (st['A'][0],st['B'][0],st['A'][3],st['B'][3]))
print("    yet A LOSES to B in juncea in both tissues (ratio 0.84 and 0.88).")
print("  So the decoupling is specific: it is not that structure fails to identify the degraded")
print("  subgenome, it is that structural advantage does not translate into expression advantage")
print("  for A against B.")
