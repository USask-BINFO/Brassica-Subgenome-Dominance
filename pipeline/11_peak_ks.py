#!/usr/bin/env python3
import glob, os, statistics, numpy as np, warnings
warnings.filterwarnings("ignore")
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV

LAB={'68107_68114':('An','rapa'),'68103_68114':('Aj','rapa'),'68108_68112':('Bj','nigra'),
     '68110_68112':('Bc','nigra'),'68109_68113':('Cn','oleracea'),'68111_68113':('Cc','oleracea')}
PUB={'An':0.037,'Aj':0.036,'Bj':0.040,'Bc':0.041,'Cn':0.022,'Cc':0.022}
ORDER=['An','Aj','Bj','Bc','Cn','Cc']

def load(p):
    recent=[]; blk=[]; sim=[]; nb=0
    for n,line in enumerate(open(p)):
        if n<3: continue
        if line[0]!='#':
            f=line.rstrip("\n").split("\t"); ks=f[0]
            sim.append(float(f[3].split('||')[8]))
            if ks not in ('NA','undef'): blk.append(float(ks))
        else:
            if blk and sim:
                nb+=1
                if 91<=statistics.mean(sim)<=100: recent+=blk
            blk=[]; sim=[]
    if blk and sim:
        nb+=1
        if 91<=statistics.mean(sim)<=100: recent+=blk
    return recent, nb

def gmm_mode(v):
    d=np.log10(np.array([x for x in v if 0.001<x<=4.1])).reshape(-1,1)
    g=GridSearchCV(GaussianMixture(random_state=0),
        param_grid={"n_components":range(1,5),"covariance_type":["spherical"]},
        scoring=lambda e,X:-e.bic(X))
    g.fit(d); b=g.best_estimator_
    return 10**max(zip(b.weights_, b.means_.ravel()))[1]

D={}
for f in sorted(glob.glob("00_inputs/synmap/*.ks.txt")):
    k=os.path.basename(f).split('.')[0]
    if k in LAB: D[LAB[k][0]]=(LAB[k][1],)+load(f)

print("Peak Ks per Allo-subgenome vs its diploid progenitor")
print("recent allotetraploidization blocks, mean per-pair identity 91-100%\n")
print("%-4s %-9s %7s %8s | %9s %9s | %9s %7s" %
      ("sub","progen","blocks","n Ks<=1","GMM mode","ratio","median","ratio"))
print("-"*78)
rows=[]
for s in ORDER:
    prog,recent,nb=D[s]
    k1=[x for x in recent if 0<x<=1.0]
    med=statistics.median(k1); gm=gmm_mode(recent)
    rows.append((s,prog,nb,len(k1),gm,med))
    print("%-4s %-9s %7d %8d | %9.4f %9.2f | %9.4f %7.3f" %
          (s,prog,nb,len(k1),gm,gm/PUB[s],med,med/PUB[s]))
print("-"*78)
gr=[r[4]/PUB[r[0]] for r in rows]; mr=[r[5]/PUB[r[0]] for r in rows]
print("%-4s %-9s %7s %8s | %9s %9.2f | %9s %7.3f" %
      ("","mean ratio vs Table 3","","","",statistics.mean(gr),"",statistics.mean(mr)))
print("%-4s %-9s %7s %8s | %9s %9.3f | %9s %7.3f" %
      ("","sd of ratio","","","",statistics.stdev(gr),"",statistics.stdev(mr)))
print("\n%-4s %10s %10s %10s" % ("sub","Table 3","median","diff"))
for s,prog,nb,n1,gm,med in rows:
    print("%-4s %10.3f %10.4f %+10.4f" % (s,PUB[s],med,med-PUB[s]))
print("\nRank order (low to high Ks):")
print("  Table 3 :", " < ".join(sorted(ORDER,key=lambda s:PUB[s])))
print("  median  :", " < ".join(r[0] for r in sorted(rows,key=lambda r:r[5])))
print("  GMM mode:", " < ".join(r[0] for r in sorted(rows,key=lambda r:r[4])))
