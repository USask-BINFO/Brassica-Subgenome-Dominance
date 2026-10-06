#!/usr/bin/env python3
import os, numpy as np, pandas as pd
from scipy.stats import binomtest
HERE=os.path.dirname(os.path.abspath(__file__))
REB=os.path.join(HERE,"..","..","expression_rebuild")
OUT=os.path.join(HERE,"supp_tables.md")
TH=[("theta0585_a05","1.5-fold"),("theta1_a05","2-fold (primary)"),
    ("span_theta1_a05","2-fold, gene span"),("theta1585_a05","3-fold"),
    ("Q10_theta1_a05","2-fold, MAPQ>=10")]
ORDER=["juncea_pollen","juncea_stigma","napus_pollen","napus_stigma",
       "carinata_pollen","carinata_stigmaE","carinata_stigmaL"]
NICE={s:s.replace("_"," ").replace("stigmaE","stigma (early)").replace("stigmaL","stigma (late)") for s in ORDER}
D={t:pd.read_csv(os.path.join(REB,"02_heb",t,"heb_summary.csv")).set_index("sample") for t,_ in TH}
L=[]

L.append("### Table S8. Dominance ratio at every analysis setting\n")
L.append("Ratio of biased pairs, subgenome 1 : subgenome 2. The primary analysis is two-fold with "
         "CDS counting; the other four columns are sensitivity checks.\n")
L.append("| sample | "+" | ".join(l for _,l in TH)+" | direction |")
L.append("|:-|"+"|".join([":-"]*len(TH))+"|:-|")
for s in ORDER:
    v=[D[t].loc[s].toward1/D[t].loc[s].toward2 for t,_ in TH]
    d="all below 1" if max(v)<1 else ("all above 1" if min(v)>1 else "CROSSES 1")
    L.append("| %s | "%NICE[s]+" | ".join("%.3f"%x for x in v)+" | %s |"%d)
L.append("\nThe direction of every contest is the same under all five settings. Only the "
         "significance of napus pollen varies, and its ratio stays between 0.890 and 0.969.\n")

def wilson(k,n,z=1.96):
    p=k/n; dd=1+z*z/n; c=(p+z*z/(2*n))/dd; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/dd
    lo,hi=c-h,c+h; return lo/(1-lo),hi/(1-hi)
SUB={"napus":("A","C"),"juncea":("A","B"),"carinata":("B","C")}
L.append("### Table S8b. The primary analysis in full\n")
L.append("| sample | contest | tested | toward sub 1 | toward sub 2 | ratio | 95% CI | p | resolved |")
L.append("|:-|:-|:-|:-|:-|:-|:-|:-|:-|")
for s in ORDER:
    r=D["theta1_a05"].loc[s]; sp=s.split("_")[0]; g1,g2=SUB[sp]
    k1,k2=int(r.toward1),int(r.toward2); n=k1+k2
    lo,hi=wilson(k1,n); p=binomtest(k1,n,0.5).pvalue
    L.append("| %s | %s vs %s | %d | %d | %d | %.3f | %.3f-%.3f | %s | %s |"%(
        NICE[s],g1,g2,int(r.tested),k1,k2,k1/k2,lo,hi,
        ("%.2g"%p if p>=1e-4 else "%.0e"%p),"yes" if not (lo<1<hi) else "no"))

NS=pd.read_csv(os.path.join(REB,"05_validation/novel_switched_split.csv"))
L.append("\n### Table S9. Biased pairs by what the diploid progenitors did\n")
L.append("Switched pairs are anti-parental by definition, so their direction is not interpretable "
         "and is shown for completeness only.\n")
L.append("| sample | maintained n (ratio) | truly novel n (ratio) | switched n (ratio) | maintained share |")
L.append("|:-|:-|:-|:-|:-|")
for s in ORDER:
    sp,ti=s.split("_",1)
    g=NS[(NS.species==sp)&(NS.tissue==ti)]
    if g.empty: continue
    tot=(g.t1+g.t2).sum()
    cells=[]
    for c in ("maintained","novel","switched"):
        r=g[g.cls==c].iloc[0]; n=int(r.t1+r.t2)
        cells.append("%d (%.3f)"%(n,r.ratio))
    m=g[g.cls=="maintained"].iloc[0]
    L.append("| %s | %s | %s | %s | %.0f%% |"%(NICE[s],cells[0],cells[1],cells[2],
             100*(m.t1+m.t2)/tot))

BR=pd.read_csv(os.path.join(REB,"05_validation/heb_by_br_subgenome.csv"))
L.append("\n### Table S11b. Dominance ratio within each Br-subgenome layer\n")
L.append("| sample | LF | MF1 | MF2 | whole sample |")
L.append("|:-|:-|:-|:-|:-|")
for s in ORDER:
    sp,ti=s.split("_",1)
    g=BR[(BR.species==sp)&(BR.tissue==ti)]
    if g.empty: continue
    def cell(l):
        r=g[g.layer==l]
        if r.empty: return "-"
        r=r.iloc[0]; star="*" if r.p<0.05 else ""
        return "%.3f%s"%(r.ratio,star)
    whole=D["theta1_a05"].loc[s]
    L.append("| %s | %s | %s | %s | %.3f |"%(NICE[s],cell("LF"),cell("MF1"),cell("MF2"),
             whole.toward1/whole.toward2))
L.append("\n\\* lean differs from 1:1 at p < 0.05. In napus and carinata the lean is carried by "
         "MF1; in napus pollen LF leans the other way, which is why the whole sample does not resolve.\n")

import glob as _glob
ELD = os.path.join(REB, "03_eld", "theta1_a05")
_cls = {}
for _s in ORDER:
    _f = os.path.join(ELD, "eld_%s.csv" % _s)
    if os.path.exists(_f):
        _cls[_s] = pd.read_csv(_f)["class"].value_counts()
_CL = pd.DataFrame(_cls).fillna(0).astype(int).T
_pat = [c for c in _CL.columns if c != "unclassified"]
def _rank(c):
    return (0 if c.startswith("PED") else 1,
            0 if "ELD_P1" in c else 1 if "ELD_P2" in c else 2 if "additive" in c
            else 3 if "unchanged" in c else 4)
_pat = sorted(_pat, key=_rank)
L.append("\n### Table S9b. Inheritance pattern counts per sample\n")
L.append("PED, the progenitors differ; PEC, they do not. Within each, the pair's total expression "
         "either resembles one progenitor (ELD), is additive, is unchanged, or is transgressive. "
         "Percentages are of classified pairs.\n")
L.append("| sample | " + " | ".join(c.replace("_", " ") for c in _pat) + " | classified |")
L.append("|:-|" + "-|" * (len(_pat) + 1))
for _s in ORDER:
    if _s not in _CL.index: continue
    _r = _CL.loc[_s, _pat]; _tot = int(_r.sum())
    L.append("| %s | %s | %d |" % (NICE[_s],
             " | ".join("%d (%.1f%%)" % (v, 100.0*v/_tot) for v in _r), _tot))
L.append("")

open(OUT,"w").write("\n".join(L))
print("wrote tables/supp_tables.md (%d lines)"%len(L))
print("\n".join(L[:14]))
