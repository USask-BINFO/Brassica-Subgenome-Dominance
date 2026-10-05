#!/usr/bin/env python3
import pickle, collections, sys
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def load_obo(path):
    parents=collections.defaultdict(set); name={}; ns={}; cur=None; obs=set()
    for line in open(path):
        line=line.rstrip("\n")
        if line=="[Term]": cur=None; continue
        if line.startswith("id: GO:"): cur=line[4:]
        elif cur and line.startswith("name: "): name[cur]=line[6:]
        elif cur and line.startswith("namespace: "): ns[cur]={"biological_process":"P","molecular_function":"F","cellular_component":"C"}.get(line[11:],"?")
        elif cur and line.startswith("is_a: "): parents[cur].add(line[6:16])
        elif cur and line.startswith("relationship: part_of "): parents[cur].add(line[22:32])
        elif cur and line.startswith("is_obsolete: true"): obs.add(cur)
    for o in obs: parents.pop(o,None); name.pop(o,None); ns.pop(o,None)
    return parents,name,ns

def ancestors(t,parents,cache):
    if t in cache: return cache[t]
    out=set()
    for p in parents.get(t,()):
        out.add(p); out |= ancestors(p,parents,cache)
    cache[t]=out; return out

parents,name,ns = load_obo("00_inputs/go-basic.obo")
raw = pickle.load(open("../reanalysis/scripts/go_map.pkl","rb"))["g2go"]
cache={}
g2go={g:set().union(*[{t}|ancestors(t,parents,cache) for t in ts]) if ts else set() for g,ts in raw.items()}
import statistics
print("propagation: %.1f -> %.1f GO terms per gene (median %d -> %d)" % (
    statistics.mean(len(v) for v in raw.values()), statistics.mean(len(v) for v in g2go.values()),
    statistics.median([len(v) for v in raw.values()]), statistics.median([len(v) for v in g2go.values()])))
pickle.dump({"g2go":g2go,"name":name,"ns":ns}, open("00_inputs/go_propagated.pkl","wb"))

def enrich(fg, bg, label, aspect="P", topn=12, minobs=3):
    fg = {g for g in fg if g in g2go and g2go[g]}
    bg = {g for g in bg if g in g2go and g2go[g]}
    fg &= bg
    cf=collections.Counter(t for g in fg for t in g2go[g] if ns.get(t)==aspect)
    cb=collections.Counter(t for g in bg for t in g2go[g] if ns.get(t)==aspect)
    rows=[(t,k,cb[t],hypergeom.sf(k-1,len(bg),cb[t],len(fg))) for t,k in cf.items() if k>=minobs and cb[t]<len(bg)*0.5]
    print("\n=== %s === (%d fg / %d bg genes, aspect %s)" % (label,len(fg),len(bg),aspect))
    if not rows: print("   no terms with >=%d genes"%minobs); return
    q=multipletests([r[3] for r in rows],method="fdr_bh")[1]
    for (t,k,b,p),qq in sorted(zip(rows,q),key=lambda x:x[1])[:topn]:
        star="  <<<" if qq<0.05 else ""
        print("   %-11s %3d/%-5d p=%.2e FDR=%.3f  %s%s" % (t,k,b,p,qq,name.get(t,"?")[:52],star))
if __name__=="__main__": pass
