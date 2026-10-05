#!/usr/bin/env python3
import os, sys, glob, collections, csv, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
A = __import__("45_dcj_anchor_segments")
M = __import__("40_dcj_normalized")

HERE  = A.HERE
ORDER = A.ORDER
PROG  = {"An":"rapa","Aj":"rapa","Bj":"nigra","Bc":"nigra","Cn":"oleracea","Cc":"oleracea"}
MINEX = 5
MINID = 95.0

def synmap_assign(sub):
    path = [p for p in glob.glob(A.SYN+"/*.ks.txt")
            if A.LAB[os.path.basename(p).split('.')[0]][0] == sub][0]
    best = {}
    for ca, oa, cb, ob, na, nb in A.anchors(path):
        best.setdefault(na, (ca, oa, cb, ob))
    return [(g,)+v for g, v in best.items()]

def mcscanx_assign(sub):
    wd = os.path.join(HERE, "10_mcscanx", sub)
    G  = M.genes(wd, sub)
    best = {}
    for line in open(os.path.join(wd, sub + ".blast")):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12: continue
        g1, g2 = f[0], f[1]
        if g1 not in G or g2 not in G: continue
        c1 = G[g1][0]; c2 = G[g2][0]
        if not (c1.startswith("aa") and c2.startswith("bb")): continue
        s = float(f[11])
        if g1 not in best or s > best[g1][0]:
            best[g1] = (s, c1[2:], G[g1][1], c2[2:], float(f[2]), G[g2][1])
    return [(g, v[1], v[2], v[3], v[5]) for g, v in best.items() if v[4] >= MINID]

def exchange_runs(assign, minex):
    tot = collections.Counter(); bychr = collections.defaultdict(list)
    prank = collections.defaultdict(list)
    for g, ac, rank, dc, drank in assign:
        tot[(ac, dc)] += 1; bychr[ac].append((rank, g, dc, drank)); prank[dc].append(drank)
    pmin = {c: min(v) for c, v in prank.items()}; pmax = {c: max(v) for c, v in prank.items()}
    partner = {ac: max((n, dc) for (a, dc), n in tot.items() if a == ac)[1] for ac in bychr}
    runs = []
    nchr = {ac: len(v) for ac, v in bychr.items()}
    for ac, items in bychr.items():
        items.sort(); cur = None
        for i, (rank, g, dc, drank) in enumerate(items):
            if cur and cur["dc"] == dc:
                cur["genes"].append(g); cur["i1"] = i
                cur["d0"] = min(cur["d0"], drank); cur["d1"] = max(cur["d1"], drank); continue
            if cur and cur["dc"] != partner[ac] and len(cur["genes"]) >= minex: runs.append(cur)
            cur = dict(ac=ac, dc=dc, partner=partner[ac], genes=[g], i0=i, i1=i,
                       d0=drank, d1=drank)
        if cur: cur["i1"] = i
        if cur and cur["dc"] != partner[ac] and len(cur["genes"]) >= minex: runs.append(cur)
    for r in runs:
        lo, hi = pmin[r["dc"]], pmax[r["dc"]]
        span = max(hi - lo, 1)
        r["pd0"] = (r["d0"] - lo) / span; r["pd1"] = (r["d1"] - lo) / span
    return runs, partner, tot, nchr

rows = []; summary = {}; matrices = {}
for sub in ORDER:
    got = {}
    for setname, assign in (("DAGChainer", synmap_assign(sub)), ("MCScanX", mcscanx_assign(sub))):
        runs, partner, tot, nchr = exchange_runs(assign, MINEX)
        got[setname] = runs; matrices[(sub, setname)] = (tot, partner, nchr)
        for r in runs:
            n = len(r["genes"])
            rows.append(dict(subgenome=sub, progenitor=PROG[sub], block_set=setname,
                             allo_chr=r["ac"], partner_chr=r["partner"], exchange_chr=r["dc"],
                             genes=n, chr_genes=nchr[r["ac"]],
                             start_frac=round(r["i0"]/max(1,nchr[r["ac"]]-1), 4),
                             end_frac=round(r["i1"]/max(1,nchr[r["ac"]]-1), 4),
                             first_gene=r["genes"][0], last_gene=r["genes"][-1],
                             donor_start_frac=round(r["pd0"], 4),
                             donor_end_frac=round(r["pd1"], 4),
                             gene_list=";".join(r["genes"])))
    gd = {g for r in got["DAGChainer"] for g in r["genes"]}
    gm = {g for r in got["MCScanX"]    for g in r["genes"]}
    summary[sub] = (len(got["DAGChainer"]), len(got["MCScanX"]), len(gd), len(gm), len(gd & gm))

os.makedirs("08_structure", exist_ok=True)
with open("08_structure/exchanges.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
import pickle
with open("08_structure/exchange_matrices.pkl","wb") as fh: pickle.dump(matrices, fh)

crows = []
for (sub, setname), (tot, partner, nchr) in matrices.items():
    for ac in sorted(nchr):
        off = sum(n for (a, dc), n in tot.items() if a == ac and dc != partner[ac])
        crows.append(dict(subgenome=sub, block_set=setname, allo_chr=ac,
                          partner_chr=partner[ac], genes=nchr[ac], genes_off_partner=off))
with open("08_structure/exchange_chromosomes.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(crows[0].keys())); w.writeheader(); w.writerows(crows)
print("wrote 08_structure/exchange_chromosomes.csv (%d rows)" % len(crows))

print("exchanges: uninterrupted chains of >= %d genes, partner identity >= %.0f%%\n" % (MINEX, MINID))
print("%-4s %-9s %11s %9s %11s %11s %13s" % ("sub","progenitor","DAGChainer","MCScanX",
      "DAG genes","MCS genes","genes shared"))
print("-"*74)
for sub in ORDER:
    nd, nm, gd, gm, sh = summary[sub]
    print("%-4s %-9s %11d %9d %11d %11d %13d" % (sub, PROG[sub], nd, nm, gd, gm, sh))
print("\nwrote 08_structure/exchanges.csv  (%d runs) and exchange_matrices.pkl" % len(rows))
