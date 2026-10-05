#!/usr/bin/env python3
import csv, re, sys, collections
GFF, OUT, SP = sys.argv[1], sys.argv[2], sys.argv[3]
SYN = "00_inputs/Subgenomes_Brassica.txt"
ABSENT = ("", "x", "-", "NA", "#N/A", "0", "#REF!", "#VALUE!", "NULL")

rows = list(csv.reader(open(SYN), delimiter="\t"))
hdr = {h.strip(): i for i, h in enumerate(rows[0])}
wanted = set()
for k, i in hdr.items():
    if k.startswith(SP + "_"):
        for r in rows[1:]:
            v = r[i].strip()
            if v not in ABSENT: wanted.add(v)
sys.stderr.write("%s: %d transcript ids referenced by the synteny table\n" % (SP, len(wanted)))

def attr(s, key):
    m = re.search(r'(?:^|;)%s=([^;]+)' % key, s)
    return m.group(1) if m else None

tx_gene, tx_cds = {}, collections.defaultdict(list)
for line in open(GFF):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9: continue
    if f[2] == "mRNA":
        tid, gid = attr(f[8], "ID"), attr(f[8], "Parent")
        if tid: tx_gene[tid] = gid or tid
    elif f[2] == "CDS":
        p = attr(f[8], "Parent")
        if p: tx_cds[p].append((f[0], int(f[3]), int(f[4]), f[6]))
for t in tx_cds:
    tx_gene.setdefault(t, t)

by_gene = collections.defaultdict(list)
for t in tx_cds: by_gene[tx_gene[t]].append(t)

def cdslen(t): return sum(e - s + 1 for _, s, e, _ in tx_cds[t])
chosen = []
for g, txs in by_gene.items():
    pick = next((t for t in txs if t in wanted), None)
    if pick is None: pick = max(txs, key=cdslen)
    chosen.append(pick)

hit = len(wanted & set(chosen))
rate = 100.0 * hit / len(wanted) if wanted else 0.0
sys.stderr.write("%s: %d genes, %d transcripts chosen; synteny ids covered %d/%d = %.2f%%\n"
                 % (SP, len(by_gene), len(chosen), hit, len(wanted), rate))
if rate < 95.0:
    sys.exit("ABORT: only %.2f%% of synteny transcript ids are in the SAF for %s" % (rate, SP))

n = 0
with open(OUT, "w") as fh:
    fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for t in chosen:
        iv = sorted(tx_cds[t], key=lambda x: (x[0], x[1]))
        merged = []
        for c, s, e, st in iv:
            if merged and merged[-1][0] == c and s <= merged[-1][2] + 1:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([c, s, e, st])
        for c, s, e, st in merged:
            fh.write("%s\t%s\t%d\t%d\t%s\n" % (t, c, s, e, st)); n += 1
sys.stderr.write("%s: %d SAF intervals written\n\n" % (SP, n))
