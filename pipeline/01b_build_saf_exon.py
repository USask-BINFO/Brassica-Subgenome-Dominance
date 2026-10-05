#!/usr/bin/env python3
import csv, re, sys, collections
ANN, OUT, SP, FEAT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
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

def tid_of(a):
    m = re.search(r'transcript_id "([^"]+)"', a)
    if m: return m.group(1)
    m = re.search(r'(?:^|;)Parent=([^;,]+)', a)
    return m.group(1) if m else None
def gid_of(a, tid):
    m = re.search(r'gene_id "([^"]+)"', a)
    if m: return m.group(1)
    return None

tx_iv = collections.defaultdict(list)
for line in open(ANN):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != FEAT: continue
    t = tid_of(f[8])
    if t: tx_iv[t].append((f[0], int(f[3]), int(f[4]), f[6]))

def clean(iv):
    u = sorted(set(iv), key=lambda x: (x[0], x[1], x[2]))
    keep = []
    for a in u:
        if any(b is not a and b[0] == a[0] and a[1] <= b[1] and b[2] <= a[2]
               and (a[1] < b[1] or b[2] < a[2]) for b in u):
            continue
        keep.append(a)
    return keep or u
tx_gene = {}
for line in open(ANN):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("mRNA", "transcript"): continue
    m = re.search(r'ID=([^;]+)', f[8]) or re.search(r'transcript_id "([^"]+)"', f[8])
    p = re.search(r'Parent=([^;,]+)', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
    if m: tx_gene[m.group(1)] = p.group(1) if p else m.group(1)
for t in tx_iv: tx_gene.setdefault(t, t)

by_gene = collections.defaultdict(list)
for t in tx_iv: by_gene[tx_gene[t]].append(t)
tlen = lambda t: sum(e - s + 1 for _, s, e, _ in clean(tx_iv[t]))
chosen = [next((t for t in txs if t in wanted), None) or max(txs, key=tlen) for txs in by_gene.values()]

hit = len(wanted & set(chosen)); rate = 100.0 * hit / len(wanted) if wanted else 0
sys.stderr.write("%-4s %-5s genes %6d  transcripts %6d  synteny covered %6d/%-6d = %6.2f%%\n"
                 % (SP, FEAT, len(by_gene), len(chosen), hit, len(wanted), rate))
if rate < 95.0: sys.exit("ABORT: synteny coverage %.2f%% for %s" % (rate, SP))

n = 0
with open(OUT, "w") as fh:
    fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for t in chosen:
        merged = []
        for c, s, e, st in sorted(clean(tx_iv[t]), key=lambda x: (x[0], x[1])):
            if merged and merged[-1][0] == c and s <= merged[-1][2] + 1:
                merged[-1][2] = max(merged[-1][2], e)
            else: merged.append([c, s, e, st])
        for c, s, e, st in merged:
            fh.write("%s\t%s\t%d\t%d\t%s\n" % (t, c, s, e, st)); n += 1
sys.stderr.write("     %d intervals\n" % n)
