#!/usr/bin/env python3
import csv, collections, sys
SYN = "00_inputs/Subgenomes_Brassica.txt"
ABSENT = ("", "x", "-", "NA", "#N/A", "0", "#REF!", "#VALUE!", "NULL")
SP_GROUPS = {"Bra": ["A"], "Bol": ["C"], "Bni": ["B"],
             "Bna": ["A", "C"], "Bju": ["A", "B"], "Bca": ["B", "C"]}

counts, libs = {}, {}
for sp in SP_GROUPS:
    path = "01_counts/%s_Q10.txt" % sp
    with open(path) as fh:
        fh.readline()
        hdr = fh.readline().rstrip("\n").split("\t")
        names = [h.split("/")[-1].replace(".bam", "") for h in hdr[6:]]
        libs[sp] = names
        d = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[0]] = [int(x) for x in f[6:]]
        counts[sp] = d
    sys.stderr.write("%s: %d genes, %d libraries %s\n" % (sp, len(counts[sp]), len(names), names))

rows = list(csv.reader(open(SYN), delimiter="\t"))
hdr = {h.strip(): i for i, h in enumerate(rows[0])}
AT = hdr["AT_geneid"]
def col(sp, g, s):
    for p in ("%s_subgenome%d", "%s_subgenome_%d", "%s_%s_subgenome%d", "%s_%s_subgenome_%d"):
        k = p % ((sp, s) if p.count("%s") == 1 else (sp, g, s))
        if k in hdr: return hdr[k]
    return None

out_cols, meta = [], []
for sp in ("Bra", "Bol", "Bni", "Bna", "Bju", "Bca"):
    for g in SP_GROUPS[sp]:
        for lib in libs[sp]:
            out_cols.append("%s_%s|%s" % (sp, g, lib)); meta.append((sp, g, lib))

n_written = 0
present = collections.Counter()
with open("01_counts/gref_counts_Q10.tsv", "w") as fh:
    fh.write("AT_gene\tbr_sub\t" + "\t".join(out_cols) + "\n")
    for r in rows[1:]:
        at = r[AT]
        for s in (1, 2, 3):
            vals, any_present = [], False
            for sp, g, lib in meta:
                c = col(sp, g, s)
                gid = r[c].strip() if c is not None else ""
                if gid in ABSENT or gid not in counts[sp]:
                    vals.append("NA")
                else:
                    j = libs[sp].index(lib)
                    vals.append(str(counts[sp][gid][j])); any_present = True
                    present[(sp, g)] += 1
            if any_present:
                fh.write("%s\t%d\t%s\n" % (at, s, "\t".join(vals))); n_written += 1
sys.stderr.write("\nrows written: %d  columns: %d\n" % (n_written, len(out_cols)))
