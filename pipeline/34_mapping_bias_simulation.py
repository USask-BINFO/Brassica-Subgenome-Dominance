#!/usr/bin/env python3
import os, random, subprocess, collections, sys
import pysam

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALN  = os.path.join(HERE, "envs", "aln", "bin")
OUT  = os.path.join(HERE, "09_mapbias")
os.makedirs(OUT, exist_ok=True)
READLEN  = 151
PER_SUB  = int(sys.argv[1]) if len(sys.argv) > 1 else 500000
THREADS  = 8
random.seed(0)

B = "/scratch2/data/Brassica/transcriptome_analysis_brassica_newest/transcriptome_analysis_brassica_new"
SP = {
 "napus":    dict(saf="Bna", ref=f"{B}/bnapus/Bna_genome_v3.1.fa",
                  subs={"A":[f"N{i}" for i in range(1,11)], "C":[f"N{i}" for i in range(11,20)]}),
 "juncea":   dict(saf="Bju", ref=f"{B}/bjuncea/Bjuncea.v2.genome.fasta",
                  subs={"A":[f"A{i:02d}" for i in range(1,11)], "B":[f"B{i:02d}" for i in range(1,9)]}),
 "carinata": dict(saf="Bca", ref=f"{B}/bcarinata/Bcarinata.v2.genome.fasta",
                  subs={"B":[f"B{i}" for i in range(1,9)], "C":[f"C{i}" for i in range(1,10)]}),
}

def genes(saf, keep):
    out = collections.defaultdict(list)
    seen = {}
    with open(os.path.join(HERE, "00_inputs", saf + ".saf")) as fh:
        next(fh)
        for line in fh:
            g, c, s, e, _ = line.rstrip("\n").split("\t")
            if c not in keep: continue
            s, e = int(s), int(e)
            if g in seen:
                p = seen[g]; p[1] = min(p[1], s); p[2] = max(p[2], e)
            else:
                seen[g] = [c, s, e]
    for g, (c, s, e) in seen.items():
        if e - s + 1 >= READLEN:
            out[c].append((g, s, e))
    return out

rows = []
for sp, cfg in SP.items():
    chrom2sub = {c: k for k, v in cfg["subs"].items() for c in v}
    fa = pysam.FastaFile(cfg["ref"])
    fq = os.path.join(OUT, f"{sp}.fq")
    counts = collections.Counter()
    with open(fq, "w") as out:
        for sub, chroms in cfg["subs"].items():
            g = genes(cfg["saf"], set(chroms))
            flat = [(c, a, b) for c, lst in g.items() for (_, a, b) in lst]
            if not flat:
                print("  no genes for", sp, sub); continue
            n = 0
            while n < PER_SUB:
                c, a, b = flat[random.randrange(len(flat))]
                if b - a + 1 < READLEN: continue
                st = random.randint(a, b - READLEN + 1)
                seq = fa.fetch(c, st - 1, st - 1 + READLEN).upper()
                if seq.count("N") > READLEN // 10: continue
                out.write("@%s_%s_%d\n%s\n+\n%s\n" % (sub, c, n, seq, "I" * READLEN))
                n += 1
            counts[sub] = n
    print("%-9s simulated %s" % (sp, dict(counts)), flush=True)

    sam = os.path.join(OUT, f"{sp}.sam")
    with open(sam, "w") as so:
        subprocess.run([os.path.join(ALN, "bwa"), "mem", "-t", str(THREADS),
                        cfg["ref"], fq], stdout=so,
                       stderr=open(os.path.join(OUT, f"{sp}.bwa.log"), "w"), check=True)

    tally = collections.Counter(); mapped = collections.Counter()
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"): continue
            f = line.split("\t", 4)
            flag = int(f[1])
            if flag & 0x100 or flag & 0x800: continue
            origin = f[0].split("_")[0]
            if flag & 0x4:
                tally[(origin, "unmapped")] += 1; continue
            dest = chrom2sub.get(f[2], "other")
            tally[(origin, dest)] += 1; mapped[origin] += 1
    subs = list(cfg["subs"])
    s1, s2 = subs
    x = 100.0 * tally[(s1, s2)] / max(mapped[s1], 1)
    y = 100.0 * tally[(s2, s1)] / max(mapped[s2], 1)
    rows.append((sp, f"{s1} <-> {s2}", x, y, x - y,
                 tally[(s1,"unmapped")], tally[(s2,"unmapped")]))
    print("  %-9s %s->%s %.3f%%   %s->%s %.3f%%   asymmetry %+.3f pp"
          % (sp, s1, s2, x, s2, s1, y, x - y), flush=True)
    os.remove(sam); os.remove(fq)

print("\n%-10s %-10s %10s %10s %12s" % ("species","contrast","copy1->copy2","copy2->copy1","asymmetry"))
print("-"*58)
for sp, con, x, y, d, u1, u2 in rows:
    print("%-10s %-10s %9.3f%% %9.3f%% %+10.3f pp" % (sp, con, x, y, d))
import csv
with open(os.path.join(OUT, "mapping_bias.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["species","contrast","fwd_pct","rev_pct","asymmetry_pp",
                                    "unmapped_sub1","unmapped_sub2","reads_per_subgenome"])
    for r in rows: w.writerow(list(r) + [PER_SUB])
print("\nwrote 09_mapbias/mapping_bias.csv")
