#!/usr/bin/env python3
import os, sys, multiprocessing as mp
import pandas as pd, pysam

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BAMS = os.path.join(HERE, "00_inputs", "bams.tsv")
OUT  = os.path.join(HERE, "08_structure", "table_s1_alignment.csv")

SPECIES = {"Bra": "rapa", "Bni": "nigra", "Bol": "oleracea",
           "Bna": "napus", "Bju": "juncea", "Bca": "carinata"}
TISSUE  = {"pollen": "pollen", "stigma": "stigma",
           "stigma_E": "stigma (early)", "stigma_L": "stigma (late)"}

def tally(row):
    sp, tissue, rep, bam = row
    qc = mapped = paired = 0
    af = pysam.AlignmentFile(bam, "rb", check_sq=False)
    nref = af.header.nreferences
    for a in af.fetch(until_eof=True):
        if a.flag & 0x200:
            continue
        qc += 1
        if not (a.flag & 0x4):
            mapped += 1
            if a.flag & 0x2:
                paired += 1
    af.close()
    return dict(species=SPECIES.get(sp, sp), tissue=TISSUE.get(tissue, tissue), replicate=rep,
                library=os.path.basename(bam).replace(".bam", ""),
                references=nref, qc_passed=qc, mapped=mapped,
                properly_paired=paired,
                alignment_rate=100.0 * mapped / qc if qc else float("nan"))

if __name__ == "__main__":
    d = pd.read_csv(BAMS, sep="\t", header=None, names=["sp", "tissue", "rep", "bam"])
    rows = list(d.itertuples(index=False, name=None))
    print("streaming %d BAMs (%.0f GB)" % (len(rows), sum(os.path.getsize(r[3]) for r in rows)/1e9),
          flush=True)
    with mp.Pool(6) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(tally, rows), 1):
            out.append(r)
            print("  [%2d/%d] %-18s %12d records  %.2f%% mapped"
                  % (i, len(rows), r["library"], r["qc_passed"], r["alignment_rate"]), flush=True)
    T = pd.DataFrame(out)
    ORDER = ["rapa", "nigra", "oleracea", "napus", "juncea", "carinata"]
    T["_s"] = pd.Categorical(T.species, ORDER)
    T = T.sort_values(["_s", "tissue", "replicate"]).drop(columns="_s")
    T.to_csv(OUT, index=False)
    print("\nalignment rate: %.1f%% to %.1f%% over %d libraries"
          % (T.alignment_rate.min(), T.alignment_rate.max(), len(T)))
    print("total QC-passed records: %s" % format(int(T.qc_passed.sum()), ","))
    print("wrote", OUT)
