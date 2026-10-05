#!/usr/bin/env python3
import glob, os, statistics, collections, gzip
import numpy as np, pandas as pd

IN  = "00_inputs/synmap"
OUT = "08_structure"
os.makedirs(OUT, exist_ok=True)

LAB = {'68107_68114':('An','napus','A','rapa'),   '68103_68114':('Aj','juncea','A','rapa'),
       '68108_68112':('Bj','juncea','B','nigra'), '68110_68112':('Bc','carinata','B','nigra'),
       '68109_68113':('Cn','napus','C','oleracea'),'68111_68113':('Cc','carinata','C','oleracea')}
ORDER = ['An','Aj','Bj','Bc','Cn','Cc']

def parse(path):
    pairs, blocks = [], []
    blk = []
    for n, line in enumerate(open(path)):
        if n < 3:
            continue
        if line.startswith('#'):
            if blk and any(np.isfinite(r['ks']) for r in blk):
                blocks.append(blk)
            blk = []
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8:
            continue
        ks, kn = f[0], f[1]
        i1, i2 = f[3].split('||'), f[7].split('||')
        try:
            pid = float(i1[8])
        except (IndexError, ValueError):
            continue
        rec = dict(chr_allo=i1[0], chr_prog=i2[0],
                   ks=(float(ks) if ks not in ('NA', 'undef') else np.nan),
                   kn=(float(kn) if kn not in ('NA', 'undef') else np.nan),
                   pid=pid)
        blk.append(rec)
    if blk and any(np.isfinite(r['ks']) for r in blk):
        blocks.append(blk)
    for b in blocks:
        bm = statistics.mean(r['pid'] for r in b)
        for r in b:
            r['block_mean_pid'] = bm
            r['recent'] = 91 <= bm <= 100
            pairs.append(r)
    return pd.DataFrame(pairs), len(blocks)

allp, summ, coll = [], [], []
for p in sorted(glob.glob(f"{IN}/*.ks.txt")):
    key = os.path.basename(p).split('.')[0]
    sub, sp, gen, prog = LAB[key]
    d, nblk = parse(p)
    d['subg'] = sub; d['species'] = sp; d['genome'] = gen; d['progenitor'] = prog
    allp.append(d)

    rec = d[d.recent]
    ks_rec = rec.ks.dropna()
    ks_rec = ks_rec[(ks_rec > 0) & (ks_rec <= 1.0)]
    summ.append(dict(subg=sub, species=sp, genome=gen, progenitor=prog,
                     blocks=nblk, pairs=len(d), recent_pairs=int(rec.shape[0]),
                     median_ks=float(ks_rec.median()), n_ks=int(len(ks_rec)),
                     median_pid=float(rec.pid.median()),
                     mean_pid=float(rec.pid.mean())))

    rc = rec[~rec.chr_allo.str.contains('caffold') & ~rec.chr_prog.str.contains('caffold')]
    for ca, g in rc.groupby('chr_allo'):
        if len(g) < 50:
            continue
        vc = g.chr_prog.value_counts()
        coll.append(dict(subg=sub, species=sp, genome=gen, chr_allo=ca,
                         pairs=len(g), top_prog=vc.index[0], top_share=float(vc.iloc[0]/len(g)),
                         n_prog_chr=int((vc >= 0.05*len(g)).sum())))

D = pd.concat(allp, ignore_index=True)
D.to_csv(f"{OUT}/divergence_pairs.csv.gz", index=False, compression="gzip")
S = pd.DataFrame(summ).set_index('subg').loc[ORDER].reset_index()
S.to_csv(f"{OUT}/divergence_summary.csv", index=False)
C = pd.DataFrame(coll)
C.to_csv(f"{OUT}/collinearity_by_chr.csv", index=False)

print("=== divergence over recent allotetraploidization blocks ===")
print(S[['subg','species','genome','progenitor','blocks','recent_pairs',
         'median_ks','median_pid']].round(4).to_string(index=False))

print("\n=== collinearity, averaged over chromosomes ===")
cs = C.groupby(['subg','species','genome']).agg(
        chroms=('chr_allo','size'), mean_top_share=('top_share','mean'),
        min_top_share=('top_share','min'),
        chr_with_split=('n_prog_chr', lambda s: int((s > 1).sum()))).reset_index()
cs['subg'] = pd.Categorical(cs['subg'], ORDER)
cs = cs.sort_values('subg')
cs.to_csv(f"{OUT}/collinearity_summary.csv", index=False)
print(cs.round(4).to_string(index=False))

print("\n=== per-genome contrast within each allotetraploid ===")
for sp, sub1, sub2 in [("napus","An","Cn"), ("juncea","Aj","Bj"), ("carinata","Bc","Cc")]:
    a = cs[cs.subg == sub1].iloc[0]; b = cs[cs.subg == sub2].iloc[0]
    sa = S[S.subg == sub1].iloc[0];  sb = S[S.subg == sub2].iloc[0]
    print("  %-9s %s collinearity %.4f vs %s %.4f   |   median identity %.2f vs %.2f   |   median Ks %.4f vs %.4f"
          % (sp, sub1, a.mean_top_share, sub2, b.mean_top_share,
             sa.median_pid, sb.median_pid, sa.median_ks, sb.median_ks))
print("\nwrote %s/divergence_pairs.csv.gz, divergence_summary.csv, collinearity_by_chr.csv, collinearity_summary.csv" % OUT)
