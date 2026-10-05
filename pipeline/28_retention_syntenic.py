#!/usr/bin/env python3
import glob, os, re, statistics
import numpy as np, pandas as pd

SYN="00_inputs/synmap"; OUT="08_structure"
LAB={'68107_68114':('An','napus','A','Bra'), '68103_68114':('Aj','juncea','A','Bra'),
     '68108_68112':('Bj','juncea','B','Bni'), '68110_68112':('Bc','carinata','B','Bni'),
     '68109_68113':('Cn','napus','C','Bol'),  '68111_68113':('Cc','carinata','C','Bol')}
ORDER=['An','Aj','Bj','Bc','Cn','Cc']

def genes(saf):
    d=pd.read_csv(saf, sep="\t")
    g=d.groupby("GeneID").agg(chrom=("Chr","first"), start=("Start","min")).reset_index()
    return g

PROG={p: genes(f"00_inputs/{p}.saf") for p in ("Bra","Bol","Bni")}
for p,g in PROG.items(): print("%s: %d genes" % (p, len(g)))

def partners(path):
    allg=set(); recent=set(); blk=[]; sim=[]
    for n,line in enumerate(open(path)):
        if n<3: continue
        if line.startswith('#'):
            if blk and sim and 91<=statistics.mean(sim)<=100: recent.update(blk)
            blk=[]; sim=[]; continue
        f=line.rstrip("\n").split("\t")
        if len(f)<8: continue
        i1=f[3].split('||'); i2=f[7].split('||')
        try: pid=float(i1[8])
        except (IndexError,ValueError): continue
        allg.add(i2[3]); blk.append(i2[3]); sim.append(pid)
    if blk and sim and 91<=statistics.mean(sim)<=100: recent.update(blk)
    return allg, recent

rows=[]
for path in sorted(glob.glob(f"{SYN}/*.ks.txt")):
    key=os.path.basename(path).split('.')[0]
    sub,sp,gen,prog=LAB[key]
    allg,rec=partners(path)
    G=PROG[prog].copy()
    G["kept_any"]=G.GeneID.isin(allg)
    G["kept_recent"]=G.GeneID.isin(rec)
    meds={}
    for col in ("kept_any","kept_recent"):
        fr=[]
        for ch,g in G.sort_values(["chrom","start"]).groupby("chrom"):
            v=g[col].to_numpy()
            for i in range(0,len(v)-49,50):
                fr.append(v[i:i+50].mean())
        meds[col]=(float(np.median(fr)), float(np.percentile(fr,25)),
                   float(np.percentile(fr,75)), len(fr))
    rows.append(dict(sub=sub, species=sp, genome=gen, progenitor=prog,
                     prog_genes=len(G), kept_any=int(G.kept_any.sum()),
                     frac_any=G.kept_any.mean(),
                     med_bin_any=meds["kept_any"][0], q1=meds["kept_any"][1],
                     q3=meds["kept_any"][2], bins=meds["kept_any"][3],
                     med_bin_recent=meds["kept_recent"][0]))
    print("  %-3s vs %-4s  %6d progenitor genes, %6d retained (%.3f)  median 50-gene bin %.3f [%.2f-%.2f]"
          % (sub, prog, len(G), G.kept_any.sum(), G.kept_any.mean(),
             meds["kept_any"][0], meds["kept_any"][1], meds["kept_any"][2]))

R=pd.DataFrame(rows).set_index("sub").loc[ORDER].reset_index()
R.to_csv(f"{OUT}/retention_syntenic.csv", index=False)

print("\n=== within each allotetraploid: which subgenome retains more of its progenitor? ===")
PUB={"napus":("An",0.860,"Cn",0.760), "juncea":("Aj",0.860,"Bj",0.800),
     "carinata":("Bc",0.800,"Cc",0.780)}
for sp,(a,pa,b,pb) in PUB.items():
    ra=R[R["sub"]==a].iloc[0]; rb=R[R["sub"]==b].iloc[0]
    win_new = a if ra.med_bin_any>rb.med_bin_any else b
    win_pub = a if pa>pb else b
    print("  %-9s manuscript %s %.3f vs %s %.3f  ->  %s"
          % (sp, a, pa, b, pb, win_pub))
    print("  %-9s rerun      %s %.3f vs %s %.3f  ->  %s   %s"
          % ("", a, ra.med_bin_any, b, rb.med_bin_any, win_new,
             "AGREES" if win_new==win_pub else "*** DISAGREES ***"))
print("\nwrote %s/retention_syntenic.csv" % OUT)
