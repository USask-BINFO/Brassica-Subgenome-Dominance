#!/usr/bin/env python3
import glob, os, statistics, sys

LAB={'68107_68114':'An','68103_68114':'Aj','68108_68112':'Bj',
     '68110_68112':'Bc','68109_68113':'Cn','68111_68113':'Cc'}
PUB_KS   ={'An':0.037,'Aj':0.036,'Bj':0.040,'Bc':0.041,'Cn':0.022,'Cc':0.022}
PUB_SIM  ={'Cn':99.179,'Cc':99.187}
PUB_PAIRS={'An':334.56,'Cn':193.48,'Bc':169.03,'Cc':186.52}
ORDER=['An','Aj','Bj','Bc','Cn','Cc']

def model(path, merge_bool, e1=(80,95), e2=(95,100)):
    hits=open(path).readlines(); flen=len(hits)-1
    E1L,E1U=e1; E2L,E2U=e2
    recent=[]; gamma=[]
    bseqA=[]; bKs=[]; perc=[]; b_merge=False; b_hits=0
    chromA=chromB=None
    for num,hit in enumerate(hits):
        if num in (0,1,2) or hit[0:3]=='#Ks': continue
        if hit[0]!='#' and b_merge and chromA is not None and (
              chromA!=hit.split('\t')[3].split('||')[0] or
              chromB!=hit.split('\t')[7].split('||')[0]):
            bKs=[]; perc=[]; b_merge=False; b_hits=0
        if hit[0]=='#':
            pass
        else:
            b_hits+=1
            f=hit.split('\t')
            perc.append(float(f[3].split('||')[8]))
            chromA=f[3].split('||')[0]; chromB=f[7].split('||')[0]
            Ks_=f[0]
            if Ks_ not in ('NA','undef'):
                v=float(Ks_)
                if 0.001<=v<=4: bKs.append(v)
        if hit[0]=='#' or num==flen:
            if chromA is None: continue
            if chromA.title().startswith('Scaffold') or chromB.title().startswith('Scaffold'):
                bKs=[]; perc=[]; b_merge=False; b_hits=0; continue
            avgKs = 'NA' if not bKs else statistics.mean(bKs)
            if not perc: continue
            avgPerc = statistics.mean(perc)
            if merge_bool:
                pass
            else:
                if avgPerc < E1L:
                    b_merge=True; continue
            rec=[chromA,avgKs,avgPerc,float(b_hits),b_merge]
            if E1L<=avgPerc<E1U: gamma.append(rec)
            if E2L<=avgPerc<=E2U: recent.append(rec)
            bKs=[]; perc=[]; b_merge=False; b_hits=0
    return recent,gamma

def nums(l,i): return [r[i] for r in l if isinstance(r[i],float)]

CONFIGS=[("merge ON  (Merge_bool=False), 80-95/95-100", False,(80,95),(95,100)),
         ("merge OFF (Merge_bool=True),  80-95/95-100", True, (80,95),(95,100)),
         ("merge ON  (Merge_bool=False), 66-78/91-100", False,(66,78),(91,100)),
         ("merge OFF (Merge_bool=True),  66-78/91-100", True, (66,78),(91,100))]
FILES={}
for f in sorted(glob.glob("00_inputs/synmap/*.ks.txt")):
    k=os.path.basename(f).split('.')[0]
    if k in LAB: FILES[LAB[k]]=f

for label,mb,e1,e2 in CONFIGS:
    print("="*94); print(label); print("="*94)
    print("%-4s %8s %8s | %9s %9s | %9s %9s | %7s" %
          ("sub","blocks","merged","avg Ks","pub Ks","avg ident","pub sim","pairs/blk"))
    print("-"*94)
    for s in ORDER:
        rec,gam = model(FILES[s], mb, e1, e2)
        if not rec: print("%-4s   no recent blocks"%s); continue
        aKs=statistics.mean(nums(rec,1)); aId=statistics.mean(nums(rec,2))
        aPr=statistics.mean(nums(rec,3)); nm=sum(1 for r in rec if r[4])
        print("%-4s %8d %8d | %9.4f %9.3f | %9.3f %9s | %7.2f" %
              (s,len(rec),nm,aKs,PUB_KS[s],aId,
               ("%.3f"%PUB_SIM[s]) if s in PUB_SIM else "-", aPr))
    print()
