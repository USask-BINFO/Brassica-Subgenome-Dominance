#!/usr/bin/env python3
import os, sys, re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
ST = os.path.join(HERE, "..", "..", "expression_rebuild", "08_structure")

E = pd.read_csv(os.path.join(ST, "exchanges.csv"))
C = pd.read_csv(os.path.join(ST, "exchange_chromosomes.csv"))
ALL = E.copy()
E = E[E.block_set == "DAGChainer"]; C = C[C.block_set == "DAGChainer"]
_mcs = [set(str(g).split(";")) for g in ALL.query("block_set=='MCScanX'").gene_list]
both = {r.first_gene for _, r in E.iterrows()
        if any(set(str(r.gene_list).split(";")) & y for y in _mcs)}

GEN  = {"An":"A","Aj":"A","Bj":"B","Bc":"B","Cn":"C","Cc":"C"}
SPEC = [("napus", "An", "Cn"), ("juncea", "Aj", "Bj"), ("carinata", "Bc", "Cc")]
MINW = 0.011
LABEL_MIN = 10

def disp(sub, c):
    c = str(c)
    if sub == "Cn":
        m = re.match(r'^C1(\d)$', c)
        if m: return "C%s" % m.group(1)
    m = re.match(r'^([A-Za-z]+)0*(\d+)$', c)
    return "%s%s" % (m.group(1), m.group(2)) if m else c

def nat(c):
    m = re.match(r'^([A-Za-z]*)(\d+)$', str(c))
    return (m.group(1), int(m.group(2))) if m else (str(c), 0)

fig, axes = plt.subplots(3, 2, figsize=(7.1, 7.4), layout="constrained")
fig.get_layout_engine().set(w_pad=0.07, h_pad=0.07, wspace=0.10, hspace=0.10)

for row, (sp, s1, s2) in enumerate(SPEC):
    for col, sub in enumerate((s1, s2)):
        ax  = axes[row, col]
        cc  = C[C.subgenome == sub].sort_values("allo_chr", key=lambda s: s.map(nat))
        chroms = list(cc.allo_chr)
        d   = E[E.subgenome == sub]
        colg = style.GENOME[GEN[sub]]
        for i, c in enumerate(chroms):
            ax.add_patch(Rectangle((0, i-0.24), 1, 0.48, facecolor="#EAEAEA",
                                   edgecolor="none", zorder=1))
        placed = {}
        for _, r in d.sort_values("genes", ascending=False).iterrows():
            if r.allo_chr not in chroms: continue
            i = chroms.index(r.allo_chr)
            w = max(r.end_frac - r.start_frac, MINW)
            x = min(r.start_frac, 1 - w)
            conf = r.first_gene in both
            ax.add_patch(Rectangle((x, i-0.24), w, 0.48,
                                   facecolor=colg if conf else "white",
                                   edgecolor=colg, linewidth=0.9, zorder=3))
            if r.genes >= LABEL_MIN and all(abs(x - p) > 0.11 for p in placed.get(i, [])):
                ax.text(min(x + w/2, 0.985), i - 0.30, disp(sub, r.exchange_chr),
                        ha="center", va="bottom", fontsize=5.8, color=style.ANNOT, zorder=4)
                placed.setdefault(i, []).append(x)
        ax.set_yticks(range(len(chroms)))
        ax.set_yticklabels([disp(sub, c) for c in chroms], fontsize=6)
        ax.set_ylim(len(chroms) - 0.45, -0.95); ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"], fontsize=6)
        ax.set_title("%s   %d exchanges, %d genes" % (sub, len(d), int(d.genes.sum())), fontsize=7.5)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        if row == 2: ax.set_xlabel("position along the chromosome", fontsize=7)
        if col == 0: ax.set_ylabel(style.ital(sp), fontsize=8)

h = [Rectangle((0,0),1,1, facecolor=style.GREY, edgecolor=style.GREY, label="recovered from both block sets"),
     Rectangle((0,0),1,1, facecolor="white",   edgecolor=style.GREY, label="recovered from one block set")]
fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=6.6,
           bbox_to_anchor=(0.5, -0.022))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "FigS9_exchanges.%s" % ext), bbox_inches="tight")
print("wrote FigS9_exchanges")
