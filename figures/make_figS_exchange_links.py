#!/usr/bin/env python3
import os, sys, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
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
PROG = {"An":"rapa","Aj":"rapa","Bj":"nigra","Bc":"nigra","Cn":"oleracea","Cc":"oleracea"}
SPEC = [("napus","An","Cn"), ("juncea","Aj","Bj"), ("carinata","Bc","Cc")]
GAP, TOP, BOT, H = 0.14, 1.0, 0.0, 0.085

def disp(sub, c):
    c = str(c)
    if sub == "Cn":
        m = re.match(r'^C1(\d)$', c)
        if m: return "C%s" % m.group(1)
    m = re.match(r'^([A-Za-z]+)0*(\d+)$', c)
    return m.group(1) + m.group(2) if m else c

def nat(c):
    m = re.match(r'^([A-Za-z]*)(\d+)$', str(c))
    return (m.group(1), int(m.group(2))) if m else (str(c), 0)

def layout(names):
    n = len(names); w = (1.0 - GAP) / n
    return {c: (i * (w + GAP / max(n - 1, 1)), w) for i, c in enumerate(names)}

def ribbon(ax, x0, x1, y0, y1, w0, w1, color, alpha, lw=0.0):
    my = (y0 + y1) / 2
    verts = [(x0 - w0/2, y0), (x0 - w0/2, my), (x1 - w1/2, my), (x1 - w1/2, y1),
             (x1 + w1/2, y1), (x1 + w1/2, my), (x0 + w0/2, my), (x0 + w0/2, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    ax.add_patch(PathPatch(Path(verts + [(x0 - w0/2, y0)], codes + [Path.CLOSEPOLY]),
                           facecolor=color, edgecolor="none", alpha=alpha, zorder=3, lw=lw))

fig, axes = plt.subplots(3, 2, figsize=(7.3, 7.0), layout="constrained")
fig.get_layout_engine().set(w_pad=0.07, h_pad=0.07, wspace=0.08, hspace=0.10)

for row, (sp, s1, s2) in enumerate(SPEC):
    for col, sub in enumerate((s1, s2)):
        ax  = axes[row, col]
        cc  = C[C.subgenome == sub].sort_values("allo_chr", key=lambda s: s.map(nat))
        allo = list(cc.allo_chr)
        prog = sorted(set(cc.partner_chr) | set(E[E.subgenome == sub].exchange_chr), key=nat)
        LA, LP = layout(allo), layout(prog)
        colg = style.GENOME[GEN[sub]]
        for c, (x, w) in LP.items():
            ax.add_patch(Rectangle((x, TOP), w, H, facecolor="#C9C9C9", edgecolor="none", zorder=4))
            ax.text(x + w/2, TOP + H + 0.035, disp(sub, c), ha="center", va="bottom", fontsize=5.2)
        for c, (x, w) in LA.items():
            ax.add_patch(Rectangle((x, BOT - H), w, H, facecolor="#C9C9C9", edgecolor="none", zorder=4))
            ax.text(x + w/2, BOT - H - 0.035, disp(sub, c), ha="center", va="top", fontsize=5.2)
        for _, r in cc.iterrows():
            if r.partner_chr not in LP: continue
            xp, wp = LP[r.partner_chr]; xa, wa = LA[r.allo_chr]
            ribbon(ax, xp + wp/2, xa + wa/2, TOP, BOT, wp*0.80, wa*0.80, "#E2E2E2", 0.85)
        for _, r in E[E.subgenome == sub].sort_values("genes").iterrows():
            if r.exchange_chr not in LP or r.allo_chr not in LA: continue
            xp, wp = LP[r.exchange_chr]; xa, wa = LA[r.allo_chr]
            x0 = xp + wp * (r.donor_start_frac + r.donor_end_frac) / 2
            x1 = xa + wa * (r.start_frac + r.end_frac) / 2
            tw = min(0.004 + r.genes / 1600.0, 0.030)
            ribbon(ax, x0, x1, TOP, BOT, tw, tw, colg, 0.95 if r.first_gene in both else 0.45)
        ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.30, 1.30); ax.axis("off")
        ax.set_title("%s   %s above, %s below" % (sub, style.ital(PROG[sub]), sub), fontsize=7.3)

fig.text(0.5, -0.012, "ribbon width is the number of genes; a solid ribbon was recovered from "
         "both block sets, a pale one from one", ha="center", fontsize=6.3, color=style.ANNOT)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "FigS10_exchange_links.%s" % ext), bbox_inches="tight")
print("wrote FigS10_exchange_links")
