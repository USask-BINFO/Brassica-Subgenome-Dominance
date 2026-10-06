#!/usr/bin/env python3
import os, sys, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); import style; style.use()
RB = os.path.join(HERE, "..", "..", "expression_rebuild")

SAMPLES = ["napus_pollen", "napus_stigma", "carinata_pollen", "carinata_stigmaE",
           "carinata_stigmaL", "juncea_pollen", "juncea_stigma"]
PRETTY = {"napus_pollen": "napus pollen", "napus_stigma": "napus stigma",
          "carinata_pollen": "carinata pollen", "carinata_stigmaE": "carinata stigma (early)",
          "carinata_stigmaL": "carinata stigma (late)", "juncea_pollen": "juncea pollen",
          "juncea_stigma": "juncea stigma"}
SUBS = {"napus": ("A", "C"), "juncea": ("A", "B"), "carinata": ("B", "C")}

def letters(fig, axes, dx=0.055, dy=0.05):
    fig.canvas.draw()
    for ax, L in zip(axes, "abcdef"):
        p = ax.get_position()
        fig.text(p.x0 - dx, p.y1 + dy, L, fontsize=10, fontweight="bold", va="top", ha="left")

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches=None)
    print("wrote", name)

ELD = os.path.join(RB, "03_eld", "theta1_a05")
summ = pd.read_csv(os.path.join(ELD, "eld_summary.csv"))
cls = {}
for s in SAMPLES:
    d = pd.read_csv(os.path.join(ELD, f"eld_{s}.csv"))
    cls[s] = d["class"].value_counts()
CL = pd.DataFrame(cls).fillna(0).T[SAMPLES[0:0] or None] if False else pd.DataFrame(cls).fillna(0).T
CL = CL.loc[SAMPLES]
keep = [c for c in CL.columns if c != "unclassified"]
def rank(c):
    return (0 if c.startswith("PED") else 1,
            0 if "ELD_P1" in c else 1 if "ELD_P2" in c else 2 if "additive" in c
            else 3 if "unchanged" in c else 4)
keep = sorted(keep, key=rank)
PCOL = plt.get_cmap("tab20")(np.linspace(0, 1, 20))

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), layout="constrained")
fig.get_layout_engine().set(w_pad=0.08, wspace=0.12)

ax = axes[0]
bot = np.zeros(len(SAMPLES))
for i, c in enumerate(keep):
    v = 100 * CL[c].to_numpy() / CL[keep].sum(axis=1).to_numpy()
    ax.bar(range(len(SAMPLES)), v, bottom=bot, width=0.7, color=PCOL[i % 20],
           edgecolor="white", linewidth=0.4, label=c.replace("_", " "))
    bot += v
ax.set_xticks(range(len(SAMPLES)))
ax.set_xticklabels([style.ital(PRETTY[s]) for s in SAMPLES], rotation=45, ha="right")
ax.set_ylabel("classified pairs (%)")
ax.set_title("Inheritance patterns", fontsize=8)
ax.legend(fontsize=5.6, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5),
          handlelength=1.1, labelspacing=0.3)

ax = axes[1]
x = np.arange(len(SAMPLES)); w = 0.36
sm = summ.set_index("sample").loc[SAMPLES]
for j, col in enumerate(["ELD_P1_pct", "ELD_P2_pct"]):
    gl = [SUBS[s.split("_")[0]][j] for s in SAMPLES]
    pos = x + (j - 0.5) * w
    ax.bar(pos, sm[col], width=w, color=[style.GENOME[g] for g in gl],
           edgecolor="white", linewidth=0.4)
    for xi, v, g in zip(pos, sm[col], gl):
        ax.text(xi, v + 0.8, g, ha="center", va="bottom", fontsize=6.5, color=style.ANNOT)
ax.set_xticks(x)
ax.set_xticklabels([style.ital(PRETTY[s]) for s in SAMPLES], rotation=45, ha="right")
ax.set_ylim(0, 58)
ax.set_ylabel("share of PED pairs (%)")
ax.set_title("Which progenitor the allotetraploid resembles", fontsize=8)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=style.GENOME[g]) for g in "ABC"],
          labels=["A", "B", "C"], title="genome", fontsize=6.5, title_fontsize=6.5,
          loc="upper left", ncol=3, handlelength=1.1, columnspacing=0.9)
letters(fig, axes)
save(fig, "FigS11_eld")

H = pd.read_csv(os.path.join(RB, "05_validation", "heb_by_br_subgenome.csv"))
H["sample"] = H.species + "_" + H.tissue
LAY = ["LF", "MF1", "MF2"]
fig, ax = plt.subplots(figsize=(6.8, 3.1), layout="constrained")
x = np.arange(len(SAMPLES)); w = 0.26
MK = {"LF": "o", "MF1": "s", "MF2": "^"}
for j, lay in enumerate(LAY):
    for i, s in enumerate(SAMPLES):
        r = H[(H["sample"] == s) & (H.layer == lay)]
        if r.empty:
            continue
        r = r.iloc[0]
        c = style.GENOME[r.sub1] if r.ratio > 1 else style.GENOME[r.sub2]
        ax.plot(i + (j - 1) * w, r.ratio, MK[lay], ms=5.5, color=c,
                mec="black" if r.p < 0.05 else "none", mew=0.8, zorder=4)
ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=1)
ax.set_xticks(x)
ax.set_xticklabels([style.ital(PRETTY[s]) for s in SAMPLES], rotation=45, ha="right")
ax.set_ylabel("dominance ratio within the layer")
ax.set_title("The contest inside each Br-subgenome layer", fontsize=8)
leg1 = ax.legend(handles=[Line2D([], [], marker=MK[l], ls="none", ms=5.5, color=style.GREY,
                                 label=l) for l in LAY] +
                         [Line2D([], [], marker="o", ls="none", ms=5.5, mfc=style.GREY,
                                 mec="black", mew=0.8, label="p < 0.05")],
                 fontsize=6.5, loc="upper left", ncol=2, handlelength=1.2,
                 title="Br-subgenome layer", title_fontsize=6.5)
ax.add_artist(leg1)
ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.5, color=style.GENOME[g],
                          label=g) for g in "ABC"],
          fontsize=6.5, loc="upper right", ncol=3, handlelength=1.0, columnspacing=0.8,
          title="leans toward", title_fontsize=6.5)
save(fig, "FigS15_brsubgenome")

E = pd.read_csv(os.path.join(RB, "05_validation", "enrichment_novel_switched.csv"))
E = E[E.fdr < 0.05]
COLS = ["An", "Aj", "Bj", "Bc", "Cn", "Cc"]

def wrap(t, n=44):
    if len(t) <= n:
        return t
    cut = t[:n].rsplit(" ", 1)[0]
    return cut + "\u2026"

tissues = [t for t in ["pollen", "stigma"] if t in set(E.tissue)]
mats = []
for tis in tissues:
    sub = E[E.tissue == tis]
    piv = -np.log10(sub.pivot_table(index="term_name", columns="col",
                                    values="fdr", aggfunc="min"))
    piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).index[:28]]
    for c in COLS:
        if c not in piv.columns:
            piv[c] = np.nan
    mats.append((tis, piv[COLS]))

vmax = max(np.nanmax(m.to_numpy()) for _, m in mats)
nrows = [len(m) for _, m in mats]
fig, axes = plt.subplots(len(mats), 1, figsize=(5.6, 0.175 * sum(nrows) + 1.0),
                         layout="constrained", height_ratios=nrows, squeeze=False)
axes = [a[0] for a in axes]
for ax, (tis, m) in zip(axes, mats):
    im = ax.imshow(m.to_numpy(), aspect="auto", cmap="viridis", vmin=1.3, vmax=vmax)
    ax.set_xticks(range(len(COLS))); ax.set_xticklabels(COLS, fontsize=7)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels([wrap(t) for t in m.index], fontsize=5.8)
    ax.set_title(tis, fontsize=8, loc="left")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
cb.set_label(r"$-\log_{10}$ FDR", fontsize=7)
cb.ax.tick_params(labelsize=6)
save(fig, "FigS20_enrichment_full")
