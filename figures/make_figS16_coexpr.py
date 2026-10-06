import numpy as np, pandas as pd, matplotlib.pyplot as plt
import style
style.use()
np.random.seed(0)

RB = "../../expression_rebuild"
df = pd.read_csv(f"{RB}/01_counts/gref_counts.tsv", sep="\t")
libs = [c for c in df.columns if "|" in c]
tracks = ["Bra_A", "Bni_B", "Bol_C", "Bna_A", "Bna_C", "Bju_A", "Bju_B", "Bca_B", "Bca_C"]
LABEL = {"Bra_A": "rapa A", "Bni_B": "nigra B", "Bol_C": "oleracea C",
         "Bna_A": "napus A", "Bna_C": "napus C", "Bju_A": "juncea A",
         "Bju_B": "juncea B", "Bca_B": "carinata B", "Bca_C": "carinata C"}

M = df[libs].to_numpy(float)
M = M / np.nansum(M, 0, keepdims=True) * 1e6
LD = pd.DataFrame(np.log2(M + 1), columns=libs, index=df["AT_gene"].values)

def tissue_of(lib):
    s = lib.split("|")[1]
    if "Pol" in s: return "pollen"
    if "E_Stig" in s: return "stigma (early)"
    if "L_Stig" in s: return "stigma (late)"
    return "stigma"

res, rdist = {}, {}
for tr in tracks:
    cols = [c for c in libs if c.startswith(tr + "|")]
    s = LD[cols]
    s = s.loc[s.notna().all(axis=1) & (s > 1).all(axis=1)]
    X = s.to_numpy()[np.random.choice(len(s), min(2500, len(s)), replace=False)]
    Xc = X - X.mean(1, keepdims=True)
    sd = Xc.std(1, keepdims=True); sd[sd == 0] = np.nan
    Z = Xc / sd
    R = (Z @ Z.T) / X.shape[1]
    r = R[np.triu_indices(len(R), 1)]; r = r[np.isfinite(r)]
    ts = np.array([tissue_of(c) for c in cols]); uniq = sorted(set(ts))
    grand = X.mean(1, keepdims=True)
    between = sum((ts == u).sum() * (X[:, ts == u].mean(1, keepdims=True) - grand) ** 2
                  for u in uniq).ravel()
    total = ((X - grand) ** 2).sum(1)
    frac = np.nanmedian(np.where(total > 0, between / np.where(total > 0, total, np.nan), np.nan))
    rdist[tr] = r
    res[tr] = dict(n_libs=len(cols), n_tissues=len(uniq), genes=len(s),
                   frac_gt9=float((np.abs(r) > 0.9).mean()),
                   med=float(np.median(np.abs(r))), tissue_frac=float(frac))

out = pd.DataFrame(res).T
out.index.name = "track"
out.to_csv(f"{RB}/07_coexpr/identifiability.csv")

fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.9), layout="constrained")
fig.get_layout_engine().set(w_pad=0.06, wspace=0.07)

ax = axes[0]
grid = np.linspace(0, 1, 201)
for tr in tracks:
    a = np.sort(np.abs(rdist[tr]))
    surv = 1.0 - np.searchsorted(a, grid, side="left") / len(a)
    ax.plot(grid, 100 * surv, lw=1.0, color=style.GENOME[style.SUB[tr]], alpha=0.85)
ax.axvline(0.9, color="black", lw=0.9, ls="--")
ax.annotate("conventional\nedge threshold", xy=(0.9, 45), xytext=(0.30, 22),
            fontsize=6.5, color=style.ANNOT, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="black",
                            connectionstyle="arc3,rad=-0.15"))
ax.set_xlim(0, 1); ax.set_ylim(0, 100)
ax.set_xlabel("|r| threshold")
ax.set_ylabel("random pairs retained (%)")
ax.set_title("Random pairs\nalready correlate", fontsize=8)

ax = axes[1]
y = np.arange(len(tracks))
vals = [100 * res[t]["frac_gt9"] for t in tracks]
ax.barh(y, vals, color=[style.GENOME[style.SUB[t]] for t in tracks], height=0.68)
ax.set_yticks(y); ax.set_yticklabels([style.ital(LABEL[t]) for t in tracks])
ax.invert_yaxis()
for yi, v in zip(y, vals):
    ax.text(v + 2, yi, "%.0f%%" % v, va="center", fontsize=6.5, color=style.ANNOT)
ax.set_xlim(0, 108); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("retained at |r| > 0.9 (%)")
ax.set_title("The threshold excludes\nalmost nothing", fontsize=8)

ax = axes[2]
vals = [100 * res[t]["tissue_frac"] for t in tracks]
ax.barh(y, vals, color=[style.GENOME[style.SUB[t]] for t in tracks], height=0.68)
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis()
for yi, v in zip(y, vals):
    ax.text(v - 2, yi, "%.0f%%" % v, va="center", ha="right", fontsize=6.5, color="white")
ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("variance that is tissue (%)")
ax.set_title("One axis of variation,\nnot many", fontsize=8)

fig.canvas.draw()
for ax, letter in zip(axes, "abc"):
    p = ax.get_position()
    fig.text(p.x0 - 0.052, p.y1 + 0.055, letter, fontsize=10, fontweight="bold",
             va="top", ha="left")

for ext in ("png", "pdf"):
    fig.savefig("FigS16_coexpr_identifiability." + ext, bbox_inches=None)
print(out.to_string())
print("\nwrote FigS16_coexpr_identifiability.png / .pdf")
