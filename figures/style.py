import matplotlib as mpl
import matplotlib.pyplot as plt

A = "#0072B2"
B = "#D55E00"
C = "#009E73"
GREY = "#555555"; LIGHT = "#BBBBBB"
ANNOT = "black"
SPNAMES = ("rapa","nigra","oleracea","napus","juncea","carinata")
GENOME = {"A": A, "B": B, "C": C}
SPECIES = {"rapa": A, "nigra": B, "oleracea": C,
           "napus": "#56189C", "juncea": "#B8860B", "carinata": "#CC79A7"}
SUB = {"Bra_A":"A","Bni_B":"B","Bol_C":"C","Bna_A":"A","Bna_C":"C",
       "Bju_A":"A","Bju_B":"B","Bca_B":"B","Bca_C":"C"}
TISSUE_MARKER = {"pollen":"o","stigma":"s","stigmaE":"^","stigmaL":"D"}

def use():
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 400,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans","Helvetica","Arial"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "legend.frameon": False,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def panel(ax, letter, dx=-0.17, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")

def tidy(ax):
    ax.tick_params(direction="out")
    return ax

def ital(txt):
    out=txt
    for sp in SPNAMES:
        if sp in out:
            out=out.replace(sp, r"$\mathit{%s}$" % sp)
    return out
