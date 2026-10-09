"""Fig. 4 - device placement of the societal and the utility optimum in the IEEE 33-bus system.
Run:  python make_placement.py        (needs ../evaluator and paper_data/front_ieee33.csv)
Output: figures/fig_placement.eps, .pdf, .png

Device symbols are those of the Bawean diagram (Fig. 2): white bar = manual switch, black bar =
remote-controlled switch, hexagon = recloser, dashed line with an ellipse = tie switch; all bars are drawn
across their line.  Tie lines run along straight segments with right-angle corners."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "evaluator"))
from systems import System
import factor_guard; factor_guard.check()

D, O = os.path.join(HERE, "paper_data"), os.path.join(HERE, "figures")
os.makedirs(O, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
                     "font.size": 8, "legend.fontsize": 8, "savefig.dpi": 600, "ps.fonttype": 42, "pdf.fonttype": 42})
MM, PT = 1 / 25.4, 72 / 25.4                      # inch per mm, point per mm

# symbol sizes in mm, identical to Fig. 2
BAR_L, BAR_T, BAR_S, HEX_R, ELL_A, ELL_B, ELL_S = 3.6, 1.2, 0.25, 1.7, 2.1, 1.05, 0.3
r = BAR_T / BAR_L
VBAR = Path([(-r, -1), (r, -1), (r, 1), (-r, 1), (-r, -1)], closed=True)     # bar across a horizontal line
HBAR = Path([(-1, -r), (1, -r), (1, r), (-1, r), (-1, -r)], closed=True)     # bar across a vertical line
c = Path.unit_circle()
VELL = Path(c.vertices * [ELL_B / ELL_A, 1.0], c.codes)                      # ellipse across a horizontal line
HELL = Path(c.vertices * [1.0, ELL_B / ELL_A], c.codes)
def device_style(kind, horizontal_line):
    if kind == 3:
        return dict(marker="h", ms=2 * HEX_R * PT, mfc="black", mec="black", mew=0)
    return dict(marker=VBAR if horizontal_line else HBAR, ms=BAR_L * PT, mfc="white" if kind == 1 else "black",
                mec="black", mew=BAR_S * PT)
TIE_STYLE = dict(ms=2 * ELL_A * PT, mfc="white", mec="black", mew=ELL_S * PT)

# bus positions: main feeder on y = 0, laterals on the rows above and below
P = {i: (i - 1.0, 0.0) for i in range(1, 19)}
P.update({i: (i - 18.0, -1.1) for i in range(19, 23)})
P.update({i: (i - 21.0, 1.1) for i in range(23, 26)})
P.update({i: (i - 21.0, 1.1) for i in range(26, 34)})
ROUTE = {(8, 21): -0.55, (12, 22): -0.78, (9, 15): -0.38, (18, 33): 0.55, (25, 29): 0.78}   # level of each tie line
LABELLED = (1, 6, 9, 12, 15, 18, 19, 22, 23, 25, 26, 29, 33)

S = System("ieee33")
df = pd.read_csv(os.path.join(D, "front_ieee33.csv"))
F = df[["f1", "f2", "f3"]].to_numpy(); G = df[[c for c in df.columns if c.startswith("g")]].to_numpy().astype(int)
plans = {"(a)": int(np.argmin(F[:, 0] + F[:, 1])), "(b)": int(np.argmin(F[:, 0] + F[:, 2]))}

def midpoint(pts):
    seg = [np.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts[:-1], pts[1:])]
    half, acc = sum(seg) / 2, 0.0
    for (a, b), L in zip(zip(pts[:-1], pts[1:]), seg):
        if acc + L >= half:
            t = (half - acc) / L
            return a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), abs(a[1] - b[1]) < 1e-9
        acc += L

fig, axs = plt.subplots(2, 1, figsize=(174 * MM, 80 * MM))
used = set()                                                             # symbols that occur, for the legend
for ax, (title, k) in zip(axs, plans.items()):
    x, w = S.decode(G[k])
    used |= {("dev", v) for v in x if v} | {("tie", v) for v in w}
    for c, t in enumerate(S.net.ties):                                   # ties first, so that the feeder is on top
        a, b = t["a"] + 1, t["b"] + 1
        yr = ROUTE[(a, b)] if (a, b) in ROUTE else ROUTE[(b, a)]
        if (a, b) not in ROUTE: a, b = b, a
        pts = [P[a], (P[a][0], yr), (P[b][0], yr), P[b]]
        if w[c]:
            ax.plot(*zip(*pts), color="black", lw=1.0, ls="--", zorder=2)
            mx, my, horizontal = midpoint(pts)
            ax.plot([mx], [my], ls="", marker=VELL if horizontal else HELL, zorder=5,
                    **dict(TIE_STYLE, mfc="white" if w[c] == 1 else "black"))
        else:
            ax.plot(*zip(*pts), color="0.6", lw=0.7, ls=":", zorder=1)
    for s, sec in enumerate(S.net.sec):
        a, b = P[sec["frm"] + 1], P[sec["to"] + 1]
        pts = [a, b] if a[0] == b[0] or a[1] == b[1] else [a, (a[0], b[1]), b]
        ax.plot(*zip(*pts), color="black", lw=1.0, zorder=3, solid_capstyle="butt")
        if x[s]:
            q = pts[1]
            ax.plot([a[0] + 0.45 * (q[0] - a[0])], [a[1] + 0.45 * (q[1] - a[1])], ls="", zorder=6,
                    **device_style(x[s], horizontal_line=abs(a[1] - q[1]) < 1e-9))
    ax.plot(*zip(*P.values()), ls="", marker="o", ms=2.6, color="black", zorder=4)
    for i in LABELLED:
        px, py = P[i]
        dy = 0.22 if py > 0 else (-0.42 if py < 0 else (-0.40 if i == 6 else 0.16))
        ax.text(px, py + dy, str(i), ha="center", va="bottom")
    ax.set_title(title, loc="left", fontsize=8)
    ax.set_xlim(-0.6, 17.6); ax.set_ylim(-1.75, 1.6); ax.axis("off")
    o = S.evaluate(G[k])
    print(title, S.describe(G[k]), f"f1 {o['f1'] / 1e3:.1f}  f2 {o['f2'] / 1e3:.1f}  f3 {o['f3'] / 1e3:.1f} kUSD/yr")
entries = [(("dev", 1), Line2D([], [], ls="", label="MS", **device_style(1, True))),
           (("dev", 2), Line2D([], [], ls="", label="RCS", **device_style(2, True))),
           (("dev", 3), Line2D([], [], ls="", label="Recloser", **device_style(3, True))),
           (("tie", 1), Line2D([], [], color="black", lw=1.0, ls="--", marker=VELL, label="Tie built (manual)", **TIE_STYLE)),
           (("tie", 2), Line2D([], [], color="black", lw=1.0, ls="--", marker=VELL, label="Tie built (remote-controlled)",
                               **dict(TIE_STYLE, mfc="black"))),
           (("tie", 0), Line2D([], [], color="0.6", lw=0.7, ls=":", label="Candidate tie not built"))]
handles = [h for key, h in entries if key in used]
fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, handlelength=3.4, columnspacing=2.5)
fig.tight_layout(rect=(0, 0.07, 1, 1))
for ext in ("eps", "pdf", "png"):
    fig.savefig(os.path.join(O, f"fig_placement.{ext}"), bbox_inches="tight", dpi=600 if ext != "png" else 200)
print("selesai: figures/fig_placement.*")
