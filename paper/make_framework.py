"""Fig. 1 - planning framework (flow of one evaluation and of the whole study).
Run:  python make_framework.py
Output: figures/fig_framework.eps, .pdf, .png   (174 mm x 104 mm, all coordinates below are in mm)"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(HERE, "figures"); os.makedirs(O, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
                     "font.size": 8, "savefig.dpi": 600, "ps.fonttype": 42, "pdf.fonttype": 42})
MM = 1 / 25.4
W, H = 174.0, 104.0
fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

def box(x, y, w, h, text, fill="white", dashed=False, title=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5", lw=0.8, ec="black",
                                fc="none" if dashed else fill, ls="--" if dashed else "-"))
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center")
    if title:
        ax.text(x + w / 2, y + h - 3.25, title, ha="center", va="center", fontweight="bold")

def arrow(pts, label=None, lpos=None):
    """Polyline with an arrowhead at its last point; an optional label is set on a white patch over the line."""
    for a, b in zip(pts[:-2], pts[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color="black", lw=0.8, solid_capstyle="butt")
    ax.annotate("", xy=pts[-1], xytext=pts[-2], arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=9))
    if label:
        ax.text(*lpos, label, ha="center", va="center", bbox=dict(facecolor="white", edgecolor="white"))

# inputs
box(58, 86, 36, 15, "Network data\ntopology, impedances,\ncustomers, sectors")
box(97, 86, 36, 15, "Reliability and cost data\nfailure rates, repair and\nswitching times, costs")
box(136, 86, 36, 15, "Regulation\nthreshold, multipliers,\ntariffs, minimum bill")
# evaluator
box(56, 3, 118, 78, None, dashed=True, title="Exact evaluator, applied to every plan")
box(60, 61, 110, 12, "1  Fault response for every section $k$ and load level $\\omega$\n(failure mode and effect analysis)")
box(60, 45, 110, 12, "2  Two-stage restoration: remote devices at $\\tau^{\\mathrm{RCS}}$, all devices at $\\tau^{\\mathrm{MS}}$;\n"
                     "customers restored in stage 1 are not interrupted again")
box(60, 27, 110, 13, "3  AC power flow at four load levels (voltage, current and source limits);\n"
                     "exact search over tie and block combinations, maximum restored demand")
box(60, 7, 52, 15, "4a  Interruption durations $r^{\\omega}_{jk}$\n$\\rightarrow$ ECOST $f_2$, SAIFI, SAIDI, ENS")
box(118, 7, 52, 15, "4b  Panjer recursion $\\rightarrow\\ \\mathbb{E}[m(D_j)]$\n$\\rightarrow$ compensation $f_3$")
# search and analyses
box(3, 56, 42, 20, "NSGA-II\npopulation of plans\n$\\mathbf{y}=(\\mathbf{x},\\mathbf{w},\\mathbf{z})$\n($f_1$ computed directly)", fill="0.9")
box(3, 31, 42, 16, "Pareto set\nannual cost $f_1$, ECOST $f_2$,\ncompensation $f_3$")
box(3, 3, 42, 20, "Analyses\nrepresentative plans, budget\ncomparison, load-level study,\nsensitivity")
# flows (each arrow stops 0.3 mm before the box it enters)
ax.plot([56, 50.5], [14.5, 14.5], color="black", lw=0.8, solid_capstyle="butt")
ax.plot([50.5, 50.5], [14.5, 62], color="black", lw=0.8, solid_capstyle="butt")
for x in (76, 115, 154):
    arrow([(x, 86), (x, 81.3)])
arrow([(115, 61), (115, 57.3)]); arrow([(115, 45), (115, 40.3)])
arrow([(86, 27), (86, 22.3)]); arrow([(144, 27), (144, 22.3)])
arrow([(45, 70), (59.7, 70)], "plan", (52, 72.5))
arrow([(50.5, 62), (45.3, 62)], "$f_2$, $f_3$", (50.5, 38))
arrow([(24, 56), (24, 47.3)]); arrow([(24, 31), (24, 23.3)])
for ext in ("eps", "pdf", "png"):
    fig.savefig(os.path.join(O, f"fig_framework.{ext}"), dpi=600 if ext != "png" else 200)
print("selesai: figures/fig_framework.*")
