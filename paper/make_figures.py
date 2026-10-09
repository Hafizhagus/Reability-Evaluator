"""Gambar final paper (dapat dijalankan ulang): python make_figures.py
Data: paper_data/ (front E2, konvergensi E2, E3 gabungan, E4). Keluaran: figures/*.eps, *.pdf, *.png
Gaya: huruf Arial/Helvetica (cadangan Liberation Sans) 8 pt, terbaca hitam-putih."""
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D, O = "paper_data", "figures"
os.makedirs(O, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
                     "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 8, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "savefig.dpi": 600, "ps.fonttype": 42, "pdf.fonttype": 42})
MM = 1 / 25.4
SYS = [("ieee33", "IEEE 33-bus"), ("ieee69", "IEEE 69-bus"), ("rbts2", "RBTS Bus 2"), ("bawean", "Bawean")]
STY = {"ieee33": dict(marker="o", ls="-"), "ieee69": dict(marker="s", ls="--"), "rbts2": dict(marker="^", ls=":"),
       "bawean": dict(marker="v", ls="-.")}
MS_ANN = (0.08 * 1.08 ** 15 / (1.08 ** 15 - 1) + 0.046) * 22000.0      # annualised MS cost, $/yr
BUDGETS = [3, 5, 8, 10, 12]

def save(fig, name, pad=0.1):
    """pad = white margin (inch) around the drawing; the saved width must stay within 174 mm."""
    for ext in ("eps", "pdf", "png"):
        fig.savefig(os.path.join(O, f"{name}.{ext}"), bbox_inches="tight", pad_inches=pad,
                    dpi=600 if ext != "png" else 200)
    plt.close(fig)

def front(s):
    df = pd.read_csv(os.path.join(D, f"front_{s}.csv"))
    return df[["f1", "f2", "f3"]].to_numpy() / 1e3                      # kUSD/yr

def reps(F):
    span = F.max(0) - F.min(0); nrm = (F - F.min(0)) / np.where(span > 0, span, 1)
    return {"social": int(np.argmin(F[:, 0] + F[:, 1])), "utility": int(np.argmin(F[:, 0] + F[:, 2])),
            "compromise": int(np.argmin(np.linalg.norm(nrm, axis=1)))}

MK = {"social": ("^", "Societal optimum (min $f_1+f_2$)"), "utility": ("s", "Utility optimum (min $f_1+f_3$)"),
      "compromise": ("D", "Compromise")}

# ---------------------------------------------------------------- Fig. Pareto fronts
fig, ax = plt.subplots(2, 4, figsize=(174 * MM, 90 * MM))
for c, (s, title) in enumerate(SYS):
    F = front(s); R = reps(F)
    for r, (j, lab) in enumerate(((1, "ECOST $f_2$ (kUSD/yr)"), (2, "Compensation $f_3$ (kUSD/yr)"))):
        a = ax[r, c]
        a.scatter(F[:, 0], F[:, j], s=3, c="0.55", lw=0, label="Union Pareto front")
        for k, (m, name) in MK.items():
            a.scatter(F[R[k], 0], F[R[k], j], s=28, marker=m, facecolor="white", edgecolor="black", lw=0.9, label=name, zorder=3)
        a.set_xlabel("Annual cost $f_1$ (kUSD/yr)"); a.set_ylabel(lab)
        if r == 0: a.set_title(f"({'abcd'[c]}) {title}", loc="left", fontsize=8)
        a.grid(lw=0.3, color="0.85")
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(rect=(0, 0.08, 1, 1)); save(fig, "fig_pareto", pad=0.04)

# ---------------------------------------------------------------- Fig. E3 objective comparison vs budget
def e3_union(s):
    rows = json.load(open(os.path.join(D, f"{s}_t60_h1.json")))
    F = pd.read_csv(os.path.join(D, f"front_{s}.csv"))[["f1", "f2", "f3"]].to_numpy()
    out = []
    for n in BUDGETS:
        B = n * MS_ANN + 1e-6
        cand = [(r["f2"], r["f3"], r["share_above"]) for r in rows if float(r["budget"]) == n]
        cand += [(F[i, 1], F[i, 2], np.nan) for i in np.where(F[:, 0] <= B)[0]]
        xE = min(cand, key=lambda c: c[0]); xK = min(cand, key=lambda c: c[1])
        out.append(((xE[1] - xK[1]) / xK[1], (xK[0] - xE[0]) / xE[0]))
    return np.array(out)

fig, ax = plt.subplots(1, 2, figsize=(174 * MM, 58 * MM))
for s, title in SYS:
    U = 100 * e3_union(s)
    ax[0].plot(BUDGETS, U[:, 0], color="black", label=title, **STY[s], ms=4, mfc="white")
    ax[1].plot(BUDGETS, U[:, 1], color="black", label=title, **STY[s], ms=4, mfc="white")
ax[0].set_ylabel("Compensation regret (%)")
ax[1].set_ylabel("Societal price (%)")
for a, t in zip(ax, ("(a)", "(b)")):
    a.set_xlabel("Budget (equivalent number of manual switches)"); a.set_xticks(BUDGETS)
    a.grid(lw=0.3, color="0.85"); a.set_title(t, loc="left", fontsize=8)
ax[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.06, 1.0))
fig.tight_layout(); save(fig, "fig_regret_budget")

# ---------------------------------------------------------------- Fig. E4 load level: valuation error vs decision regret
A = json.load(open(os.path.join(D, "e4_analysis.json")))
TR = [10, 25, 50, 100, 200]
fig, ax = plt.subplots(1, 2, figsize=(174 * MM, 58 * MM))
for lm, st, lab in (("peak", dict(marker="o", ls="-"), "Peak snapshot"), ("avg", dict(marker="s", ls="--"), "Average snapshot")):
    v = [100 * (A[f"fixed|{t}"][lm] - A[f"fixed|{t}"]["ldc4"]) / A[f"fixed|{t}"]["ldc4"] for t in TR]
    ax[0].plot(TR, v, color="black", label=lab, ms=4, mfc="white", **st)
for ob, st, lab in (("social", dict(marker="o", ls="-"), "Societal cost $f_1+f_2$"),
                    ("becost:5", dict(marker="^", ls="--"), "ECOST, budget of 5 manual switches")):
    v = [100 * A[f"{t}|{ob}|peak"]["regret"] for t in TR]
    ax[1].plot(TR, v, color="black", label=lab, ms=4, mfc="white", **st)
ax[0].axhline(0, color="0.5", lw=0.5); ax[1].set_ylim(bottom=0)
ax[0].set_ylabel("ECOST valuation error (%)")
ax[1].set_ylabel("Decision regret, peak snapshot (%)")
for a, t in zip(ax, ("(a)", "(b)")):
    a.set_xscale("log"); a.set_xticks(TR); a.set_xticklabels([str(t) for t in TR])
    a.set_xlabel("Distribution transformer restoration time (h)")
    a.grid(lw=0.3, color="0.85"); a.legend(frameon=False, loc="best"); a.set_title(t, loc="left", fontsize=8)
fig.tight_layout(); save(fig, "fig_loadlevel")

# ---------------------------------------------------------------- Fig. E6 convergence
fig, ax = plt.subplots(1, 1, figsize=(84 * MM, 60 * MM))
for s, title in SYS:
    C = np.loadtxt(os.path.join(D, f"conv_{s}.csv"), delimiter=",", ndmin=2)
    C = C / C[:, -1:].mean()                                           # relative to the mean final HV
    g = np.arange(1, C.shape[1] + 1)
    ax.plot(g, C.mean(0), color="black", ls=STY[s]["ls"], label=title)
    ax.fill_between(g, C.min(0), C.max(0), color="0.8", lw=0)
ax.set_xlabel("Generation"); ax.set_ylabel("Hypervolume / mean final hypervolume")
ax.set_ylim(0.9, 1.005); ax.grid(lw=0.3, color="0.85"); ax.legend(frameon=False, loc="lower right", fontsize=8)
fig.tight_layout(); save(fig, "fig_convergence")
print("gambar selesai:", sorted(os.listdir(O)))
