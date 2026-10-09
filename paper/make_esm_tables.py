"""Isi tabel Online Resource 1 (ESM) dari data: python make_esm_tables.py  -> tables/esm_*.tex
esm_scdf (Tabel S3 fungsi kerugian sektoral), esm_budget (S5 perbandingan anggaran), esm_loadlevel (S6 studi
tingkat beban), esm_social (S7 sensitivitas optimum sosial), esm_utility (S8 sensitivitas optimum utilitas).
Baris yang dihasilkan disalin ke manuscript/ESM_1.tex; tabel S1, S2 dan S4 diisi dari data model (bukan hasil).
Membutuhkan ../evaluator dan paper_data/ (front E2, E3, e4_analysis.json, e5.json, e5bw.json)."""
import json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "evaluator"))
import ecost
from systems import System
import factor_guard; factor_guard.check()
D, O = os.path.join(HERE, "paper_data"), os.path.join(HERE, "tables"); os.makedirs(O, exist_ok=True)
SYS = [("ieee33", "IEEE 33-bus"), ("ieee69", "IEEE 69-bus"), ("rbts2", "RBTS Bus 2"), ("bawean", "Bawean")]
MS_ANN = (0.08 * 1.08 ** 15 / (1.08 ** 15 - 1) + 0.046) * 22000.0
def w(name, lines): open(os.path.join(O, name), "w").write("\n".join(lines) + "\n"); print("tulis", name)

# --- Table S3: sector customer damage functions
f = ecost.conversion_factor()
rows = [f"{s.capitalize()} & " + " & ".join(f"{f * v:.3f}" for v in ecost.SCDF_RBTS[s]) + r" \\"
        for s in ("residential", "commercial", "industrial", "government")]
w("esm_scdf.tex", [f"% faktor konversi {f:.4f}"] + rows)

# --- Table S5: budget comparison (same best-known rule as make_tables.py)
lines = []
for s, title in SYS:
    rows = json.load(open(os.path.join(D, f"{s}_t60_h1.json")))
    F = pd.read_csv(os.path.join(D, f"front_{s}.csv"))[["f1", "f2", "f3"]].to_numpy()
    for k, n in enumerate((3, 5, 8, 10, 12)):
        cand = [(r["f2"], r["f3"]) for r in rows if float(r["budget"]) == n]
        cand += [(F[i, 1], F[i, 2]) for i in np.where(F[:, 0] <= n * MS_ANN + 1e-6)[0]]
        xE = min(cand, key=lambda c: c[0]); xK = min(cand, key=lambda c: c[1])
        R = [r["f3"] for r in rows if float(r["budget"]) == n and r["obj"] == "comp" and r["seed"] != "exact"]
        g = [100 * (v - xK[1]) / xK[1] for v in R]
        lines.append(f"{title if k == 0 else ''} & {n} & {len(R)} & {100*(xE[1]-xK[1])/xK[1]:.1f} & "
                     f"{100*(xK[0]-xE[0])/xE[0]:.1f} & {np.median(g):.1f} [{min(g):.1f}, {max(g):.1f}] \\\\")
    lines.append(r"\midrule" if s != "bawean" else r"\bottomrule\end{tabular}\end{table}")
w("esm_budget.tex", lines)

# --- Table S6: load-level study
A = json.load(open(os.path.join(D, "e4_analysis.json")))
lines = []
for t in (10, 25, 50, 100, 200):
    fx = A[f"fixed|{t}"]; e = lambda lm: 100 * (fx[lm] - fx["ldc4"]) / fx["ldc4"]
    r = lambda ob, lm: 100 * A[f"{t}|{ob}|{lm}"]["regret"]
    lines.append(f"{t} & {e('peak'):+.1f} & {e('avg'):+.2f} & {r('social','peak'):.2f} & {r('social','avg'):.2f} & "
                 f"{r('becost:5','peak'):.2f} & {r('becost:5','avg'):.2f} \\\\")
w("esm_loadlevel.tex", lines)

# --- Tables S7 and S8: sensitivity of the societal and utility optima
E5 = json.load(open(os.path.join(D, "e5.json"))) + json.load(open(os.path.join(D, "e5bw.json")))
SY = {s: System(s) for s, _ in SYS}
def best(s, var, ob):
    cfg = dict(name=s, **var)
    R = [r for r in E5 if r["cfg"] == cfg and r["obj"] == ob]
    key = (lambda r: r["own"]["f1"] + r["own"]["f2"]) if ob == "social" else (lambda r: r["own"]["f1"] + r["own"]["f3"])
    return min(R, key=key)
def devices(s, g):
    sy = SY[s]; x = np.array(g[:sy.nsec]); t = np.array(g[sy.nsec:])
    parts = [f"{int((x == k).sum())} {lab}" for k, lab in ((1, "MS"), (2, "RCS"), (3, "REC")) if (x == k).any()]
    if s == "bawean":
        if (t > 0).any(): parts.append("tie upgrade")
        return " + ".join(parts) if parts else "existing"
    if (t > 0).any(): parts.append(f"{int((t > 0).sum())} tie")
    return " + ".join(parts) if parts else "none"
V5 = [("Base", {}), (r"$\tau^{\mathrm{MS}}=45$ min", dict(tau_ms_min=45.0)), (r"$\tau^{\mathrm{MS}}=90$ min", dict(tau_ms_min=90.0)),
      ("RCS/MS 1.3", dict(cost_ratio=1.3)), ("RCS/MS 2.0", dict(cost_ratio=2.0)),
      (r"CLPU $a=0.6$", dict(clpu_a=0.6)), (r"CLPU $a=1.0$", dict(clpu_a=1.0))]
lines = []
for lab, var in V5:
    cells = []
    for s, _ in SYS:
        b = best(s, var, "social"); cells.append(f"{devices(s, b['g'])}; {(b['own']['f1'] + b['own']['f2']) / 1e3:.1f}")
    lines.append(f"{lab} & " + " & ".join(cells) + r" \\")
w("esm_social.tex", lines)
V6 = [("Base", {}), (r"$H=2$ h", dict(thr_h=2.0)), (r"$H=6$ h", dict(thr_h=6.0)),
      (r"Residential tariff $\times2$", dict(tariff_scale={"residential": 2.0})),
      (r"Residential tariff $\times0.5$", dict(tariff_scale={"residential": 0.5})),
      (r"$\tau^{\mathrm{MS}}=90$ min", dict(tau_ms_min=90.0)), ("RCS/MS cost ratio 1.3", dict(cost_ratio=1.3)),
      (r"CLPU, $a=1.0$", dict(clpu_a=1.0))]
lines = []
for lab, var in V6:
    cells = []
    for s, _ in SYS:
        b = best(s, var, "utility"); cells.append(f"{devices(s, b['g'])}; {100 * b['own']['share_above']:.0f}\\%")
    lines.append(f"{lab} & " + " & ".join(cells) + r" \\")
w("esm_utility.tex", lines)

# --- numbers quoted in the text of Sections 5.2 and 5.5
print("\nSelisih refinemen optimum sosial terhadap front (%):")
for s, title in SYS:
    F = pd.read_csv(os.path.join(D, f"front_{s}.csv"))[["f1", "f2", "f3"]].to_numpy()
    fs, fu = (F[:, 0] + F[:, 1]).min(), (F[:, 0] + F[:, 2]).min()
    bs, bu = best(s, {}, "social")["own"], best(s, {}, "utility")["own"]
    print(f"  {title:12s} sosial {100 * (fs - bs['f1'] - bs['f2']) / fs:.3f} | utilitas {100 * (fu - bu['f1'] - bu['f3']) / fu:.3f}")
print("Perubahan biaya optimum sosial terhadap kasus dasar (%):")
for lab, var in V5[1:]:
    ch = []
    for s, _ in SYS:
        b0 = best(s, {}, "social")["own"]; b = best(s, var, "social")["own"]
        ch.append(100 * ((b["f1"] + b["f2"]) / (b0["f1"] + b0["f2"]) - 1))
    print(f"  {lab:32s} " + "  ".join(f"{v:+6.2f}" for v in ch))
