"""Tabel final paper dalam LaTeX (booktabs): python make_tables.py  -> tables/*.tex
Membutuhkan folder evaluator di ../evaluator (untuk evaluasi ulang solusi representatif)."""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "evaluator"))
from systems import System
import factor_guard; factor_guard.check()
HERE = os.path.dirname(os.path.abspath(__file__))
D, O = os.path.join(HERE, "paper_data"), os.path.join(HERE, "tables"); os.makedirs(O, exist_ok=True)
SYS = [("ieee33", "IEEE 33-bus"), ("ieee69", "IEEE 69-bus"), ("rbts2", "RBTS Bus 2"), ("bawean", "Bawean")]
MS_ANN = (0.08 * 1.08 ** 15 / (1.08 ** 15 - 1) + 0.046) * 22000.0
def w(name, body): open(os.path.join(O, name), "w").write(body); print("tulis", name)

# --- representative solutions (E2)
lines = [r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}llrrrrrrl@{}}", r"\toprule",
         r"System & Solution & $f_1$ & $f_2$ & $f_3$ & SAIFI & SAIDI & $>H$ (\%) & Devices (MS/RCS/REC/tie) \\",
         r" & & \multicolumn{3}{c}{(kUSD/yr)} & (1/yr) & (h/yr) & & \\", r"\midrule"]
for s, title in SYS:
    df = pd.read_csv(os.path.join(D, f"front_{s}.csv")); F = df[["f1","f2","f3"]].to_numpy(); X = df.iloc[:, 3:].to_numpy().astype(int)
    span = F.max(0) - F.min(0); nrm = (F - F.min(0)) / np.where(span > 0, span, 1)
    picks = [("Societal", int(np.argmin(F[:,0]+F[:,1]))), ("Utility", int(np.argmin(F[:,0]+F[:,2]))),
             ("Compromise", int(np.argmin(np.linalg.norm(nrm, axis=1))))]
    sy = System(s)
    for k, (lab, i) in enumerate(picks):
        o = sy.evaluate(X[i]); x = X[i][:sy.nsec]; t = X[i][sy.nsec:]
        dev = f"{(x==1).sum()}/{(x==2).sum()}/{(x==3).sum()}/{(t>0).sum()}"
        lines.append(f"{title if k == 0 else ''} & {lab} & {o['f1']/1e3:.1f} & {o['f2']/1e3:.1f} & {o['f3']/1e3:.2f} & "
                     f"{o['SAIFI']:.3f} & {o['SAIDI']:.2f} & {100*o['share_above']:.0f} & {dev} \\\\")
    lines.append(r"\midrule" if s != "bawean" else r"\bottomrule")
lines.append(r"\end{tabular*}"); w("tab_representative.tex", "\n".join(lines))

# --- E3 objective comparison (best-known plans from E2 fronts and E3 runs)
lines = [r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{IEEE 33-bus} & \multicolumn{2}{c}{IEEE 69-bus} & \multicolumn{2}{c}{RBTS Bus 2} & \multicolumn{2}{c}{Bawean} \\",
         r"Budget (MS equiv.) & Regret & Price & Regret & Price & Regret & Price & Regret & Price \\", r"\midrule"]
vals = {}
for s, _ in SYS:
    rows = json.load(open(os.path.join(D, f"{s}_t60_h1.json")))
    F = pd.read_csv(os.path.join(D, f"front_{s}.csv"))[["f1","f2","f3"]].to_numpy()
    for n in (3, 5, 8, 10, 12):
        cand = [(r["f2"], r["f3"]) for r in rows if float(r["budget"]) == n]
        cand += [(F[i,1], F[i,2]) for i in np.where(F[:,0] <= n*MS_ANN + 1e-6)[0]]
        xE = min(cand, key=lambda c: c[0]); xK = min(cand, key=lambda c: c[1])
        vals[(s, n)] = (100*(xE[1]-xK[1])/xK[1], 100*(xK[0]-xE[0])/xE[0])
for n in (3, 5, 8, 10, 12):
    cells = " & ".join(f"{vals[(s,n)][0]:.1f} & {vals[(s,n)][1]:.1f}" for s, _ in SYS)
    lines.append(f"{n} & {cells} \\\\")
lines += [r"\bottomrule", r"\end{tabular*}"]; w("tab_objectives_budget.tex", "\n".join(lines))

# --- E7 computation time (Intel Core i7-13620H, 16 logical cores, 32 GB RAM; final code, idle machine)
BS = chr(92)                                   # one backslash, written explicitly
TS = BS + ","                                  # LaTeX thin space
rows_e7 = [("IEEE 33-bus", "24.5", "26.1", "80" + TS + "000", "1" + TS + "547$^{a}$"),
           ("IEEE 69-bus", "54.5", "87.8", "80" + TS + "000", "6" + TS + "721$^{a}$"),
           ("RBTS Bus 2", "31.3", "19.2", "50" + TS + "000", "461$^{a}$"),
           ("Bawean", "124.6", "168.6", "160" + TS + "000", "5" + TS + "403$^{b}$")]
EOL = " " + BS + BS
lines = [BS + "begin{tabular}{@{}lrrrr@{}}", BS + "toprule",
         "System & " + BS + "multicolumn{2}{c}{Time per evaluation (ms)} & Evaluations & Wall time of one run" + EOL,
         " & Connectivity & AC, four levels & per run & (s, 16 logical cores)" + EOL, BS + "midrule"]
lines += [" & ".join(r) + EOL for r in rows_e7]
lines += [BS + "bottomrule",
          BS + "multicolumn{5}{@{}l}{" + BS + "footnotesize $^{a}$ fastest of ten seeds, measured during the "
          "experiments; $^{b}$ one seed run with the final code.}" + EOL,
          BS + "end{tabular}"]
w("tab_time.tex", "\n".join(lines))
