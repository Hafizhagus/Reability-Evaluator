"""Analisis E2: hypervolume antar seed, front gabungan, dan solusi representatif.
  python analyze_front.py results/ieee33_t60_h1
Normalisasi HV memakai titik ideal/nadir gabungan semua seed; titik acuan 1,1."""
import glob, json, os, sys
import numpy as np
from pymoo.indicators.hv import HV
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from systems import System

d = sys.argv[1]
files = sorted(glob.glob(os.path.join(d, "seed*.npz")))
runs = [np.load(f, allow_pickle=True) for f in files]
cfg = json.loads(str(runs[0]["cfg"]))
allF = np.vstack([r["F"] for r in runs]); allX = np.vstack([r["X"] for r in runs])
nd = NonDominatedSorting().do(allF, only_non_dominated_front=True)
Fg, Xg = allF[nd], allX[nd]
_, uniq = np.unique(Xg, axis=0, return_index=True); Fg, Xg = Fg[uniq], Xg[uniq]
ideal, nadir = Fg.min(0), Fg.max(0); span = np.where(nadir > ideal, nadir - ideal, 1.0)
hv = HV(ref_point=np.full(3, 1.1))
norm = lambda F: (np.asarray(F) - ideal) / span
hvs = [hv(norm(r["F"])) for r in runs]
print(f"{cfg['name']}: {len(runs)} seed | HV = {np.mean(hvs):.4f} +- {np.std(hvs, ddof=1) if len(hvs) > 1 else 0:.4f} "
      f"| front gabungan {len(Fg)} solusi | waktu rata-rata {np.mean([float(r['secs']) for r in runs]):.0f} s/seed")
# convergence (HV of the running front per generation), one row per seed
conv = []
for r in runs:
    hist = json.loads(str(r["hist"]))
    conv.append([hv(norm(F)) if len(F) else 0.0 for F in hist])
np.savetxt(os.path.join(d, "convergence_hv.csv"), np.array(conv, dtype=float), delimiter=",")
sy = System(**cfg)
picks = {"min f1+f2 (sosial)": int(np.argmin(Fg[:, 0] + Fg[:, 1])),
         "min f1+f3 (utilitas)": int(np.argmin(Fg[:, 0] + Fg[:, 2])),
         "kompromi (terdekat ke ideal, ternormalisasi)": int(np.argmin(np.linalg.norm(norm(Fg), axis=1)))}
rows = []
for lab, i in picks.items():
    o = sy.evaluate(Xg[i])
    comp = sum(1 for v in Xg[i][:sy.nsec] if v)
    print(f"\n{lab}: f1={o['f1']/1e3:.1f} k$/th, ECOST={o['f2']/1e3:.1f} k$/th, KOMP={o['f3']/1e3:.2f} k$/th | "
          f"SAIFI={o['SAIFI']:.3f} SAIDI={o['SAIDI']:.3f} ENS={o['ENS']:.2f} | >{cfg.get('thr_h',1):g} jam: {100*o['share_above']:.0f}% pelanggan")
    print("   perangkat:", sy.describe(Xg[i]))
    rows.append(dict(label=lab, plan=sy.describe(Xg[i]), **{k: float(v) for k, v in o.items()}))
json.dump(dict(hv=hvs, representatives=rows), open(os.path.join(d, "summary.json"), "w"), indent=1)
np.savetxt(os.path.join(d, "front_global.csv"), np.hstack([Fg, Xg]), delimiter=",",
           header="f1,f2,f3," + ",".join(f"g{i}" for i in range(Xg.shape[1])), comments="")
