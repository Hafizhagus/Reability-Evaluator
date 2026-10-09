"""Seberapa peka hasil IEEE 69 terhadap reaktansi tiga saluran tie?
Kode memakai X = R (15-46: 1+j1, 50-59: 2+j2, 27-65: 1+j1 ohm). Data yang umum dipakai di
literatur rekonfigurasi: 15-46: 1+j0,5; 50-59: 2+j1; 27-65: 1+j0,5 ohm.
Rencana front IEEE 69 dengan f1 <= 250 kUSD/th dievaluasi dengan kedua data."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "evaluator"))
import ieee69_data
from systems import System
HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "paper_data", "front_ieee69.csv"))
X = df.iloc[:, 3:].to_numpy().astype(int); F = df[["f1", "f2", "f3"]].to_numpy()
# rencana termahal (lebih dari 25 perangkat) butuh puluhan detik per evaluasi dan jauh di atas
# optimum sosial maupun utilitas; cek dibatasi pada rencana dengan f1 <= 250 kUSD/th
keep = F[:, 0] <= 250e3
X, F = X[keep], F[keep]
res = {}
for label, ties in (("kode (X = R)", list(ieee69_data.TIES_IEEE)),
                    ("literatur (X = R/2 pada tiga tie)", [(11, 43, 0.5, 0.5), (13, 21, 0.5, 0.5), (15, 46, 1.0, 0.5),
                                                          (50, 59, 2.0, 1.0), (27, 65, 1.0, 0.5)])):
    ieee69_data.TIES_IEEE = ties
    sy = System("ieee69")
    O = [sy.evaluate(g) for g in X]
    res[label] = np.array([[o["f1"], o["f2"], o["f3"], o["SAIDI"], o["share_above"]] for o in O])
    print(label, "selesai", flush=True)
a, b = res["kode (X = R)"], res["literatur (X = R/2 pada tiga tie)"]
print("cek: f2 kode == front:", float(np.abs(a[:, 1] - F[:, 1]).max()))
d2 = (b[:, 1] - a[:, 1]) / a[:, 1] * 100; d3 = (b[:, 2] - a[:, 2]) / np.maximum(a[:, 2], 1e-9) * 100
print(f"rencana: {len(X)}; f2 berubah pada {(np.abs(d2) > 1e-9).sum()} rencana, maks {np.abs(d2).max():.3f}%; "
      f"f3 berubah pada {(np.abs(d3) > 1e-9).sum()} rencana, maks {np.abs(d3).max():.3f}%")
for nm, col in (("sosial", 1), ("utilitas", 2)):
    ia, ib = int(np.argmin(a[:, 0] + a[:, col])), int(np.argmin(b[:, 0] + b[:, col]))
    print(f"optimum {nm}: indeks {ia} -> {ib}; nilai {(a[ia,0]+a[ia,col])/1e3:.3f} -> {(b[ib,0]+b[ib,col])/1e3:.3f} kUSD/th")
json.dump(dict(n=len(X), f2_changed=int((np.abs(d2) > 1e-9).sum()), f2_max_pct=float(np.abs(d2).max()),
               f3_changed=int((np.abs(d3) > 1e-9).sum()), f3_max_pct=float(np.abs(d3).max())),
          open(os.path.join(HERE, "check_ieee69_ties.json"), "w"), indent=1)
