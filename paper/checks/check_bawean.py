"""Bawean claims of Section 5.2: tie upgrade, single additions, cheapest plan below the threshold."""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "evaluator"))
from systems import System
S = System("bawean"); n = S.nvar
z = np.zeros(n, dtype=int); base = S.evaluate(z)
g = z.copy(); g[-1] = 1; up = S.evaluate(g)
print("gen:", n, "| MS tahunan", round(S.ms_ann, 2))
print("dasar      :", {k: round(float(v), 4) for k, v in base.items()})
print("upgrade tie:", {k: round(float(v), 4) for k, v in up.items()}, S.describe(g))
print(f"  biaya {up['f1']/1e3:.3f} kUSD/th | SAIDI {base['SAIDI']:.3f} -> {up['SAIDI']:.3f} | ECOST turun {(base['f2']-up['f2'])/1e3:.3f} | kompensasi turun {(base['f3']-up['f3'])/1e3:.3f} | share {base['share_above']:.4f} -> {up['share_above']:.4f}")
print(f"  f1+f2: dasar {(base['f1']+base['f2'])/1e3:.3f} vs upgrade {(up['f1']+up['f2'])/1e3:.3f} | f1+f3: dasar {(base['f1']+base['f3'])/1e3:.3f} vs upgrade {(up['f1']+up['f3'])/1e3:.3f}")
# all single additions: one device (MS/RCS/REC) on one candidate section
t0 = time.time(); best = {"soc": (1e18, None), "util": (1e18, None)}; zero_share = []; n_bad = 0
for i in range(S.nsec):
    for v in (1, 2, 3):
        g = z.copy(); g[i] = v; o = S.evaluate(g)
        if not o["rec_ok"]: n_bad += 1; continue
        for key, val in (("soc", o["f1"] + o["f2"]), ("util", o["f1"] + o["f3"])):
            if val < best[key][0]: best[key] = (val, S.describe(g), o["f1"])
        if o["share_above"] < 1e-12: zero_share.append((o["f1"], S.describe(g)))
print(f"penambahan tunggal: {3 * S.nsec} rencana dievaluasi ({n_bad} melanggar batas recloser seri), {time.time()-t0:.0f} s")
print(f"  f1+f2 terbaik {best['soc'][0]/1e3:.3f} ({best['soc'][1]}) vs dasar {(base['f1']+base['f2'])/1e3:.3f} -> {'TIDAK ada yang lebih baik' if best['soc'][0] >= base['f1']+base['f2'] else 'ADA yang lebih baik'}")
print(f"  f1+f3 terbaik {best['util'][0]/1e3:.3f} ({best['util'][1]}) vs dasar {(base['f1']+base['f3'])/1e3:.3f} -> {'TIDAK ada yang lebih baik' if best['util'][0] >= base['f1']+base['f3'] else 'ADA yang lebih baik'}")
print("  penambahan tunggal yang membuat share = 0:", sorted(zero_share)[:5], "(jumlah", len(zero_share), ")")
