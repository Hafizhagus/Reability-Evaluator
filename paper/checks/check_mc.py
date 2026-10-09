"""Cek ulang Tabel 1: kompensasi analitik (rekursi Panjer) vs Monte Carlo.
  python check_mc.py [bulan]        (default 3e6 bulan per titik beban)
Untuk tiap sistem dan dua rencana (tanpa perangkat baru; campuran MS/RCS + semua tie):
pengali E[m(D_j)] SEMUA titik beban disimulasikan, lalu f3 sistem dibandingkan."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "evaluator"))
from systems import System
from compensation import expected_multiplier, BRACKETS_H, MULT

MONTHS = int(float(sys.argv[1])) if len(sys.argv) > 1 else 3_000_000
HERE = os.path.dirname(os.path.abspath(__file__))


def mc(sev_one, thr_h, months, seed):
    rng = np.random.default_rng(seed)
    ds = np.array(list(sev_one.keys()), dtype=float); rates = np.array(list(sev_one.values())) / 12.0
    if len(ds) == 0:
        return 0.0, 0.0
    tot = 0.0; tot2 = 0.0; done = 0
    while done < months:
        m_ = min(500_000, months - done)
        counts = rng.poisson(rates[None, :], size=(m_, len(ds)))
        D = (counts * ds[None, :]).sum(axis=1) / 60.0
        ex = D - thr_h; m = np.zeros(m_)
        for lo, v in zip(BRACKETS_H, MULT):
            m[ex > lo + 1e-12] = v
        tot += m.sum(); tot2 += (m * m).sum(); done += m_
    mean = tot / months; var = max(tot2 / months - mean * mean, 0.0) / months
    return mean, var


def plans(sy):
    yield "tanpa perangkat baru", np.zeros(sy.nvar, int)
    g = np.zeros(sy.nvar, int)
    g[:sy.nsec:4] = 1; g[2:sy.nsec:9] = 2; g[sy.nsec:] = 1      # campuran MS/RCS + semua tie
    yield "MS/RCS campuran + tie", g


rows = []; worst_sys = 0.0; worst_pt = 0.0; worst_pt_big = 0.0; worst_sig = 0.0
for si, name in enumerate(("ieee33", "ieee69", "rbts2", "bawean")):
    sy = System(name)
    N = np.array([sy.net.N[j] for j in sy.loads], float)
    for pi, (label, g) in enumerate(plans(sy)):
        x, w = sy.decode(g)
        o = sy.ev.evaluate(x, w, severity=True)
        sev = [o["sev"][j] for j in sy.loads]
        Em, ED, _ = expected_multiplier(sev, sy.thr_h)
        m = np.zeros(len(sev)); v = np.zeros(len(sev))
        for i in range(len(sev)):
            m[i], v[i] = mc(sev[i], sy.thr_h, MONTHS, seed=1_000_000 * si + 100_000 * pi + i)
        wgt = 12.0 * N * sy.base
        f3_an = float((wgt * Em).sum()); f3_mc = float((wgt * m).sum()); se = float(np.sqrt((wgt ** 2 * v).sum()))
        rel = abs(f3_mc - f3_an) / f3_an
        pos = Em > 0
        relp = np.abs(m[pos] - Em[pos]) / Em[pos]
        sig = np.abs(m[pos] - Em[pos]) / np.sqrt(np.maximum(v[pos], 1e-300))
        big = Em[pos] >= 0.01
        worst_sys = max(worst_sys, rel); worst_pt = max(worst_pt, float(relp.max()))
        worst_pt_big = max(worst_pt_big, float(relp[big].max()) if big.any() else 0.0)
        worst_sig = max(worst_sig, float(sig.max()))
        rows.append(dict(system=name, plan=label, f3_panjer=f3_an, f3_mc=f3_mc, f3_mc_se=se, rel_diff=rel,
                         n_points=int(pos.sum()), max_point_rel=float(relp.max()),
                         max_point_rel_mult_ge_0p01=float(relp[big].max()) if big.any() else None,
                         max_point_sigma=float(sig.max())))
        print(f"{name:7s} {label:22s}: f3 Panjer {f3_an:11.2f}  MC {f3_mc:11.2f} +- {se:7.2f}  selisih {100*rel:.3f}% "
              f"({abs(f3_mc-f3_an)/se:.1f} sigma) | per titik: maks {100*relp.max():.2f}%, "
              f"maks (pengali >= 0,01) {100*(relp[big].max() if big.any() else 0):.2f}%, maks {sig.max():.1f} sigma, {int(pos.sum())} titik",
              flush=True)
print(f"\nf3 sistem: selisih terbesar {100*worst_sys:.3f}% | per titik beban: terbesar {100*worst_pt:.2f}%, "
      f"terbesar untuk pengali >= 0,01: {100*worst_pt_big:.2f}%, terbesar {worst_sig:.1f} sigma ({MONTHS:.0e} bulan per titik)")
json.dump(dict(months=MONTHS, worst_system_rel=worst_sys, worst_point_rel=worst_pt,
               worst_point_rel_mult_ge_0p01=worst_pt_big, worst_point_sigma=worst_sig, rows=rows),
          open(os.path.join(HERE, "check_mc.json"), "w"), indent=1)
