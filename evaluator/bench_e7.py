"""E7: waktu per evaluasi yang bersih (satu core, tanpa beban lain).
Jalankan di folder evaluator saat PC tidak menjalankan apa pun:  py bench_e7.py
Mengukur waktu rata-rata per evaluasi rencana acak berkepadatan khas front Pareto,
untuk mode konektivitas dan mode AC empat kelompok beban (formulasi v4)."""
import os, time, platform
import numpy as np
from systems import System

def plans(sy, n, rng):
    out = []
    for _ in range(n):
        g = np.zeros(sy.nvar, int)
        dens = rng.uniform(0.05, 0.35)
        on = rng.random(sy.nvar) < dens
        g[on] = [rng.integers(1, sy.xu[j] + 1) for j in np.where(on)[0]]
        out.append(g)
    return out

if __name__ == "__main__":
    print(f"CPU: {platform.processor()} | core logis: {os.cpu_count()} | Python {platform.python_version()}")
    rng = np.random.default_rng(2026)
    for name in ("ieee33", "ieee69", "rbts2", "bawean"):
        P = plans(System(name), 60, rng)
        for mode in ("connectivity", "loadflow"):
            sy = System(name, mode=mode)
            sy.evaluate(P[0])                      # warm-up
            t0 = time.perf_counter()
            for g in P:
                sy.evaluate(g)
            dt = (time.perf_counter() - t0) / len(P)
            print(f"{name:7s} {mode:12s}: {1000*dt:7.1f} ms per evaluasi (f1, f2, f3 dan indeks)")
