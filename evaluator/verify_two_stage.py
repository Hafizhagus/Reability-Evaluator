"""Verifikasi independen restorasi DUA TAHAP (campuran MS/RCS, tie manual/remote).
Tahap 1 (5 menit): oracle brute-force dengan partisi perangkat remote saja dan tie
remote saja. Tahap 2 (60 menit): oracle brute-force dengan partisi semua perangkat dan
semua tie terbangun, DENGAN SYARAT pelanggan yang pulih di tahap 1 tetap bertegangan
dan tie yang ditutup di tahap 1 tetap tertutup. Nilai = (kW pulih, pelanggan pulih)."""
import itertools, random, sys
import numpy as np
from network import Network, NONE, MS, RCS, bfs_loadflow
from reliability import Evaluator, default_params
from verify_oracle import random_feeder

def brute(net, xpart, allowed_ties, k, p, must_on=frozenset(), must_closed=()):
    blk = {net.root: 0}; head = {0: None}; nb = 1
    for v in net.order[1:]:
        s = net.sec_of_bus[v]
        if xpart[s] != NONE:
            blk[v] = nb; head[nb] = s; nb += 1
        else:
            blk[v] = blk[net.parent[v]]
    par = {b: (blk[net.sec[head[b]]["frm"]] if head[b] is not None else -1) for b in range(nb)}
    F = blk[net.sec[k]["to"]]
    def desc(b):
        while b != -1:
            if b == F: return True
            b = par[b]
        return False
    base = {v for v in range(net.n) if not desc(blk[v])}
    if net.root not in base:
        return (0.0, 0.0), base, set()
    dead = sorted({blk[v] for v in range(net.n) if desc(blk[v]) and blk[v] != F})
    buses = {b: {v for v in range(net.n) if blk[v] == b} for b in dead}
    forced = [b for b in dead if buses[b] & must_on]
    free = [b for b in dead if b not in forced]
    opt_ties = [c for c in allowed_ties if c not in must_closed]
    best = (-1.0, -1.0)
    for tmask in itertools.product([0, 1], repeat=len(opt_ties)):
        closed = list(must_closed) + [c for c, m in zip(opt_ties, tmask) if m]
        for bmask in itertools.product([0, 1], repeat=len(free)):
            en = set(base)
            for b in forced: en |= buses[b]
            for b, m in zip(free, bmask):
                if m: en |= buses[b]
            if any(net.ties[c]["a"] not in en or net.ties[c]["b"] not in en for c in closed):
                continue
            edges = [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"]) for d in net.sec
                     if d["frm"] in en and d["to"] in en]
            edges += [(net.ties[c]["a"], net.ties[c]["b"], net.ties[c]["R"], net.ties[c]["X"],
                       net.ties[c]["ampacity_a"]) for c in closed]
            res = bfs_loadflow(net, en, edges)
            if not (res["radial"] and res["converged"]): continue
            V = res["V"][list(en)]
            if V.min() < p["v_min"] - 1e-12 or V.max() > p["v_max"] + 1e-12: continue
            if any(a > edges[e][4] for e, a in res["I_a"].items()): continue
            rest = {v for b in dead for v in buses[b]} & en
            best = max(best, (round(sum(net.P[v] for v in rest), 6), float(sum(net.N[v] for v in rest))))
    return best, base, {v for b in dead for v in buses[b]}

if __name__ == "__main__":
    cases = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    tlim = float(sys.argv[2]) if len(sys.argv) > 2 else 250
    import time; t0 = time.perf_counter()
    rng = random.Random(77); checked = mism1 = mism2 = bind = 0
    for case in range(cases):
        n = rng.randint(10, 13)
        net = random_feeder(rng, n)
        lf = bfs_loadflow(net, range(n), [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"]) for d in net.sec])
        p = default_params(); p.update(v_min=float(np.nanmin(lf["V"])) - 0.004, tau_rcs=5/60, tau_ms=1.0)
        ev = Evaluator(net, p, mode="loadflow")
        x = [rng.choice([NONE, MS, RCS]) for _ in net.sec]
        w = [rng.choice([1, 2]) for _ in net.ties]
        xr = [d if d == RCS else NONE for d in x]
        rem_ties = [c for c in range(len(w)) if w[c] == 2]; all_ties = list(range(len(w)))
        for k in range(len(net.sec)):
            r, _ = ev.fault_durations(x, w, k)
            s1_rest, s1_closed = ev._dbg_stage["remote"]
            b1, base1, dead1 = brute(net, xr, rem_ties, k, p)
            got1 = {v for v in dead1 if abs(r[v] - p["tau_rcs"]) < 1e-9}
            v1 = (round(sum(net.P[v] for v in got1), 6), float(sum(net.N[v] for v in got1)))
            b2, base2, dead2 = brute(net, x, all_ties, k, p, must_on=frozenset(s1_rest), must_closed=tuple(s1_closed))
            got2 = {v for v in dead2 if r[v] <= p["tau_ms"] + 1e-9}
            v2 = (round(sum(net.P[v] for v in got2), 6), float(sum(net.N[v] for v in got2)))
            checked += 1
            if b2[0] >= 0 and b2[0] < round(sum(net.P[v] for v in dead2), 6) - 1e-6: bind += 1
            if b1[0] >= 0 and v1 != b1:
                mism1 += 1; print(f"TAHAP1 beda case {case} fault {k}: eval {v1} oracle {b1}")
            if b2[0] >= 0 and v2 != b2:
                mism2 += 1; print(f"TAHAP2 beda case {case} fault {k}: eval {v2} oracle {b2}")
        if time.perf_counter() - t0 > tlim: break
    print(f"diperiksa {checked} gangguan ({case+1} jaringan acak); tahap 2 dengan batas mengikat: {bind}; "
          f"tidak cocok tahap 1: {mism1}, tahap 2: {mism2}")
