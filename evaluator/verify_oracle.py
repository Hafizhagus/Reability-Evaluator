"""Verifikasi independen: restorasi eksak di evaluator vs oracle brute-force.

Oracle ini TIDAK memakai logika pencarian evaluator. Untuk setiap gangguan ia
menghitung sendiri blok, zona gangguan, dan daerah padam, lalu mencoba SEMUA kombinasi
(tie ditutup/dibuka) x (blok padam dinyalakan/tidak). Kombinasi valid bila jaringan
berenergi radial dan terhubung ke sumber, setiap tie tertutup menghubungkan dua bagian
berenergi, dan load flow AC memenuhi batas tegangan dan arus. Nilai terbaik =
(kW pulih, pelanggan pulih). Hasilnya harus sama persis dengan evaluator.
Jaringan uji: semua perangkat dan tie remote, sehingga hanya ada satu tahap efektif.
"""
import itertools, random, sys
from collections import deque
import numpy as np
from network import Network, NONE, RCS, bfs_loadflow
from reliability import Evaluator, default_params


def random_feeder(rng, n):
    sections = []
    for v in range(1, n):
        u = rng.randrange(max(0, v - 3), v)          # mostly chain-like with laterals
        sections.append(dict(frm=u, to=v, R=rng.uniform(0.3, 1.2), X=rng.uniform(0.3, 1.0),
                             length_km=1.0, ampacity_a=rng.choice([150.0, 200.0, 400.0])))
    parent = {s["to"]: s["frm"] for s in sections}
    ties, tries = [], 0
    while len(ties) < rng.randint(2, 3) and tries < 100:
        tries += 1
        a, b = rng.sample(range(1, n), 2)
        if parent.get(a) == b or parent.get(b) == a or any({a, b} == {t["a"], t["b"]} for t in ties):
            continue
        ties.append(dict(a=a, b=b, R=0.6, X=0.6, length_km=1.0, ampacity_a=400.0))
    P = [0.0] + [rng.uniform(60, 260) for _ in range(n - 1)]
    Q = [0.3 * p for p in P]
    N = [0] + [rng.randint(5, 60) for _ in range(n - 1)]
    return Network(n, sections, ties, P, Q, N, root=0, V_kv=11.0, sectors=["residential"] * n)


def _zone(net, x, k):
    """Buses of the faulted block (same block as the receiving end of section k)."""
    blk = {net.root: 0}; nb = 1
    for v in net.order[1:]:
        s = net.sec_of_bus[v]
        if x[s] != NONE:
            blk[v] = nb; nb += 1
        else:
            blk[v] = blk[net.parent[v]]
    F = blk[net.sec[k]["to"]]
    return {v for v in range(net.n) if blk[v] == F}


def oracle(net, x, w, k, p):
    """Brute-force best (kW, customers) restorable in the dead region of fault k."""
    blk = {net.root: 0}; head = {0: None}; nb = 1
    for v in net.order[1:]:
        s = net.sec_of_bus[v]
        if x[s] != NONE:
            blk[v] = nb; head[nb] = s; nb += 1
        else:
            blk[v] = blk[net.parent[v]]
    par = {b: (blk[net.sec[head[b]]["frm"]] if head[b] is not None else -1) for b in range(nb)}
    F = blk[net.sec[k]["to"]]
    def is_desc(b):
        while b != -1:
            if b == F:
                return True
            b = par[b]
        return False
    base = {v for v in range(net.n) if not is_desc(blk[v])}
    dead = sorted({blk[v] for v in range(net.n) if is_desc(blk[v]) and blk[v] != F})
    buses = {b: [v for v in range(net.n) if blk[v] == b] for b in dead}
    ties = [c for c in range(len(net.ties)) if w[c] > 0]
    best = (0.0, 0.0)
    if net.root not in base:                  # fault in the source block: nothing to feed from
        return best, dead, base
    for tmask in itertools.product([0, 1], repeat=len(ties)):
        closed = [c for c, m in zip(ties, tmask) if m]
        for bmask in itertools.product([0, 1], repeat=len(dead)):
            en = set(base)
            for b, m in zip(dead, bmask):
                if m:
                    en |= set(buses[b])
            if any(net.ties[c]["a"] not in en or net.ties[c]["b"] not in en for c in closed):
                continue
            edges = [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"]) for d in net.sec
                     if d["frm"] in en and d["to"] in en]
            edges += [(net.ties[c]["a"], net.ties[c]["b"], net.ties[c]["R"], net.ties[c]["X"],
                       net.ties[c]["ampacity_a"]) for c in closed]
            res = bfs_loadflow(net, en, edges)
            if not (res["radial"] and res["converged"]):
                continue
            V = res["V"][list(en)]
            if V.min() < p["v_min"] - 1e-12 or V.max() > p["v_max"] + 1e-12:
                continue
            if any(a > edges[e][4] for e, a in res["I_a"].items()):
                continue
            rest = en - base
            val = (round(sum(net.P[v] for v in rest), 6), float(sum(net.N[v] for v in rest)))
            best = max(best, val)
    return best, dead, base


if __name__ == "__main__":
    cases = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    rng = random.Random(2026)
    mism = checked = shed_cases = nonempty = 0
    for case in range(cases):
        n = rng.randint(10, 14)
        net = random_feeder(rng, n)
        base_lf = bfs_loadflow(net, range(n), [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"])
                                               for d in net.sec])
        p = default_params()
        p.update(v_min=float(np.nanmin(base_lf["V"])) - 0.004, tau_rcs=5 / 60, tau_ms=1.0)
        ev = Evaluator(net, p, mode="loadflow")
        x = [RCS if rng.random() < 0.5 else NONE for _ in net.sec]
        w = [2] * len(net.ties)
        for k in range(len(net.sec)):
            best, dead, base = oracle(net, x, w, k, p)
            if len(dead) > 11:
                continue
            r, _ = ev.fault_durations(x, w, k)
            dead_bus = {v for v in range(n) if v not in base}
            got = {v for v in dead_bus if abs(r[v] - p["tau_rcs"]) < 1e-9}
            val = (round(sum(net.P[v] for v in got), 6), float(sum(net.N[v] for v in got)))
            checked += 1
            blkF = None
            deadblk_bus = set()
            # buses of the dead blocks only (the faulted block itself cannot be restored)
            for v in dead_bus:
                pass
            deadblk_bus = {v for v in dead_bus if v not in _zone(net, x, k)}
            dead_kw = round(sum(net.P[v] for v in deadblk_bus), 6)
            if dead_kw > 0:
                nonempty += 1
                if best[0] < dead_kw - 1e-6:
                    shed_cases += 1                 # constraints actually bind here
            if val != best:
                mism += 1
                print(f"BEDA case {case} fault {k}: evaluator {val} vs oracle {best}")
    print(f"diperiksa {checked} kasus gangguan pada {cases} jaringan acak; "
          f"{nonempty} kasus dengan daerah padam yang bisa dipulihkan; "
          f"{shed_cases} di antaranya batas tegangan/arus mengikat (tidak semua bisa pulih); "
          f"tidak cocok: {mism}")
