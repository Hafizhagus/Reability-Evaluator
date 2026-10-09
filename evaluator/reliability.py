"""Exact reliability evaluator (steps 2b, 2d, 2e) -- stage-based restoration.

For a permanent fault on section k:
  1. Clearing: the nearest recloser above the fault (or the feeder CB) trips;
     every customer below it is interrupted.
  2. Restoration proceeds in stages. Stage 1 (time tau_rcs): only remote devices
     (RCS, REC) and remote ties can be operated. Stage 2 (time tau_ms): every
     installed device and every built tie can be operated.
     In each stage the operable devices cut the network into blocks; the block
     containing the fault (F) is isolated; interrupted customers outside the
     subtree of F are re-energized from the source; the de-energized subtree
     below F is restored through operable ties, split into radial sub-islands
     (one tie each) and checked by load flow in mode "loadflow".
  3. A customer's interruption duration is the time of the first stage in which
     it is restored, otherwise the repair time.
mode = "connectivity" (no load flow) or "loadflow" (restorations checked by AC LF).
"""
import itertools
from collections import deque
import numpy as np

from network import NONE, MS, RCS, REC, bfs_loadflow
import ecost as _ecost


def default_params():
    return dict(
        lambda_per_km=0.2,        # f/km/yr   (SPLN 59:1985)
        r_rep=3.0,                # h         (SPLN 59:1985)
        tau_ms=1.0,               # h         (manual switching)
        tau_rcs=5.0 / 60.0,       # h         (remote switching)
        v_min=0.90, v_max=1.05,   # p.u.      (SPLN 1:1995)
        load_level=1.0,           # multiplier on peak for LF checks and ENS
        cost=dict(MS=22000.0, RCS=61000.0, REC=105000.0, line_per_km=0.0),
        delta=0.08, h_sw=15, h_line=35, rho_om=0.046,
        m_rec=2,                  # max reclosers in series (constraint)
        order_search_max_ties=5,  # try all tie orders when demand is left dark
        scdf=_ecost.SCDF_RBTS,    # sector damage functions (USD via cost_factor)
        cost_factor=None,         # None -> ecost.conversion_factor()
        load_levels=None,         # [(multiplier of peak, probability)]; None -> [(load_level, 1)]
        clpu=None,                # None, or dict(kmax={sector: k}, tc_h=build-up time constant,
                                  #   tdec_h=decay time constant after pickup, step_h=pickup step)
    )


def crf(delta, h):
    return delta * (1 + delta) ** h / ((1 + delta) ** h - 1)


# ---------------------------------------------------------------------------
class Blocks:
    """Block partition induced by a device vector (only non-NONE entries cut)."""

    def __init__(self, net, x):
        self.blk = [-1] * net.n
        self.head = {0: None}
        self.blk[net.root] = 0
        nb = 1
        for v in net.order[1:]:
            s = net.sec_of_bus[v]
            if x[s] != NONE:
                self.blk[v] = nb
                self.head[nb] = s
                nb += 1
            else:
                self.blk[v] = self.blk[net.parent[v]]
        self.nb = nb
        self.par = {0: -1}
        self.kids = {b: [] for b in range(nb)}
        for b in range(1, nb):
            p = self.blk[net.sec[self.head[b]]["frm"]]
            self.par[b] = p
            self.kids[p].append(b)
        self.buses = {b: [] for b in range(nb)}
        for v in range(net.n):
            self.buses[self.blk[v]].append(v)
        self._desc = {}

    def desc(self, b):
        if b not in self._desc:
            out, q = {b}, deque([b])
            while q:
                u = q.popleft()
                for c in self.kids[u]:
                    out.add(c)
                    q.append(c)
            self._desc[b] = out
        return self._desc[b]


# ---------------------------------------------------------------------------
class Evaluator:
    def __init__(self, net, params=None, mode="connectivity"):
        self.net = net
        self.p = default_params() if params is None else params
        self.mode = mode
        self.lam = np.array([d.get("lam", self.p["lambda_per_km"] * d["length_km"])
                             for d in net.sec])
        self.rrep = np.array([d.get("r_rep", self.p["r_rep"]) for d in net.sec])
        self.n_lf = 0
        self._mult = self.p["load_level"]
        self._ctx = dict(interrupted=np.zeros(net.n, bool), restored_at={}, t_stage=0.0)

    # -- feasibility ---------------------------------------------------------
    def _clpu_factor(self, sector, d):
        """Undiversified pickup factor after an outage of d hours (delayed-exponential
        build-up towards kmax)."""
        c = self.p["clpu"]
        if not c or d <= 0:
            return 1.0
        kmax = c["kmax"].get(sector, 1.0)
        return 1.0 + (kmax - 1.0) * (1.0 - np.exp(-d / c["tc_h"]))

    def _clpu_now(self, sector, t_pick, t_now):
        """Factor at time t_now of a load re-energized at t_pick (after an outage of
        t_pick hours); decays exponentially towards 1 after pickup."""
        k0 = self._clpu_factor(sector, t_pick)
        if k0 == 1.0 or t_now <= t_pick:
            return k0
        return 1.0 + (k0 - 1.0) * np.exp(-(t_now - t_pick) / self.p["clpu"]["tdec_h"])

    def _scale(self, energized):
        """Per-bus load multiplier: level x CLPU factor (built up during the outage,
        decaying after pickup). Buses not yet picked up are taken at t_now."""
        ctx, m = self._ctx, self._mult
        t_now, pick = ctx["t_stage"], ctx["restored_at"]
        out = {}
        for v in energized:
            if ctx["interrupted"][v]:
                tp = pick.get(v, t_now)
                out[v] = m * self._clpu_now(self.net.sector[v], tp, t_now)
            else:
                out[v] = m
        return out

    def _feasible(self, energized, closed_ties):
        if self.mode == "connectivity":
            return True
        net = self.net
        scale = self._scale(energized)
        if not self.p["clpu"]:
            # without CLPU every energised load has the same multiplier: a compact key
            key = (frozenset(energized), tuple(sorted(closed_ties)), round(self._mult, 6))
        else:
            key = (frozenset(energized), tuple(sorted(closed_ties)),
                   tuple(round(scale[v], 4) for v in sorted(energized)))
        cache = self.__dict__.setdefault("_lfcache", {})
        if key in cache:
            return cache[key]
        edges = [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"])
                 for d in net.sec if d["frm"] in energized and d["to"] in energized]
        for c in closed_ties:
            t = net.ties[c]
            edges.append((t["a"], t["b"], t["R"], t["X"], t["ampacity_a"]))
        self.n_lf += 1
        res = bfs_loadflow(net, energized, edges, load_scale=scale)
        ok = bool(res["radial"] and res["converged"])
        if ok:
            V = res["V"][list(energized)]
            ok = V.min() >= self.p["v_min"] - 1e-12 and V.max() <= self.p["v_max"] + 1e-12
        if ok:
            ok = all(amps <= edges[e][4] for e, amps in res["I_a"].items())
        if ok and net.s_max_kva is not None:
            ok = res["S_root_kva"] <= net.s_max_kva
        if len(cache) > 20000:            # bounded memory per evaluator
            cache.clear()
        cache[key] = ok
        return ok

    def _restore_source(self, B, nonint, cand):
        """Re-supply from the normal source; with CLPU the pickup may violate limits,
        so the farthest leaf blocks are shed until the re-energized network is feasible."""
        if self.mode == "connectivity" or not cand:
            return set(cand)
        if not self.p["clpu"] and self._mult <= 1.0 + 1e-12:
            return set(cand)            # subset of the normal network at <= peak demand
        blocks = {B.blk[v] for v in cand}
        dist = self._distances(self.net.root, set(range(self.net.n)))
        while blocks:
            buses = {v for b in blocks for v in B.buses[b]} & cand
            if self._feasible(nonint | buses, []):
                return buses
            leaves = [b for b in blocks if not any(c in blocks for c in B.kids[b])]
            blocks.discard(max(leaves, key=lambda b: min(dist[v] for v in B.buses[b])))
        return set()

    # -- sub-island restoration with block shedding --------------------------
    def _distances(self, src, allowed):
        net, dist, q, adj = self.net, {src: 0.0}, deque([src]), {}
        for s, d in enumerate(net.sec):
            if d["frm"] in allowed and d["to"] in allowed:
                z = net.zabs(s)
                adj.setdefault(d["frm"], []).append((d["to"], z))
                adj.setdefault(d["to"], []).append((d["frm"], z))
        while q:
            u = q.popleft()
            for v, z in adj.get(u, []):
                if v not in dist:
                    dist[v] = dist[u] + z
                    q.append(v)
        return dist

    def _component(self, B, blocks, start):
        comp, q = {start}, deque([start])
        while q:
            b = q.popleft()
            for n in B.kids[b] + [B.par[b]]:
                if n in blocks and n not in comp:
                    comp.add(n)
                    q.append(n)
        return comp

    def _restorable(self, B, comp, t_blk, t_bus, energized, closed, c, topk=None):
        """Largest restorable set of blocks through tie c (exact): among the connected
        block sets that contain the tie block, maximise restored demand (kW), then the
        number of customers, subject to AC feasibility. Feasibility is monotone in the
        set (adding blocks only adds load), so infeasible sets are not extended and a
        demand upper bound prunes the search."""
        comp = set(comp)
        net = self.net
        bus = lambda S: {v for b in S for v in B.buses[b]}
        cl = closed + [c]
        if self._feasible(energized | bus(comp), cl):
            return [comp] if topk else comp
        if not self._feasible(energized | bus({t_blk}), cl):
            return [] if topk else set()
        adj = {b: [n for n in B.kids[b] + [B.par[b]] if n in comp] for b in comp}
        found = {}
        kw = {b: float(sum(net.P[v] for v in B.buses[b])) for b in comp}
        nc = {b: float(sum(net.N[v] for v in B.buses[b])) for b in comp}
        val = lambda S: (round(sum(kw[b] for b in S), 6), sum(nc[b] for b in S))
        best = [{t_blk}, val({t_blk})]
        limit = self.p.get("shed_search_limit", 20000)
        checks = [0]

        def rec(S, frontier, excluded):
            avail = comp - S - excluded
            if sum(kw[b] for b in S) + sum(kw[b] for b in avail) < best[1][0] - 1e-9:
                return
            ex = set(excluded)
            for cb in sorted(frontier, key=lambda b: -kw[b]):
                if cb in ex:
                    continue
                if checks[0] >= limit:
                    self.n_trunc = getattr(self, "n_trunc", 0) + 1
                    return
                S2 = S | {cb}
                checks[0] += 1
                if self._feasible(energized | bus(S2), cl):
                    v = val(S2)
                    if topk:
                        found[frozenset(S2)] = v
                    if v > best[1]:
                        best[0], best[1] = set(S2), v
                    nf = [n for n in frontier if n != cb and n not in ex] + \
                         [n for n in adj[cb] if n not in S2 and n not in frontier]
                    rec(S2, nf, ex)
                ex.add(cb)

        rec({t_blk}, list(adj[t_blk]), set())
        if topk:
            found[frozenset({t_blk})] = val({t_blk})
            ranked = sorted(found.items(), key=lambda kv: kv[1], reverse=True)
            return [set(k) for k, _ in ranked[:topk]]
        return best[0]

    def _candidate(self, B, remaining, energized, closed, c):
        t = self.net.ties[c]
        rem_bus = {v for b in remaining for v in B.buses[b]}
        for ti, far in ((t["a"], t["b"]), (t["b"], t["a"])):
            if ti in rem_bus and far in energized:
                comp = self._component(B, remaining, B.blk[ti])
                R = self._restorable(B, comp, B.blk[ti], ti, energized, closed, c)
                if R:
                    return sum(self.net.P[v] for b in R for v in B.buses[b]), R
        return None

    def _all_feasible_sets(self, B, comp, t_blk, energized, closed, c, limit=50000):
        """Every AC-feasible connected block set of comp that contains the tie block.
        Feasibility is monotone in the set, so only feasible sets are extended."""
        net = self.net
        bus = lambda S: {v for b in S for v in B.buses[b]}
        cl = closed + ([c] if c is not None else [])
        if not self._feasible(energized | bus({t_blk}), cl):
            return []
        adj = {b: [n for n in B.kids[b] + [B.par[b]] if n in comp] for b in comp}
        out = [frozenset({t_blk})]
        checks = [0]

        def rec(S, frontier, excluded):
            ex = set(excluded)
            for cb in list(frontier):
                if cb in ex:
                    continue
                if checks[0] >= limit:
                    self.n_trunc = getattr(self, "n_trunc", 0) + 1
                    return
                S2 = S | {cb}
                checks[0] += 1
                if self._feasible(energized | bus(S2), cl):
                    out.append(frozenset(S2))
                    nf = [n for n in frontier if n != cb and n not in ex] + \
                         [n for n in adj[cb] if n not in S2 and n not in frontier]
                    rec(S2, nf, ex)
                ex.add(cb)

        rec({t_blk}, list(adj[t_blk]), set())
        return out

    def _reachable(self, B, F, dead, base, ties, prev_restored):
        """Dead blocks that could be energised in this stage if network limits were
        ignored: components of the dead region holding a terminal of an operable tie whose
        far end is energised (or reachable, so chains count), and, in the second stage,
        components adjacent to a block that is already energised."""
        net = self.net
        reach, en = set(), set(base)
        changed = True
        while changed:
            changed = False
            left = set(dead) - reach
            for c in ties:
                t = net.ties[c]
                for ti, far in ((t["a"], t["b"]), (t["b"], t["a"])):
                    bt = B.blk[ti]
                    if bt in left and far in en:
                        comp = self._component(B, left, bt)
                        reach |= comp; left -= comp
                        en |= {v for b in comp for v in B.buses[b]}; changed = True
            if prev_restored:
                for b in list(left):
                    if b in left and any(n != -1 and n != F and B.buses.get(n) and B.buses[n][0] in en
                                         for n in B.kids[b] + [B.par[b]]):
                        comp = self._component(B, left, b)
                        reach |= comp; left -= comp
                        en |= {v for bb in comp for v in B.buses[bb]}; changed = True
        return reach

    def _restore_exact(self, B, dead, base, ties, budget=300000, closed0=()):
        """Exact restoration of the dead region: over all orders of operable ties and all
        feasible block sets per tie (chains allowed), maximise restored kW, then
        customers. Memoised on (restored blocks, used ties); pruned by a kW bound."""
        net = self.net
        kw = {b: float(sum(net.P[v] for v in B.buses[b])) for b in dead}
        nc = {b: float(sum(net.N[v] for v in B.buses[b])) for b in dead}
        best = [frozenset(), (0.0, 0.0), []]
        memo = set()
        nodes = [0]

        def rec(energized, closed, remaining, used, rblocks, val):
            key = (rblocks, frozenset(used))
            if key in memo:
                return
            memo.add(key)
            nodes[0] += 1
            if nodes[0] > budget:
                self.n_trunc = getattr(self, "n_trunc", 0) + 1
                return
            if val > best[1]:
                best[0], best[1], best[2] = rblocks, val, list(closed)
            if val[0] + sum(kw[b] for b in remaining) < best[1][0] - 1e-9:
                return
            rem_bus = {v for b in remaining for v in B.buses[b]}
            # (a) extension through a device boundary from an energized neighbour block
            seen = set()
            for b in remaining:
                if b in seen:
                    continue
                if any(n != -1 and n not in remaining and B.buses.get(n) and
                       B.buses[n][0] in energized for n in B.kids[b] + [B.par[b]]):
                    comp = self._component(B, remaining, b)
                    seen |= comp
                    for S in self._all_feasible_sets(B, comp, b, energized, closed, None):
                        nb = {v for bb in S for v in B.buses[bb]}
                        rec(energized | nb, list(closed), remaining - S, used,
                            rblocks | S, (round(val[0] + sum(kw[bb] for bb in S), 6),
                                          val[1] + sum(nc[bb] for bb in S)))
            # (b) through an operable tie whose far end is energized
            for c in ties:
                if c in used:
                    continue
                t = net.ties[c]
                for ti, far in ((t["a"], t["b"]), (t["b"], t["a"])):
                    if ti in rem_bus and far in energized:
                        comp = self._component(B, remaining, B.blk[ti])
                        for S in self._all_feasible_sets(B, comp, B.blk[ti], energized,
                                                         closed, c):
                            nb = {v for b in S for v in B.buses[b]}
                            rec(energized | nb, closed + [c], remaining - S, used | {c},
                                rblocks | S, (round(val[0] + sum(kw[b] for b in S), 6),
                                              val[1] + sum(nc[b] for b in S)))
        rec(set(base), list(closed0), frozenset(dead), frozenset(), frozenset(), (0.0, 0.0))
        restored = {v for b in best[0] for v in B.buses[b]}
        return restored, best[1], best[2]

    def _restore_joint(self, B, dead, base, ties, k=5, budget=4000):
        """Joint restoration over ties: each tie may take any of its k best feasible
        block sets; the combination restoring the most demand is kept. Used when the
        greedy/order search leaves demand dark."""
        net = self.net
        best = [set(), 0.0]
        calls = [0]

        def rec(energized, closed, remaining, used, restored, total):
            if total > best[1] + 1e-9:
                best[0], best[1] = set(restored), total
            if calls[0] > budget or not remaining:
                return
            rem_kw = sum(net.P[v] for b in remaining for v in B.buses[b])
            if total + rem_kw <= best[1] + 1e-9:
                return
            rem_bus = {v for b in remaining for v in B.buses[b]}
            for c in ties:
                if c in used:
                    continue
                t = net.ties[c]
                for ti, far in ((t["a"], t["b"]), (t["b"], t["a"])):
                    if ti in rem_bus and far in energized:
                        comp = self._component(B, remaining, B.blk[ti])
                        calls[0] += 1
                        for R in self._restorable(B, comp, B.blk[ti], ti, energized,
                                                  closed, c, topk=k):
                            nb = {v for b in R for v in B.buses[b]}
                            rec(energized | nb, closed + [c], remaining - R, used | {c},
                                restored | nb, total + sum(net.P[v] for v in nb))
        rec(set(base), [], set(dead), set(), set(), 0.0)
        return best[0], best[1]

    def _restore(self, B, dead, base_energized, ties, greedy, closed0=()):
        """Restore dead blocks through ties; returns (restored buses, demand)."""
        energized, closed, remaining = set(base_energized), list(closed0), set(dead)
        restored, total, pool = set(), 0.0, list(ties)
        while pool and remaining:
            if greedy:
                best = None
                for c in pool:
                    cand = self._candidate(B, remaining, energized, closed, c)
                    if cand and (best is None or cand[0] > best[0]):
                        best = (cand[0], c, cand[1])
                if best is None:
                    break
                dem, c, R = best
            else:
                c = pool[0]
                cand = self._candidate(B, remaining, energized, closed, c)
                if cand is None:
                    pool.pop(0)
                    continue
                dem, R = cand
            pool.remove(c)
            new = {v for b in R for v in B.buses[b]}
            restored |= new
            energized |= new
            closed.append(c)
            remaining -= R
            total += dem
        self._last_closed = list(closed)
        return restored, total

    # -- one fault -----------------------------------------------------------
    def fault_durations(self, x, w, k, blocks=None):
        net, p = self.net, self.p
        if blocks is None:
            blocks = self._stage_blocks(x)
        B_all = blocks["all"]
        # 1) clearing: nearest recloser at/above the fault block, else the CB
        P = B_all.blk[net.sec[k]["to"]]
        while P != 0 and x[B_all.head[P]] != REC:
            P = B_all.par[P]
        intr_blocks = B_all.desc(P)
        interrupted = np.zeros(net.n, bool)
        for b in intr_blocks:
            interrupted[B_all.buses[b]] = True
        r = np.where(interrupted, self.rrep[k], 0.0)
        # 2) restoration stages
        stages = [("remote", p["tau_rcs"], [c for c in range(len(w)) if w[c] == 2], p["tau_ms"]),
                  ("all", p["tau_ms"], [c for c in range(len(w)) if w[c] > 0], self.rrep[k])]
        restored_at = {}
        nonint = {v for v in range(net.n) if not interrupted[v]}
        prev_restored, prev_closed = set(), []
        self._dbg_stage = {}
        for key, t_stage, ties_all, t_end in stages:
            self._ctx = dict(interrupted=interrupted, restored_at=restored_at, t_stage=t_stage)
            B = blocks[key]
            F = B.blk[net.sec[k]["to"]]
            descF = B.desc(F)
            # customers restored in an earlier stage stay energized (no re-interruption)
            # and the ties closed then stay closed; this stage only adds to that state
            dead = {b for b in descF - {F}
                    if not any(v in prev_restored for v in B.buses[b])}
            dead_bus = {v for b in dead for v in B.buses[b]}
            ties = [c for c in ties_all if c not in prev_closed and
                    (net.ties[c]["a"] in dead_bus or net.ties[c]["b"] in dead_bus)]
            cand = {v for v in range(net.n) if interrupted[v] and B.blk[v] not in descF}
            src = self._restore_source(B, nonint, cand)
            base = nonint | src | prev_restored
            restored = set(src) | prev_restored
            closed = list(prev_closed)
            if dead and (ties or prev_restored):
                got, total = self._restore(B, dead, base, ties, greedy=True,
                                           closed0=prev_closed)
                closed = self._last_closed
                if self.mode == "loadflow":
                    # Blocks that no tie or operable device can ever connect to an energised
                    # part cannot be restored in any configuration; dropping them keeps the
                    # search exact. If greedy already restores every reachable block, it is
                    # optimal and the exact search is skipped.
                    reach = self._reachable(B, F, dead, base, ties, prev_restored)
                    reach_kw = sum(net.P[v] for b in reach for v in B.buses[b])
                    kw_g = sum(net.P[v] for v in got); nc_g = sum(net.N[v] for v in got)
                    if kw_g < reach_kw - 1e-9:
                        g, val, cl = self._restore_exact(B, reach, base, ties, closed0=prev_closed)
                        if val > (round(kw_g, 6), nc_g):
                            got, total, closed = g, val[0], cl
                if self.mode == "connectivity":
                    # without limits, every dead block reachable from the energized region
                    # through operable device boundaries is restored as well
                    en = base | got
                    grew = True
                    while grew:
                        grew = False
                        for b in dead:
                            if B.buses[b][0] in en:
                                continue
                            if any(n != -1 and B.buses.get(n) and B.buses[n][0] in en and
                                   n != F for n in B.kids[b] + [B.par[b]]):
                                en |= set(B.buses[b]); got |= set(B.buses[b]); grew = True
                restored |= got
            for v in restored:
                r[v] = min(r[v], t_stage)
                restored_at.setdefault(v, t_stage)
            prev_restored, prev_closed = set(restored), list(closed)
            self._dbg_stage[key] = (set(restored), list(closed))
            # 3) staged pickup: blocks left dark only because of CLPU are picked up
            #    later, one step at a time, as the CLPU of restored load decays
            if self.p["clpu"] and self.mode == "loadflow":
                self._staged_pickup(B, F, descF, cand, dead, set(nonint | restored),
                                    list(closed), ties, interrupted, restored_at, r,
                                    t_stage, t_end)
        return r, interrupted

    def _damage(self, sector, r):
        key = (sector, round(r, 9))
        cache = self.__dict__.setdefault("_dcache", {})
        if key not in cache:
            cache[key] = _ecost.damage(sector, r, table=self.p["scdf"],
                                       factor=self.p["cost_factor"])
        return cache[key]

    def _staged_pickup(self, B, F, descF, cand, dead, energized, closed, ties,
                       interrupted, restored_at, r, t_stage, t_end):
        net, step = self.net, self.p["clpu"]["step_h"]
        pend = {B.blk[v] for v in cand if v not in energized} | \
               {b for b in dead if not any(v in energized for v in B.buses[b])}
        pend.discard(F)
        if not pend:
            return
        dist = self._distances(net.root, set(range(net.n)))
        t = t_stage + step
        while pend and t < t_end - 1e-9:
            self._ctx = dict(interrupted=interrupted, restored_at=restored_at, t_stage=t)
            progress = True
            while progress and pend:
                progress = False
                for b in sorted(pend, key=lambda b: min(dist.get(v, 1e9) for v in B.buses[b])):
                    bb = set(B.buses[b])
                    nbrs = B.kids[b] + [B.par[b]]
                    options = []
                    if any(n in B.buses and any(v in energized for v in B.buses[n])
                           for n in nbrs if n != -1 and n != F):
                        options.append(None)
                    for c in ties:
                        if c in closed:
                            continue
                        a, z = net.ties[c]["a"], net.ties[c]["b"]
                        if (a in bb and z in energized) or (z in bb and a in energized):
                            options.append(c)
                    for opt in options:
                        new_en = energized | bb
                        new_cl = closed + ([opt] if opt is not None else [])
                        if self._feasible(new_en, new_cl):
                            energized, closed = new_en, new_cl
                            for v in bb:
                                restored_at.setdefault(v, t)
                                r[v] = min(r[v], t)
                            pend.discard(b)
                            progress = True
                            break
                    if progress:
                        break
            t += step

    def _full(self, x):
        if not self.net.fixed:
            return list(x)
        xf = list(x)
        for sec, dev in self.net.fixed.items():
            xf[sec] = dev
        return xf

    def _stage_blocks(self, x):
        remote = [d if d in (RCS, REC) else NONE for d in x]
        return {"remote": Blocks(self.net, remote), "all": Blocks(self.net, x)}

    # -- full evaluation -----------------------------------------------------
    def evaluate(self, x, w, severity=False):
        """severity=True also returns, per load point, {duration_min: annual rate}
        (weighted by load-level probability) for the compensation objective f3."""
        net, p = self.net, self.p
        x_dec = list(x)
        x = self._full(x)
        blocks = self._stage_blocks(x)
        levels = p["load_levels"] or [(p["load_level"], 1.0)]
        lam_j = np.zeros(net.n)
        U_j = np.zeros(net.n)
        ens = 0.0
        ecost_total = 0.0
        Lpk = net.P                                    # SCDF are per kW of peak demand
        use_cost = all(net.sector[j] is not None for j in net.load_buses)
        first = True
        sev = {j: {} for j in net.load_buses} if severity else None
        for mult, prob in levels:
            self._mult = mult
            U_lv = np.zeros(net.n)
            for k in range(len(net.sec)):
                r, intr = self.fault_durations(x, w, k, blocks)
                if severity:
                    for j in net.load_buses:
                        if r[j] > 0:
                            dm = int(round(r[j] * 60))
                            sev[j][dm] = sev[j].get(dm, 0.0) + self.lam[k] * prob
                if first:
                    lam_j += self.lam[k] * intr
                U_lv += self.lam[k] * r
                if use_cost:
                    for j in net.load_buses:
                        if r[j] > 0:
                            ecost_total += prob * self.lam[k] * Lpk[j] * \
                                self._damage(net.sector[j], r[j])
            first = False
            U_j += prob * U_lv
            ens += prob * float((net.P * mult * U_lv)[net.load_buses].sum()) / 1000.0
        loads = net.load_buses
        N = net.N[loads]
        saifi = float((lam_j[loads] * N).sum() / N.sum())
        saidi = float((U_j[loads] * N).sum() / N.sum())
        return dict(SAIFI=saifi, SAIDI=saidi, CAIDI=saidi / saifi if saifi else 0.0,
                    ENS=ens, ASAI=1 - saidi / 8760.0, f1=self.cost(x_dec, w),
                    ECOST=ecost_total if use_cost else None,
                    f2=ecost_total if use_cost else None,
                    rec_ok=self.rec_series_ok(x_dec), lam_j=lam_j, U_j=U_j, sev=sev)

    # -- objective f1 and constraint -----------------------------------------
    def cost(self, x, w):
        p, c = self.p, self.p["cost"]
        name = {MS: "MS", RCS: "RCS", REC: "REC"}
        ann_sw = crf(p["delta"], p["h_sw"]) + p["rho_om"]
        cap = sum(c[name[d]] for d in x if d != NONE)
        cap += sum(c["MS"] if wc == 1 else c["RCS"] for wc in w if wc > 0)
        line = sum(c["line_per_km"] * self.net.ties[i]["length_km"]
                   for i, wc in enumerate(w) if wc > 0)
        return ann_sw * cap + crf(p["delta"], p["h_line"]) * line

    def rec_series_ok(self, x):
        net, cnt = self.net, [0] * self.net.n
        for v in net.order[1:]:
            cnt[v] = cnt[net.parent[v]] + (1 if x[net.sec_of_bus[v]] == REC else 0)
        return max(cnt) <= self.p["m_rec"]
