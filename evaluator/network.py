"""Radial network topology (step 2a) and backward/forward sweep load flow (step 2c)."""
from collections import deque
import math
import numpy as np

NONE, MS, RCS, REC = 0, 1, 2, 3
DEVICE_NAME = {NONE: "-", MS: "MS", RCS: "RCS", REC: "REC"}


class Network:
    """Radial distribution network rooted at the source bus.

    Section s connects parent bus frm(s) to child bus to(s). A device installed at
    section s sits at its sending (parent) end, so its 'child side' is the section
    line itself plus the subtree of to(s).
    """

    def __init__(self, n_bus, sections, ties, P_kw, Q_kvar, customers, root=0,
                 V_kv=12.66, sectors=None, fixed=None, s_max_kva=None):
        self.n, self.root, self.V_kv = n_bus, root, V_kv
        self.sec, self.ties = sections, ties
        self.P = np.asarray(P_kw, float)
        self.Q = np.asarray(Q_kvar, float)
        self.N = np.asarray(customers, float)
        self.sector = list(sectors) if sectors is not None else [None] * n_bus
        self.fixed = dict(fixed) if fixed else {}      # section -> device (not a decision)
        self.s_max_kva = s_max_kva                     # source capacity (None = unlimited)
        self.parent = [-1] * n_bus
        self.sec_of_bus = [-1] * n_bus
        self.children = [[] for _ in range(n_bus)]
        for s, d in enumerate(sections):
            if self.parent[d["to"]] != -1:
                raise ValueError(f"bus {d['to']} has two parents: not radial")
            self.parent[d["to"]] = d["frm"]
            self.sec_of_bus[d["to"]] = s
            self.children[d["frm"]].append(d["to"])
        # topological (BFS) order from the root
        self.order, q = [], deque([root])
        while q:
            v = q.popleft()
            self.order.append(v)
            q.extend(self.children[v])
        if len(self.order) != n_bus:
            raise ValueError("network is not connected to the root")
        self.load_buses = [i for i in range(n_bus) if self.P[i] > 0 or self.Q[i] > 0]

    # -- helpers -------------------------------------------------------------
    def zabs(self, s):
        d = self.sec[s]
        return math.hypot(d["R"], d["X"])


# ---------------------------------------------------------------------------
# Load flow
# ---------------------------------------------------------------------------
def bfs_loadflow(net, energized, edges, s_base_kva=1000.0, tol=1e-9, max_iter=100,
                 load_scale=1.0):
    """Backward/forward sweep for a radial energized sub-network.

    energized : iterable of bus indices that are energized (must include root)
    edges     : list of (u, v, R_ohm, X_ohm, ampacity_a) closed branches among them
    Returns dict(V=|V| p.u. per bus (nan if dead), I_a={edge_index: amps},
                 loss_kw, converged, radial).
    """
    energized = set(energized)
    z_base = net.V_kv ** 2 / (s_base_kva / 1000.0)            # ohm
    i_base = s_base_kva / (math.sqrt(3) * net.V_kv)            # A
    adj = {b: [] for b in energized}
    for e, (u, v, R, X, _) in enumerate(edges):
        adj[u].append((v, e))
        adj[v].append((u, e))
    # orient the tree from the root
    par, par_edge, order = {net.root: None}, {net.root: None}, []
    q = deque([net.root])
    while q:
        u = q.popleft()
        order.append(u)
        for v, e in adj[u]:
            if v not in par:
                par[v], par_edge[v] = u, e
                q.append(v)
    radial = (len(order) == len(energized)) and (len(edges) == len(energized) - 1)
    if not radial:
        return dict(V=None, I_a=None, loss_kw=None, converged=False, radial=False)
    if isinstance(load_scale, dict):
        s_load = {b: load_scale.get(b, 1.0) * (net.P[b] + 1j * net.Q[b]) / s_base_kva
                  for b in order}
    else:
        s_load = {b: load_scale * (net.P[b] + 1j * net.Q[b]) / s_base_kva for b in order}
    z = {e: (R + 1j * X) / z_base for e, (_, _, R, X, _) in enumerate(edges)}
    V = {b: 1.0 + 0j for b in order}
    converged = False
    for _ in range(max_iter):
        I_br = {}
        inj = {b: np.conj(s_load[b] / V[b]) for b in order}
        acc = dict(inj)
        for b in reversed(order[1:]):                 # backward sweep
            e = par_edge[b]
            I_br[e] = acc[b]
            acc[par[b]] += acc[b]
        V_new = {net.root: 1.0 + 0j}
        for b in order[1:]:                            # forward sweep
            V_new[b] = V_new[par[b]] - z[par_edge[b]] * I_br[par_edge[b]]
        dv = max(abs(V_new[b] - V[b]) for b in order)
        V = V_new
        if dv < tol:
            converged = True
            break
    Vmag = np.full(net.n, np.nan)
    for b in order:
        Vmag[b] = abs(V[b])
    I_a = {e: abs(I_br[e]) * i_base for e in I_br}
    loss_kw = sum(abs(I_br[e]) ** 2 * z[e].real for e in I_br) * s_base_kva
    I_root = sum(I_br[par_edge[b]] for b in order[1:] if par[b] == net.root)
    S_root_kva = abs(np.conj(I_root)) * s_base_kva if order[1:] else 0.0
    return dict(V=Vmag, I_a=I_a, loss_kw=loss_kw, converged=converged, radial=True,
                S_root_kva=S_root_kva)


def normal_edges(net):
    """All sections closed, all ties open."""
    return [(d["frm"], d["to"], d["R"], d["X"], d["ampacity_a"]) for d in net.sec]
