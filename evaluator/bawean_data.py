"""Sistem Bawean (sistem uji ke-4) dari ekspor PowerFactory (bawean_export.xlsx).

Pemodelan untuk evaluator radial:
  * busbar 20 kV PLTD (Bus2, Bus3, Bus9, Bus23, digabung kopel tertutup) = satu akar;
    pembangkit (PLTD/PLTMG) dianggap sebagai sumber di akar;
  * saluran 20 kV -> seksi (panjang, R1, X1, Inom dari ElmLne);
  * kopel/LBS/recloser eksisting -> seksi panjang nol dengan perangkat TETAP:
    CB penyulang -> pemutus kepala penyulang, REC -> recloser, LBSM -> RCS (motorised),
    LBS -> MS; LBS KOPEL TAMBAK (normal terbuka) -> tie eksisting (manual);
  * beban 0,4 kV dipindah ke terminal 20 kV trafo distribusinya (trafo lossless di model PF).
Pelanggan: 1 pelanggan R-1 per 2,2 kVA beban (konfirmasi pengguna); data keandalan: default PLN (SPLN 59:1985).
Tie eksisting: LBS KOPEL TAMBAK di ujung Line193 (Tambak) ke Bus371 (Sangkapura). Line80 tidak tersambung -> diabaikan.
"""
import collections, os
import numpy as np, pandas as pd
_HERE = os.path.dirname(os.path.abspath(__file__))

ROOT_BUSES = {"Bus2", "Bus3", "Bus9", "Bus23"}
FEEDER_CB = {"CB KOTA", "CB TAMBAK", "CB SANGKAPURA"}


def _sheet(x, n):
    return x.parse(n).iloc[1:].reset_index(drop=True)


def load_raw(path=None):
    x = pd.ExcelFile(path or os.path.join(_HERE, "bawean_export.xlsx"))
    return dict(topo=x.parse("Topology"), lne=_sheet(x, "ElmLne"), lod=_sheet(x, "ElmLod"),
                nodes=x.parse("LF_Nodes"), el=x.parse("LF_Elements"), typtr=_sheet(x, "TypTr2"))


def build(path=None, kw_per_customer=None):
    raw = load_raw(path); topo, lne, lod = raw["topo"], raw["lne"], raw["lod"]
    for c in ("dline", "R1", "X1", "Inom"):
        lne[c] = pd.to_numeric(lne[c], errors="coerce")
    L = lne.set_index("loc_name")
    # MV graph
    mvkv = {}
    for _, r in topo.iterrows():
        for s in (1, 2):
            t, kv = r.get(f"terminal_{s}"), r.get(f"terminal_{s}_kV")
            if isinstance(t, str) and pd.notna(kv):
                mvkv[t] = float(kv)
    is_mv = lambda t: mvkv.get(t, 0) > 1.0
    alias = {t: "ROOT" if t in ROOT_BUSES else t for t in mvkv}
    adj = collections.defaultdict(list); ties = []
    skipped = []
    for _, r in topo.iterrows():
        cl, a, b = r["class"], r["terminal_1"], r["terminal_2"]
        if cl in ("ElmLne", "ElmCoup") and not (isinstance(a, str) and isinstance(b, str)):
            skipped.append((cl, r.element, a, b, r.get("n_connections"))); continue
        if cl == "ElmLne" and r["outserv"] != 1:
            ln = L.loc[r.element]
            adj[alias[a]].append((alias[b], dict(kind="line", name=r.element, R=float(ln.R1), X=float(ln.X1),
                                                 length=float(ln.dline), amp=1000 * float(ln.Inom))))
            adj[alias[b]].append((alias[a], adj[alias[a]][-1][1]))
        elif cl == "ElmCoup" and is_mv(a) and is_mv(b):
            if alias[a] == alias[b]:
                continue                                      # bus coupler inside the PLTD busbar
            if r.on_off == 0:
                ties.append((alias[a], alias[b], r.element)); continue
            name = r.element
            dev = ("REC" if name in FEEDER_CB else "REC" if name.startswith("REC") else
                   "RCS" if name.startswith("LBSM") else "MS")
            e = dict(kind="device", name=name, R=1e-5, X=1e-5, length=0.0, amp=1e4, dev=dev)
            adj[alias[a]].append((alias[b], e)); adj[alias[b]].append((alias[a], e))
    # orient from ROOT (BFS)
    order, par, pedge = ["ROOT"], {"ROOT": None}, {}
    q = collections.deque(["ROOT"])
    while q:
        u = q.popleft()
        for v, e in adj[u]:
            if v not in par:
                par[v] = u; pedge[v] = e; order.append(v); q.append(v)
    idx = {t: i for i, t in enumerate(order)}
    # loads: LV terminal -> MV terminal of its distribution transformer
    lv2mv = {}
    for _, r in topo[topo["class"] == "ElmTr2"].iterrows():
        a, b = r["terminal_1"], r["terminal_2"]
        if is_mv(a) and not is_mv(b): lv2mv[b] = alias[a]
        elif is_mv(b) and not is_mv(a): lv2mv[a] = alias[b]
    P = np.zeros(len(order)); Q = np.zeros(len(order)); nload = collections.Counter()
    S = np.zeros(len(order))
    lodp = lod.set_index("loc_name")
    for _, r in topo[topo["class"] == "ElmLod"].iterrows():
        mv = lv2mv.get(r["terminal_1"], alias.get(r["terminal_1"]))
        j = idx[mv]
        pl, ql = 1000 * float(lodp.loc[r.element, "plini"]), 1000 * float(lodp.loc[r.element, "qlini"])
        P[j] += pl; Q[j] += ql; S[j] += float(np.hypot(pl, ql)); nload[mv] += 1
    # reactive losses of the distribution transformers (lossless in R, as in the PF model)
    typ = raw["typtr"].set_index("loc_name")
    load_S_lv = collections.defaultdict(float)
    for _, r in topo[topo["class"] == "ElmLod"].iterrows():
        load_S_lv[r["terminal_1"]] += 1000 * float(np.hypot(float(lodp.loc[r.element, "plini"]),
                                                              float(lodp.loc[r.element, "qlini"])))
    for _, r in topo[topo["class"] == "ElmTr2"].iterrows():
        tname = str(r["typ_id"]).split("\\")[-1].replace(".TypTr2", "")
        if tname.startswith("GT") or tname not in typ.index:
            continue
        lv = r["terminal_2"] if not is_mv(r["terminal_2"]) else r["terminal_1"]
        mv = lv2mv.get(lv)
        if mv is None:
            continue
        srat = 1000 * float(typ.loc[tname, "strn"]); uk = float(typ.loc[tname, "uktr"]) / 100
        Q[idx[mv]] += uk * load_S_lv[lv] ** 2 / srat
    sections, fixed, labels, heads = [], {}, {}, []
    for v in order[1:]:
        e = pedge[v]
        sections.append(dict(frm=idx[par[v]], to=idx[v], R=e["R"], X=e["X"], length_km=e["length"],
                             ampacity_a=e["amp"]))
        labels[len(sections) - 1] = e["name"]
        if e["kind"] == "device":
            fixed[len(sections) - 1] = e["dev"]
            if e["name"] in FEEDER_CB:
                heads.append(len(sections) - 1)
    tie_list = [dict(a=idx[a], b=idx[b], R=1e-5, X=1e-5, length_km=0.0, ampacity_a=1e4) for a, b, _ in ties]
    # one residential customer (R-1, 2.2 kVA) per 2.2 kVA of load (confirmed by the user)
    customers = [max(1, round(sv / 2.2)) if sv > 0 else 0 for sv in S]
    sectors = ["residential" if sv > 0 else None for sv in S]
    net = dict(n_bus=len(order), sections=sections, ties=tie_list, P_kw=list(P), Q_kvar=list(Q),
               customers=customers, root=0, V_kv=20.0, sectors=sectors)
    unreached = [t for t in adj if t not in par]
    info = dict(order=order, fixed=fixed, labels=labels, ties=ties, raw=raw, skipped=skipped,
                unreached=unreached, heads=heads, contract_kva=2.2)
    return net, info
