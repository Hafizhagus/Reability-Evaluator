#!/usr/bin/env python3
"""Independent check of the Bawean diagram (figures/fig_bawean.svg, Fig. 2 of the paper) against the model.

The network that the drawing shows is rebuilt from the geometry of the file alone (which line ends on which
symbol, which symbols are squares, dots, bars, hexagons, triangles) and compared with the model of
evaluator/bawean_data.py: tree structure of every feeder, load buses, branching buses, type and position of
every device, bus names, tie ends, symbol spacing and lettering size.  It shares no code with
make_bawean_sld.py.

Usage:  python3 check_bawean_sld.py [figures/fig_bawean.svg] [--evaluator ../evaluator]
"""
import argparse, collections, math, os, re, sys
import xml.etree.ElementTree as ET

ap = argparse.ArgumentParser()
here = os.path.dirname(os.path.abspath(__file__))
ap.add_argument("svg", nargs="?", default=os.path.join(here, "figures", "fig_bawean.svg"))
ap.add_argument("--evaluator", default=os.path.join(here, "..", "evaluator"))
args = ap.parse_args()
sys.path.insert(0, os.path.abspath(args.evaluator))
import bawean_data as B

ok = True
def check(cond, msg):
    global ok
    print(("  ok    " if cond else "  GAGAL ") + msg); ok = ok and bool(cond)


# ----------------------------------------------------------------------------- contracted tree (shared form)
def contract(root, children, load_of, dev_in, name_of):
    """Keep the root, load buses and branching buses; drop dead ends without load; devices stay on their edge."""
    def rec(u, devs):
        kids = [k for k in (rec(c, list(dev_in.get(c, []))) for c in children(u)) if k is not None]
        ld = load_of(u)
        if ld == 0 and not kids and u != root:
            return None
        if ld == 0 and len(kids) == 1 and u != root:
            kids[0]["devs"] = devs + kids[0]["devs"]; return kids[0]
        return dict(name=name_of(u), load=ld, devs=devs, kids=kids)
    return rec(root, [])

def canon(t):
    return "(" + "".join(d[0] + d[-1] for d in t["devs"]) + ("L" if t["load"] else "J") + "".join(sorted(canon(k) for k in t["kids"])) + ")"

def count(t, what):
    return what(t) + sum(count(k, what) for k in t["kids"])


# ----------------------------------------------------------------------------- model
kw, info = B.build()
secs, order, fixed, labels, topo = kw["sections"], info["order"], info["fixed"], info["labels"], info["raw"]["topo"]
idx = {t: i for i, t in enumerate(order)}
kv = {}
for _, r in topo.iterrows():
    for s in (1, 2):
        t, v = r.get(f"terminal_{s}"), r.get(f"terminal_{s}_kV")
        if isinstance(t, str) and v == v: kv[t] = float(v)
alias = lambda t: "ROOT" if t in B.ROOT_BUSES else t
lv2mv = {}
for _, r in topo[topo["class"] == "ElmTr2"].iterrows():
    a, b = r["terminal_1"], r["terminal_2"]
    if kv.get(a, 0) > 1 and kv.get(b, 0) < 1: lv2mv[b] = alias(a)
    elif kv.get(b, 0) > 1 and kv.get(a, 0) < 1: lv2mv[a] = alias(b)
nload = collections.Counter()
for _, r in topo[topo["class"] == "ElmLod"].iterrows():
    nload[idx[lv2mv.get(r["terminal_1"], alias(r["terminal_1"]))]] += 1
mch, mpar, sec_to = collections.defaultdict(list), {}, {}
for k, s in enumerate(secs):
    mch[s["frm"]].append(s["to"]); mpar[s["to"]] = s["frm"]; sec_to[s["to"]] = k
heads = {labels[h]: secs[h]["to"] for h in info["heads"]}
head_secs = set(info["heads"])
kind_of = lambda nm: "REC" if nm.startswith("REC") else "LBSM" if nm.startswith("LBSM") else "LBS"
mdev = {v: [kind_of(labels[k])] for v, k in sec_to.items() if k in fixed and k not in head_secs}
model = {nm: contract(h, lambda u: mch[u], lambda u: nload[u], mdev, lambda u: order[u]) for nm, h in heads.items()}
(ta, tb, tname), = info["ties"]

# ----------------------------------------------------------------------------- drawing
NS = "{http://www.w3.org/2000/svg}"
root = ET.fromstring(open(args.svg, encoding="utf-8").read())
mm_w, mm_h = float(root.get("width").replace("mm", "")), float(root.get("height").replace("mm", ""))
S = float(root.get("viewBox").split()[2]) / mm_w
num = lambda v: float(v) / S
lines, squares, dots, devs, texts, tie_line, tie_sw, bars = [], [], [], [], [], [], [], []
def walk(el, gid, in_legend):
    tag = el.tag.replace(NS, "")
    if tag in ("g", "svg"):
        for c in el: walk(c, el.get("id") or gid, in_legend or el.get("id") == "legend")
        return
    if tag == "text":
        texts.append(num(el.get("font-size"))); return
    if in_legend or el.get("id") == "background":
        return
    if tag == "path":
        pts = [(num(a), num(b)) for a, b in re.findall(r"(-?[\d.]+) (-?[\d.]+)", el.get("d"))]
        (tie_line if el.get("stroke-dasharray") else lines).append(dict(pts=pts, color=el.get("stroke"), w=num(el.get("stroke-width"))))
    elif tag == "rect":
        x, y, w, h = (num(el.get(k)) for k in ("x", "y", "width", "height"))
        it = dict(p=(x + w / 2, y + h / 2), gid=gid)
        if w > 100: bars.append(it)
        elif abs(w - h) < 1e-6: squares.append(it)
        else: devs.append(dict(it, kind="LBSM" if el.get("fill").upper() == "#000000" else "LBS"))
    elif tag == "circle":
        dots.append(dict(p=(num(el.get("cx")), num(el.get("cy"))), gid=gid))
    elif tag == "polygon":
        P = [tuple(num(v) for v in q.split(",")) for q in el.get("points").split()]
        devs.append(dict(p=(sum(p[0] for p in P) / len(P), (min(p[1] for p in P) + max(p[1] for p in P)) / 2),
                         kind={6: "REC", 3: "CB"}[len(P)], gid=gid))
    elif tag == "ellipse":
        tie_sw.append((num(el.get("cx")), num(el.get("cy"))))
walk(root, None, False)

eps = 0.05
mangle = lambda s: s.replace("(", "_").replace(")", "").replace(" ", "_")
unmangle = {mangle(t): t for t in order}; unmangle.update({mangle(l): l for l in labels.values()})
nodes = []
for d in dots:
    check(not any(math.dist(d["p"], n["p"]) < eps for n in nodes), "titik bus ganda di " + str(d["p"])) if any(math.dist(d["p"], n["p"]) < eps for n in nodes) else None
    nodes.append(dict(p=d["p"], load=any(math.dist(q["p"], d["p"]) < eps for q in squares), name=unmangle.get(d["gid"], d["gid"])))
check(len(squares) == sum(n["load"] for n in nodes), f"setiap kotak beban berisi tepat satu titik bus ({len(squares)} kotak)")
roots = {}
for ln in lines:
    if abs(ln["pts"][0][1] - bars[0]["p"][1]) < eps:
        nodes.append(dict(p=ln["pts"][0], load=False, name="ROOT")); roots[ln["color"]] = len(nodes) - 1
check(len(roots) == 3, f"tiga penyulang berangkat dari busbar sumber ({len(roots)})")

def on_seg(a, b, p):
    L = math.dist(a, b)
    t = ((p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1])) / L
    d = abs((p[0] - a[0]) * (b[1] - a[1]) - (p[1] - a[1]) * (b[0] - a[0])) / L
    return t if d < eps and -eps <= t <= L + eps else None

adj, edge_dev, used, clean, pitch = collections.defaultdict(set), collections.defaultdict(list), set(), True, collections.Counter()
for ln in lines:
    hits, s0 = set(), 0.0
    for a, b in zip(ln["pts"][:-1], ln["pts"][1:]):
        horizontal = abs(a[1] - b[1]) < eps
        for k, n in enumerate(nodes):
            t = on_seg(a, b, n["p"])
            if t is not None: hits.add((round(s0 + t, 3), "n", k, horizontal))
        for j, d in enumerate(devs):
            t = on_seg(a, b, d["p"])
            if t is not None: hits.add((round(s0 + t, 3), "d", j, horizontal))
        s0 += math.dist(a, b)
        if (round(s0, 3), "n") not in {(h[0], h[1]) for h in hits} and (round(s0, 3), "d") not in {(h[0], h[1]) for h in hits}:
            hits.add((round(s0, 3), "c", -1, horizontal))                      # corner of the line
    hits = sorted(hits)
    ks = [h for h in hits if h[1] == "n"]
    if not (len(ks) == 2 and ks[0][0] < eps and abs(ks[1][0] - s0) < eps):
        clean = False; continue
    u, v = ks[0][2], ks[1][2]
    adj[u].add(v); adj[v].add(u)
    for h in hits:
        if h[1] == "d":
            edge_dev[frozenset((u, v))].append((h[0], devs[h[2]]["kind"], unmangle.get(devs[h[2]]["gid"], devs[h[2]]["gid"]), u)); used.add(h[2])
    for h1, h2 in zip(hits[:-1], hits[1:]):                                    # spacing between neighbouring symbols/corners
        if nodes[u]["name"] != "ROOT" or h1[0] > eps:
            pitch[(round(h2[0] - h1[0], 2), "mendatar" if h2[3] else "tegak")] += 1
    for k in (u, v):
        nodes[k].setdefault("colors", set()).add(ln["color"])
check(clean, "setiap garis berawal dan berakhir tepat pada simbol bus, tanpa bus lain di tengahnya")
check(len(used) == len(devs), f"semua {len(devs)} simbol perangkat berada tepat di atas sebuah garis")
check(all(len(n.get("colors", ())) == 1 for n in nodes), "tidak ada bus yang tersambung ke garis dua penyulang")
par = {}
for r in roots.values():
    par[r] = None; st = [r]
    while st:
        u = st.pop()
        for v in adj[u]:
            if v not in par: par[v] = u; st.append(v)
check(len(par) == len(nodes), f"semua {len(nodes) - 3} simbol bus tersambung ke sumber")
check(sum(len(v) for v in adj.values()) // 2 == len(nodes) - 3, "jaringan yang digambar radial (tanpa loop)")
ch = collections.defaultdict(list)
for v, u in par.items():
    if u is not None: ch[u].append(v)
dev_in, dev_names = {}, {}
for v, u in par.items():
    if u is None: continue
    ds = sorted(edge_dev[frozenset((u, v))], key=lambda d: d[0] if d[3] == u else -d[0])
    dev_in[v] = [d[1] for d in ds if d[1] != "CB"]; dev_names[v] = [d[2] for d in ds if d[1] != "CB"]
drawn = [contract(r, lambda u: ch[u], lambda u: int(nodes[u]["load"]), dev_in, lambda u: nodes[u]["name"]) for r in roots.values()]
dcanon = {canon(t) for t in drawn}
print("\nStruktur per penyulang (bus beban, percabangan, dan perangkat pada tiap ruas):")
for nm, t in model.items():
    check(canon(t) in dcanon, f"{nm}: gambar identik dengan model ({count(t, lambda x: int(x['load'] > 0))} bus beban, "
                              f"{count(t, lambda x: x['load'])} beban, {count(t, lambda x: len(x['devs']))} perangkat)")
# names: the bus drawn under a name has the same parent and the same devices as that bus in the model
name2node = {n["name"]: k for k, n in enumerate(nodes)}
bad = []
for nm, k in name2node.items():
    if nm == "ROOT": continue
    v, dv = idx[nm], []
    while True:
        kk = sec_to[v]
        if kk in fixed and kk not in head_secs: dv.append(labels[kk])
        v = mpar[v]
        if order[v] in name2node or v in heads.values(): break
    exp = "ROOT" if v in heads.values() else order[v]
    if exp != nodes[par[k]]["name"] or dv[::-1] != dev_names[k] or (nload[idx[nm]] > 0) != nodes[k]["load"]:
        bad.append(nm)
check(not bad, "nama: setiap bus pada gambar berinduk, berperangkat, dan berstatus beban sama dengan model " + (str(bad[:6]) if bad else ""))
cnt = collections.Counter(d["kind"] for d in devs)
check(cnt == {"CB": 3, "REC": 3, "LBSM": 4, "LBS": 7}, f"perangkat: {cnt['CB']} CB, {cnt['REC']} recloser, {cnt['LBSM']} LBSM, {cnt['LBS']} LBS")
check(len(squares) == len(nload), f"kotak beban {len(squares)} = bus berbeban pada model {len(nload)} (beban PowerFactory {sum(nload.values())})")
def up(v):
    while order[v] not in name2node: v = mpar[v]
    return order[v]
tp = tie_line[0]["pts"]
at = lambda p: next((n["name"] for n in nodes if math.dist(n["p"], p) < eps), None)
check({at(tp[0]), at(tp[-1])} == {up(idx[ta]), up(idx[tb])}, f"tie {tname}: {at(tp[0])} <-> {at(tp[-1])} (model: {tb} <-> {ta})")
check(len(tie_sw) == 1 and any(on_seg(a, b, tie_sw[0]) is not None for a, b in zip(tp[:-1], tp[1:])), "simbol tie berada pada garis tie")
P = [n["p"] for n in nodes if n["name"] != "ROOT"] + [d["p"] for d in devs]
dmin = min(math.dist(a, b) for i, a in enumerate(P) for b in P[i + 1:])
print("\nTata letak:")
check(len(pitch) == 2, "jarak antar-simbol seragam: " + ", ".join(f"{d:g} mm {o} ({n}x)" for (d, o), n in sorted(pitch.items())))
check(dmin > 4.0, f"tidak ada simbol bertumpuk (jarak pusat terkecil {dmin:.2f} mm)")
check(all(abs(f * 72 / 25.4 - 8) < 0.01 for f in texts), f"semua tulisan 8 pt pada ukuran akhir {mm_w:g} x {mm_h:.1f} mm")
check(min(l["w"] for l in lines + tie_line) >= 0.1, f"tebal garis {min(l['w'] for l in lines + tie_line):.2f}-{max(l['w'] for l in lines):.2f} mm (minimum jurnal 0,1 mm)")
print("\nHASIL:", "SEMUA PEMERIKSAAN LOLOS" if ok else "ADA YANG GAGAL")
sys.exit(0 if ok else 1)
