#!/usr/bin/env python3
"""Fig. 2 - single-line diagram of the Bawean 20 kV network.

The diagram is generated from the network model (evaluator/bawean_data.py, built from the PowerFactory export),
so the drawn topology cannot differ from the data: every load bus and every branching bus of the model must
have a grid cell in POS, the lines and the devices on them are taken from the model, and the script stops if a
bus is missing, two symbols share a cell, or two lines touch.  Layout and symbols follow the hand-drawn
diagram of the first author (source bar on top, Tambak on the left and bottom, Kota in the middle,
Sangkapura on the right).

Output (in --out): fig_bawean.svg (editable), fig_bawean.pdf, fig_bawean.eps (fonts embedded), fig_bawean.png.
Usage:  python3 make_bawean_sld.py [--evaluator ../evaluator] [--out figures]
Check:  python3 check_bawean_sld.py   (rebuilds the drawn network from the SVG and compares it with the model)
Needs:  pandas, openpyxl (for the model), cairosvg (for PDF/EPS/PNG), pillow (text widths).
"""
import argparse, collections, math, os, sys

# ----------------------------------------------------------------------------- geometry (mm, final size)
W = 174.0                    # figure width (two-column width of the journal)
NCOL, NROW = 25, 25          # last grid column / last grid row
MX = 2.0                     # centre of column 0
DX = (W - 2 * MX) / NCOL     # column pitch
DY = 4.3                     # row pitch
YBAR = 4.7                   # centre of the source bar (row -1); row 0 = feeder breaker, rows 1.. = buses
H = YBAR + (NROW + 1) * DY + 2.2
SQ, RDOT, LW = 3.1, 1.05, 0.45              # load square, bus dot radius, line width
BAR_L, BAR_T, BAR_S = 3.6, 1.2, 0.25        # switch bar: length, thickness, outline
RHEX = 1.7                                   # recloser hexagon (circumradius, pointy top)
TRI_W, TRI_H = 3.4, 3.1                      # feeder breaker
SRC_T = 1.4                                  # thickness of the source bar
TIE_LW, TIE_DASH = 0.35, (1.0, 0.8)
ELL = (2.1, 1.05, 0.3)                       # tie switch: rx, ry, outline
PT = 25.4 / 72
FS = 8 * PT                                  # 8 pt lettering
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"
INK = "#000000"
C_LOAD_FILL, C_LOAD_EDGE, LOAD_S = "#5FA3AE", "#0E5661", 0.3   # square: fill, outline, outline width
COL = {"CB TAMBAK": "#0064FA", "CB KOTA": "#1FA637", "CB SANGKAPURA": "#580ACE"}
FEEDER_LABEL = {"CB TAMBAK": "Feeder 1 (Tambak)", "CB KOTA": "Feeder 2 (Kota)", "CB SANGKAPURA": "Feeder 3 (Sangkapura)"}
SRC_LABEL = "Power plant, 20 kV busbar"
X = lambda c: MX + c * DX
Y = lambda r: YBAR + (r + 1) * DY

# ----------------------------------------------------------------------------- layout: bus -> (column, row)
HEAD_COL = {"CB TAMBAK": 2, "CB KOTA": 11, "CB SANGKAPURA": 19}
POS = {
    # Feeder 1 (Tambak): down the left side, then to the right along the bottom
    "Bus184": (2, 1), "Bus14": (2, 2), "Bus229": (1, 2),
    "Bus28": (2, 4), "Bus36": (2, 5), "Bus231": (1, 5), "Bus234": (2, 6),
    "Bus235": (2, 7), "Bus74": (1, 7), "Bus78": (0, 7), "Bus238": (2, 8),
    "Bus39": (2, 10), "Bus82": (2, 11), "Bus240": (2, 12),
    "Bus91": (3, 12), "Bus93": (4, 12), "Bus241": (4, 11),
    "Bus244": (5, 12), "Bus247": (6, 12), "Bus250": (7, 12), "Bus254": (8, 12), "Bus257": (9, 12),
    "Bus84": (2, 14), "Bus86": (2, 15), "Bus262": (1, 15), "Bus88": (2, 16),
    "Bus95": (2, 18), "Bus97": (2, 19), "Bus264": (2, 20), "Terminal(16)": (1, 20),
    "Bus99": (2, 21), "Bus101": (2, 22), "Bus270": (2, 23),
    "Bus105": (3, 23), "Terminal(17)": (3, 22), "Bus276": (4, 23), "Bus279": (5, 23),
    "Bus103": (2, 24),
    "Bus107": (4, 25), "Bus110": (5, 25), "Bus112": (6, 25), "Bus280": (6, 24),
    "Bus284": (7, 25), "Bus286": (8, 25),
    "Bus289": (10, 25), "Bus115": (10, 24), "Bus117": (10, 23), "Bus290": (9, 23), "Bus119": (10, 22),
    "Bus121": (11, 23), "Bus122": (12, 23), "Bus294": (13, 23), "Bus298": (14, 23),
    "Bus124": (11, 25), "Bus126": (13, 25), "Bus128": (14, 25), "Bus130": (15, 25),
    # Feeder 2 (Kota)
    "Bus25": (11, 1), "Bus188": (12, 1), "Bus192": (13, 1), "Bus27": (11, 2), "Bus31": (11, 3),
    "Bus194": (11, 5), "Bus41": (11, 6), "Bus198": (11, 7), "Bus43": (11, 8),
    "Bus49": (12, 8), "Bus50": (13, 8), "Bus202": (14, 8), "Bus205": (15, 8), "Bus201": (13, 9),
    "Bus63": (11, 9), "Bus206": (11, 11), "Bus71": (11, 12), "Bus209": (11, 13), "Bus210": (11, 14),
    "Bus216": (12, 14), "Bus219": (13, 14), "Bus220": (14, 14), "Bus223": (15, 14), "Bus227": (14, 15),
    "Bus213": (11, 15),
    # Feeder 3 (Sangkapura): down, to the right, down, then to the left
    "Bus311": (19, 1), "Bus48": (19, 2), "Bus150": (18, 2), "Bus314": (20, 2),
    "Bus155": (19, 4), "Bus316": (20, 4), "Bus319": (21, 4), "Bus157": (19, 5),
    "Bus159": (19, 6), "Bus321": (20, 6), "Bus324": (21, 6), "Bus328": (22, 6), "Bus161": (19, 7),
    "Bus163": (19, 8), "Bus165": (20, 8), "Bus167": (21, 8), "Bus169": (22, 8),
    "Terminal(31)": (18, 8), "Terminal(32)": (17, 8),
    "Bus171": (19, 10), "Bus173": (18, 10), "Bus175": (20, 10), "Bus177": (21, 10), "Bus335": (21, 11),
    "Bus179": (23, 11), "Bus181": (23, 12),
    "Bus183": (23, 14), "Bus132": (24, 14), "Bus338": (25, 14),
    "Bus136": (23, 15), "Bus139": (23, 16), "Bus141": (23, 17), "Bus143": (23, 18),
    "Bus145": (21, 19), "Bus147": (20, 19), "Bus365": (19, 19), "Bus371": (18, 19),
    "Bus370": (17, 19), "Bus369": (16, 19),
    "Bus343": (20, 20), "Terminal(39)": (20, 21), "Terminal(40)": (20, 22), "Terminal(41)": (21, 21),
}
# corner cells of the lines that are not straight (edge = line that ends at the named bus)
VIA = {"Bus107": [(2, 25)], "Bus179": [(23, 10)], "Bus145": [(23, 19)]}
TIE_VIA = [(18, 25)]                         # corner of the tie route (from the Tambak end to the Sangkapura end)
TIE_SWITCH = (18, 22)
LEGEND_BETWEEN, LEGEND_ROW0 = (2, 11), 14    # the legend is centred between these two columns; its first row
LEGEND = [("load", "Load bus"), ("bus", "Bus without load"), ("CB", "Feeder circuit breaker"),
          ("REC", "Recloser"), ("LBSM", "Motorized load-break switch"), ("LBS", "Manual load-break switch"),
          ("tie", "Normally open tie switch")]


# ----------------------------------------------------------------------------- model -> drawn tree
def model_tree(evaluator_dir):
    """Load buses and branching buses of the model, their parent among these, and the devices in between."""
    sys.path.insert(0, evaluator_dir)
    import bawean_data as B
    kw, info = B.build()
    secs, order, fixed, labels, topo = kw["sections"], info["order"], info["fixed"], info["labels"], info["raw"]["topo"]
    idx = {t: i for i, t in enumerate(order)}
    kv = {}
    for _, r in topo.iterrows():
        for s in (1, 2):
            t, v = r.get(f"terminal_{s}"), r.get(f"terminal_{s}_kV")
            if isinstance(t, str) and v == v:
                kv[t] = float(v)
    alias = lambda t: "ROOT" if t in B.ROOT_BUSES else t
    lv2mv = {}
    for _, r in topo[topo["class"] == "ElmTr2"].iterrows():
        a, b = r["terminal_1"], r["terminal_2"]
        if kv.get(a, 0) > 1 and kv.get(b, 0) < 1: lv2mv[b] = alias(a)
        elif kv.get(b, 0) > 1 and kv.get(a, 0) < 1: lv2mv[a] = alias(b)
    nload = collections.Counter()
    for _, r in topo[topo["class"] == "ElmLod"].iterrows():
        nload[idx[lv2mv.get(r["terminal_1"], alias(r["terminal_1"]))]] += 1
    ch, par, sec_to = collections.defaultdict(list), {}, {}
    for k, s in enumerate(secs):
        ch[s["frm"]].append(s["to"]); par[s["to"]] = s["frm"]; sec_to[s["to"]] = k
    (ta, tb, tname), = info["ties"]
    tie_ends = {idx[ta], idx[tb]}
    feeders = {}
    for h in info["heads"]:
        head = secs[h]["to"]
        sub = [head]; st = [head]
        while st:
            u = st.pop()
            for c in ch[u]:
                sub.append(c); st.append(c)
        alive = {}
        for u in reversed(sub):                                  # children come after their parent in `sub`
            alive[u] = nload[u] > 0 or any(alive[c] for c in ch[u])
        live = lambda u: [c for c in ch[u] if alive[c]]
        # a tie end without load that is not a branching point is drawn at its nearest drawn bus upstream,
        # unless it lies between two drawn buses (then it is a real attachment point on the line)
        keep = {u for u in sub if nload[u] > 0 or len(live(u)) >= 2 or (u in tie_ends and alive[u])}
        keep.add(head)
        nodes = {}
        for u in sub:
            if u not in keep or u == head:
                continue
            devs, v = [], u
            while True:
                k = sec_to[v]
                if k in fixed:
                    devs.append(labels[k])
                v = par[v]
                if v in keep:
                    break
            nodes[order[u]] = dict(parent=order[v], devs=devs[::-1], load=nload[u], kids=len(live(u)))
        feeders[labels[h]] = dict(head=order[head], nodes=nodes, nload=sum(nload[u] for u in sub))
    def drawn(u):                                               # nearest drawn bus at or above a tie end
        names = {n for f in feeders.values() for n in f["nodes"]} | {f["head"] for f in feeders.values()}
        while order[u] not in names:
            u = par[u]
        return order[u]
    return feeders, (drawn(idx[tb]), drawn(idx[ta]), tname)


def dev_kind(name):
    return "CB" if name.startswith("CB") else "REC" if name.startswith("REC") else "LBSM" if name.startswith("LBSM") else "LBS"


# ----------------------------------------------------------------------------- scene
class Scene:
    def __init__(self, w, h, scale=20):
        self.w, self.h, self.s, self.items = w, h, scale, []
    def u(self, v):
        return f"{v * self.s:.2f}".rstrip("0").rstrip(".")
    def group(self, name, body):
        self.items.append(f'<g id="{name}">\n' + "\n".join("  " + b for b in body) + "\n</g>")
    def line(self, pts, color, lw, dash=None):
        d = "M" + " L".join(f"{self.u(x)} {self.u(y)}" for x, y in pts)
        extra = f' stroke-dasharray="{self.u(dash[0])} {self.u(dash[1])}"' if dash else ""
        return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{self.u(lw)}" stroke-linejoin="miter"{extra}/>'
    def rect(self, cx, cy, w, h, fill, stroke=None, sw=0.0, cls=None):
        st = f' stroke="{stroke}" stroke-width="{self.u(sw)}"' if stroke else ""
        c = f' class="{cls}"' if cls else ""
        return f'<rect{c} x="{self.u(cx - w / 2)}" y="{self.u(cy - h / 2)}" width="{self.u(w)}" height="{self.u(h)}" fill="{fill}"{st}/>'
    def circle(self, cx, cy, r, fill=INK, cls=None):
        c = f' class="{cls}"' if cls else ""
        return f'<circle{c} cx="{self.u(cx)}" cy="{self.u(cy)}" r="{self.u(r)}" fill="{fill}"/>'
    def poly(self, pts, fill=INK, cls=None):
        c = f' class="{cls}"' if cls else ""
        return f'<polygon{c} points="' + " ".join(f"{self.u(x)},{self.u(y)}" for x, y in pts) + f'" fill="{fill}"/>'
    def ellipse(self, cx, cy, rx, ry, sw, cls=None):
        c = f' class="{cls}"' if cls else ""
        return (f'<ellipse{c} cx="{self.u(cx)}" cy="{self.u(cy)}" rx="{self.u(rx)}" ry="{self.u(ry)}" fill="#FFFFFF" '
                f'stroke="{INK}" stroke-width="{self.u(sw)}"/>')
    def text(self, x, yc, s, anchor="start"):
        """Text whose capital letters are centred on yc."""
        return (f'<text x="{self.u(x)}" y="{self.u(yc + 0.358 * FS)}" font-family="{FONT}" font-size="{self.u(FS)}" '
                f'text-anchor="{anchor}" fill="{INK}">{s}</text>')
    def svg(self):
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w:g}mm" height="{self.h:g}mm" '
                f'viewBox="0 0 {self.u(self.w)} {self.u(self.h)}">\n'
                f'<rect id="background" x="0" y="0" width="{self.u(self.w)}" height="{self.u(self.h)}" fill="#FFFFFF"/>\n'
                + "\n".join(self.items) + "\n</svg>\n")

    # symbols (shared by the network and the legend)
    def sym(self, kind, x, y, horizontal_line=False):
        if kind == "load":
            return [self.rect(x, y, SQ - LOAD_S, SQ - LOAD_S, C_LOAD_FILL, C_LOAD_EDGE, LOAD_S, cls="load"),
                    self.circle(x, y, RDOT, cls="bus")]
        if kind == "bus":
            return [self.circle(x, y, RDOT, cls="bus")]
        if kind in ("LBSM", "LBS"):
            w, h = (BAR_T, BAR_L) if horizontal_line else (BAR_L, BAR_T)
            if kind == "LBSM":
                return [self.rect(x, y, w, h, INK, cls="LBSM")]
            return [self.rect(x, y, w - BAR_S, h - BAR_S, "#FFFFFF", INK, BAR_S, cls="LBS")]
        if kind == "REC":
            return [self.poly([(x + RHEX * math.sin(math.radians(a)), y - RHEX * math.cos(math.radians(a)))
                               for a in range(0, 360, 60)], cls="REC")]
        if kind == "CB":
            return [self.poly([(x - TRI_W / 2, y - TRI_H / 2), (x + TRI_W / 2, y - TRI_H / 2), (x, y + TRI_H / 2)], cls="CB")]
        if kind == "tie":
            return [self.ellipse(x, y, *ELL, cls="tie")]
        raise ValueError(kind)


def text_width(s):
    from PIL import ImageFont
    for p in ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "C:/Windows/Fonts/arial.ttf",
              "/Library/Fonts/Arial.ttf", "/System/Library/Fonts/Supplemental/Arial.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, 1000).getlength(s) / 1000 * FS
    return 0.5 * FS * len(s)


# ----------------------------------------------------------------------------- build
def build(evaluator_dir):
    feeders, tie = model_tree(evaluator_dir)
    cells, segs, problems = {}, [], []
    def occupy(cell, what):
        if cell in cells:
            problems.append(f"sel {cell} dipakai dua kali: {cells[cell]} dan {what}")
        cells[cell] = what
    def straight(a, b):
        return a[0] == b[0] or a[1] == b[1]
    def cells_between(a, b):
        n = max(abs(a[0] - b[0]), abs(a[1] - b[1]))
        sx, sy = (b[0] > a[0]) - (b[0] < a[0]), (b[1] > a[1]) - (b[1] < a[1])
        return [(a[0] + sx * i, a[1] + sy * i) for i in range(1, n)]
    sc = Scene(W, H)
    used = set()
    plan = []                                                    # (feeder, lines, symbols)
    for cb, f in feeders.items():
        hc = HEAD_COL[cb]
        pos = dict(POS); pos[f["head"]] = (hc, -1)
        lines, symbols = [], []
        occupy((hc, 0), cb); symbols.append(("CB", cb, (hc, 0), False))
        for name, n in f["nodes"].items():
            if name not in POS:
                problems.append(f"bus {name} ({cb}) belum punya posisi"); continue
            used.add(name)
            occupy(POS[name], name)
            symbols.append(("load" if n["load"] else "bus", name, POS[name], False))
        for name, n in f["nodes"].items():
            if name not in POS or n["parent"] not in pos:
                continue
            path = [pos[n["parent"]]] + VIA.get(name, []) + [POS[name]]
            inner = []
            for a, b in zip(path[:-1], path[1:]):
                if not straight(a, b):
                    problems.append(f"garis {n['parent']} -> {name} tidak lurus: {a} -> {b} (tambahkan sudut di VIA)")
                inner += cells_between(a, b) + ([b] if b != path[-1] else [])
            corners = set(VIA.get(name, []))
            free = [c for c in inner if c not in corners and not (n["parent"] == f["head"] and c == (hc, 0))]
            if len(free) != len(n["devs"]):
                problems.append(f"garis {n['parent']} -> {name}: {len(free)} sel kosong untuk {len(n['devs'])} perangkat {n['devs']}")
            for c in corners:
                occupy(c, f"sudut garis ke {name}")
            for c, d in zip(free, n["devs"]):
                occupy(c, d)
                # a switch bar is drawn across its line: is the line horizontal at this cell?
                seg_h = any(a[1] == b[1] == c[1] and min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) for a, b in zip(path[:-1], path[1:]))
                symbols.append((dev_kind(d), d, c, seg_h))
            lines.append((name, path))
            segs += [(a, b, f"{n['parent']} -> {name}") for a, b in zip(path[:-1], path[1:])]
        plan.append((cb, f, lines, symbols))
    for name in POS:
        if name not in used:
            problems.append(f"posisi untuk {name} tidak dipakai (bukan bus beban/percabangan pada model)")
    # tie route
    tpath = [POS[tie[0]]] + TIE_VIA + [POS[tie[1]]]
    for a, b in zip(tpath[:-1], tpath[1:]):
        if not straight(a, b):
            problems.append(f"rute tie tidak lurus: {a} -> {b}")
        for c in cells_between(a, b) + ([b] if b != tpath[-1] else []):
            occupy(c, "tie")
        segs.append((a, b, "tie"))
    # no two lines may overlap or cross away from their common end
    def pts(a, b):
        return set([a] + cells_between(a, b) + [b])
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            common = pts(*segs[i][:2]) & pts(*segs[j][:2])
            ends = {segs[i][0], segs[i][1]} & {segs[j][0], segs[j][1]}
            if common - ends:
                problems.append(f"garis bertumpuk/bersilang: {segs[i][2]} dan {segs[j][2]} di {sorted(common - ends)}")
    if problems:
        raise SystemExit("TATA LETAK BERMASALAH:\n  " + "\n  ".join(problems))

    P = lambda c: (X(c[0]), Y(c[1]))
    for cb, f, lines, symbols in plan:
        body = [sc.line([P(c) for c in path], COL[cb], LW) for _, path in lines]
        for kind, name, c, seg_h in symbols:
            tag = name.replace("(", "_").replace(")", "").replace(" ", "_")
            body.append(f'<g id="{tag}">' + "".join(sc.sym(kind, *P(c), horizontal_line=seg_h)) + "</g>")
        body.append(sc.text(X(HEAD_COL[cb]) + TRI_W / 2 + 1.0, Y(0), FEEDER_LABEL[cb]))
        sc.group(FEEDER_LABEL[cb].replace(" ", "_").replace("(", "").replace(")", ""), body)
    sc.group("tie_" + tie[2].replace(" ", "_"), [sc.line([P(c) for c in tpath], INK, TIE_LW, TIE_DASH)] + sc.sym("tie", *P(TIE_SWITCH)))
    # source bar (drawn after the feeders, whose lines end under it) and its label
    sc.group("source", [sc.rect((X(0) + X(NCOL)) / 2, YBAR, X(NCOL) - X(0) + SQ, SRC_T, INK, cls="source"),
                        sc.text((X(0) + X(NCOL)) / 2, YBAR - SRC_T / 2 - 0.9 - 0.358 * FS, SRC_LABEL, "middle")])
    # legend
    wmax = max(text_width(t) for _, t in LEGEND)
    width = 1.4 + BAR_L + 1.6 + wmax + 1.6
    x0 = (X(LEGEND_BETWEEN[0]) + X(LEGEND_BETWEEN[1]) - width) / 2
    x1, lx = x0 + width, x0 + 1.4 + BAR_L / 2
    tx = lx + BAR_L / 2 + 1.6
    stub = ELL[1] + TIE_DASH[0] + 0.6                               # half length of the dashed stub of the tie entry
    ys = [Y(LEGEND_ROW0 + i) + (stub - 2.15 if kind == "tie" else 0.0) for i, (kind, _) in enumerate(LEGEND)]
    y0, y1 = ys[0] - 3.0, ys[-1] + (stub + 0.9 if LEGEND[-1][0] == "tie" else 3.0)
    body = [sc.rect((x0 + x1) / 2, (y0 + y1) / 2, x1 - x0, y1 - y0, "#FFFFFF", INK, 0.2, cls="frame")]
    for yy, (kind, label) in zip(ys, LEGEND):
        if kind == "tie":                                           # one dash on each side, drawn towards the switch
            body += [sc.line([(lx, yy - stub), (lx, yy)], INK, TIE_LW, TIE_DASH),
                     sc.line([(lx, yy + stub), (lx, yy)], INK, TIE_LW, TIE_DASH)]
        body += sc.sym(kind, lx, yy)
        body.append(sc.text(tx, yy, label))
    sc.group("legend", body)
    stats = {cb: dict(load_buses=sum(1 for n in f["nodes"].values() if n["load"]), loads=f["nload"],
                      junctions=sum(1 for n in f["nodes"].values() if not n["load"]),
                      devices=[d for n in f["nodes"].values() for d in n["devs"]]) for cb, f, _, _ in plan}
    return sc, stats, tie


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--evaluator", default=os.path.join(here, "..", "evaluator"))
    ap.add_argument("--out", default=os.path.join(here, "figures"))
    ap.add_argument("--name", default="fig_bawean")
    a = ap.parse_args()
    sc, stats, tie = build(os.path.abspath(a.evaluator))
    svg = sc.svg()
    os.makedirs(a.out, exist_ok=True)
    base = os.path.join(a.out, a.name)
    open(base + ".svg", "w", encoding="utf-8").write(svg)
    import cairosvg
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=base + ".pdf")
    cairosvg.svg2eps(bytestring=svg.encode(), write_to=base + ".eps")
    cairosvg.svg2png(bytestring=svg.encode(), write_to=base + ".png", output_width=2055)       # 300 dpi at 174 mm
    print(f"ukuran {W:g} x {H:.1f} mm | jarak kolom {DX:.2f} mm, jarak baris {DY:.2f} mm | huruf {FS / PT:g} pt")
    for cb, s in stats.items():
        print(f"{FEEDER_LABEL[cb]:24s} bus beban {s['load_buses']:3d} (beban PowerFactory {s['loads']:3d}) | bus percabangan tanpa beban {s['junctions']:2d} | perangkat {len(s['devices']) + 1}")
    print("tie:", tie)


if __name__ == "__main__":
    main()
