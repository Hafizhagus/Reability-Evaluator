"""RBTS Bus 2 distribution system (Allan et al., IEEE Trans. Power Syst. 6(2), 1991).

Source data: Fig. 2 (topology), Table 2 (feeder section lengths by type), Table 3
(customers and loads per load point), Table 5 (reliability data).
Topology reconstructed from Fig. 2: every feeder has main sections with laterals
tapped at the main-feeder nodes; each lateral ends at one load point (LP).

ASSUMPTIONS (to be stated in the paper):
  * 11 kV line impedance per km is not given in the source. Assumed overhead
    conductor AAAC 150 mm2: z' = 0.22 + j0.35 ohm/km, ampacity 412 A (ESB Networks MV
    line data) -- same conductor assumption as the IEEE 33 impedance-based lengths.
  * Power factor 0.98 lagging (Q = 20 % of P), the value suggested by the RBTS basic
    data paper (Billinton et al. 1989) when reactive power matters.
  * LV transformers (11/0.415 kV) are modelled as zero-impedance branches between the
    lateral end and the load point; small users (LP8, LP9) have no utility
    transformer (Allan et al. 1991, comment (a)).
  * Normally open points: F1-F2 (between the ends of laterals to LP7 and LP9) and
    F3-F4 (between the ends of the main feeders at LP15 and LP21), read from Fig. 2.
"""
V_KV = 11.0
Z_PER_KM = complex(0.22, 0.35)
AMPACITY_A = 412.0
LINE_LAMBDA = 0.065          # f/km/yr, 11 kV overhead lines (Table 5)
LINE_REPAIR = 5.0            # h (Table 5)
TR_LAMBDA = 0.015            # f/yr, 11/0.415 kV transformer (Table 5)
TR_REPAIR = 200.0            # h (Table 5)
SWITCHING = 1.0              # h, "line" system switching time (Table 5)
SP_CAPACITY_KVA = 2 * 16000  # two 16 MVA 33/11 kV transformers (Table 5)

# Table 2: length (km) by feeder section number
TYPE_LEN = {1: 0.60, 2: 0.75, 3: 0.80}
SEC_TYPE = {}
for s in (2, 6, 10, 14, 17, 21, 25, 28, 30, 34):
    SEC_TYPE[s] = 1
for s in (1, 4, 7, 9, 12, 16, 19, 22, 24, 27, 29, 32, 35):
    SEC_TYPE[s] = 2
for s in (3, 5, 8, 11, 13, 15, 18, 20, 23, 26, 31, 33, 36):
    SEC_TYPE[s] = 3

# Feeders: ordered main sections; laterals tapped at the node AFTER each main section.
# (feeder, [(main_section, [(lateral_section, LP), ...]), ...])
FEEDERS = [
    ("F1", [(1, [(2, 1), (3, 2)]), (4, [(5, 3), (6, 4)]), (7, [(8, 5), (9, 6)]),
            (10, [(11, 7)])]),
    ("F2", [(12, [(13, 8)]), (14, [(15, 9)])]),
    ("F3", [(16, [(17, 10)]), (18, [(19, 11), (20, 12)]), (21, [(22, 13), (23, 14)]),
            (24, [(25, 15)])]),
    ("F4", [(26, [(27, 16), (28, 17)]), (29, [(30, 18), (31, 19)]), (32, [(33, 20)]),
            (34, [(35, 21), (36, 22)])]),
]

# Table 3: (customer type, average MW, peak MW, customers)
LP_DATA = {}
for lp in (1, 2, 3, 10, 11):
    LP_DATA[lp] = ("residential", 0.535, 0.8668, 210)
for lp in (12, 17, 18, 19):
    LP_DATA[lp] = ("residential", 0.450, 0.7291, 200)
LP_DATA[8] = ("small user", 1.00, 1.6279, 1)
LP_DATA[9] = ("small user", 1.15, 1.8721, 1)
for lp in (4, 5, 13, 14, 20, 21):
    LP_DATA[lp] = ("govt/inst", 0.566, 0.9167, 1)
for lp in (6, 7, 15, 16, 22):
    LP_DATA[lp] = ("commercial", 0.454, 0.7500, 10)

# customer type -> damage-function sector (Billinton et al. 1989, Table XV)
SECTOR = {"residential": "residential", "commercial": "commercial",
          "govt/inst": "government", "small user": "industrial"}

TIE_LPS = [(7, 9), (15, 21)]          # normally open points between feeder ends


def build(load="peak"):
    """Return (network kwargs, info). Buses: 0 = 11 kV SP busbar."""
    buses = {"SP": 0}
    sections, P, Q, N, sector = [], [0.0], [0.0], [0], [None]
    info = dict(main=[], lateral=[], transformer=[], feeder_head=[], lp_bus={},
                lp_node={}, sec_no={}, feeder_of={})

    def new_bus(name, p=0.0, n=0, sec=None):
        buses[name] = len(P)
        P.append(p); Q.append(0.2 * p); N.append(n); sector.append(sec)
        return buses[name]

    def add_sec(frm, to, length, kind, number=None, lam=None, rrep=None, feeder=None):
        z = Z_PER_KM * length if length > 0 else complex(1e-5, 1e-5)
        d = dict(frm=frm, to=to, R=z.real, X=z.imag, length_km=length,
                 ampacity_a=AMPACITY_A)
        if lam is not None:
            d["lam"] = lam
        if rrep is not None:
            d["r_rep"] = rrep
        sections.append(d)
        idx = len(sections) - 1
        info[kind].append(idx)
        info["feeder_of"][idx] = feeder
        if number is not None:
            info["sec_no"][idx] = number
        return idx

    for fname, mains in FEEDERS:
        prev = 0
        for mi, (ms, laterals) in enumerate(mains):
            node = new_bus(f"{fname}_n{mi}")
            idx = add_sec(prev, node, TYPE_LEN[SEC_TYPE[ms]], "main", ms,
                          lam=LINE_LAMBDA * TYPE_LEN[SEC_TYPE[ms]], feeder=fname)
            if mi == 0:
                info["feeder_head"].append(idx)
            for ls, lp in laterals:
                ctype, avg, pk, ncust = LP_DATA[lp]
                mw = pk if load == "peak" else avg
                lp_node = new_bus(f"LP{lp}_hv")
                add_sec(node, lp_node, TYPE_LEN[SEC_TYPE[ls]], "lateral", ls,
                        lam=LINE_LAMBDA * TYPE_LEN[SEC_TYPE[ls]], feeder=fname)
                info["lp_node"][lp] = lp_node
                if ctype == "small user":
                    # no utility transformer: load sits at the lateral end
                    P[lp_node] = mw * 1000; Q[lp_node] = 0.2 * mw * 1000
                    N[lp_node] = ncust; sector[lp_node] = SECTOR[ctype]
                    info["lp_bus"][lp] = lp_node
                else:
                    b = new_bus(f"LP{lp}", mw * 1000, ncust, SECTOR[ctype])
                    add_sec(lp_node, b, 0.0, "transformer", lam=TR_LAMBDA,
                            rrep=TR_REPAIR, feeder=fname)
                    info["lp_bus"][lp] = b
            prev = node
    ties = []
    for a, b in TIE_LPS:
        z = Z_PER_KM * 0.75
        ties.append(dict(a=info["lp_node"][a], b=info["lp_node"][b], R=z.real, X=z.imag,
                         length_km=0.75, ampacity_a=AMPACITY_A))
    net = dict(n_bus=len(P), sections=sections, ties=ties, P_kw=P, Q_kvar=Q,
               customers=N, root=0, V_kv=V_KV, sectors=sector)
    return net, info
