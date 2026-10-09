"""Satu pintu untuk semua sistem uji: jaringan, parameter, gen keputusan, basis kompensasi.
Formulasi v4: mode load flow, 4 kelompok beban (LDC RBTS), restorasi dua tahap eksak."""
import numpy as np
from network import Network, NONE, REC
from reliability import Evaluator, default_params, crf
from compensation import monthly_base_usd, f3 as comp_f3
from go_nogo_rbts import ldc_levels
import ieee33_data, ieee69_data, rbts2_data as R
import bawean_data as BW

LEVELS = ldc_levels(4)


class System:
    def __init__(self, name, tau_ms_min=60.0, thr_h=1.0, tariff_scale=None, cost_ratio=None,
                 mode="loadflow", tr_repair=None, load_model="ldc4", clpu_a=None):
        self.name, self.thr_h = name, thr_h
        p = default_params()
        if name in ("ieee33", "ieee69"):
            mod = ieee33_data if name == "ieee33" else ieee69_data
            self.net = Network(**mod.build())
            self.cand = list(range(len(self.net.sec)))
            self.labels = [f"{s['frm']+1}-{s['to']+1}" for s in self.net.sec]
        elif name in ("rbts2", "rbts2_tr10"):
            kw, info = R.build("peak")
            tr = 10.0 if name == "rbts2_tr10" else tr_repair
            if tr is not None:
                for s in info["transformer"]:
                    kw["sections"][s]["r_rep"] = float(tr)
            self.net = Network(**kw, fixed={s: REC for s in info["feeder_head"]},
                               s_max_kva=R.SP_CAPACITY_KVA)
            p.update(r_rep=R.LINE_REPAIR)
            self.cand = [s for s in info["main"] if s not in info["feeder_head"]] + list(info["lateral"])
            self.labels = [f"s{info['sec_no'][s]}" for s in range(len(self.net.sec)) if s in info["sec_no"]]
            self.labels = {s: f"s{info['sec_no'][s]}" for s in info["sec_no"]}
        elif name == "bawean":
            kw, info = BW.build()
            fixed = {k: {"MS": 1, "RCS": 2, "REC": 3}[d] for k, d in info["fixed"].items()}
            self.net = Network(**kw, fixed=fixed)
            self.cand = [k for k in range(len(self.net.sec)) if k not in fixed]   # line sections
            self.labels = {k: info["labels"][k] for k in range(len(self.net.sec))}
            self.bw_heads = set(info["heads"]); self.existing_ties = list(range(len(self.net.ties)))
            self.contract_kva = info["contract_kva"]
        else:
            raise ValueError(f"sistem tidak dikenal: {name}")
        if load_model == "ldc4":
            levels = LEVELS
        elif load_model == "peak":
            levels = [(1.0, 1.0)]
        elif load_model == "avg":                    # single snapshot at the LDC mean
            levels = [(sum(m * q for m, q in LEVELS), 1.0)]
        else:
            raise ValueError(load_model)
        p.update(tau_ms=tau_ms_min / 60.0, load_levels=levels)
        if clpu_a is not None:
            a = float(clpu_a)
            p["clpu"] = dict(kmax={"residential": 1 + a, "commercial": 1 + a, "government": 1 + a,
                                   "industrial": 1 + 0.2 * a}, tc_h=1.0, tdec_h=0.5, step_h=10 / 60)
        if cost_ratio is not None:
            p["cost"] = dict(p["cost"]); p["cost"]["RCS"] = cost_ratio * p["cost"]["MS"]
        self.ev = Evaluator(self.net, p, mode=mode)
        self.p = p
        self.nsec, self.ntie = len(self.cand), len(self.net.ties)
        self.nvar = self.nsec + self.ntie
        self.loads = [j for j in self.net.load_buses if self.net.N[j] > 0]
        self.base, self.classes = monthly_base_usd(self.net, self.loads, tariff_scale=tariff_scale,
                                                   contract_kva=getattr(self, "contract_kva", None))
        self.ms_ann = (crf(p["delta"], p["h_sw"]) + p["rho_om"]) * p["cost"]["MS"]
        # existing ties (Bawean): gene 0 = keep manual, 1 = upgrade to remote
        tie_ub = [1 if c in getattr(self, "existing_ties", []) else 2 for c in range(self.ntie)]
        self.xu = np.array([3] * self.nsec + tie_ub)

    def decode(self, g):
        x = [NONE] * len(self.net.sec)
        for i, s in enumerate(self.cand):
            x[s] = int(g[i])
        w = [int(v) for v in g[self.nsec:]]
        for c in getattr(self, "existing_ties", []):
            w[c] = 1 + w[c]                       # existing tie is always there (manual or remote)
        return x, w

    def evaluate(self, g):
        x, w = self.decode(g)
        o = self.ev.evaluate(x, w, severity=True)
        if getattr(self, "existing_ties", None):
            # sunk cost of existing manual ties; an upgrade costs (RCS - MS)
            o["f1"] -= len(self.existing_ties) * self.ms_ann
            # recloser-series limit counts existing reclosers too (feeder breakers excluded)
            xf = self.ev._full(x)
            for h in self.bw_heads:
                xf[h] = NONE
            o["rec_ok"] = self.ev.rec_series_ok(xf)
        sev = [o["sev"][j] for j in self.loads]
        comp, share = comp_f3(self.net, self.loads, sev, self.thr_h, base=self.base)
        return dict(f1=o["f1"], f2=o["ECOST"], f3=comp, SAIFI=o["SAIFI"], SAIDI=o["SAIDI"],
                    ENS=o["ENS"], share_above=share, rec_ok=o["rec_ok"])

    def describe(self, g):
        out = []
        names = {0: None, 1: "MS", 2: "RCS", 3: "REC"}
        for i, s in enumerate(self.cand):
            if g[i]:
                lab = self.labels[s] if isinstance(self.labels, dict) else self.labels[s]
                out.append(f"{names[int(g[i])]}@{lab}")
        for c in range(self.ntie):
            v = int(g[self.nsec + c])
            if v:
                if c in getattr(self, "existing_ties", []):
                    out.append("tie KOPEL TAMBAK: upgrade ke remote")   # gene 1 = upgrade
                    continue
                if self.name.startswith("rbts2"):
                    lab = ["F1-F2", "F3-F4"][c]            # normally open points (Allan 1991, Fig. 2)
                else:
                    t = self.net.ties[c]; lab = f"{t['a']+1}-{t['b']+1}"
                out.append(f"tie{lab}:{'M' if v == 1 else 'R'}")
        return out
