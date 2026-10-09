# Switch, recloser and tie placement with AC-feasible two-stage restoration

Code and data for the article

> H. A. Yudiansyah, R. S. Wibowo, H. A. Wardana, "Value-Based and Regulation-Aware Placement of
> Sectionalizing Switches, Reclosers, and Tie Switches in Distribution Networks with AC-Feasible
> Two-Stage Restoration".

Department of Electrical Engineering, Institut Teknologi Sepuluh Nopember, Surabaya, Indonesia.

## Contents

* **Evaluator** of a switching-device plan: failure mode and effect analysis of every permanent section
  fault, two-stage restoration (remote devices after 5 min, all devices after the manual switching time
  of 60 min), an exact search over tie switches and restorable blocks, and AC feasibility of every
  restoration by backward-forward sweep load flow at four load levels. A plan is scored with three
  objectives, all in USD per year:
  * f1: annualized device cost,
  * f2: expected customer interruption cost (ECOST),
  * f3: expected guaranteed-service compensation under Indonesian regulation (Regulation of the Minister
    of Energy and Mineral Resources No. 2 of 2025).
* **Test systems**: IEEE 33-bus, IEEE 69-bus, RBTS Bus 2 and the 20 kV island network of Bawean, Indonesia.
* **Experiments**: NSGA-II (pymoo) for the three-objective fronts, and a single-objective genetic
  algorithm with local search for the budget, load-level and sensitivity studies.
* **Results** used in the article (`paper/paper_data/`) and the scripts that reproduce its figures and tables.
* **Verification scripts** for Table 1 of the article.

Comments and console messages in the code are partly in Indonesian; the expected outputs are given below.

## Repository structure

| Path | Content |
| --- | --- |
| `evaluator/network.py` | radial topology and backward-forward sweep load flow |
| `evaluator/reliability.py` | fault response, two-stage restoration, reliability indices, f1 and f2 |
| `evaluator/compensation.py` | f3: monthly interruption duration (Panjer recursion), tariffs, compensation |
| `evaluator/ecost.py` | sector customer damage functions and their conversion to 2025 USD |
| `evaluator/systems.py` | test systems, parameters and encoding of a plan |
| `evaluator/ieee33_data.py`, `ieee69_data.py`, `ieee69_raw.py`, `rbts2_data.py` | benchmark system data |
| `evaluator/bawean_data.py`, `bawean_export.xlsx` | Bawean network, built from a PowerFactory export |
| `evaluator/go_nogo_rbts.py` | RBTS load-duration data (four load levels) used by `systems.py` |
| `evaluator/run_nsga2.py`, `analyze_front.py` | Pareto fronts (E2) |
| `evaluator/e3_budget.py`, `e3_analyze.py` | objective comparison under a budget (E3, exact enumeration E6) |
| `evaluator/scalar_opt.py`, `scalar_analyze.py` | load-level study (E4) and sensitivity study (E5) |
| `evaluator/make_warm_start.py`, `warm_start.json` | best-known plans added to the initial populations of E5 |
| `evaluator/bench_e7.py` | time per evaluation (E7) |
| `evaluator/regression.py`, `test_evaluator.py`, `verify_oracle.py`, `verify_two_stage.py` | validation |
| `paper/paper_data/` | results used in the article |
| `paper/make_*.py`, `check_bawean_sld.py`, `factor_guard.py` | figures and tables of the article |
| `paper/checks/` | additional checks of numbers reported in the article |

## Installation

Python 3.10 or newer.

    pip install -r requirements.txt

`cairosvg` is needed only to convert Fig. 2 to PDF, EPS and PNG, and it requires the Cairo library.
`pandapower` is needed only for `paper/checks/check_loadflow_nr.py`.

Quick check (about 15 s), inside `evaluator/`:

    python regression.py         # last line: SEMUA REGRESI LOLOS (all regression checks passed)
    python test_evaluator.py     # every test prints PASS

## Reproducing the article

The experiments run inside `evaluator/`. `--workers` sets the number of parallel processes (default:
all logical cores). The budget, load-level and sensitivity scripts can be stopped and restarted with the
same command; tasks already stored in the output file are skipped.

### Validation (Section 3.4, Table 1)

| Check | Command | Expected result |
| --- | --- | --- |
| IEEE 33-bus reliability (9 configurations) and base-case load flow vs PowerFactory; RBTS Bus 2 cases B, D, F vs published indices; IEEE 69-bus and Bawean base-case load flow | `python regression.py` | `SEMUA REGRESI LOLOS` |
| Hand-calculated cases and monotonicity tests | `python test_evaluator.py` | `PASS` for every test |
| Restoration search vs brute-force enumeration, single stage | `python verify_oracle.py 40` | 444 fault cases, `tidak cocok: 0` (no mismatch) |
| Same, two stages | `python verify_two_stage.py 30` | 307 fault cases, `tidak cocok tahap 1: 0, tahap 2: 0` |
| Load flow vs Newton-Raphson (pandapower) | `python ../paper/checks/check_loadflow_nr.py` | voltage differences below 1e-10 pu |
| f3 vs Monte Carlo | `python ../paper/checks/check_mc.py` | f3 of every tested plan within 0.1% of the Monte Carlo estimate |

The PowerFactory reference values are stored in `regression.py` and `test_evaluator.py`.

### Pareto fronts (Section 5.2: Figs. 3 and 4, Table 3; Fig. S1)

    python run_nsga2.py --system ieee33 --seeds 0-9 --pop 200 --gen 400 --out results_e2
    python run_nsga2.py --system ieee69 --seeds 0-9 --pop 200 --gen 400 --out results_e2
    python run_nsga2.py --system rbts2  --seeds 0-9 --pop 200 --gen 250 --out results_e2
    python run_nsga2.py --system bawean --seeds 0-9 --pop 200 --gen 800 --out results_e2
    python analyze_front.py results_e2/ieee33_t60_h1      # likewise for the other systems

`analyze_front.py` writes `front_global.csv` (union front of the ten runs) and `convergence_hv.csv`,
which correspond to `paper/paper_data/front_<system>.csv` and `conv_<system>.csv`. One run took
26 min (IEEE 33-bus), 1.9 h (IEEE 69-bus), 8 min (RBTS Bus 2) and 1.5 h (Bawean) on the 16 threads
of an Intel Core i7-13620H.

### Objective comparison under a budget (Section 5.3: Fig. 5, Table 4, Table S5)

    python e3_budget.py --system ieee33 --budgets 3,5,8,10,12 --seeds 5
    python e3_budget.py --system ieee33 --budgets 8,10,12 --seeds 10 --seed-offset 5 --no-exact
    # the same two commands for ieee69 and rbts2
    python e3_budget.py --system bawean --budgets 3,5,8,10,12 --seeds 5 --no-exact
    python e3_analyze.py e3_results/ieee33_t60_h1.json   # likewise for the other systems

Budgets are expressed as the equivalent number of manual switches (MS). The first command also
enumerates all plans at the smallest budget. The output `e3_results/<system>_t60_h1.json` corresponds to
`paper/paper_data/<system>_t60_h1.json`.

### Load level used for feasibility checks (Section 5.4: Fig. 6, Table S6)

    python scalar_opt.py --exp e4 --out e45_results
    python scalar_analyze.py e45_results/e4.json         # writes e45_results/e4_analysis.json

### Sensitivity analysis (Section 5.5: Tables S7 and S8)

    python scalar_opt.py --exp e5   --out e45_results    # IEEE 33-bus, IEEE 69-bus, RBTS Bus 2
    python scalar_opt.py --exp e5bw --out e45_results    # Bawean

The initial populations contain the best-known plans of each case from `warm_start.json`, which
`make_warm_start.py` builds from the union fronts and earlier single-objective runs.

### Computational performance (Section 5.6)

    python bench_e7.py     # time per evaluation on one thread

### Figures and tables

Run inside `paper/`. The scripts read `paper_data/` and the evaluator in `../evaluator`, and write to
`figures/` and `tables/`.

| Command | Output |
| --- | --- |
| `python make_framework.py` | Fig. 1 |
| `python make_bawean_sld.py`, then `python check_bawean_sld.py` | Fig. 2 and an automatic check of the drawing against the network model |
| `python make_figures.py` | Figs. 3, 5, 6 and S1 |
| `python make_placement.py` | Fig. 4 |
| `python make_tables.py` | Tables 3 and 4 |
| `python make_esm_tables.py` | Tables S3 and S5 to S8 of Online Resource 1 |

With the data in `paper_data/`, these scripts reproduce the figures and tables of the article.

### Additional checks (`paper/checks/`)

| Script | Checks |
| --- | --- |
| `check_52.py` | system data, plans without new devices, representative plans and ratios in Section 5.2 |
| `check_bawean.py`, `check_bawean_single.py` | Bawean results in Section 5.2: tie upgrade and every single-device addition |
| `check_mc.py` | f3 from the Panjer recursion vs Monte Carlo, 3 million months per load point (Table 1) |
| `check_loadflow_nr.py` | backward-forward sweep vs Newton-Raphson of pandapower (Table 1) |
| `check_ieee69_ties.py` | sensitivity of the IEEE 69-bus results to the reactance of the tie lines |

## Results (`paper/paper_data/`)

| File | Content |
| --- | --- |
| `front_<system>.csv` | union Pareto front of the ten NSGA-II runs: `f1`, `f2`, `f3` and the plan `g0, g1, ...` |
| `conv_<system>.csv` | hypervolume of the current front, one row per run and one column per generation |
| `<system>_t60_h1.json` | budget study: one record per budget, objective (`ecost` or `comp`) and run (`seed`, or `exact` for complete enumeration), with the plan `g`, the objectives, SAIFI, SAIDI, ENS and the share of customers above the monthly threshold |
| `e4.json`, `e4_analysis.json` | load-level study (RBTS Bus 2, transformer repair times of 10 to 200 h): plans optimized with peak, average and four-level load, evaluated with their own load model (`own`) and with the four-level reference (`atref`) |
| `e5.json`, `e5bw.json` | sensitivity study (benchmark systems; Bawean) |
| `FAKTOR.txt` | note on the ECOST conversion factor (see below) |

`t60_h1` denotes a manual switching time of 60 min and a monthly threshold of 1 h. A plan has one gene
per candidate section (0 none, 1 MS, 2 RCS, 3 REC) followed by one gene per candidate tie (0 none,
1 manual, 2 remote) or, for Bawean, one gene for the existing tie (0 manual, 1 upgraded to remote).
`System(name).describe(g)` in `evaluator/systems.py` lists the devices of a plan. The per-run NSGA-II
files (`seed<k>.npz`) are not included because of their size.

## ECOST conversion factor

The sector customer damage functions of Billinton et al. (1989), in 1987 Canadian dollars, are converted
to 2025 US dollars with the factor 164.2 / 68.5 / 1.3978 = 1.7149 (`evaluator/ecost.py`, Table S3).
The NSGA-II fronts and the budget study were computed with an earlier factor of 1.6543, which used a
2025 consumer price index of 158.4 instead of 164.2. Non-dominated sorting, crowding distance, the
normalized hypervolume and the budget-constrained problems do not change when f2 is multiplied by a
constant, so the stored f2 values of these results were multiplied by 164.2 / 158.4 = 1.036616. The
load-level and sensitivity studies, whose societal objective f1 + f2 depends on the factor, were rerun
with the factor 1.7149.

## Data sources

* IEEE 33-bus: Baran and Wu (1989); section lengths, customers and sectors as described in the article.
* IEEE 69-bus: MATPOWER `case69` (Zimmerman et al. 2011), from Baran and Wu (1989).
* RBTS Bus 2: Allan et al. (1991) and Billinton et al. (1989).
* Bawean: the 20 kV network data were provided by the distribution network operator and are shared with
  its permission.

All parameters and their sources are listed in Table 2 of the article and in Online Resource 1.

## License

The code is released under the MIT License (`LICENSE`). The Bawean network data and the results in
`paper/paper_data/` are released under the Creative Commons Attribution 4.0 International license
(CC BY 4.0). The benchmark system data are taken from the publications cited above.

## Citation

If you use this code or data, please cite the article (reference to be added on publication) and this
repository (see `CITATION.cff`).
