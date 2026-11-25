"""Python CSI-GEP vs R CSI-GEP on the same NMF outputs (tutorial demo, 2026.11 run).

Reference: tests/reference/demo/ (export_demo_reference.R). Python output: tests/output/demo/ (run_demo_python.bsub).
Checks, all must hold:
  1. Jaccard step: same significant pairs (same order) and community number for every (K, JL); p-values at K = 22
     equal to 1e-9 relative.
  2. Selection: per JL the same elbow rank and X_selected; slopes equal to 1e-9 relative; same selected (K, JL).
  3. Final GEPs (rescue 0 and 1): same parameters, same GEP membership; GEP scores equal to 1e-12 relative or 1e-15
     absolute (R's file reader is not exact in the last bit, ~1e-16), matched by membership (GEP numbering is reported separately: same / permuted).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

T = os.path.dirname(__file__)
REF = os.environ.get("REF", os.path.join(T, "reference", "demo"))
OUT = os.environ.get("OUT", os.path.join(T, "output", "demo"))
NAME = os.environ.get("NAME", "demo")
CHECK_JACCARD = os.environ.get("CHECK_JACCARD", "1") == "1"   # 0 when Python used R's converted Jaccard results
PV_K = os.environ.get("PV_K", "22")
log = []


def say(m):
    print(m, flush=True); log.append(m)


def check_jaccard():
    sp = pd.read_csv(os.path.join(REF, "sig_pairs.tsv"), sep="\t")
    cm = pd.read_csv(os.path.join(REF, "communities.tsv"), sep="\t")
    bad = []
    for r in cm.itertuples():
        z = np.load(os.path.join(OUT, "results", f"jaccardtest_results_K{r.K}_JL{r.JL}.npz"))
        ref_sig = sp[(sp.K == r.K) & (sp.JL == r.JL)].sort_values("ord")[["row", "col"]].values
        if int(z["community_number"]) != r.community_number or not np.array_equal(z["sig"].reshape(-1, 2), ref_sig.reshape(-1, 2)):
            bad.append((r.K, r.JL, r.community_number, int(z["community_number"]), len(ref_sig), len(z["sig"])))
    say(f"1. Jaccard: {len(cm)} (K, JL) compared; mismatches {len(bad)} {bad[:5]}")
    pv = pd.read_csv(os.path.join(REF, f"pvalues_K{PV_K}.tsv"), sep="\t")
    worst = 0.0
    for jl, sub in pv.groupby("JL"):
        z = np.load(os.path.join(OUT, "results", f"jaccardtest_results_K{PV_K}_JL{jl}.npz"))["pvalues"]
        py = z[sub.i.values - 1, sub.j.values - 1]
        worst = max(worst, float(np.max(np.abs(py - sub.p.values) / np.maximum(np.abs(sub.p.values), 1e-300))))
    say(f"   p-values K{PV_K}: max relative difference {worst:.2e}")
    return not bad and worst < 1e-9


def check_selection():
    ref = pd.read_csv(os.path.join(REF, "selection.tsv"), sep="\t")
    py = pd.read_csv(os.path.join(OUT, "GEP_results_rescue0", f"{NAME}_selection.tsv"), sep="\t")
    m = ref.merge(py, on="JL", suffixes=("_r", "_py"))
    ok_rank = (m.Rank_r == m.Rank_py).all()
    ok_x = np.allclose(m.X_selected_r, m.X_selected_py, rtol=1e-12)
    rel = np.abs(m.Slope_r - m.Slope_py) / np.abs(m.Slope_r)
    say(f"2. Selection: ranks equal {ok_rank}, X_selected equal {ok_x}, slope max relative difference {rel.max():.2e}")
    if not (ok_rank and ok_x):
        say(m.to_string())
    return ok_rank and ok_x and rel.max() < 1e-9


def check_final(rs):
    d_ref = os.path.join(REF, f"rescue{rs}"); d_py = os.path.join(OUT, f"GEP_results_rescue{rs}")
    pr = pd.read_csv(os.path.join(d_ref, "params.tsv"), sep="\t").iloc[0]
    pp = pd.read_csv(os.path.join(d_py, f"{NAME}_optimized_parameters.tsv"), sep="\t").iloc[0]
    same_params = all(int(pr[c]) == int(pp[c]) for c in ["Community_Number", "Rank", "Jaccard_Length"])
    gi = pd.read_csv(os.path.join(d_ref, "gepindex.tsv"), sep="\t")
    ref_sets = {g: (tuple(sorted(s[s.split == 0].program)), tuple(sorted(s[s.split == 1].program))) for g, s in gi.groupby("GEP")}
    py_idx = json.load(open(os.path.join(d_py, f"{NAME}_gepindex.json")))
    py_sets = {i + 1: (tuple(sorted(g["Split0"])), tuple(sorted(g["Split1"]))) for i, g in enumerate(py_idx)}
    same_members = sorted(ref_sets.values()) == sorted(py_sets.values())
    same_numbering = ref_sets == py_sets
    same_order = all(list(gi[(gi.GEP == i + 1) & (gi.split == s)].sort_values("pos").program) == g[f"Split{s}"]
                     for i, g in enumerate(py_idx) for s in (0, 1)) if same_numbering else False
    gr = pd.read_csv(os.path.join(d_ref, "community_geps.tsv"), sep="\t", index_col=0)
    gp = pd.read_csv(os.path.join(d_py, f"{NAME}_community_geps.tsv"), sep="\t", index_col=0)
    inv = {v: k for k, v in py_sets.items()}
    worst = np.inf
    if same_members and list(gr.index) == list(gp.index):
        worst = 0.0
        for g, s in ref_sets.items():
            a = gr.iloc[:, g - 1].values; b = gp.iloc[:, inv[s] - 1].values
            worst = max(worst, float(np.max(np.abs(a - b) / (1e-12 * np.abs(a) + 1e-15))))
    say(f"3. Final (rescue {rs}): params {dict(pr)} same {same_params}; membership same {same_members}; numbering same "
        f"{same_numbering}; member order same {same_order}; genes {gr.shape[0]} vs {gp.shape[0]}; GEP score max diff / (1e-12 |x| + 1e-15) = {worst:.2f} (pass < 1)")
    return same_params and same_members and worst < 1


if __name__ == "__main__":
    ok = [check_jaccard() if CHECK_JACCARD else True, check_selection(), check_final(0), check_final(1)]
    status = "OK" if all(ok) else "FAIL"
    say(f"STATUS: {status} | jaccard {ok[0]} | selection {ok[1]} | final rescue0 {ok[2]} rescue1 {ok[3]}")
    open(os.path.join(OUT, "..", f"{os.path.basename(OUT.rstrip('/'))}_equivalence.summary.txt"), "w").write("\n".join(log) + "\n")
    sys.exit(0 if all(ok) else 1)
