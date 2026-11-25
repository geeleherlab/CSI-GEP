"""GEP analysis step of CSI-GEP (Python port of CSI-GEP/GEP_analysis.R).

Selects (K, JL) from the community-number curves, optionally rescues unmatched programs at the next larger JL, and
averages each community's programs into a final GEP.

Selection, as in the R script: for each JL, loess (span 0.3, degree 1) of community number vs K, evaluated at 80
points; elbow of the smoothed curve (elbow package, ahasverus/elbow@cf58784); K = scanned rank closest to the elbow;
slope of a linear fit to the curve's 'benefits' after the elbow. The JL with the most negative slope is selected.

Usage:  python -m csigep.gep_analysis <output_dir/> <name> <gep_dir/> <K list, comma-separated> <rescue 0/1>
Input:  <output_dir>/results/jaccardtest_results_K<K>_JL<JL>.npz (csigep.jaccard), gene_spectra_score files
Output: <gep_dir>/<name>_community_geps.tsv (genes x GEPs), <name>_gepindex.json, <name>_labels.tsv,
        <name>_optimized_parameters.tsv, <name>_jaccardtest_communities.tsv, <name>_selection.tsv, <name>_Curves.pdf,
        output_summary.txt
"""
import argparse
import json
import math
import os

import numpy as np
import pandas as pd

from .jaccard import JL_DEFAULT, calc_split_stats, communities, read_spectra


def r_round(x, digits=0):
    """R (>= 4.0.0) round(): nearest of the two decimal candidates, ties to even (src/nmath/fround.c)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    for n, v in np.ndenumerate(x):
        if not math.isfinite(v) or v == 0.0:
            out[n] = v; continue
        sgn = -1.0 if v < 0 else 1.0; a = abs(v)
        if digits == 0:
            out[n] = sgn * float(np.rint(a)); continue
        l10x = math.log10(2) * (0.5 + math.floor(math.log2(a)))
        if l10x + digits > 15:
            out[n] = sgn * a; continue
        p10 = 10.0 ** digits
        x10 = p10 * a; i10 = math.floor(x10)
        xd = i10 / p10; xu = math.ceil(x10) / p10
        du = xu - a; dd = a - xd
        out[n] = sgn * (xu if (du < dd or (du == dd and math.fmod(i10, 2.0) == 1)) else xd)
    return out if out.ndim else float(out)


def loess_curve(x, y, span=0.3, degree=1, n=80):
    """ggplot2 geom_smooth(method = 'loess'): stats::loess fit, predicted at n points over range(x)."""
    from skmisc.loess import loess
    lo = loess(np.asarray(x, float), np.asarray(y, float), span=span, degree=degree)
    lo.fit()
    lo_x, hi_x = float(np.min(x)), float(np.max(x))
    xs = lo_x + np.arange(n) * ((hi_x - lo_x) / (n - 1))   # R seq(from, to, length.out = n)
    return xs, lo.predict(xs, stderror=False).values


def elbow(X, Y):
    """elbow::elbow(data, plot = FALSE): returns (X_selected, data frame with X, Y, constant, benefits)."""
    d = pd.DataFrame({"X": X, "Y": Y}).sort_values("X", kind="stable").reset_index(drop=True)
    (x1, y1), (x2, y2) = d.iloc[0][["X", "Y"]], d.iloc[-1][["X", "Y"]]
    coef = np.linalg.lstsq(np.array([[1.0, x1], [1.0, x2]]), np.array([y1, y2]), rcond=None)[0]
    d["constant"] = r_round(coef[0] + coef[1] * d["X"].values, 3)
    pos = int(r_round(len(d) / 2)) - 1
    if d.loc[pos, "constant"] < d.loc[pos, "Y"]:                  # concave down
        d["benefits"] = d["Y"].min() + r_round(d["Y"].values - d["constant"].values, 3)
        sel = int(np.argmax(d["benefits"].values))
    else:                                                         # concave up
        d["benefits"] = d["Y"].max() - r_round(d["constant"].values - d["Y"].values, 3)
        sel = int(np.argmin(d["benefits"].values))
    return float(d.loc[sel, "X"]), d


def select_k_jl(com, ranks, jls):
    """Per JL: elbow rank and post-elbow slope; returns (selection table, chosen rows of com)."""
    rows = []
    for jl in jls:
        sub = com[com.Jaccard_Length == jl]
        xs, ys = loess_curve(sub.Rank.values, sub.Community_Number.values)
        xsel, d = elbow(xs, ys)
        opt = ranks[int(np.argmin(np.abs(np.asarray(ranks) - xsel)))]
        after = d[d.X >= xsel]
        slope = np.polyfit(after.X.values, after.benefits.values, 1)[0] if len(after) > 1 else np.nan
        rows.append({"JL": jl, "Rank": opt, "Slope": slope, "X_selected": xsel})
    sel = pd.DataFrame(rows)
    best = sel[sel.JL.isin(sel.JL[sel.Slope == np.nanmin(sel.Slope)])]
    return sel, best


def community_geps(sig_names, split0, split1):
    """Edge-betweenness communities of the significant pairs -> GEP index and averaged GEP scores (genes x GEPs)."""
    left = [a for a, _ in sig_names]; right = [b for _, b in sig_names]
    names = list(dict.fromkeys(left + right))
    idx = {n: i for i, n in enumerate(names)}
    import igraph as ig
    g = ig.Graph(n=len(names), edges=[(idx[a], idx[b]) for a, b in sig_names], directed=False)
    memb = np.asarray(g.community_edge_betweenness(directed=False).as_clustering().membership) + 1
    order = np.argsort(memb, kind="stable")                      # R: sort(membership(ceb))
    gep_index = []
    for c in range(1, memb.max() + 1):
        s0, s1 = [], []
        for i in order[memb[order] == c]:
            num = int(names[i][10:])
            if names[i].startswith("Split0"):
                s0.append(num)
            else:
                s1.insert(0, num)                                # R prepends split1 members
        gep_index.append({"Split0": s0, "Split1": s1})
    S0 = split0.values; S1 = split1.values
    mean = lambda M: np.nanmean(M.astype(np.longdouble), axis=1).astype(float) if M.shape[1] > 1 else M[:, 0]
    vecs = []
    for g_ in gep_index:
        s0 = [i - 1 for i in g_["Split0"]]; s1 = [i - 1 for i in g_["Split1"]]
        if not s0 or not s1:
            vecs.append(mean(S0[:, s0]) if len(s0) > len(s1) else mean(S1[:, s1]))
        else:
            vecs.append(mean(np.column_stack([mean(S0[:, s0]), mean(S1[:, s1])])))
    mat = pd.DataFrame(np.column_stack(vecs), index=split0.index, columns=[f"GEP{i}" for i in range(1, len(vecs) + 1)])
    return gep_index, mat.dropna()


def plot_curves(com, params, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 7))
    for jl in sorted(com.Jaccard_Length.unique()):
        sub = com[com.Jaccard_Length == jl]
        xs, ys = loess_curve(sub.Rank.values, sub.Community_Number.values, span=0.2, degree=2)
        ax.plot(xs, ys, label=str(jl))
    ax.scatter([params["Rank"]], [params["Community_Number"]], color="black", zorder=5)
    ax.set_xlabel("Rank"); ax.set_ylabel("Community_Number"); ax.legend(title="Jaccard_Length", fontsize=8)
    ax.set_title("Community curves", fontsize=10)
    fig.savefig(path); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("output_dir"); ap.add_argument("name"); ap.add_argument("gep_dir")
    ap.add_argument("ranks"); ap.add_argument("rescue", type=int)
    ap.add_argument("--jl", default=",".join(map(str, JL_DEFAULT)))
    a = ap.parse_args(argv)
    ranks = [int(float(x)) for x in a.ranks.split(",")]; jls = [int(x) for x in a.jl.split(",")]
    os.makedirs(a.gep_dir, exist_ok=True)
    out = lambda f: os.path.join(a.gep_dir, f"{a.name}_{f}")
    res = {}
    for k in ranks:
        for jl in jls:
            with np.load(os.path.join(a.output_dir, "results", f"jaccardtest_results_K{k}_JL{jl}.npz")) as z:
                res[(k, jl)] = {"community_number": int(z["community_number"]), "sig": z["sig"].reshape(-1, 2)}
    com = pd.DataFrame([{"Community_Number": int(res[(k, jl)]["community_number"]), "Rank": k, "Jaccard_Length": jl}
                        for k in ranks for jl in jls])
    com.to_csv(out("jaccardtest_communities.tsv"), sep="\t", index=False)
    sel, best = select_k_jl(com, ranks, jls)
    sel.to_csv(out("selection.tsv"), sep="\t", index=False)
    b = best.iloc[0]
    K, JL = int(b.Rank), int(b.JL)
    n_com = int(com[(com.Rank == K) & (com.Jaccard_Length == JL)].Community_Number.iloc[0])
    summary = [f"Community_Number: {n_com}", f"Rank: {K}", f"Jaccard_Length: {JL}"]
    rescue = bool(a.rescue) and len(best) == 1
    sig = res[(K, JL)]["sig"]
    sig_names = [(f"Split0_GEP{r}", f"Split1_GEP{c}") for r, c in sig]
    split0 = read_spectra(a.output_dir, 0, K); split1 = read_spectra(a.output_dir, 1, K)
    if rescue:
        keep0 = [i for i in range(1, K + 1) if i not in set(sig[:, 0])]
        keep1 = [i for i in range(1, K + 1) if i not in set(sig[:, 1])]
        larger = [j for j in jls if j > JL]
        rescue = False
        if keep0 and keep1 and larger:
            r = calc_split_stats(split0.values[:, [i - 1 for i in keep0]], split1.values[:, [i - 1 for i in keep1]], min(larger))
            if len(r["sig"]):
                sig_names += [(f"Split0_GEP{keep0[i - 1]}", f"Split1_GEP{keep1[j - 1]}") for i, j in r["sig"]]
                rescue = True
    summary.append(f"Rescued more programs: {'TRUE' if rescue else 'FALSE'}")
    pd.DataFrame(sig_names, columns=["row", "col"]).to_csv(out("labels.tsv"), sep="\t", index=False)
    gep_index, mat = community_geps(sig_names, split0, split1)
    json.dump(gep_index, open(out("gepindex.json"), "w"), indent=1)
    pd.DataFrame([{"Community_Number": len(gep_index), "Rank": K, "Jaccard_Length": JL}]).to_csv(
        out("optimized_parameters.tsv"), sep="\t", index=False)
    mat.to_csv(out("community_geps.tsv"), sep="\t")
    plot_curves(com, {"Rank": K, "Community_Number": n_com}, out("Curves.pdf"))
    if rescue:
        summary.append(f"Number of programs rescued: {len(gep_index) - n_com}")
    summary += [f"Automatically selected rank k = {K}", "CSIGEP finished!"]
    open(os.path.join(a.gep_dir, "output_summary.txt"), "w").write("\n".join(summary) + "\n")
    print(f"STATUS: OK | K={K} JL={JL} | GEPs={len(gep_index)} | rescue={rescue} | out={a.gep_dir}")


if __name__ == "__main__":
    main()
