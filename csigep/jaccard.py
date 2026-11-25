"""Jaccard step of CSI-GEP (Python port of CSI-GEP/Jaccard.R).

For one rank K, every split0 program is compared with every split1 program by the overlap of their top genes
(Jaccard test, MCA p-value), Bonferroni-corrected over all pairs, at Jaccard lengths (JL) 10, 20, ..., 100.
Significant pairs form a graph whose edge-betweenness communities are the programs found in both halves.

Usage:  python -m csigep.jaccard <K> <output_dir/> [--jl 10,20,...] [--n-jobs 16]
Input:  <output_dir>/cNMF_split{0,1}/cNMF_split{0,1}.gene_spectra_score.k_<K>.dt_0_1.txt
Output: <output_dir>/results/jaccardtest_results_K<K>_JL<JL>.npz, one per JL (pvalues, statistics: split0 x split1;
        sig: significant (split0, split1) pairs, 1-based, in R's order; community_number)
"""
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import igraph as ig
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .jaccard_mca import jaccard_test_mca

JL_DEFAULT = list(range(10, 101, 10))


def read_spectra(output_dir, split, k):
    """genes x programs matrix of gene_spectra_score (as t(vroom(...)) in R; columns = programs 1..K)."""
    f = os.path.join(output_dir, f"cNMF_split{split}", f"cNMF_split{split}.gene_spectra_score.k_{k}.dt_0_1.txt")
    return pd.read_csv(f, sep="\t", index_col=0, float_precision="round_trip").T


def top_genes(spectra, jval):
    """Binary genes x programs matrix: rank(-score) < jval (R rank, ties averaged)."""
    return np.column_stack([rankdata(-spectra[:, i], method="average") < jval for i in range(spectra.shape[1])])


def _pair_row(args):
    a, B = args
    out = [jaccard_test_mca(a, B[:, j]) for j in range(B.shape[1])]
    return [o["pvalue"] for o in out], [o["statistics"] for o in out]


def communities(sig):
    """Edge-betweenness communities of the significant pairs, built as R graph_from_data_frame + cluster_edge_betweenness.

    sig: (n, 2) array of 1-based (split0, split1) program indices in R's order. Returns (vertex names, membership 1-based).
    """
    left = [f"Split0_GEP{r}" for r in sig[:, 0]]
    right = [f"Split1_GEP{c}" for c in sig[:, 1]]
    names = list(dict.fromkeys(left + right))          # R: unique(c(d[, 1], d[, 2]))
    idx = {n: i for i, n in enumerate(names)}
    g = ig.Graph(n=len(names), edges=[(idx[a], idx[b]) for a, b in zip(left, right)], directed=False)
    g.vs["name"] = names
    cl = g.community_edge_betweenness(directed=False).as_clustering()
    return names, np.asarray(cl.membership) + 1


def calc_split_stats(split0, split1, jval, n_jobs=1):
    """Jaccard p-values for all split0 x split1 program pairs at one JL; Bonferroni; significant pairs; community number."""
    A = top_genes(np.asarray(split0, float), jval)
    B = top_genes(np.asarray(split1, float), jval)
    tasks = [(A[:, i], B) for i in range(A.shape[1])]
    if n_jobs > 1:
        with ProcessPoolExecutor(n_jobs) as ex:
            rows = list(ex.map(_pair_row, tasks))
    else:
        rows = [_pair_row(t) for t in tasks]
    pv = np.array([r[0] for r in rows]); st = np.array([r[1] for r in rows])
    fwer = np.minimum(pv * pv.size, 1.0)                 # p.adjust(method = "bonferroni")
    ii, jj = np.nonzero(fwer < 0.05)
    order = np.lexsort((ii, jj))                          # R which(arr.ind = TRUE): column-major order
    sig = np.column_stack([ii[order] + 1, jj[order] + 1]).astype(int)
    n_com = int(communities(sig)[1].max()) if len(sig) else 0
    return {"pvalues": pv, "statistics": st, "sig": sig, "community_number": n_com}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("k", type=int); ap.add_argument("output_dir")
    ap.add_argument("--jl", default=",".join(map(str, JL_DEFAULT)), help="Jaccard lengths (default 10,20,...,100)")
    ap.add_argument("--n-jobs", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    s0 = read_spectra(a.output_dir, 0, a.k); s1 = read_spectra(a.output_dir, 1, a.k)
    os.makedirs(os.path.join(a.output_dir, "results"), exist_ok=True)
    for jl in [int(x) for x in a.jl.split(",")]:
        r = calc_split_stats(s0.values, s1.values, jl, a.n_jobs)
        np.savez(os.path.join(a.output_dir, "results", f"jaccardtest_results_K{a.k}_JL{jl}.npz"), **r)
        print(f"K={a.k} JL={jl}: significant pairs {len(r['sig'])}, communities {r['community_number']}", flush=True)
    print(f"STATUS: OK | K={a.k} | JL={a.jl} | out={os.path.join(a.output_dir, 'results')}")


if __name__ == "__main__":
    main()
