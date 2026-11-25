"""Convert an exported R Jaccard result (sig_pairs.tsv, communities.tsv from export_demo_reference.R) to the
Python npz layout, so csigep.gep_analysis can be tested on runs whose Jaccard step was done in R.
Usage: python -m tests.r_jaccard_to_npz <ref_dir/> <out_dir/>   -> <out_dir>/results/jaccardtest_results_K<K>_JL<JL>.npz
"""
import os
import sys
import numpy as np
import pandas as pd

ref, out = sys.argv[1], sys.argv[2]
sp = pd.read_csv(os.path.join(ref, "sig_pairs.tsv"), sep="\t")
cm = pd.read_csv(os.path.join(ref, "communities.tsv"), sep="\t")
os.makedirs(os.path.join(out, "results"), exist_ok=True)
for r in cm.itertuples():
    s = sp[(sp.K == r.K) & (sp.JL == r.JL)].sort_values("ord")[["row", "col"]].values.astype(int).reshape(-1, 2)
    np.savez(os.path.join(out, "results", f"jaccardtest_results_K{r.K}_JL{r.JL}.npz"), sig=s, community_number=r.community_number)
print(f"STATUS: OK | converted {len(cm)} (K, JL) | out={out}results")
