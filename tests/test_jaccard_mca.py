"""jaccard_test_mca must reproduce R jaccard.test.mca (tests/reference/mca_reference.tsv, from make_mca_reference.R)."""
import os
import numpy as np
import pandas as pd
try:
    import pytest
    _skip = pytest.mark.skipif
except ImportError:                       # runnable as a plain script too
    def _skip(cond, reason):
        return lambda f: f
from csigep.jaccard_mca import jaccard_test_mca

REF = os.path.join(os.path.dirname(__file__), "reference", "mca_reference.tsv")


def _xy(m, nx, ny, c):
    x = np.zeros(m, int); y = np.zeros(m, int)
    x[:nx] = 1; y[:c] = 1; y[nx:nx + ny - c] = 1
    return x, y


@_skip(not os.path.exists(REF), reason="R reference not generated")
def test_matches_r():
    ref = pd.read_csv(REF, sep="\t")
    bad = []
    for r in ref.itertuples():
        out = jaccard_test_mca(*_xy(int(r.m), int(r.nx), int(r.ny), int(r.c)))
        ok = (np.isclose(out["statistics"], r.statistics, rtol=1e-12, atol=1e-15)
              and np.isclose(out["pvalue"], r.pvalue, rtol=1e-9, atol=1e-300)
              and np.isclose(out["expectation"], r.expectation, rtol=1e-14))
        if not ok:
            bad.append((r.m, r.nx, r.ny, r.c, r.pvalue, out["pvalue"]))
    print(f"compared {len(ref)} cases; mismatches {len(bad)}")
    assert not bad, f"{len(bad)}/{len(ref)} mismatches, first: {bad[:5]}"


def test_degenerate():
    assert jaccard_test_mca(np.zeros(10), np.ones(10))["pvalue"] == 1.0


if __name__ == "__main__":     # runnable without pytest: python -m tests.test_jaccard_mca
    test_degenerate(); test_matches_r(); print("STATUS: OK | jaccard_test_mca matches R reference")
