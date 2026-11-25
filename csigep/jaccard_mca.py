"""Jaccard similarity test by the measure concentration algorithm (MCA).

Python port of jaccard::jaccard.test.mca (ncchung/jaccard@ad1d8fc; Chung et al. 2019), which CSI-GEP uses to test
the overlap of two programs' top-gene lists. The C++ enumeration (jaccard_mca_rcpp.cpp, adapted from IsoSpec) is
reproduced step by step: same mode search, same layer thresholds (-0.8 per layer), same neighbour order, same
Kahan summation and stopping rule. The p-value is therefore the same approximation R computes, not an exact one.
"""
import ctypes
import ctypes.util
import math

import numpy as np
from numba import njit

_libm = ctypes.CDLL(ctypes.util.find_library("m") or "libm.so.6")
_libm.lgamma.restype = ctypes.c_double
_libm.lgamma.argtypes = [ctypes.c_double]


def _lfactorials(m):
    """log(i!) for i = 0..m with the C library's lgamma, as the C++ code does."""
    return np.array([_libm.lgamma(i + 1.0) for i in range(m + 1)], dtype=np.float64)


_LFACT_CACHE = {}


@njit(cache=True)
def _logprob(c0, c1, c2, c3, lp, lf, nom):
    res = 0.0
    res -= lf[c0]; res += c0 * lp[0]
    res -= lf[c1]; res += c1 * lp[1]
    res -= lf[c2]; res += c2 * lp[2]
    res -= lf[c3]; res += c3 * lp[3]
    return res + nom


@njit(cache=True)
def _mca(px, py, m, t_obs, accuracy, known_p, lf):
    probs = np.array([px * py, px * (1 - py), py * (1 - px), (1 - px) * (1 - py)])
    lp = np.log(probs)
    nom = lf[m]
    # --- mode (Marginal::setupInitialConf) ---
    res = np.empty(4, dtype=np.int64)
    for i in range(4):
        res[i] = np.int64(m * probs[i]) + 1
    diff = m - res.sum()
    if diff > 0:
        res[0] += diff
    if diff < 0:
        diff = -diff
        i = 0
        while diff > 0:
            cd = res[i] - diff
            if cd >= 0:
                res[i] -= diff
                diff = 0
            else:
                res[i] = 0
                i += 1
                diff = abs(cd)
    LP = _logprob(res[0], res[1], res[2], res[3], lp, lf, nom)
    modified = True
    while modified:
        modified = False
        for ii in range(4):
            for jj in range(4):
                if ii != jj and res[ii] > 0:
                    res[ii] -= 1
                    res[jj] += 1
                    NLP = _logprob(res[0], res[1], res[2], res[3], lp, lf, nom)
                    if NLP > LP or (NLP == LP and ii > jj):
                        modified = True
                        LP = NLP
                    else:
                        res[ii] += 1
                        res[jj] -= 1
    M1 = m + 1
    # configurations are int64 keys (c0, c1, c2); c3 = m - c0 - c1 - c2
    mode_key = (res[0] * M1 + res[1]) * M1 + res[2]
    fringe = [mode_key]; fringe.pop()            # typed empty list
    new_fringe = [mode_key]
    visited = {mode_key}; visited.clear()
    expectation = px * py / (px + py - px * py)
    cur_thr = 0.0
    new_thr = 1.0
    s_sum = 0.0; s_c = 0.0                       # Kahan summators
    o_sum = 0.0; o_c = 0.0
    layer_thr = -0.8
    conf = np.empty(4, dtype=np.int64)
    while s_sum < 1 - accuracy and not (len(fringe) == 0 and len(new_fringe) == 0):
        # LayeredMarginal::extend
        if len(new_fringe) > 0:
            cur_thr = new_thr
            new_thr = layer_thr
            visited.clear()
            for k in new_fringe:
                visited.add(k)
            tmp = fringe
            fringe = new_fringe
            new_fringe = tmp
            new_fringe.clear()
        while s_sum < 1 - accuracy:
            # LayeredMarginal::next
            found = False
            opc = 0.0
            while len(fringe) > 0:
                key = fringe.pop()
                conf[2] = key % M1
                conf[1] = (key // M1) % M1
                conf[0] = key // (M1 * M1)
                conf[3] = m - conf[0] - conf[1] - conf[2]
                opc = _logprob(conf[0], conf[1], conf[2], conf[3], lp, lf, nom)
                if opc < new_thr:
                    new_fringe.append(key)
                else:
                    found = True
                    break
            if not found:
                break
            for ii in range(4):
                for jj in range(4):
                    if ii != jj and conf[jj] > 0:
                        conf[ii] += 1
                        conf[jj] -= 1
                        lpc = _logprob(conf[0], conf[1], conf[2], conf[3], lp, lf, nom)
                        nk = (conf[0] * M1 + conf[1]) * M1 + conf[2]
                        if (nk not in visited) and lpc < cur_thr and (opc > lpc or (opc == lpc and ii > jj)):
                            visited.add(nk)
                            if lpc >= new_thr:
                                fringe.append(nk)
                            else:
                                new_fringe.append(nk)
                        conf[ii] -= 1
                        conf[jj] += 1
            # statistic of the current configuration
            if conf[3] == m:
                stat = 0.0
            elif known_p:
                stat = conf[0] / (m - conf[3]) - expectation
            else:
                cpx = (conf[0] + conf[1]) / m
                cpy = (conf[0] + conf[2]) / m
                cpxpy = cpx * cpy
                stat = conf[0] / (m - conf[3]) - cpxpy / (cpx + cpy - cpxpy)
            p = math.exp(opc)
            if abs(stat) >= abs(t_obs):
                y = p - o_c; t = o_sum + y; o_c = (t - o_sum) - y; o_sum = t
            y = p - s_c; t = s_sum + y; s_c = (t - s_sum) - y; s_sum = t
        layer_thr -= 0.8
    return o_sum, s_sum


def jaccard_test_mca(x, y, px=None, py=None, accuracy=1e-05, error_type="average"):
    """Same arguments and result as R jaccard::jaccard.test.mca: dict(statistics, pvalue, expectation, accuracy)."""
    x = np.asarray(x).astype(bool); y = np.asarray(y).astype(bool)
    if x.shape != y.shape:
        raise ValueError("Length mismatch")
    m = x.size
    known_p = not (px is None or py is None)
    if not known_p:
        px = x.sum() / m; py = y.sum() / m
    expectation = px * py / (px + py - px * py)
    union = np.logical_or(x, y).sum()
    j_obs = np.logical_and(x, y).sum() / union - expectation if union > 0 else math.nan
    if px == 1 or py == 1 or x.sum() == m or y.sum() == m or px == 0 or py == 0 or x.sum() == 0 or y.sum() == 0:
        return {"statistics": 0.0, "pvalue": 1.0, "expectation": expectation}
    if m not in _LFACT_CACHE:
        _LFACT_CACHE[m] = _lfactorials(m)
    pv, acc = _mca(float(px), float(py), m, float(j_obs), float(accuracy), known_p, _LFACT_CACHE[m])
    pvalue = {"lower": pv, "average": pv / acc, "upper": pv + 1 - acc}[error_type]
    return {"statistics": float(j_obs), "pvalue": pvalue, "expectation": float(expectation), "accuracy": 1 - acc}
