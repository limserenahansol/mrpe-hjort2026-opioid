"""The six learning models of Hjort et al. 2026 (Nature 655:174, Fig. 1k-m), ported line by line from the authors'
notebook BehaviorModels.ipynb (github.com/stuberlab/Hjort-et-al.-2026-PFC-and-reversal-learning, commit e06ae22).

Same equations, same starting values (Q = [0.15, 0.85] for an H->L cue), same softmax choice rule, same bounds and
same grid of starting guesses for scipy.optimize.minimize (L-BFGS-B), same BIC = k ln(T) + 2 negLL.
Only change: the trial loops are compiled with numba so all 40 mice fit in minutes. Quirks kept on purpose:
 - every model except PH/eMack updates only the value of the option that was chosen (lick or no lick);
 - PH: the first parameter is Be although the notebook's start grid calls it gamma (order of unpacking kept);
 - eMack: the learning rate can only grow (alpha <- alpha + theta*|...|, clipped to [0.05, 1]).
c = choice per trial (1 = anticipatory lick, 0 = none), r = reward per trial (1/0).
"""
import itertools
import warnings

import numpy as np
from numba import njit
from scipy.optimize import minimize

warnings.filterwarnings("ignore")


@njit(cache=True)
def _sm(beta, q0, q1):
    # p(choice 0) exactly as written in the notebook (exp form; overflow-safe rewrite of the same number)
    d = beta * (q1 - q0)
    if d > 0:
        e = np.exp(-d)
        return e / (1.0 + e)
    return 1.0 / (1.0 + np.exp(d))


@njit(cache=True)
def nll_mrpe(params, c, r, q0, q1):
    alpha, beta, m, n = params[0], params[1], params[2], params[3]
    Q = np.array([q0, q1])
    md = 0.0
    me = 0.0
    s = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q[0], Q[1])
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        delta = r[t] - Q[c[t]]
        md = (1 - m) * md + m * delta          # CD: rolling average of RPE at rate d (= m)
        me = (1 - n) * me + n * delta          # CE: rolling average of RPE at rate e (= n)
        M = (-md + abs(md)) / 2 + (me + abs(me)) / 2
        Q[c[t]] = Q[c[t]] + alpha * delta * abs(M)
    return -s


@njit(cache=True)
def sim_mrpe(params, c, r, q0, q1):
    alpha, beta, m, n = params[0], params[1], params[2], params[3]
    Q = np.array([q0, q1])
    T = c.shape[0]
    out = np.zeros((T, 6))  # Q_lick before trial, delta, CD, CE, M, p_lick
    md = 0.0
    me = 0.0
    for t in range(T):
        p0 = _sm(beta, Q[0], Q[1])
        out[t, 0] = Q[1]
        out[t, 5] = 1 - p0
        delta = r[t] - Q[c[t]]
        md = (1 - m) * md + m * delta
        me = (1 - n) * me + n * delta
        M = (-md + abs(md)) / 2 + (me + abs(me)) / 2
        Q[c[t]] = Q[c[t]] + alpha * delta * abs(M)
        out[t, 1] = delta
        out[t, 2] = md
        out[t, 3] = me
        out[t, 4] = M
    return out


@njit(cache=True)
def nll_rw(params, c, r, q0, q1):
    alpha, beta = params[0], params[1]
    Q = np.array([q0, q1])
    s = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q[0], Q[1])
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        delta = r[t] - Q[c[t]]
        Q[c[t]] = Q[c[t]] + alpha * delta
    return -s


@njit(cache=True)
def nll_rw2a(params, c, r, q0, q1):
    a1, a2, beta = params[0], params[1], params[2]
    Q = np.array([q0, q1])
    s = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q[0], Q[1])
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        delta = r[t] - Q[c[t]]
        if delta <= 0:
            Q[c[t]] = Q[c[t]] + a1 * delta
        else:
            Q[c[t]] = Q[c[t]] + a2 * delta
    return -s


@njit(cache=True)
def nll_rwt(params, c, r, q0, q1):
    alpha, beta = params[0], params[1]
    Q = np.array([q0, q1])
    s = 0.0
    dprev = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q[0], Q[1])
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        delta = r[t] - Q[c[t]]
        Q[c[t]] = Q[c[t]] + alpha * delta * dprev
        dprev = abs(delta)
    return -s


@njit(cache=True)
def nll_ph(params, c, r, q0, q1):
    Be, Bi, gamma, beta = params[0], params[1], params[2], params[3]
    vy = q1
    vn = q0
    a = 0.05
    Q0, Q1 = q0, q1
    s = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q0, Q1)
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        vnet = vy - vn
        R = r[t] - (vy - vn)
        dvy = 0.0
        dvn = 0.0
        if R > 0:
            dvy = a * Be * r[t]
        elif R < 0:
            dvn = a * Bi * abs(R)
        a = gamma * abs(r[t] - vnet) + (1 - gamma) * a
        vy = vy + dvy
        vn = vn + dvn
        if vy > 1:
            vy = 1.0
        if vn < 0:
            vn = 0.0
        Q0 = vn
        Q1 = vy
    return -s


@njit(cache=True)
def nll_emack(params, c, r, q0, q1):
    Be, Bi, te, ti, beta = params[0], params[1], params[2], params[3], params[4]
    vy = q1
    vn = q0
    a = 0.1
    Q0, Q1 = q0, q1
    s = 0.0
    for t in range(c.shape[0]):
        p0 = _sm(beta, Q0, Q1)
        p = p0 if c[t] == 0 else 1.0 - p0
        s += np.log(max(p, 1e-300))
        R = r[t] - (vy - vn)
        dvy = 0.0
        dvn = 0.0
        if R > 0:
            dvy = a * Be * (1 - vy + vn) * abs(R)
        elif R < 0:
            dvn = a * Bi * (1 - vn + vy) * abs(R)
        if R > 0:
            da = -te * abs(r[t] - vy + vn)
        elif R < 0:
            da = -ti * abs(R * abs(-vn + vy))
        else:
            da = 0.0
        vy = vy + dvy
        vn = vn + dvn
        Q0 = vn
        Q1 = vy
        a = a - da
        if a < 0.05:
            a = 0.05
        if a > 1:
            a = 1.0
    return -s


# start grids and bounds exactly as in the notebook
SPECS = {
    "mRPE": dict(f=nll_mrpe, bounds=((0.01, 10), (1, 10), (0.01, 0.1), (0.01, 0.1)),
                 grid=[np.linspace(0, 1, 2), np.linspace(1, 5, 2), np.linspace(0.01, 0.1, 3), np.linspace(0.01, 0.1, 3)]),
    "RPE2a": dict(f=nll_rw2a, bounds=((0.01, 1), (0.01, 1), (0, 10)),
                  grid=[np.linspace(0, 1, 3), np.linspace(0, 1, 3), np.linspace(1, 5, 3)]),
    "PH": dict(f=nll_ph, bounds=((0.01, 1), (0.05, 1), (0.05, 1), (1, 10)),
               grid=[np.linspace(0.05, 0.9, 5), np.linspace(0.05, 0.9, 3), np.linspace(0.05, 0.9, 3), np.linspace(1, 5, 3)]),
    "eMack": dict(f=nll_emack, bounds=((0.01, 0.3), (0.009, 0.3), (0.002, 0.2), (0.002, 0.19), (1, 5)),
                  grid=[np.linspace(0.01, 0.25, 3), np.linspace(0.009, 0.23, 3), np.linspace(0.01, 0.1, 3),
                        np.linspace(0.009, 0.09, 3), np.linspace(1, 5, 3)],
                  keep=lambda g: g[0] > g[1] and g[2] > g[3]),
    "RPEt-1": dict(f=nll_rwt, bounds=((0.01, 0.5), (0, 50)), grid=[np.linspace(0, 1, 3), np.linspace(1, 5, 3)]),
    "RPE": dict(f=nll_rw, bounds=((0.01, 1), (0, 50)), grid=[np.linspace(0, 1, 3), np.linspace(1, 5, 3)]),
}
ORDER = ["mRPE", "RPE2a", "PH", "eMack", "RPEt-1", "RPE"]


def fit(name, c, r, q0=0.15, q1=0.85):
    sp = SPECS[name]
    lo = np.array([b[0] for b in sp["bounds"]], float)
    hi = np.array([b[1] for b in sp["bounds"]], float)
    best = (np.inf, None)
    for g in itertools.product(*sp["grid"]):
        if "keep" in sp and not sp["keep"](g):
            continue
        x0 = np.clip(np.array(g, float), lo, hi)
        res = minimize(lambda x: sp["f"](x, c, r, q0, q1), x0, bounds=sp["bounds"], method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x.copy())
    k = len(sp["bounds"])
    T = len(c)
    return dict(model=name, negLL=best[0], BIC=k * np.log(T) + 2 * best[0], params=best[1].tolist(), k=k, T=T)
