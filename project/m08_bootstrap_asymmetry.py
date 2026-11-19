"""Step 8: is learning slower after morphine -> water (negative errors) than after water -> morphine (positive errors)?
Group fits with shared learning rates (as in m05b), bootstrapped over mice (1000 resamples per group):
RPE2a (alpha_minus for negative RPE, alpha_plus for positive RPE). Also Active vs Passive difference in alpha_minus.
Output: data/m08_bootstrap_asymmetry.json
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m05_recovery_and_group_fits as m5  # noqa: E402


def one(args):
    g, seed = args
    rng = np.random.default_rng(seed)
    ms = list(rng.choice(m5.MICE[g], size=len(m5.MICE[g]), replace=True))
    f = m5.fit_group("RPE2a", ms)
    return dict(group=g, seed=seed, a_minus=f["params"][0], a_plus=f["params"][1])


if __name__ == "__main__":
    full = {g: m5.fit_group("RPE2a", m5.MICE[g]) for g in ("Active", "Passive")}
    with ProcessPoolExecutor(max_workers=30) as ex:
        B = pd.DataFrame(list(ex.map(one, [(g, 5000 + s) for g in ("Active", "Passive") for s in range(1000)])))
    out = {}
    for g in ("Active", "Passive"):
        b = B[B.group == g]
        out[g] = dict(full_fit=dict(a_minus=round(full[g]["params"][0], 4), a_plus=round(full[g]["params"][1], 4)),
                      a_minus_CI=b.a_minus.quantile([.025, .5, .975]).round(4).tolist(),
                      a_plus_CI=b.a_plus.quantile([.025, .5, .975]).round(4).tolist(),
                      frac_a_minus_less_than_a_plus=float((b.a_minus < b.a_plus).mean()),
                      log10_ratio_plus_over_minus_CI=np.log10(b.a_plus / b.a_minus).quantile([.025, .5, .975]).round(3).tolist())
    a = B[B.group == "Active"].a_minus.values
    p = B[B.group == "Passive"].a_minus.values
    out["Active_minus_Passive_a_minus_CI"] = np.quantile(a - p, [.025, .5, .975]).round(4).tolist()
    out["frac_Active_a_minus_lower"] = float((a < p).mean())
    json.dump(out, open(HERE / "data" / "m08_bootstrap_asymmetry.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
