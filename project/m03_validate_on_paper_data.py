"""Step 3: check the port. Fit the six models to the authors' own 40 mice (Sample_Choice_Data, H->L cue, 160 trials:
pre 0-40, reversal at trial 60, post 120-160) and compare with Fig. 1k,l of the paper
(mRPE best on negLL, P < 0.0001 vs all; on BIC still best: RPE2a P = 0.0064, PH < 0.0001, eMack < 0.0001,
RPEt-1 P = 0.0258, RPE P = 0.0339; RM one-way ANOVA + Holm-Sidak in Prism).
Output: data/m03_paper_data_fits.csv, data/m03_paper_data_summary.json
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m02_hjort_models as hm  # noqa: E402

SRC = HERE / "data" / "hjort2026_sample_choice"


def one(animal):
    c = np.load(SRC / animal / "HL_choice.npy", allow_pickle=False).astype(np.int64)
    r = np.load(SRC / animal / "HL_outcome.npy", allow_pickle=False).astype(np.float64)
    return [dict(animal=animal, **hm.fit(m, c, r)) for m in hm.ORDER]


if __name__ == "__main__":
    animals = sorted(p.name for p in SRC.iterdir() if p.is_dir())
    with ProcessPoolExecutor(max_workers=30) as ex:
        rows = [x for res in ex.map(one, animals) for x in res]
    F = pd.DataFrame(rows)
    F.to_csv(HERE / "data" / "m03_paper_data_fits.csv", index=False)
    out = {"n_animals": len(animals)}
    for metric in ("negLL", "BIC"):
        W = F.pivot(index="animal", columns="model", values=metric)[hm.ORDER]
        chg = W.sub(W.min(axis=1), axis=0)          # 'change from best model', as in Fig. 1k,l
        res = {"mean_change_from_best": chg.mean().round(2).to_dict(),
               "sem": (chg.std(ddof=1) / np.sqrt(len(chg))).round(2).to_dict(),
               "n_best": W.idxmin(axis=1).value_counts().to_dict(), "vs_mRPE": {}}
        ps = []
        for m in hm.ORDER[1:]:
            d = W[m] - W["mRPE"]
            t = stats.ttest_rel(W[m], W["mRPE"])
            w = stats.wilcoxon(W[m], W["mRPE"])
            res["vs_mRPE"][m] = dict(mean_diff=round(float(d.mean()), 2), n_mRPE_better=int((d > 0).sum()),
                                     p_t=float(t.pvalue), p_wilcoxon=float(w.pvalue))
            ps.append(t.pvalue)
        # Holm over the 5 comparisons with mRPE
        order = np.argsort(ps)
        adj = np.empty(5)
        run = 0
        for i, j in enumerate(order):
            run = max(run, min(1, (5 - i) * ps[j]))
            adj[j] = run
        for m, a in zip(hm.ORDER[1:], adj):
            res["vs_mRPE"][m]["p_t_holm"] = float(a)
        out[metric] = res
    P = F[F.model == "mRPE"].params.apply(pd.Series)
    P.columns = ["alpha", "beta", "d_CD", "e_CE"]
    out["mRPE_params_median"] = P.median().round(3).to_dict()
    out["mRPE_d_less_than_e_n"] = int((P.d_CD < P.e_CE).sum())
    json.dump(out, open(HERE / "data" / "m03_paper_data_summary.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
