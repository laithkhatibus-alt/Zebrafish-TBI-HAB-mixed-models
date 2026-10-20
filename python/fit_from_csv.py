"""Fit every model in Supplemental Fig. 8 from the long-format CSVs.
Usage: python fit_from_csv.py <data_dir> <out_csv>
Needs glmm.py (NBGLMM: negative-binomial GLMM, random intercept per larva, adaptive Gauss-Hermite quadrature)."""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from glmm import NBGLMM

def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); adj = np.empty(m); run = 0
    for k, i in enumerate(o):
        run = max(run, (m - k) * p[i]); adj[i] = min(1, run)
    return adj

def prep(d):
    tps = d.sort_values('tp_order').tp.unique().tolist()
    T = pd.get_dummies(pd.Categorical(d.tp, categories=tps)).astype(float).values
    S = pd.get_dummies(d.session).astype(float).values
    return tps, T, S

def fit(d, X):
    return NBGLMM(d['count'].values, X, np.log(d.bins.values / 120), d.larva.values).fit()

def contrast(m, idx_sets, tps):
    rows = []
    for t in range(len(tps)):
        c = np.zeros(len(m.params))
        for idx in idx_sets: c[idx[t]] = 1
        e = c @ m.params; se = np.sqrt(c @ m.cov @ c)
        rows.append((e, se))
    e, se = np.array(rows).T; p = 2 * stats.norm.sf(np.abs(e / se))
    return pd.DataFrame(dict(tp=tps, RR=np.exp(e), lo=np.exp(e - 1.96 * se), hi=np.exp(e + 1.96 * se), p=p, p_holm=holm(p)))

def injury_by_time(d):
    tps, T, S = prep(d); tbi = d.tbi.values[:, None]
    m0 = fit(d, np.c_[S, tbi]); m1 = fit(d, np.c_[S, T * tbi])
    lr = 2 * (m1.llf - m0.llf); df = len(tps) - 1
    k = S.shape[1]
    tab = contrast(m1, [list(range(k, k + len(tps)))], tps); tab['contrast'] = 'TBI vs sham'
    return tab, dict(test='injury x time', chi2=lr, df=df, p=stats.chi2.sf(lr, df), larvae=d.larva.nunique(), n=len(d))

def injury_by_drug(d):
    tps, T, S = prep(d); tbi = d.tbi.values[:, None]; dr = d.drug.values[:, None].astype(float); n = len(tps)
    m_no = fit(d, np.c_[S, T * tbi, T * dr]); m_full = fit(d, np.c_[S, T * tbi, T * dr, T * tbi * dr])
    lr = 2 * (m_full.llf - m_no.llf)
    k = S.shape[1]; I_tbi = list(range(k, k + n)); I_dr = list(range(k + n, k + 2 * n)); I_int = list(range(k + 2 * n, k + 3 * n))
    tabs = []
    for name, sets in [('TBI vs sham (vehicle)', [I_tbi]), ('Drug vs vehicle (injured)', [I_dr, I_int]), ('Drug vs vehicle (sham)', [I_dr])]:
        t = contrast(m_full, sets, tps); t['contrast'] = name; tabs.append(t)
    return pd.concat(tabs), dict(test='injury x drug (time-varying)', chi2=lr, df=n, p=stats.chi2.sf(lr, n), larvae=d.larva.nunique(), n=len(d))

if __name__ == '__main__':
    D, out = sys.argv[1], sys.argv[2]
    jobs = [('A', 'Fig7', 'fig7', injury_by_time), ('B', 'Fig8', 'fig8_pre', injury_by_time), ('B', 'Fig8', 'fig8_post', injury_by_drug),
            ('C', 'Fig9', 'fig9_pre', injury_by_time), ('C', 'Fig9', 'fig9_post', injury_by_drug),
            ('D', 'Fig10', 'fig10_pre', injury_by_time), ('D', 'Fig10', 'fig10_post', injury_by_drug)]
    tabs, tests = [], []
    for panel, exp, key, f in jobs:
        d = pd.read_csv(f'{D}/{key}_long.csv', dtype={'plate': str})
        tab, test = f(d)
        seg = 'whole series' if key == 'fig7' else key.split('_')[1] + '-intervention'
        tab.insert(0, 'segment', seg); tab.insert(0, 'experiment', exp); tab.insert(0, 'panel', panel); tabs.append(tab)
        tests.append(dict(panel=panel, experiment=exp, segment=seg, **test))
        print(f'{panel} {exp:5s} {seg:18s} {test["test"]:28s} chi2({test["df"]})={test["chi2"]:.2f} p={test["p"]:.3g} larvae={test["larvae"]} n={test["n"]}', flush=True)
    pd.concat(tabs).to_csv(out, index=False)
    pd.DataFrame(tests).to_csv(out.replace('.csv', '_tests.csv'), index=False)
