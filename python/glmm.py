"""Negative-binomial GLMM with a per-larva random intercept, fitted by
maximum likelihood using adaptive Gauss-Hermite quadrature (AGHQ).

log mu_ij = offset_ij + X_ij beta + b_i,   b_i ~ N(0, sigma^2)
y_ij | b_i ~ NB2(mu_ij, theta)            Var = mu + mu^2/theta
"""
import numpy as np
from scipy import optimize, special, stats


def load(path):
    import pandas as pd
    raw = pd.read_csv(path, sep='\t', header=None, dtype=str)
    groups = raw.iloc[0, 1:].tolist()
    tps = raw.iloc[1:, 0].tolist()
    rows = []
    for r, tp in enumerate(tps):
        for c in range(96):
            s = raw.iloc[r + 1, c + 1]
            if s is None or (isinstance(s, float)) or str(s).strip() in ('', 'nan'):
                continue
            s = str(s).strip()
            col = c + 1
            plate = 'A' if (col <= 24 or 49 <= col <= 72) else 'B'
            rows.append(dict(larva=col, group='TBI' if col > 48 else 'Sham',
                             plate=plate, tp=tp, raw=s, value=float(s)))
    df = pd.DataFrame(rows)
    # a recording block (plate x timepoint x group) is "scaled" if any value
    # in it was written with a decimal point (values are multiples of 2/3)
    df['scaled'] = df.groupby(['plate', 'tp', 'group'])['raw'].transform(
        lambda s: s.str.contains(r'\.').any())
    df['count'] = np.where(df.scaled, np.round(df.value * 1.5), np.round(df.value)).astype(int)
    df['bins'] = np.where(df.scaled, 180, 120)
    df['session'] = df.plate + ':' + df.tp
    df['tp'] = pd.Categorical(df.tp, categories=tps, ordered=True)
    return df, tps


def nb_logpmf(y, mu, theta):
    return (special.gammaln(y + theta) - special.gammaln(theta) - special.gammaln(y + 1)
            + theta * (np.log(theta) - np.log(theta + mu)) + y * (np.log(mu) - np.log(theta + mu)))


class NBGLMM:
    def __init__(self, y, X, offset, cluster, nq=15):
        self.y = np.asarray(y, float)
        self.X = np.asarray(X, float)
        self.off = np.asarray(offset, float)
        cl = np.asarray(cluster)
        self.ucl, self.cidx = np.unique(cl, return_inverse=True)
        self.K = len(self.ucl)
        self.p = self.X.shape[1]
        self.gh_x, self.gh_w = special.roots_hermite(nq)  # for exp(-x^2)

    def _cluster_sum(self, v):
        return np.bincount(self.cidx, weights=v, minlength=self.K)

    def _modes(self, eta, theta, sigma):
        # Newton for mode of h(b) = sum_j log f(y|eta+b) - b^2/(2 sigma^2), per cluster
        b = np.zeros(self.K)
        for _ in range(50):
            mu = np.exp(eta + b[self.cidx])
            # d/d eta of NB loglik: theta (y - mu) / (theta + mu)
            g1 = theta * (self.y - mu) / (theta + mu)
            g2 = -theta * mu * (theta + self.y) / (theta + mu) ** 2
            G = self._cluster_sum(g1) - b / sigma ** 2
            H = self._cluster_sum(g2) - 1 / sigma ** 2
            step = G / H
            b = b - step
            if np.max(np.abs(step)) < 1e-10:
                break
        return b, -H

    def _quad(self, params):
        beta = params[:self.p]
        theta = np.exp(params[self.p])
        sigma = np.exp(params[self.p + 1])
        eta = self.off + self.X @ beta
        bhat, prec = self._modes(eta, theta, sigma)
        s = 1 / np.sqrt(prec)
        bq = bhat[None, :] + np.sqrt(2) * s[None, :] * self.gh_x[:, None]      # (Q, K)
        terms = np.empty_like(bq)
        for q in range(len(self.gh_x)):
            mu = np.exp(eta + bq[q][self.cidx])
            terms[q] = (np.log(self.gh_w[q]) + self.gh_x[q] ** 2
                        + self._cluster_sum(nb_logpmf(self.y, mu, theta))
                        + stats.norm.logpdf(bq[q], 0, sigma))
        lse = special.logsumexp(terms, axis=0)
        return eta, theta, sigma, bq, terms, lse, s

    def loglik(self, params):
        eta, theta, sigma, bq, terms, lse, s = self._quad(params)
        return float(np.sum(lse + np.log(np.sqrt(2) * s)))

    def grad(self, params):
        eta, theta, sigma, bq, terms, lse, s = self._quad(params)
        W = np.exp(terms - lse[None, :])          # posterior weights per cluster (Q, K)
        gb = np.zeros(self.p); gth = 0.0; gsg = 0.0
        for q in range(len(self.gh_x)):
            mu = np.exp(eta + bq[q][self.cidx])
            wq = W[q][self.cidx]
            r = theta * (self.y - mu) / (theta + mu)
            gb += self.X.T @ (wq * r)
            dth = (special.digamma(self.y + theta) - special.digamma(theta) + np.log(theta) + 1
                   - np.log(theta + mu) - (theta + self.y) / (theta + mu))
            gth += np.sum(wq * dth) * theta
            gsg += np.sum(W[q] * (-1 + bq[q] ** 2 / sigma ** 2))
        return np.r_[gb, gth, gsg]

    def fit(self, start=None, cov=True):
        if start is None:
            z = np.log(self.y + 0.5) - self.off
            b0 = np.linalg.lstsq(self.X, z, rcond=None)[0]
            start = np.r_[b0, np.log(1.0), np.log(1.0)]
        f = lambda p: -self.loglik(p)
        g = lambda p: -self.grad(p)
        res = optimize.minimize(f, start, jac=g, method='BFGS', options=dict(maxiter=10000, gtol=1e-5))
        res = optimize.minimize(f, res.x, jac=g, method='L-BFGS-B', options=dict(maxiter=10000, gtol=1e-8, ftol=1e-14))
        self.res = res
        self.params = res.x
        self.llf = -res.fun
        self.gnorm = np.max(np.abs(g(res.x)))
        self.cov = self._num_cov(res.x) if cov else None
        return self

    def _num_cov(self, x, h=1e-5):
        n = len(x)
        H = np.zeros((n, n))
        for i in range(n):
            e = np.zeros(n); e[i] = h
            H[:, i] = -(self.grad(x + e) - self.grad(x - e)) / (2 * h)
        H = (H + H.T) / 2
        return np.linalg.inv(H)

    @property
    def theta(self):
        return np.exp(self.params[self.p])

    @property
    def sigma(self):
        return np.exp(self.params[self.p + 1])


def simulate(X, off, cidx, beta, theta, sigma, rng):
    K = cidx.max() + 1
    b = rng.normal(0, sigma, K)
    mu = np.exp(off + X @ beta + b[cidx])
    p = theta / (theta + mu)
    return rng.negative_binomial(theta, p)
