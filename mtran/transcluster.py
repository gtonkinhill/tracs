import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson
from numba import jit

SECONDS_IN_YEAR = 31556952


def memoize(f):
    memo = {}

    def helper(N, K, delta, lamb, beta, lgamma):
        x = (N, K, delta, lamb, beta)
        if x not in memo:
            memo[x] = f(N, K, delta, lamb, beta, lgamma)
        return memo[x]

    return helper

def memoize2(f):
    memo = {}

    def helper(N, delta, lamb, beta, max_k, lgamma):
        x = (N, delta, lamb, beta)
        if x not in memo:
            memo[x] = f(N, delta, lamb, beta, max_k, lgamma)
        return memo[x]

    return helper


def calculate_trans_prob(
    sparse_snp_dist, sample_dates, K, lamb, beta, samplenames=None, log=False
):
    npairs = len(sparse_snp_dist[0])

    # precalculate lgamma
    max_nk = 1e5 #max([t[2] for t in sparse_snp_dist]) + K
    lgamma = gammaln(np.arange(max_nk + 2))

    lprob = np.zeros(npairs, dtype=float)
    expectedk = np.zeros(npairs, dtype=float)
    datediff = np.zeros(npairs, dtype=float)
    for c, i, j, d in zip(range(npairs), sparse_snp_dist[0], sparse_snp_dist[1], sparse_snp_dist[2]):
        delta = np.abs(
            (
                sample_dates[samplenames[i]][1] - sample_dates[samplenames[j]][1]
            ).total_seconds()
            / SECONDS_IN_YEAR
        )
        lprob[c] = lprob_k_given_N(int(d), 0, delta, lamb, beta, lgamma)
        expectedk[c] = expected_k(int(d), delta, lamb, beta, 100, lgamma)
        datediff[c] = delta

    return lprob, expectedk, datediff


@memoize2
def expected_k(N, delta, lamb, beta, max_k, lgamma):
    lprob = -np.inf
    for k in range(1, max_k + 1):
        lprob = np.logaddexp(lprob,
                             lprob_k_given_N(N, k, delta, lamb, beta, lgamma) + np.log(k))
    return np.exp(lprob)

@memoize
@jit(nopython=True)
def lprob_k_given_N(N, k, delta, lamb, beta, lgamma):

    if delta > 0:
        lprob = (
            (N + 1) * np.log(lamb)
            - delta * (lamb + beta)
            + k * np.log(beta)
            - lgamma[k + 1]
        )

        # ugly poisson cdf but allows for use of numba
        pois_cdf = -np.inf
        for i in range(N + 1):
            pois_cdf = np.logaddexp(i * np.log(lamb * delta) - lgamma[i + 1], pois_cdf)
        pois_cdf -= lamb * delta
        lprob -= pois_cdf

        integral = -np.inf
        for i in range(N + k + 1):
            integral = np.logaddexp(
                lgamma[N + k + 1]
                - lgamma[i + 1]
                - lgamma[N + k - i + 1]
                + (N + k - i) * np.log(delta)
                + lgamma[i + 1]
                - (i + 1) * np.log(lamb + beta),
                integral,
            )

        integral -= lgamma[N + 1]
        lprob += integral
    else:
        lprob = (
            (N + 1) * np.log(lamb)
            + k * np.log(beta)
            + lgamma[N + k + 1]
            - lgamma[N + 1]
            - lgamma[k + 1]
            - (N + k + 1) * np.log(lamb + beta)
        )
    return lprob
