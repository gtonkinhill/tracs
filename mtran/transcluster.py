import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson
from numba import jit


def memoize(f):
    memo = {}

    def helper(N, K, delta, lamb, beta, lgamma):
        x = (N, K, delta, lamb, beta)
        if x not in memo:
            memo[x] = f(N, K, delta, lamb, beta, lgamma)
        return memo[x]

    return helper


def calculate_trans_prob(sparse_snp_dist,
                         sample_dates,
                         K,
                         lamb,
                         beta,
                         samplenames=None,
                         log=False):

    # precalculate lgamma
    max_nk = max([t[2] for t in sparse_snp_dist]) + K
    lgamma = gammaln(np.arange(max_nk + 2))

    lprob = np.zeros(len(sparse_snp_dist), dtype=float)
    for c, i, j, d in enumerate(sparse_snp_dist):
        delta = np.abs(sample_dates[samplenames[i]][1] - sample_dates[samplenames[j]][1])
        lp = lprob_transmission(d, K, delta, lamb, beta, lgamma)
        lprob[c] = lp
        # if lp >= lthreshold:
        #     row_ind.append(i)
        #     col_ind.append(j)
        #     lprob.append(lp)
        #     # write out log probabilities if requested.
        #     if outputfile is not None:
        #         outfile.write(",".join([samplenames[i], samplenames[j],
        #                       str(d), str(int(np.ceil(delta*365)))]))
        #         for k in range(K + 1):
        #             lkgN = lprob_k_given_N(d, k, delta, lamb, beta, lgamma)
        #             if not log:
        #                 lkgN = np.exp(lkgN)
        #             outfile.write(
        #                 "," +
        #                 str(lkgN))
        #         outfile.write("\n")

    return lprob


@memoize
def lprob_transmission(N, K, delta, lamb, beta, lgamma):
    lprob = -np.inf
    for k in range(K + 1):
        lprob = np.logaddexp(lprob,
                             lprob_k_given_N(N, k, delta, lamb, beta, lgamma))
    return lprob


@memoize
@jit(nopython=True)
def lprob_k_given_N(N, k, delta, lamb, beta, lgamma):

    if delta > 0:
        lprob = (N + 1) * np.log(lamb) - delta * (
            lamb + beta) + k * np.log(beta) - lgamma[k + 1]

        # ugly poisson cdf but allows for use of numba
        pois_cdf = -np.inf
        for i in range(N + 1):
            pois_cdf = np.logaddexp(i * np.log(lamb * delta) - lgamma[i + 1],
                                    pois_cdf)
        pois_cdf -= lamb * delta
        lprob -= pois_cdf

        integral = -np.inf
        for i in range(N + k + 1):
            integral = np.logaddexp(
                lgamma[N + k + 1] - lgamma[i + 1] - lgamma[N + k - i + 1] +
                (N + k - i) * np.log(delta) + lgamma[i + 1] -
                (i + 1) * np.log(lamb + beta), integral)

        integral -= lgamma[N + 1]
        lprob += integral
    else:
        lprob = (N + 1) * np.log(lamb) + k * np.log(beta) + lgamma[
            N + k +
            1] - lgamma[N + 1] - lgamma[k +
                                        1] - (N + k + 1) * np.log(lamb + beta)
    return (lprob)
