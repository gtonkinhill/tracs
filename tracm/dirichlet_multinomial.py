# Finds the paramters of the dirichlet multinomial distribution
# by maximising the leave-one-out (LOO) likelihood or fixed point iteration as described in
# Minka, 2000. https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf

import numpy as np
from scipy.special import psi


def find_dirichlet_priors(counts, max_iter=1000, tol=1e-5, method="FPI", error_filt_threshold=None):
    data =  np.ndarray.copy(counts)
    K = np.shape(data)[1]

    if error_filt_threshold is not None:
        freq = data.T/np.sum(data,1)
        data[freq.T < error_filt_threshold] = 0

    # Here we make the assumption that the allele frequencies are ordered
    # i.e. if there are multiple strains the combined frequencies will
    # always appear in the same order.
    nz = np.count_nonzero(data, 1)
    # if np.sum(nz>2) >= 100:
    #     keep_rows_thresh = 2
    # el
    
    if np.sum(nz>1) > 5:
        keep_rows_thresh = 1
    else:
        return(np.array([0,0,0,1.0]))
        print("Less than 5 polymorphic loci found!")
        print("Calculated alphas: ", alpha)

    data = data[
        nz > keep_rows_thresh,
    ]
    data.sort(axis=1)

    # Now fit data
    total_counts = np.sum(data, 1)
    alpha = np.mean(data, 0) + 0.5
    nalpha = np.zeros(K)
    if method == "LOO":
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = (
                    alpha[k]
                    * np.sum(data[:, k] / (data[:, k] - 1 + alpha[k]), 0)
                    / np.sum(total_counts / (total_counts - 1 + a0), 0)
                )
            if np.max(np.abs(nalpha - alpha)) < tol:
                alpha = np.array(nalpha)
                break
            alpha = np.array(nalpha)
    else:
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = (
                    alpha[k]
                    * np.sum(psi(data[:, k] + alpha[k]) - psi(alpha[k]))
                    / np.sum(psi(total_counts + a0) - psi(a0))
                )
            if np.sum(np.abs(nalpha - alpha)) < tol:
                alpha = np.array(nalpha)
                break
            alpha = np.array(nalpha)
            alpha[alpha<1e-16] = 1e-16

    alpha[::-1].sort()
    print("Calculated alphas: ", alpha)

    return alpha
