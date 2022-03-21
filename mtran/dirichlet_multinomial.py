# Finds the paramters of the dirichlet multinomial distribution
# by maximising the leave-one-out (LOO) likelihood or fixed point iteration as described in
# Minka, 2000. https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf

import numpy as np
from scipy.special import psi

from numba import jit


def find_dirichlet_priors(data, max_iter=1000, tol=1e-5, method='LOO'):
    K = np.shape(data)[1]

    # Here we make the assumption that the allele frequencies are ordered
    # i.e. if there are multiple strains the combined frequencies will
    # always appear in the same order. 
    data = data[np.count_nonzero(data, 1)>1, ]
    data.sort(axis=1)
    
    # Now fite data
    total_counts = np.sum(data, 1)
    alpha=np.mean(data, 0) + 0.5
    nalpha = np.zeros(K)
    if method=='LOO':
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = alpha[k] * np.sum(data[:,k] / (data[:,k] - 1 + alpha[k]), 0) / np.sum(total_counts/ (total_counts - 1 + a0), 0)
            if np.max(np.abs(nalpha - alpha)) < tol:
                alpha = np.array(nalpha)
                break
            alpha = np.array(nalpha)
    else:
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = (alpha[k] * 
                    np.sum(psi(data[:,k] + alpha[k]) - psi(alpha[k])) / 
                    np.sum(psi(total_counts + a0) - psi(a0)))
            if np.sum(np.abs(nalpha - alpha)) < tol: 
                alpha = np.array(nalpha)
                break
            alpha = np.array(nalpha)


    alpha[::-1].sort()
    print("Calculated alphas: ", alpha)

    return (alpha)

@jit(nopython=True)
def calculate_posteriors(counts, alphas, keep_all, expected_freq_threshold):
    a0 = np.sum(alphas)
    amin = alphas[1]/a0
    for i in range(counts.shape[0]):
        if keep_all:
            k = counts[i,:]>0
        
        if np.sum(counts[i,:]) == 0:
            counts[i,:] = amin
        else:
            denom = np.sum(counts[i,:]) + a0
            for a, c in enumerate(np.sort(np.unique(counts[i,:]))[::-1]):
                a = alphas[a]
                for j in range(4):
                    if counts[i,j]==c:
                        counts[i,j] = (counts[i,j] + a)/denom

        if keep_all:
            counts[i,:][(counts[i,:]<expected_freq_threshold) & k] = expected_freq_threshold

    return

