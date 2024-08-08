import numpy as np
from TRACS import trans_dist
from datetime import date

SECONDS_IN_YEAR = 31556952.0


def calculate_trans_prob(
    sparse_snp_dist,
    sample_dates,
    K,
    lamb,
    beta,
    samplenames=None,
    log=False,
    precision=0.01,
):
    i = np.array(sparse_snp_dist[0])
    j = np.array(sparse_snp_dist[1])
    d = np.array(sparse_snp_dist[2], dtype=int)

    npairs = len(sparse_snp_dist[0])
    nsamples = max(max(sparse_snp_dist[0]), max(sparse_snp_dist[1]))

    # compute difference in sampling times
    reftime = date.fromisoformat("1970-01-01")
    time_array = np.array(
        [
            (sample_dates[samplenames[s]][1] - reftime).total_seconds()
            for s in range(nsamples + 1)
        ]
    )
    time_diff = np.abs(time_array[i] - time_array[j]) / SECONDS_IN_YEAR

    # calculate transmission distances
    p0, eK = trans_dist(d, time_diff, lamb, beta, precision)

    if not log:
        p0 = np.exp(p0)

    return p0, eK, time_diff
