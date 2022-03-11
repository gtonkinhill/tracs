import os
import gzip
import numpy as np
from collections import defaultdict
from tqdm import tqdm
from mtrans.dirichlet_multinomial import find_dirichlet_priors

MAX_LEN = 1000000


def pileup_dist(pileups,
                max_dist=9999999999999,
                quiet=False):

    # load pileups generated using the bampileup function
    if not quiet:
        print('Loading pileups...')

    sample_arrays = []
    sample_names = []
    for i, f in tqdm(enumerate(pileups)):
        with gzip.open(f, 'rt') as infile:
            for line in infile:
                if line.strip()=='': continue
                line = line.strip().split(',')
                sample_names.append(line[0])
                a = np.array(line[1:]).astype(float)
                a = np.reshape(a, (int(a.size/4),4), order='C')
                sample_arrays.append(a)
    
    n_samples = len(sample_names)
    ref_length = sample_arrays[0].shape[0]

    # calculate pairwise distances
    if not quiet:
        print('calculating distances...')

    distances = []
    for i in tqdm(range(n_samples)):
        for j in range(i + 1, n_samples):
            d = ref_length - np.count_nonzero(np.sum(sample_arrays[i] * sample_arrays[j], 1))
            if d <= max_dist:
                distances.append((i, j, int(d)))

    sample_to_index = {}
    for i, s in enumerate(sample_names):
        sample_to_index[s] = i

    return (sample_names, sample_to_index, distances)
