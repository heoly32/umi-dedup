from __future__ import division
import numpy as np
import collections
import MCMC_algorithm

DEFAULT_NSAMP = 1000
DEFAULT_NTHIN = 1
DEFAULT_NBURN = 200

def deduplicate_counts (umi_counts, nsamp=DEFAULT_NSAMP, nthin=DEFAULT_NTHIN, nburn=DEFAULT_NBURN, uniform=True, total_counts = None):

    # Remove zeros from data, to shorten the vector
    data = []
    if uniform:
        for value in umi_counts.values():
            if value > 0:
                data.append(value)
    else:
        total_data = []
        for key, value in umi_counts.items():
            if value > 0:
                data.append(value)
                total_data.append(total_counts[key])

    n = len(data)
    N = sum(data)

    # Set priors for the different parameters
    pi_prior = [1., 1.]
    S_prior = [1.] * n
    if uniform:
        # The 'uniform' algorithm assumes equi-probability for all tags, before amplification
        C_prior = [1./n] * n
    else:
        # The non-uniform algorithm illicits prior from data
        C_prior = [float(total_data[j])/N for j in range(n)]

    # Run Gibbs sampler
    p_post = MCMC_algorithm.MCMC_algorithm(data, \
                                           n, N, \
                                           S_prior, C_prior, pi_prior, \
                                           nsamp, nthin, nburn, \
                                           True)

    # Compute median for each tag
    median_list = [0] * n
    for i in range(n):
        median_list[i] = computeMedian(p_post[i::n])

    # Return ordered dictionary with estimated number of true molecules
    umi_true = collections.OrderedDict()
    index = 0
    for key, value in umi_counts.items():
        if value == 0:
            umi_true[key] = value
        else:
            umi_true[key] = int(np.ceil(median_list[index] * data[index]))
            index += 1

    return umi_true

def computeMedian(list):
    list.sort()
    lens = len(list)
    if lens % 2 != 0:
        midl = (lens / 2)
        res = list[midl]
    else:
        odd = int((lens / 2) -1)
        ev = int((lens / 2))
        res = float(list[odd] + list[ev]) / float(2)
    return res
