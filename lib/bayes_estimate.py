from __future__ import division
import numpy as np
import collections
import sys
from . import MCMC_algorithm_pi, umi_data, apportion_counts

DEFAULT_NSAMP = 1000
DEFAULT_NTHIN = 1
DEFAULT_NBURN = 200
DEFAULT_ALPHA1 = 1.5
DEFAULT_ALPHA2 = 0.1

def compute_prior (umi_counts):
  denom = sum(umi_counts.nonzero_values())
  count_iter = umi_counts.nonzero_items()
  (first_umi, first_count) = next(count_iter)
  result = umi_data.UmiValues([(first_umi, first_count / denom)])
  for umi, count in count_iter: result[umi] = count / denom
  return result

def deduplicate_counts (umi_counts, nsamp=DEFAULT_NSAMP, nthin=DEFAULT_NTHIN, nburn=DEFAULT_NBURN, uniform=True, alpha2 = DEFAULT_ALPHA2, total_counts = None, prior=None, filter_counts = True):
    if max(umi_counts.nonzero_values()) == 1: return(umi_counts) # shortcut when there are no duplicates

    if filter_counts:
        # Remove zeros from data, to shorten the vector
        data = list(umi_counts.nonzero_values())
        umi_list = umi_counts.nonzero_keys()
        if uniform:
            n = len(data)
            C_prior = [1./n] * n
        else:
            C_prior = [prior[key] for key in umi_counts.nonzero_keys()]
            n = len(data)
    else:
        data = umi_counts.values()
        umi_list = umi_counts.keys()
        n = len(data)

    N = sum(data)

    # Set priors for the different parameters
    k = len(list(umi_counts.nonzero_values()))
    pi_prior = [k * alpha2, (N - k) * alpha2]
    S_prior = [1.] * n

    # Run Gibbs sampler
    pi_post = MCMC_algorithm_pi.MCMC_algorithm(data, \
                                            n, N, \
                                            S_prior, C_prior, pi_prior, \
                                            nsamp, nthin, nburn, \
                                            True)

    # Distribute counts across tags
    p = computeMedian(pi_post)
    data_dedup = apportion_counts.apportion_counts(data, round(p * sum(data)))

    # Return UmiValues with estimated number of true molecules
    return umi_data.UmiValues(zip(umi_list, data_dedup))

def computeMedian(list):
    list.sort()
    lens = len(list)
    if lens % 2 != 0:
        midl = int(lens / 2)
        res = list[midl]
    else:
        odd = int((lens / 2) -1)
        ev = int((lens / 2))
        res = float(list[odd] + list[ev]) / float(2)
    return res

