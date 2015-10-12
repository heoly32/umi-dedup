from __future__ import division
import numpy as np
import collections

DEFAULT_NSAMP = 1000
DEFAULT_NTHIN = 1
DEFAULT_NBURN = 200

def deduplicate_counts (umi_counts, nsamp=DEFAULT_NSAMP, nthin=DEFAULT_NTHIN, nburn=DEFAULT_NBURN):

    # Remove zeros from data, to shorten the vector
    data = []
    for value in umi_counts.values():
        if value > 0:
            data.append(value)

    n = len(data)
    N = sum(data)
    nits = nburn + nsamp * nthin
    # The 'uniform' algorithm assumes equi-probability for all tags, before amplification
    C_fixed = [1./n] * n
    # Results will be stored in a matrix
    p_post = [[0] * n] * nsamp

    # Set priors for the different parameters
    pi_prior = [1, 1]
    S_prior = [1] * n

    # Starting values
    pi_old = np.random.beta(pi_prior[0], pi_prior[1])
    S_alphas = S_prior
    S_old = np.random.dirichlet(S_alphas)
    Y_old = [0] * n
    Y_new = [0] * n
    p_old = [0] * n
    p_new = [0] * n

    # Start the algorithm
    for i in range(n):
        Y_old[i] = np.random.binomial(data[i], 0.5)
        p_old[i] = (pi_old * C_fixed[i]) / (pi_old * C_fixed[i] + (1 - pi_old) * S_old[i])

    index = 0
    for i in range(nits):

        #1. Update count of true molecules, using binomial distribution
        for j in range(n):
            Y_new[j] = np.random.binomial(data[j], p_old[j])
            # Accept-reject to stay in sample space
            if Y_new[j] == 0 and data[j] >0:
                Y_new[j] = Y_old[j]
            S_alphas[j] = data[j] - Y_new[j] + S_prior[j]

        #2. Update probability of being true molecule, using beta distribution
        Y = sum(Y_new)
        pi_new = np.random.beta(Y + pi_prior[0], N - Y + pi_prior[1])

        #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
        S_new = np.random.dirichlet(S_alphas)

        #4.  Probability of being true molecule given tag, using Bayes' theorem
        p_new = [(pi_new * C_fixed[j]) / (pi_new * C_fixed[j] + (1 - pi_new) * S_new[j]) for j in range(n)]

        # Should we record this draw?
        if i >nburn and i % nthin == 0:
            p_post[index] = p_new
            index += 1

        # Re-initialize vectors before next draw
        pi_old = pi_new
        S_old = S_new
        Y_old = Y_new
        p_old = p_new

    # Compute median for each tag
    p_post_flat = [item for sublist in p_post for item in sublist]

    median_list = [0] * n

    for i in range(n):
        median_list[i] = computeMedian(p_post_flat[i::n])

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

def dirichletvariate(alphas):

    sample = [random.gammavariate(alpha,1) for alpha in alphas]
    sample = [x/sum(sample) for x in sample]

    return sample

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
