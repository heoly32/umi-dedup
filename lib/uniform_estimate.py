from __future__ import division
import numpy as np
from scipy import stats
import collections

def deduplicate_counts (umi_counts, nsamp=1000, nthin=1, nburn=200):

    n = len(umi_counts)

    # Create array containing data, to simplify computations
    data = np.full((n,), 0)
    for i in range(n):
        data[i] = umi_counts.values()[i]

    N = sum(data)
    nits = nburn + nsamp * nthin
    # The 'uniform' algorithm assumes equi-probability for all tags, before amplification
    C_fixed = np.full((n,), 1./n)
    # Results will be stored in a matrix
    p_post = np.full((nsamp, n), 0)

    # Set priors for the different parameters
    pi_prior = np.array((1, 1))
    S_prior = np.full((n, ), 1)

    # Initialize random variables
    pi_rv = stats.beta(pi_prior[0], pi_prior[1])
    S_rv = stats.dirichlet(S_prior)
    Y_rv = [0] * n
    for i in range(n):
        Y_rv[i] = stats.binom(data[i], 0.5)

    # Starting values
    pi_old = pi_rv.rvs(size = 1)
    S_old = S_rv.rvs(size = 1)
    Y_old = np.full((n,), 0)
    Y_new = np.full((n,), 0)
    for i in range(n):
        Y_old[i] = Y_rv[i].rvs(size = 1)
    p_old = ((pi_old * C_fixed) / (pi_old * C_fixed + (1 - pi_old) * S_old)).flatten()

    index = 0
    for i in range(nits):

        #1. Update count of true molecules, using binomial distribution
        for j in range(n):
            Y_rv[j].args = (data[j], p_old[j])
            Y_new[j] = Y_rv[j].rvs(size = 1)
            # Accept-reject to stay in sample space
            if Y_new[j] == 0 and data[j] >0:
                Y_new[j] = Y_old[j]

        Y = sum(Y_new)

        #2. Update probability of being true molecule, using beta distribution
        pi_rv.args = (Y + pi_prior[0], N - Y + pi_prior[1])
        pi_new = pi_rv.rvs(size = 1)

        #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
        S_rv.args = data - Y_new + S_prior
        S_new = S_rv.rvs(size = 1)

        #4.  Probability of being true molecule given tag, using Bayes' theorem
        p_new = ((pi_new * C_fixed) / (pi_new * C_fixed + (1 - pi_new) * S_new)).flatten()

        # Should we record this draw?
        if i >nburn and i % nthin == 0:
            p_post[index,] = p_new
            index += 1

        # Re-initialize vectors before next draw
        pi_old = pi_new
        S_old = S_new
        Y_old = Y_new
        p_old = p_new

    # Compute median for each tag
    median = np.full((n,), 0)
    for i in range(n):
        median[i] = np.median(p_post[:,i])

    # Return ordered dictionary with estimated number of true molecules
    umi_true = collections.OrderedDict()
    for i in range(n):
        umi_true[umi_counts.keys()[i]] = int(np.ceil(median[i] * data[i]))

    return umi_true
