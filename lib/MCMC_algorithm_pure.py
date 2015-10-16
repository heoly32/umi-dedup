import numpy as np

def MCMC_algorithm(data, n, N, S_prior, C_prior, pi_prior, nsamp=1000, nthin=1, nburn=200, uniform=True):

    # Number of iterations
    nits = nburn + nsamp * nthin

    # Starting values
    pi_old = np.random.beta(pi_prior[0], pi_prior[1])
    S_alphas = S_prior
    S_old = np.random.dirichlet(S_alphas)
    if uniform:
        C_old = C_prior
        C_new = C_prior
    else:
        C_alphas = C_prior
        C_old = np.random.dirichlet(C_alphas)

    Y_old = [0] * n
    Y_new = [0] * n
    p_old = [0] * n
    p_new = [0] * n

    # Results will be stored in a list
    p_post = [0] * (n * nsamp)

    # Start the algorithm
    for i in range(n):
        Y_old[i] = np.random.binomial(data[i], 0.5)
        p_old[i] = (pi_old * C_old[i]) / (pi_old * C_old[i] + (1 - pi_old) * S_old[i])

    index = 0
    for i in range(nits):

        #1. Update count of true molecules, using binomial distribution
        for j in range(n):
            Y_new[j] = np.random.binomial(data[j], p_old[j])
            # Accept-reject to stay in sample space
            if Y_new[j] == 0 and data[j] >0:
                Y_new[j] = Y_old[j]
            S_alphas[j] = data[j] - Y_new[j] + S_prior[j]
            if not uniform:
                C_alphas[j] = Y_new[j] + C_prior[j]

        #2. Update probability of being true molecule, using beta distribution
        Y = sum(Y_new)
        pi_new = np.random.beta(Y + pi_prior[0], N - Y + pi_prior[1])

        #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
        S_new = np.random.dirichlet(S_alphas)
        if not uniform:
            C_new = np.random.dirichlet(C_alphas)

        #4.  Probability of being true molecule given tag, using Bayes' theorem
        p_new = [(pi_new * C_new[j]) / (pi_new * C_new[j] + (1 - pi_new) * S_new[j]) for j in range(n)]

        # Should we record this draw?
        if i > nburn and i % nthin == 0:
            for j in range(n):
                p_post[index*n + j] = p_new[j]
            index += 1

        # Re-initialize vectors before next draw
        pi_old = pi_new
        S_old = S_new
        C_old = C_new
        Y_old = Y_new
        p_old = p_new

    return p_post
