from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
import numpy as np

def MCMC_algorithm(data, unsigned int n, unsigned int N, S_prior, C_prior, pi_prior, unsigned int nsamp=1000, unsigned int nthin=1, unsigned int nburn=200, uniform=True):

    # Type iterating variables
    cdef int i, j

    # Number of iterations
    cdef unsigned int nits = nburn + nsamp * nthin
    cdef unsigned int index = 0

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

    # Containers -- dynamic memory allocation necessary
    cdef double *Y_old = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *Y_new = <double *>PyMem_Malloc(n * sizeof(double))
    cdef float *p_old = <float *>PyMem_Malloc(n * sizeof(float))
    cdef float *p_new = <float *>PyMem_Malloc(n * sizeof(float))

    # Results will be stored in an array
    cdef float *p_post = <float *>PyMem_Malloc(n * nsamp * sizeof(float))

    for i in range(n):
            Y_old[i] = np.random.binomial(data[i], 0.5)
            p_old[i] = (pi_old * C_old[i]) / (pi_old * C_old[i] + (1 - pi_old) * S_old[i])

    # Start the algorithm
    try:
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
            Y = 0
            for j in range(n):
                Y += Y_new[j]
            pi_new = np.random.beta(Y + pi_prior[0], N - Y + pi_prior[1])

            #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
            S_new = np.random.dirichlet(S_alphas)
            if not uniform:
                C_new = np.random.dirichlet(C_alphas)

            #4.  Probability of being true molecule given tag, using Bayes' theorem
            for j in range(n):
                p_new[j] = (pi_new * C_new[j]) / (pi_new * C_new[j] + (1 - pi_new) * S_new[j])

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

        return [p_post[j] for j in range(n * nsamp)]
    finally:
        PyMem_Free(Y_old)
        # Y_new points to Y_old at the end of the algorithm
        # PyMem_Free(Y_new)
        PyMem_Free(p_old)
        # p_new points to p_old at the end of the algorithm
        # PyMem_Free(p_new)
        PyMem_Free(p_post)
