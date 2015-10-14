from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
import numpy as np

# declaring external GSL functions to be used
# replacting numpy RNGs by GSL's

cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng_type:
        pass
    ctypedef struct gsl_rng:
        pass
    gsl_rng_type *gsl_rng_mt19937
    gsl_rng *gsl_rng_alloc(gsl_rng_type * T)

cdef extern from "gsl/gsl_randist.h":
   double beta "gsl_ran_beta"(const gsl_rng * r, double, double)
   void dirichlet "gsl_ran_dirichlet"(const gsl_rng * r, size_t K, const double alpha[], double theta[])
   unsigned int binomial "gsl_ran_binomial"(const gsl_rng * r, double p, unsigned int n)

cdef gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937)

def MCMC_algorithm(data, \
                   unsigned int n, unsigned int N, \
                   S_prior, C_prior, pi_prior, \
                   unsigned int nsamp=1000, unsigned int nthin=1, unsigned int nburn=200, \
                   uniform=True):

    # Type iterating variables
    cdef int i, j

    # Number of iterations
    cdef unsigned int nits = nburn + nsamp * nthin
    cdef unsigned int index = 0

    # Containers -- dynamic memory allocation necessary
    cdef long *Y_old = <long *>PyMem_Malloc(n * sizeof(long))
    cdef long *Y_new = <long *>PyMem_Malloc(n * sizeof(long))
    cdef double *p_old = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *p_new = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *S_old = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *S_new = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *S_alphas = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *C_old = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *C_new = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *C_alphas = <double *>PyMem_Malloc(n * sizeof(double))


    # Results will be stored in an array
    cdef double *p_post = <double *>PyMem_Malloc(n * nsamp * sizeof(double))

    # Starting values
    cdef double pi_old = beta(r, pi_prior[0], pi_prior[1])
    for i in range(n):
        S_alphas[i] = S_prior[i]
    dirichlet(r, n, S_alphas, S_old)
    if uniform:
        for i in range(n):
            C_old[i] = C_prior[i]
            C_new[i] = C_prior[i]
    else:
        for i in range(n):
            C_alphas[i] = C_prior[i]
        dirichlet(r, n, C_alphas, C_old)

    for i in range(n):
            Y_old[i] = binomial(r, 0.5, data[i])
            if Y_old[i] == 0:
                Y_old[i] = 1
            p_old[i] = (pi_old * C_old[i]) / (pi_old * C_old[i] + (1 - pi_old) * S_old[i])

    # Start the algorithm
    try:
        for i in range(nits):

            #1. Update count of true molecules, using binomial distribution
            for j in range(n):
                Y_new[j] = binomial(r, p_old[j], data[j])
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
            pi_new = beta(r, Y + pi_prior[0], N - Y + pi_prior[1])

            #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
            dirichlet(r, n, S_alphas, S_new)
            if not uniform:
                dirichlet(r, n, C_alphas, S_new)

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
        PyMem_Free(S_old)
        PyMem_Free(C_old)
