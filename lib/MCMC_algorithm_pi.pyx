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
                   unsigned int n, double N, \
                   S_prior, C_prior, pi_prior, \
                   unsigned int nsamp=1000, unsigned int nthin=1, unsigned int nburn=200, \
                   uniform=True):

    # Type iterating variables
    cdef int i, j
    cdef double Y_sum
    cdef long Y_single_draw
    cdef double pi_prior_alpha = pi_prior[0]
    cdef double pi_prior_beta = pi_prior[1]

    # Number of iterations
    cdef unsigned int nits = nburn + nsamp * nthin
    cdef unsigned int index = 0

    # Containers -- dynamic memory allocation necessary
    cdef long *Y_values = <long *>PyMem_Malloc(n * sizeof(long))
    cdef long *data_t = <long *>PyMem_Malloc(n * sizeof(long))
    cdef double *p_values = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *S_values = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *S_alphas = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *C_values = <double *>PyMem_Malloc(n * sizeof(double))
    cdef double *C_alphas = <double *>PyMem_Malloc(n * sizeof(double))

    # Results will be stored in an array
    cdef double *pi_post = <double *>PyMem_Malloc(nsamp * sizeof(double))

    # Type data array
    for i in range(n):
        data_t[i] = data[i]

    # Starting values
    cdef double pi_value = beta(r, pi_prior_alpha, pi_prior_beta)

    if uniform:
        for i in range(n):
            C_values[i] = C_prior[i]
            S_alphas[i] = S_prior[i]
    else:
        for i in range(n):
            C_alphas[i] = C_prior[i]
            S_alphas[i] = S_prior[i]
        # Draw starting values
        dirichlet(r, n, C_alphas, C_values)

    # Draw starting values
    dirichlet(r, n, S_alphas, S_values)

    for i in range(n):
        Y_single_draw = binomial(r, 0.5, data[i])
        if Y_single_draw == 0:
            Y_values[i] = 1
        else:
            Y_values[i] = Y_single_draw
        p_values[i] = (pi_value * C_values[i]) / (pi_value * C_values[i] + (1 - pi_value) * S_values[i])

    # Start the algorithm
    try:
        for i in range(nits):

            #1. Update count of true molecules, using binomial distribution
            for j in range(n):
                Y_single_draw = binomial(r, p_values[j], data[j])
                # Accept-reject to stay in sample space
                if Y_single_draw != 0 and data[j] >0:
                    Y_values[j] = Y_single_draw
                S_alphas[j] = data[j] - Y_values[j] + S_prior[j]
                if not uniform:
                    C_alphas[j] = Y_values[j] + C_prior[j]

            #2. Update probability of being true molecule, using beta distribution
            Y_sum = 0.0
            for j in range(n):
                Y_sum += Y_values[j]
            pi_value = beta(r, Y_sum + pi_prior_alpha, N - Y_sum + pi_prior_beta)

            #3. Update probability of having some tag given that it's a replicate, using dirichlet distribution
            dirichlet(r, n, S_alphas, S_values)
            if not uniform:
                dirichlet(r, n, C_alphas, C_values)

            #4.  Probability of being true molecule given tag, using Bayes' theorem
            for j in range(n):
                p_values[j] = (pi_value * C_values[j]) / (pi_value * C_values[j] + (1 - pi_value) * S_values[j])

            # Should we record this draw?
            if i > nburn and i % nthin == 0:
                for j in range(n):
                    pi_post[index] = pi_value
                index += 1

        return [pi_post[j] for j in range(nsamp)]
    finally:
        PyMem_Free(Y_values)
        PyMem_Free(p_values)
        PyMem_Free(pi_post)
        PyMem_Free(S_values)
        PyMem_Free(C_values)
