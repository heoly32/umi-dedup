from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

include_gsl_dir = "/usr/local/include/"
lib_gsl_dir = "/usr/local/lib/"

ext_modules = [Extension("lib.MCMC_algorithm", ["lib/MCMC_algorithm.pyx"],
                        include_dirs=[include_gsl_dir],
                        library_dirs=[lib_gsl_dir],
                        libraries=["gsl", "gslcblas"])]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules = ext_modules
)
