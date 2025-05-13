from distutils.core import setup
from Cython.Build import cythonize

setup(
    name        = "biquad_filter_optimised",
    version     = "1.1.0",
    description = "Biquad filter library with calculator",
    author      = "Joe Simon",
    url         = "https://github.com/ignaciodsimon/biquad_filter_optimised",
    ext_modules = cythonize('biquad_filter_optimised.pyx'))