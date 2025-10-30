from distutils.core import setup
from Cython.Build import cythonize

setup(
    name        = "biquad_filter_optimised",
    version     = "1.3.0",
    description = "Biquad filter library",
    author      = "Joe Simon",
    url         = "https://github.com/ignaciodsimon/biquad_filter_optimised",
    ext_modules = cythonize('biquad_filter_optimised.pyx'))