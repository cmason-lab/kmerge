from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(ext_modules = cythonize(
    "convert.pyx",                 # our Cython source
    sources=["/home/darryl/include/dlib/all/source.cpp"],  # additional source file(s)
    include_path=[np.get_include()],
))
