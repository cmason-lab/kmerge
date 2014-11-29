# distutils: language=c++ 
# distutils: include_dirs = /home/darryl/include
# distutils: sources = /home/darryl/include/dlib/all/source.cpp
# distutils: extra_compile_args = -DDLIB_NO_GUI_SUPPORT
 
from libcpp.vector cimport vector 
from libcpp.string cimport string
from libcpp cimport bool

cimport cython

cdef extern from "loader.h": 
     vector[cython.uint] load_data(string bin_filename, bool sum)

def create_list(string bin_filename, bool sum):
    return load_data(bin_filename, sum)