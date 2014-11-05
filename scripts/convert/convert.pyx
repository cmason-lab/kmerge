# distutils: language=c++ 
# distutils: include_dirs = /home/darryl/include
from libcpp.vector cimport vector 
from libcpp.string cimport string
from libcpp cimport bool

cimport cython.operator as co
from cython.operator cimport dereference
cimport cython.operator as co
cimport cython
cimport cpython
cimport libcpp

import numpy as np 
cimport numpy as np 

cdef extern from "convert.h": 
     #vector[np.uint32_t] load_data(string bin_filename, sum) 
     vector[cython.uint] load_data(string bin_filename, bool sum)

def create_list(string bin_filename, bool sum):
    return load_data(bin_filename, sum)

'''def create_ndarray(string bin_filename, sum): 
    cdef vector[np.uint32_t] vec = load_data(bin_filename, sum)
    cdef np.uint32_t* data = &(vec[0]) 
    cdef np.uint32_t[:] data_view = <np.uint32_t[:vec.size()]> data 
    cdef np.ndarray[np.uint32_t, ndim=1, mode='c'] ar 
    ar = np.asarray(data_view, dtype=np.uint32, order='C') 
    #print('array content in Cython: ' + repr(ar)) 
    return ar 
'''