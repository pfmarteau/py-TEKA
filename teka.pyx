# distutils: language = c++
# distutils: sources = TEKA.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
import numpy as np

cdef extern from "string" namespace "std":
    cdef cppclass string:
        char* c_str()

# c++ interface to cython
cdef extern from "TEKA.h" namespace "teka":

  cdef cppclass TEKA:
    TEKA() except +
    double SqEuclideanDistance(vector[vector[double]] x, vector[vector[double]] y) 
    vector[vector[double]] interpolate(vector[vector[double]])
    double kdtw(vector[vector[double]], vector[vector[double]], double sigma, double epsilon)
    vector[vector[double]] iTEKA(vector[vector[double]], vector[vector[vector[double]]], double sigma, double epsilon)
    vector[vector[vector[double]]] iTEKA_stdev(vector[vector[double]], vector[vector[vector[double]]], double sigma, double epsilon)
    vector[vector[double]] iTEKA2(vector[vector[double]], vector[vector[vector[double]]], double sigma, double epsilon)

# creating a cython wrapper class
cdef class PyTEKA:
    cdef TEKA *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new TEKA()
    def SqEuclideanDistance(self, sv1, sv2):
        return self.thisptr.SqEuclideanDistance(sv1, sv2)
    def interpolate(self, sv):
        return np.array(self.thisptr.interpolate(sv))
    def kdtw(self, sv1, sv2, sigma, epsilon):
        return self.thisptr.kdtw(sv1, sv2, sigma, epsilon)
    def iTEKA(self, sv1, ds, sigma, epsilon):
        return np.array(self.thisptr.iTEKA(sv1, ds, sigma, epsilon))
    def iTEKA_stdev(self, sv1, ds, sigma, epsilon):
        return np.array(self.thisptr.iTEKA_stdev(sv1, ds, sigma, epsilon),dtype=object)

