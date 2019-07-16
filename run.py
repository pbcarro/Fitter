import numpy as np
from ctypes import *

# https://stackoverflow.com/questions/46700043/using-ctypes-with-gsl-to-pass-array
# https://docs.python.org/3/library/ctypes.html
# http://doc.sagemath.org/html/en/thematic_tutorials/numerical_sage/ctypes_examples.html
# https://cvstuff.wordpress.com/2014/11/27/wraping-c-code-with-python-ctypes-memory-and-pointers/

###Structure definition for python
class Level(Structure):
    _fields_=[("Index",c_uint),("J",c_uint),("Ka",c_uint),("Kc",c_uint)]

class Transition(Structure):
    _fields_=[("Frequency",c_double),("Upper",c_uint),("Lower",c_uint),("Type",c_uint)]

class ETauStruct(Structure):
    _fields_=[("StatePoints",c_uint),("Delta",c_double),("ETVals",POINTER(c_double))]

class Triple(Structure):
    _fields_=[("TriplesCount",c_uint),("TransitionList",Transition),("TriplesList",POINTER(c_double))]

class Opt_Bundle(Structure):
    #Note this is currently not identical to the C code, TransitionGSL is a pointer rather than a finite block
    _fields_=[("ETGSL",ETauStruct),("MyDictionary",POINTER(Level)),("TransitionsGSL",POINTER(Transition))]

#class GSL_Bundle(Structure):
#    _fields_=[("T",c_uint),("Workspace",Transition),("fdf",POINTER(c_double)),("x",POINTER(c_double)),("f",POINTER(c_double))]




###Code starts here
FitterLib=CDLL("Fitter.so")     #Pull in the fitter library
Load_ETau_File = FitterLib.Load_ETau_File("J0_25_dk2.dat")       #Pull the ET File loader for a test

#Load_ETau_File (POINTER(c_char),POINTER(POINTER(c_double)),POINTER(c_double),POINTER(c_int),POINTER(c_int))

#int Load_ETau_File (char */*FileName*/, double **/*X*/, double */*FileDelta*/, int */*StatePoints*/, int */*StateCount*/);

