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

class Triple(Structure):
    _fields_=[("TriplesCount",c_uint),("TransitionList",Transition),("TriplesList",POINTER(c_double))]

#class GSL_Bundle(Structure):
#    _fields_=[("T",c_uint),("Workspace",Transition),("fdf",POINTER(c_double)),("x",POINTER(c_double)),("f",POINTER(c_double))]






###Code starts here
FitterLib=CDLL("Fitter.so")	#Pull in the fitter library
#Load_ETau_File = FitterLib.Load_ETau_File	#Pull the ET File loader for a test
Load_ETau_File2 = FitterLib.Load_ETau_File2
Load_Base_Catalog_Dictionary = FitterLib.Load_Base_Catalog_Dictionary
Load_Base_Catalog = FitterLib.Load_Base_Catalog
Get_Catalog = FitterLib.Get_Catalog

ETFileName = create_string_buffer(b"etau.dat")
DictionaryFileName = create_string_buffer(b"base_cat_dict.txt")
CatalogFileName = create_string_buffer(b"base_cat.txt")

ETARRAY = POINTER(c_double)()
FILEDELTA = c_double(0)
STATEPOINTS = c_int(0)
STATECOUNT = c_int(0)
ETSTATECOUNT = c_int(0)


Verbose = c_int(0)
ConstantsType = c_double*10


MyLevels = POINTER(Level)()
MyCatalog = POINTER(Transition)()
MyET = ETauStruct()


print ("Verbose: %d" % Verbose.value)

STATECOUNT = Load_Base_Catalog (	CatalogFileName,
									byref(MyCatalog),
									Verbose
)

print (STATECOUNT)

Load_Base_Catalog_Dictionary (	DictionaryFileName, 
								byref(MyLevels),  
								Verbose
)

print("Going ETau2")

Load_ETau_File2 (ETFileName,
 				byref(MyET),
 				byref(ETSTATECOUNT)
)

print("Making Constants")

Constants = ConstantsType (3000.0,2000.0,1000.0)

print("Getting catalog")

Get_Catalog (	MyCatalog, 		
				Constants, 			
				STATECOUNT,	
				Verbose,					
				MyET,
				MyLevels
)

for transition in MyCatalog:
	print(transition.Frequency)

#int Load_ETau_File (char */*FileName*/, double **/*X*/, double */*FileDelta*/, int */*StatePoints*/, int */*StateCount*/);

