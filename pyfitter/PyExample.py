import numpy as np
import pandas as pd
from pathlib import Path
from ctypes import c_uint, c_int, c_double, create_string_buffer, CDLL, POINTER, byref, Structure
import GSL_CTYPES

###Structure definition for python
class Level(Structure):
    _fields_ = [
        ("Index", c_uint),
        ("J", c_uint),
        ("Ka", c_uint),
        ("Kc", c_uint),
        ("Energy",c_double)
        ]

class Transition(Structure):
    _fields_ = [
        ("Frequency", c_double),
        ("Upper", c_uint),
        ("Lower", c_uint),
        ("Type", c_uint),
        ("Intensity",c_double)
        ]

class ETauStruct(Structure):
    _fields_ = [
        ("StatePoints", c_uint),
        ("Delta", c_double),
        ("ETVals", POINTER(c_double))
        ]

class Triple(Structure):
    _fields_ = [
        ("TriplesCount", c_uint),
        ("TransitionList", Transition),
        ("TriplesList", POINTER(c_double))
        ]

class Opt_Bundle(Structure):
    #Note this is currently not identical to the C code, TransitionGSL is a pointer rather than a finite block
    _fields_ = [
        ("ETGSL", ETauStruct),
        ("MyDictionary", POINTER(Level)),
        ("TransitionsGSL", POINTER(Transition))
        ]

def Score_Fit (Catalog,Scoring Parameters):
	###Score stuff
	return 1.0


##Import the library and pull out functions to call
FitterLib=CDLL("Fitter.so")	#Pull in the fitter library
Load_ETau_File2 = FitterLib.Load_ETau_File2
Load_Base_Catalog_Dictionary = FitterLib.Load_Base_Catalog_Dictionary
Load_Base_Catalog = FitterLib.Load_Base_Catalog
Get_Catalog = FitterLib.Get_Catalog
Fit_Triples_Bundle =  FitterLib.Fit_Triples_Bundle



#(struct Triple TransitionstoFit, double *Guess, double **FitResults, struct Transition **MyFittingCatalog, int CatalogLines, double *ExperimentalLines, int ExperimentalLineCount, struct GSL_Bundle *FitBundle, struct Opt_Bundle *MyOpt_Bundle, ScoreFunction TriplesScoreFunction, void *ScoringParameters)

####Define a scoring function in python
CMPFUNC = CFUNCTYPE(c_double, POINTER(Transition), c_void_p)

###Create file names for loading
ETFileName = create_string_buffer(b"J0_25_dk2.dat")
DictionaryFileName = create_string_buffer(b"Base Catalog Dictionary1.txt")
CatalogFileName = create_string_buffer(b"Base Catalog.txt")

###Create instances of the C structures
MyLevels = POINTER(Level)()
MyCatalog = POINTER(Transition)()
MyET = POINTER(ETauStruct)()

##C integer type initialized to 0
Statecount = c_int(0)
Verbose = c_int(0)	#0 for a non verbose program

Statecount = Load_Base_Catalog (	CatalogFileName,
									byref(MyCatalog),
									Verbose
)

print (Statecount)
Load_Base_Catalog_Dictionary (	DictionaryFileName, 
								byref(MyLevels),  
								Verbose
)
Load_ETau_File2 (ETFileName,
				byref(MyET),
				byref(ETSTATECOUNT)
)
Constants = ConstantsType (3000.0,2000.0,1000.0)
Get_Catalog (	MyCatalog, 		
				Constants, 			
				Statecount,	
				Verbose,					
				MyET,
				MyLevels
)
