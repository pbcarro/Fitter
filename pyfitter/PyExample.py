import os
import numpy as np
import pandas as pd
from pathlib import Path
from ctypes import c_uint, c_int, c_double, create_string_buffer, CDLL, POINTER, byref, Structure
from GSL_CTYPES import *

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
        ("TriplesCount", c_uint*3),
        ("TransitionList", Transition*3),
        ("TriplesList", POINTER(c_double))
        ]

class Opt_Bundle(Structure):
    #Note this is currently not identical to the C code, TransitionGSL is a pointer rather than a finite block
    _fields_ = [
        ("ETGSL", ETauStruct),
        ("MyDictionary", POINTER(Level)),
        ("TransitionsGSL", POINTER(Transition))
        ]

class GSL_Bundle(Structure):
    _fields_ = [
        ("T", 			POINTER(GSL_multifit_nlinear_type)),
        ("Workspace", 	POINTER(GSL_multifit_nlinear_workspace)),
        ("fdf", 		POINTER(GSL_fdf)),
        ("f", 			POINTER(GSL_vector)),
        ("fdf_params", 	POINTER(GSL_multifit_nlinear_parameters)),
        ]

def Score_Fit (Catalog,Scoring_Parameters):
	###Score stuff
	return 1.0


##Import the library and pull out functions to call
FitterLib=CDLL("Fitter.so")	#Pull in the fitter library
Load_ETau_File2 = FitterLib.Load_ETau_File2
Load_Base_Catalog_Dictionary = FitterLib.Load_Base_Catalog_Dictionary
Load_Base_Catalog = FitterLib.Load_Base_Catalog
Get_Catalog = FitterLib.Get_Catalog
Fit_Triples_Bundle =  FitterLib.Fit_Triples_Bundle
Initialize_Triples_Fitter =  FitterLib.Initialize_Triples_Fitter
Load_Exp_File = FitterLib.Load_Exp_File
Find_Triples = FitterLib.Find_Triples
Peak_Find = FitterLib.Peak_Find
Calculate_Intensities = FitterLib.Calculate_Intensities
Calculate_State_Energies = FitterLib.Calculate_State_Energies
Sort_Catalog = FitterLib.Sort_Catalog

####Define a scoring function in python
ScoringFunction = CFUNCTYPE(c_double, POINTER(Transition), c_void_p)
CurrentScoringFunction = ScoringFunction(Score_Fit)

###Create file names for loading
ETFileName = create_string_buffer(b"etau.dat")
DictionaryFileName = create_string_buffer(b"base_cat_dict.txt")
CatalogFileName = create_string_buffer(b"base_cat.txt")
ExperimentalFileName = create_string_buffer(b"../tests/ft2494.txt")

###Create instances of the C structures
MyLevels = POINTER(Level)()
MyCatalog = POINTER(Transition)()
MyET = ETauStruct()
MyGSLBundle = GSL_Bundle()
MyOptBundle = Opt_Bundle()
MyTriple = Triple()
TestType = Transition*3
FitTransitions = TestType()
ConstantsType = c_double * 3

###Create some variables
ExpX = POINTER(c_double)()
ExpY = POINTER(c_double)()
ExperimentalPoints = c_int(0)
LineList = POINTER(c_double)()
YMax = c_double(1.0E+3)	#Min and max values of the experimental data that we want to use for peak finding
YMin = c_double(0.003)
Window = c_double(100.0)	#Distance from each line center in MHz to use for taking potential lines for the triples
LineCount = c_int(0)
FitResults = POINTER(c_double)()
Scoring_Parameters = c_void_p(0)

##C integer type initialized to 0
Statecount = c_int(0)
Verbose = c_int(0)	#0 for a non verbose program

###Basic load of program files
CatalogStateCount = Load_Base_Catalog (		CatalogFileName,
											byref(MyCatalog),
											Verbose
										)

Load_Base_Catalog_Dictionary (	DictionaryFileName, 
								byref(MyLevels),  
								Verbose
							)
Load_ETau_File2 (	ETFileName,
					byref(MyET),
					byref(Statecount),
					Verbose
				)


##Example catalog calculation and sort
Constants = ConstantsType (3000.0,2000.0,1000.0)
Get_Catalog (	MyCatalog, 		
				Constants, 			
				CatalogStateCount,	
				Verbose,					
				MyET,
				MyLevels
)


Sort_Catalog (	MyCatalog,
				CatalogStateCount,
				c_int(2),		#Sort by frequency
				c_int(0)		#Use insertion sort to do it
			)

##Example catalog energy calculation
Calculate_State_Energies (MyLevels, MyCatalog, CatalogStateCount)

###Example Triples fit
Verbose.value = 0
ExperimentalPoints = Load_Exp_File (	ExperimentalFileName,
										byref(ExpX),
										byref(ExpY),
										Verbose
									)

Initialize_Triples_Fitter (byref(MyGSLBundle)) 
print (MyET.Delta)

LineCount = Peak_Find (byref(LineList), YMax,YMin,ExpX,ExpY,ExperimentalPoints,Verbose)
MyTriple.TransitionList[0] = MyCatalog[400]
MyTriple.TransitionList[1] = MyCatalog[500]
MyTriple.TransitionList[2] = MyCatalog[450]
Find_Triples (byref(MyTriple),LineList,Window,LineCount)
###Correct till here


MyOptBundle.ETGSL = MyET
MyOptBundle.MyDictionary = MyLevels
MyOptBundle.TransitionsGSL = FitTransitions
Constants = ConstantsType (3000.0,2000.0,1000.0)
Fit_Triples_Bundle 	(	MyTriple,
 						Constants,
 						FitResults,
 						byref(MyCatalog),
 						CatalogStateCount,
 						MyGSLBundle,
 						MyOptBundle,
 						#CurrentScoringFunction,
# 						#Scoring_Parameters
					)
					
					
					
					