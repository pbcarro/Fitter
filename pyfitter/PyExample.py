import os
import numpy as np
import pandas as pd
from pathlib import Path
from ctypes import c_uint, c_int, c_double, create_string_buffer, CDLL, POINTER, byref, Structure, c_size_t, CFUNCTYPE, c_void_p, c_char
#from GSL_CTYPES import *


###Level 0 structures
class GSL_block(Structure):
	_fields_=	[	
		("size",c_size_t),
		("data",POINTER(c_double)),
		]
				
class GSL_vector(Structure):
	_fields_=	[	
		("size",	c_size_t),
		("stride",	c_size_t),
		("data",	POINTER(c_double)),
		("block",	POINTER(GSL_block)),
		("owner",	c_int)
		]

class GSL_matrix(Structure):
	_fields_=	[	
		("size1",	c_size_t),
		("size2",	c_size_t),
		("tda",		c_size_t),
		("data",	POINTER(c_double)),
		("block",	POINTER(GSL_block)),
		("owner",	c_int)
		]

###Level 1 structures
###Function pointer types for Level 1 gsl structures
free_type = CFUNCTYPE(	c_void_p,
						c_void_p
					)

f_type = CFUNCTYPE(	c_int,
					POINTER(GSL_vector),
					c_void_p,
					POINTER(GSL_vector)
					)

df_type = CFUNCTYPE(c_int,
					POINTER(GSL_vector),
					c_void_p,
					POINTER(GSL_matrix)
					)
					
fvv_type = CFUNCTYPE(	c_int, 
						POINTER(GSL_vector), 
						POINTER(GSL_vector), 
						POINTER(c_void_p),
						POINTER(GSL_vector)
					)

trs_alloc_type = CFUNCTYPE	(	c_void_p,
							c_void_p,
							c_size_t,
							c_size_t
						)

trs_init_type = CFUNCTYPE(	c_int,
							c_void_p,
							c_void_p
						)

preloop_type = CFUNCTYPE(	c_int,
							c_void_p,
							c_void_p
						)

step_type = CFUNCTYPE	(	c_int,
							c_void_p,
							c_double,POINTER(GSL_vector),
							c_void_p
						)

preduction_type = CFUNCTYPE(c_int,
							c_void_p,
							POINTER(GSL_vector),
							c_double,
							c_void_p
							)

scale_init_type = CFUNCTYPE(c_int,
							POINTER(GSL_matrix),
							POINTER(GSL_vector)
							)

scale_update_type = CFUNCTYPE(	c_int,
								POINTER(GSL_matrix),
								POINTER(GSL_vector)
							)

solver_alloc_type = CFUNCTYPE	(	c_void_p,
							c_void_p,
							c_size_t,
							c_size_t
						)

solver_init_type = CFUNCTYPE(	c_int,
								c_void_p,
								c_void_p
							)

solver_presolve_type = CFUNCTYPE(	c_int,
									c_double,
									c_void_p,
									c_void_p
								)

solver_rcond_type = CFUNCTYPE(	c_int,
								c_double,
								c_void_p
							)

class GSL_fdf(Structure):
	_fields_=	[	
		("f",		f_type),
		("df",		df_type),
		("fvv",		fvv_type),
		("n",		c_size_t),
		("p",		c_size_t),
		("params",	c_void_p),
		("nevalf",	c_size_t),
		("nevaldf",	c_size_t),
		("nevalfvv",c_size_t)
 		]

class GSL_nlinear_trs(Structure):
	_fields_=	[	
		("name",		POINTER(c_char)),
		("alloc",		trs_alloc_type),
		("init",		trs_init_type),
		("preloop",		preloop_type),
		("step",		step_type),
		("preduction",	preduction_type),
		("free",		free_type)
 		]

class GSL_nlinear_scale(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("init",	scale_init_type),
		("update",	scale_update_type)
 		]

class GSL_nlinear_solver(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("alloc",	solver_alloc_type),
		("init",	solver_init_type),
		("presolve",solver_presolve_type),
		("rcond",	solver_rcond_type),
		("free",	free_type),
 		]

###Level 2 structures
class GSL_multifit_nlinear_parameters(Structure):
	_fields_=	[	
		("trs",			POINTER(GSL_nlinear_trs)),
		("scale",		POINTER(GSL_nlinear_scale)),
		("solver",		POINTER(GSL_nlinear_solver)),
		("fdtype",		c_int),
		("factor_up",	c_double),
		("factor_down",	c_double),
		("avmax",		c_double),
		("h_df",		c_double),
		("h_fvv",		c_double),
 		]

###Function pointer types for Level 2 gsl structures
nlinear_alloc_type = CFUNCTYPE	(	c_void_p,
									POINTER(GSL_multifit_nlinear_parameters),
									c_size_t,
									c_size_t
								)

nlinear_init_type = CFUNCTYPE	(	c_int,
									c_void_p,
									POINTER(GSL_vector),
									POINTER(GSL_fdf),
									POINTER(GSL_vector),
									POINTER(GSL_vector),
									POINTER(GSL_matrix),
									POINTER(GSL_vector)
								)

nlinear_iterate_type = CFUNCTYPE(	c_int,
									c_void_p,
									POINTER(GSL_vector),
									POINTER(GSL_fdf),
									POINTER(GSL_vector),
									POINTER(GSL_vector),
									POINTER(GSL_matrix),
									POINTER(GSL_vector),
									POINTER(GSL_vector)
								)

nlinear_type_rcond_type = CFUNCTYPE(c_int,
									c_double,
									c_void_p
									)

nlinear_avratio_type = CFUNCTYPE(	c_double,
									c_void_p
								)

class GSL_multifit_nlinear_type(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("alloc",	nlinear_alloc_type),
		("init",	nlinear_init_type),
		("iterate",	nlinear_iterate_type),
		("rcond",	nlinear_type_rcond_type),
		("avratio",	nlinear_avratio_type),
		("free",	free_type)
 		]


class GSL_multifit_nlinear_workspace(Structure):
	_fields_=	[	
		("type",			POINTER(GSL_multifit_nlinear_type)),
		("fdf",				POINTER(GSL_fdf)),
		("x",				POINTER(GSL_vector)),
		("f",				POINTER(GSL_vector)),
		("dx",				POINTER(GSL_vector)),
		("g",				POINTER(GSL_vector)),
		("J",				POINTER(GSL_matrix)),
		("sqrt_wts_work",	POINTER(GSL_vector)),
		("sqrt_wts",		POINTER(GSL_vector)),
		("niter",			c_size_t),
		("params",			GSL_multifit_nlinear_parameters),
		("state",			c_void_p)
 		]	

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
        ("fdf", 		GSL_fdf),
        ("f", 			POINTER(GSL_vector)),
        ("fdf_params", 	GSL_multifit_nlinear_parameters),
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



LineCount = Peak_Find (byref(LineList), YMax,YMin,ExpX,ExpY,ExperimentalPoints,Verbose)
MyTriple.TransitionList[0] = MyCatalog[400]
MyTriple.TransitionList[1] = MyCatalog[500]
MyTriple.TransitionList[2] = MyCatalog[450]
##Currently redundantly passing the transitions to fit twice, need to address this in the code
FitTransitions[0] = MyTriple.TransitionList[0]
FitTransitions[1] = MyTriple.TransitionList[1]
FitTransitions[2] = MyTriple.TransitionList[2]
Find_Triples (byref(MyTriple),LineList,Window,LineCount)


###Correct till here


MyOptBundle.ETGSL = MyET
MyOptBundle.MyDictionary = MyLevels
MyOptBundle.TransitionsGSL = FitTransitions

##Checking objects are set correctly going into the triples fit
# print ("%f %f %f" %(Constants[0], Constants[1], Constants[2]))	#Constants correctly passed to C
#print ("%f %d" % (MyOptBundle.ETGSL.Delta,MyOptBundle.ETGSL.StatePoints))
# print ("%d %d %d" % (MyOptBundle.MyDictionary[400].J,MyOptBundle.MyDictionary[500].Ka,MyOptBundle.MyDictionary[450].Kc))
# print ("%d" % (MyOptBundle.TransitionsGSL[2].Upper))
# 
# print ("%d %d %d" % (MyTriple.TriplesCount[0],MyTriple.TriplesCount[1],MyTriple.TriplesCount[2]))
# print ("%f %f %f" % (MyTriple.TriplesList[10],MyTriple.TriplesList[100],MyTriple.TriplesList[67]))

# print (CatalogStateCount) #CatalogStateCount


Constants = ConstantsType (3000.0,2000.0,1000.0)
Fit_Triples_Bundle 	(	MyTriple,
 						Constants,
 						byref(FitResults),
 						byref(MyCatalog),
 						CatalogStateCount,
 						byref(MyGSLBundle),
 						MyOptBundle,
 						#CurrentScoringFunction,
# 						#Scoring_Parameters
					)
					
					
					
					