import os
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from ctypes import c_uint, c_int, c_double, create_string_buffer, CDLL, POINTER, byref, Structure, c_size_t, CFUNCTYPE, c_void_p, c_char
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
        ("Intensity",c_double),
        ("Map",c_int),
        ("Error",c_double)
        ]

class ETauStruct(Structure):
    _fields_ = [
        ("StatePoints", c_int),
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
        ("TransitionsGSL", POINTER(Transition)),
        ("TransitionCount",c_uint)
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


#Note these special _Alloc versions need to be called for python routines for the first initialization
#This is due to an issue getting C to recognize the ctypes null pointer. I still dont know why it isn't working
#If calling the initialize function after the first call, use the normal version without _Alloc. 
#Failure to do so will result in a minor memory leak, but I can't guarantee that is all it will do.
Initialize_Triples_Fitter =  FitterLib.Initialize_Triples_Fitter_Alloc	
Initialize_SBFIT = FitterLib.Initialize_SBFIT_Alloc	

Load_Exp_File = FitterLib.Load_Exp_File
Find_Triples = FitterLib.Find_Triples
Peak_Find = FitterLib.Peak_Find
Calculate_Intensities = FitterLib.Calculate_Intensities
Calculate_State_Energies = FitterLib.Calculate_State_Energies
Sort_Catalog = FitterLib.Sort_Catalog
print_Transition = FitterLib.print_Transition
Fill_Catalog_Restricted_J = FitterLib.Fill_Catalog_Restricted_J
Fill_Catalog_Restricted_Ka = FitterLib.Fill_Catalog_Restricted_Ka
SBFIT = FitterLib.SBFIT
Timing_Test_Triples = FitterLib.Timing_Test_Triples
Fill_Catalog_Restricted_J2 = FitterLib.Fill_Catalog_Restricted_J2
Get_Catalog2 = FitterLib.Get_Catalog2
Make_Default_Level_List = FitterLib.Make_Default_Level_List
Get_Catalog2_DJ = FitterLib.Get_Catalog2_DJ
Get_Catalog_DJ = FitterLib.Get_Catalog_DJ
Load_DJ_File = FitterLib.Load_DJ_File

####Define a scoring function in python
ScoringFunction = CFUNCTYPE(c_double, POINTER(Transition), c_void_p)
CurrentScoringFunction = ScoringFunction(Score_Fit)

###Create file names for loading
ETFileName = create_string_buffer(b"etau.dat")
DictionaryFileName = create_string_buffer(b"base_cat_dict.txt")
CatalogFileName = create_string_buffer(b"base_cat.txt")
ExperimentalFileName = create_string_buffer(b"../tests/ft2494.txt")
DJFileName = create_string_buffer(b"DJCoeffs.txt")

###Create instances of the C structures
MyLevels = POINTER(Level)()
MyCatalog = POINTER(Transition)()
MyCatalog_JRestricted = POINTER(Transition)()
MyCatalog_KaRestricted = POINTER(Transition)()
MyCatalog_IntensityRestricted = POINTER(Transition)()
MyET = ETauStruct()
MyGSLBundle = GSL_Bundle()
MyOptBundle = Opt_Bundle()
MyTriple = Triple()
TestType = Transition*3
FitTransitions = TestType()
ConstantsType = c_double * 3
ConstantsType_DJ = c_double * 4


FittingFrequencyType = c_double*10
TestWorkspace = POINTER(GSL_multifit_nlinear_workspace)()

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
ChiSqr = c_double(-1.0E-10)	#Initialize the chi squared to a nonsense number
SBFITResults = POINTER(c_double)()


##C integer type initialized to 0
Statecount = c_int(0)
Verbose = c_int(0)	#0 for a non verbose program

#####################################
###Basic load of program files
###These are the foundational set of files needed by essentially all functions in the code
###

CatalogStateCount = Load_Base_Catalog (		CatalogFileName,
											byref(MyCatalog),
											Verbose
										)

CatalogLevelCount = Load_Base_Catalog_Dictionary (	DictionaryFileName, 
													byref(MyLevels),  
													Verbose
												)
Load_ETau_File2 (	ETFileName,
					byref(MyET),
					byref(Statecount),
					Verbose
				)


#####################################
###Flags for what examples to run
CatalogandSort = 0		#Basic example of calculating frequencies and then sorting them
EnergyCalculation = 0	#Calculate the energies of a catalog
FitConstants = 0		#Example of fitting constants to lines, and using custom catalogs for extra speed
TriplesExample = 0	#Example of running a triples fit
GetCatalog_V2 = 0
Distort_DJ = 1



#####################################
###Example catalog calculation and sort
if (CatalogandSort):
	Verbose.value = 1
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
	

#####################################
###Example catalog energy calculation
if (EnergyCalculation):
	Verbose.value = 0
	Constants = ConstantsType (3000.0,2000.0,1000.0)
	#An un initialized catalog will have random energy values, an energy calculation requires a frequency calculation first
	#Consider using Get_Catalog V2 (example below) instead, it's faster and doesn't require separate energy calculations.
	#This is included for legacy code and just general interest.
	Get_Catalog (	MyCatalog, 		
					Constants, 			
					CatalogStateCount,	
					Verbose,					
					MyET,
					MyLevels
	)
	Calculate_State_Energies (MyLevels, MyCatalog, CatalogStateCount, CatalogLevelCount, Verbose)


#####################################
###Example fit of constants
###This also demonstrates how to build/use trimmed catalogs
###In general if you know you dont need transitions (e.g. high J/Ka lines) you can gain considerable speed by dropping them

if (FitConstants):
	RealConstants = ConstantsType (3333.0,2222.0,1111.0)
	GuessConstants = ConstantsType (3000.0,2000.0,1000.0)
	SBFITResults =  ConstantsType(1.0,2.0,3.0)
	Get_Catalog (	MyCatalog, 		
					RealConstants, 			
					CatalogStateCount,	
					Verbose,					
					MyET,
					MyLevels
	)
	
	MyOptBundle.ETGSL = MyET
	MyOptBundle.MyDictionary = MyLevels
	MyOptBundle.TransitionCount = 10
	Initialize_SBFIT (byref(MyGSLBundle),byref(MyOptBundle))
	FittingFrequencies = FittingFrequencyType()
	
	JRestrictedLines = Fill_Catalog_Restricted_J (	MyCatalog,
													byref(MyCatalog_JRestricted),
													RealConstants,
													CatalogStateCount,
													0,	#JMin
													10,	#JMax
													Verbose,	
													MyLevels
												)
	#This takes the J restricted catalog as an input and further reduces its size
	KaRestrictedLines = Fill_Catalog_Restricted_Ka (MyCatalog_JRestricted,
													byref(MyCatalog_KaRestricted),
													RealConstants,
													JRestrictedLines,
													0,	#KaMin
													7,	#KaMax
													Verbose,	
													MyLevels
												)
	print ("Initial Catalog Size: %d. JRestricted Catalog Size: %d. Ka Restricted Catalog Size: %d" % (CatalogStateCount,JRestrictedLines,KaRestrictedLines))
	i = 0
	Count = 0
	while ((i<KaRestrictedLines) and (Count < 10)): 	#if (MyCatalog_KaRestricted[i].Frequency > 6000.0):
		if ((MyCatalog_KaRestricted[i].Frequency > 6000.0) and (MyCatalog_KaRestricted[i].Frequency < 18000.0)):
			FittingFrequencies[Count] = MyCatalog_KaRestricted[i].Frequency
			print ("Adding line at %f to fitting frequencies" % FittingFrequencies[Count])
			MyOptBundle.TransitionsGSL[Count] = MyCatalog_KaRestricted[i]
			Count +=1
		i+=1
	
	
	SBFIT (	GuessConstants,
			byref(ChiSqr),
			byref(MyGSLBundle),
			MyOptBundle,
			FittingFrequencies,
			SBFITResults
			
			)
	
	print ("===Fit Complete===")
	print ("Input Constants: A: %f B: %f C: %f" % (RealConstants[0],RealConstants[1],RealConstants[2]))
	print ("Guess Constants: A: %f B: %f C: %f" % (GuessConstants[0],GuessConstants[1],GuessConstants[2]))
	print ("Fit Constants: A: %f B: %f C: %f" % (SBFITResults[0],SBFITResults[1],SBFITResults[2]))
	print ("Chi Squared: %f" % ChiSqr.value)


#####################################
###Example Triples fit
if (TriplesExample):
	Verbose.value = 1
	ExperimentalPoints = Load_Exp_File (	ExperimentalFileName,
											byref(ExpX),
											byref(ExpY),
											Verbose
										)
	
	MyOptBundle.ETGSL = MyET
	MyOptBundle.MyDictionary = MyLevels
	MyOptBundle.TransitionCount = 3
	
	#Picking three random trasitions for a demo
	MyTriple.TransitionList[0] = MyCatalog[400]
	MyTriple.TransitionList[1] = MyCatalog[475]
	MyTriple.TransitionList[2] = MyCatalog[450]
	
	#Setting the 
	FitTransitions[0] = MyTriple.TransitionList[0]
	FitTransitions[1] = MyTriple.TransitionList[1]
	FitTransitions[2] = MyTriple.TransitionList[2]
	
	
	LineCount = Peak_Find (byref(LineList), YMax,YMin,ExpX,ExpY,ExperimentalPoints,Verbose)
	Find_Triples (byref(MyTriple),LineList,Window,LineCount,Verbose)
	Initialize_Triples_Fitter (byref(MyGSLBundle),byref(MyOptBundle)) 
	MyOptBundle.TransitionsGSL = FitTransitions
	
	Constants = ConstantsType (5002.0,3004.0,1998.0)
	Fit_Triples_Bundle	(			MyTriple,
	 								Constants,
	 								byref(FitResults),
	 								byref(MyCatalog),
	 								CatalogStateCount,
	 								byref(MyGSLBundle),
	 								MyOptBundle,
	 								CurrentScoringFunction,
	 								0
						)
	print ("First fit: A: %f B: %f C: %f" % (FitResults[0],FitResults[1],FitResults[2])) 					
						
	#diditwork = Timing_Test_Triples()								
					
					
#####################################
###Example using Get_Catalog_V2.  
###This version using the same structures in and out to hold results, but uses a different method for calculating the frequencies
###This also results in energies automatically being calculated
###In testing this is typically ~4x faster than the V1 and it is generally recommended that it is used.
###The only caveat here is that the function requires a list of the energy levels that need to be calculated, failure to correctly supply this will cause inaccurate results.
###The list should be given as a 1D int array, with the values in the list being the indices of the levels in the catalog's base dictionary.
###C examples can be found in the Timing_Test function

if (GetCatalog_V2):
	
	print ("Running the Get_Catalog V2 Example:")
	
	Constants = ConstantsType (5002.0,3004.0,1998.0)
	JDictionaryCount = c_int(0)		#Integer set to 0
	LevelList = POINTER(c_int)()	#Uninitialized pointer to int
		
	###If you want to use the full catalog, you can call this. Alternatively you can just make an int array of the same length as the dictionary and fill it with 0,1,2,...N-1 
	#LevelsCount = Make_Default_Level_List (CatalogLevelCount,byref(LevelList))

	JRestrictedLines = Fill_Catalog_Restricted_J2 (	MyCatalog, 						#Input Catalog
													byref(MyCatalog_JRestricted), 	#Output Catalog
													Constants, 						#Constants for the calculation
													CatalogStateCount, 				#Number of transitions in the catalog
													0, 								#J Min
													10, 							#J Max
													0, 								#Verbose Flag
													MyLevels, 						#Catalog Dictionary
													byref(LevelList), 					#List of energy levels in the catalog - see above for more detail
													byref(JDictionaryCount)						#Number of energy levels in LevelList
												);
	

	Get_Catalog2 (			Constants, 			#Input Constants
							MyET, 				#E_tau structure
							MyLevels, 			#Catalog Dictionary
							LevelList, 			#List of levels 
							JDictionaryCount,	#Length of LevelList 
							MyCatalog_JRestricted, 		#Catalog to use
							JRestrictedLines,	#Number of transitions in the catalog
							1					#Verbose
			);

#####################################
###Example using Get_Catalog_V2.  
###
###
### 
###
###
###

if (Distort_DJ):
	
	print ("Running the DJ Example:")
	
	Constants = ConstantsType_DJ (5002.0,3004.0,1998.0,0.001)
	JDictionaryCount = c_int(0)		#Integer set to 0
	LevelList = POINTER(c_int)()	#Uninitialized pointer to int
	DJSlopes = POINTER(c_double)()
	
	Test = Load_DJ_File (DJFileName,byref(DJSlopes),1)
	
		
	###If you want to use the full catalog, you can call this. Alternatively you can just make an int array of the same length as the dictionary and fill it with 0,1,2,...N-1 
	#LevelsCount = Make_Default_Level_List (CatalogLevelCount,byref(LevelList))

	JRestrictedLines = Fill_Catalog_Restricted_J2 (	MyCatalog, 						#Input Catalog
													byref(MyCatalog_JRestricted), 	#Output Catalog
													Constants, 						#Constants for the calculation
													CatalogStateCount, 				#Number of transitions in the catalog
													0, 								#J Min
													10, 							#J Max
													0, 								#Verbose Flag
													MyLevels, 						#Catalog Dictionary
													byref(LevelList), 					#List of energy levels in the catalog - see above for more detail
													byref(JDictionaryCount)						#Number of energy levels in LevelList
												);
	

	Get_Catalog2_DJ (			Constants, 			#Input Constants
								MyET, 				#E_tau structure
								MyLevels, 			#Catalog Dictionary
								LevelList, 			#List of levels 
								JDictionaryCount,	#Length of LevelList 
								MyCatalog_JRestricted, 		#Catalog to use
								JRestrictedLines,	#Number of transitions in the catalog
								0,					#Verbose
								DJSlopes
			)


			