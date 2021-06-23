#ifndef __FITTER_H__
#define __FITTER_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_linalg.h>

#define MAXLINESIZE 50000000	// A hard limit on the load buffer size, can cause issues on low RAM systems	

//=============Structures==============
struct Level
{
	unsigned int Index;
	unsigned int J;
	unsigned int Ka;
	unsigned int Kc;
	double Energy;
};

struct Transition 
{
	double Frequency;
	unsigned int Upper;
	unsigned int Lower;
	unsigned int Type;
	double Intensity;
	int Map;
	double Error;
}; 

struct ETauStruct 
{
	int StatePoints;
	double Delta;
	double *ETVals;
};

struct Triple 
{
	unsigned int TriplesCount[3];
	struct Transition TransitionList[3];
	double *TriplesList;
};

struct Opt_Bundle 
{
	//New struct to hold all the parameters we pass to the GSL optimizer function
	struct ETauStruct ETGSL;
	struct Level *MyDictionary;
	struct Transition *TransitionsGSL;
	unsigned int TransitionCount;
};

struct GSL_Bundle 
{
	const gsl_multifit_nlinear_type *T;
	gsl_multifit_nlinear_workspace *Workspace;
	gsl_multifit_nlinear_fdf fdf;
	gsl_vector *f;
	gsl_multifit_nlinear_parameters fdf_params;
};

typedef double (*ScoreFunction)(struct Transition *, void *);	//Generic function pointer for the scoring function used in the triples fitter

//=============Function Prototypes==============

//Program setup functions
int Initialize_Stuff (double ** /*ETArray*/, int * /*CatTransitions*/, int * /*DictTransitions*/, double * /*FileDelta*/, int * /*StatePoints*/, struct Level ** /*DictionaryIn*/, struct Transition ** /*CatalogIn*/);
int Load_ETau_File (char * /*FileName*/, double ** /*X*/, double * /*FileDelta*/, int * /*StatePoints*/, int * /*StateCount*/);
int Load_ETau_File2 (char * /*FileName*/, struct ETauStruct * /*StructToLoad*/, int * /*StateCount*/, int /*Verbose*/);
int Load_Base_Catalog (char * /*FileName*/, struct Transition ** /*BaseCatalog*/,  int /*Verbose*/);
int Load_Base_Catalog_Dictionary (char * /*FileName*/, struct Level ** /*DictIn*/,  int /*Verbose*/);
int Load_Exp_File  (char * /*FileName*/, double ** /*X*/, double ** /*Y*/, int /*Verbose*/);
int Load_Str_File (char * /*FileName*/, double *** /*Data*/, int /*Verbose*/);

//Frequency predicting functions
double Get_Kappa (double /*A*/, double /*B*/, double /*C*/);  
double Get_Delta (double /*A*/, double /*B*/, double /*C*/);
int Get_J (int /*TransitionIndex*/, struct Level */*MyDictionary*/);
double Partition_Function (double */*Constants*/, double /*Temperature*/);
double E_tau (int /*TransitionIndex*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
double Rigid_Rotor (double /*A*/, double /*C*/, int /*J*/, int /*Index*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
double Rigid_Rotor_Error (double /*A*/, double /*C*/, int /*J*/, int /*Index*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/, double * /*Hkappa*/);
double Get_Frequency (int /*J_Up*/, int /*J_Low*/, int /*IndexUp*/, int /*IndexLow*/, double */*Constants*/, struct ETauStruct /*ETStruct*/);
double Get_Frequency_Error (int /*J_Up*/, int /*J_Low*/, int /*IndexUp*/, int /*IndexLow*/, double * /*Constants*/, struct ETauStruct /*ETStruct*/, double * /*TransitionError*/, double * /*ConstantsError*/);
int Get_Catalog (struct Transition * restrict/*CatalogtoFill*/, double * restrict/*Constants*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level */*MyDictionary*/);
int Get_Catalog_Error (struct Transition *restrict /*CatalogtoFill*/, double *restrict /*Constants*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level * /*MyDictionary*/, double * /*ConstantsError*/);
double Get_Str (double /*Kappa*/, int /*Transition*/, int /*PointsPerSate*/, double ** /*StrData*/);

//General catalog functions
void print_Transition (struct Transition /*TransitionToPrint*/, struct Level */*MyDictionary*/);
void Sort_Catalog (struct Transition */*Catalog*/, int /*TransitionCount*/, int /*SortMethod*/, int /*SortType*/);
int Catalog_Comparator_Frequency (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Intensity (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Upper (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Lower (const void */*a*/, const void */*b*/);
void insertionSort(struct Transition */*CatalogtoSort*/, int /*TransitionCount*/);
int Fill_Catalog_Restricted_Frequency (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*FrequencyLow*/, double /*FrequencyHigh*/, int /*Verbose*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_J (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*JMin*/, int /*JMax*/, int /*Verbose*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_Intensity (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*Percentile*/, int /*Verbose*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_Intensity_Count (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*CatLinesOut*/, int /*Verbose*/, struct Level */*MyDictionary*/);
void Calculate_State_Energies (struct Level */*MyCatalog*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/, int /*DictionaryLevels*/, int /*Verbose*/);
void Calculate_Intensities (struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/, struct Level */*MyDictionary*/, double /*T*/, double */*Dipoles*/, int /*Verbose*/);
void Calculate_Intensities_Sij (struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, struct Level * /*MyDictionary*/, double /*T*/, double * /*Dipoles*/, double ** /*StrData*/, double /*Kappa*/, double * /*Constants*/, int /*Verbose*/);
int Standardize_Catalog (struct Transition * /*CatalogtoStandardize*/, int /*CatLines*/, int /*Verbose*/);

//Fitting functions
int Peak_Find (double **/*LineList*/, double /*Max*/, double /*Min*/, double */*X*/, double */*Y*/, int /*ArraySize*/, int /*Verbose*/);
int Find_Triples (struct Triple */*TripletoFit*/, double */*LineFrequencies*/, double /*Window*/, int /*LineCount*/, int /*Verbose*/);
double Fit_Triples (double */*Guess*/, int /*Verbose*/, struct Opt_Bundle /*GSLOptBundle*/);
void Initialize_Triples_Fitter (struct GSL_Bundle * /*FitBundle*/, struct Opt_Bundle * /*MyOpt_Bundle*/);
void Initialize_Triples_Fitter_Alloc (struct GSL_Bundle * /*FitBundle*/, struct Opt_Bundle * /*MyOpt_Bundle*/);
int OptFunc_gsl (const gsl_vector */*x*/, void */*params*/, gsl_vector */*f*/);
int Fit_Triples_Bundle (struct Triple /*TransitionstoFit*/, double */*Guess*/, double **/*FitResults*/, struct Transition **/*Catalog*/, int /*CatalogLines*/, struct GSL_Bundle */*FitBundle*/, struct Opt_Bundle /*MyOpt_Bundle*/, ScoreFunction /*TriplesScoreFunction*/, void */*ScoringParameters*/);
void callback (const size_t /*iter*/, void */*params*/, const gsl_multifit_nlinear_workspace */*w*/);
int Initialize_SBFIT (struct GSL_Bundle * /*FitBundle*/, struct Opt_Bundle * /*MyOpt_Bundle*/);
int Initialize_SBFIT_Alloc (struct GSL_Bundle * /*FitBundle*/, struct Opt_Bundle * /*MyOpt_Bundle*/);
int SBFIT (double */*Guess*/, double */*ChiSq*/, struct GSL_Bundle */*FitBundle*/, struct Opt_Bundle /*MyOpt_Bundle*/, double */*LineFrequencies*/, double [3]/*FinalConstants*/);
int Get_SBFIT_Error (struct GSL_Bundle * /*FitBundle*/, double [3] /*Errors*/);
int SBFIT_OptFunc_gsl (const gsl_vector * /*x*/, void * /*params*/, gsl_vector * /*f*/);

//New Functions
int Search_DR_Hits (int /*DRPairs*/, double /*ConstStart*/, double /*ConstStop*/, double /*Step*/, double */*DRFrequency*/, double /*Tolerance*/, int /*ExtraLineCount*/, double */*ExtraLines*/, int **/*DRLinks*/, int /*LinkCount*/, struct Transition */*CatalogtoFill*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level * /*MyDictionary*/, char * /*FileName*/);
int Match_Levels (int /*Match1*/, int /*Match2*/, struct Transition * /*MatchCatalog*/);

//Test Functions
int Timing_Test(void);
int Accuracy_Test(void);
void Test_SBFIT (void);
void Test_Triples (char *, struct Transition *, struct Level *, int /*CatalogLines*/, struct ETauStruct /*FittingETStruct*/);
double DummyFunction (struct Transition */*MyCatalog*/, void */*Data*/);
int Timing_Test_Triples (void);

//Functions

int Initialize_Stuff (double **ETArray, int *CatTransitions, int *DictLevels, double *FileDelta, int *StatePoints, struct Level **DictionaryIn, struct Transition **CatalogIn) 
{
//Ease of use function for loading stuff all at once, outdated and will be removed at some point

int TempPoints;
	TempPoints = 0;
	srand(time(NULL));								//Initializing the time calls elsewhere	
	*DictLevels = Load_Base_Catalog_Dictionary (	"base_cat_dict.txt", 
													DictionaryIn,//&BaseDict, 
													0 //Verbose
	);	
	printf ("Dictionary Loaded\n");
	*CatTransitions = Load_Base_Catalog (	"base_cat.txt",
											CatalogIn, 		//List of transitions in the lower state
											0				//Verbose
	);	
	printf ("Catalog Loaded\n");

	//Load all of our ET values into the array and check that it worked, using static external files cause I'm lazy
	if (Load_ETau_File ("etau.dat", ETArray,FileDelta,StatePoints,&TempPoints)) printf ("Tables Loaded\n");
	else goto Error;
	if (*DictLevels != TempPoints) {
		printf ("Warning: mismatch between number of dictionary states(%d) and number of states in the eigenvalue file (%d)\n",*DictLevels,TempPoints);
	}
	return 1;
Error:
	printf ("Error initializing program\n");
	return 0;
}

int Load_ETau_File (char *FileName, double **X, double *FileDelta, int *StatePoints, int *StateCount) 
{
//This is a deprecated version of the function. Please use Load_ETau_File2 instead unless you have a very good reason to use this. It is here mostly for compatibility/historical reasons. 

//This loads the E_Tau file in by reading a single value from the file until there are no more values found in the file
//This unwraps the 2D file into a single long array for contiguousness, files should have a state's ET values in a single ROW, not column
//Figuring out where each state starts and stops is done elsewhere, so the global variables need to be set properly to deal with this
//ET is preallocated based on the global variables, there are no checks to see if the file matches until later
int i,StateLimit;
FILE *FileHandle;
char *TempString;	
int CharLimit = 500000000;	
	TempString = malloc (CharLimit*sizeof(char));
	StateLimit = 100000;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;	
	*StateCount = 0;
	while (fgets(TempString, CharLimit, FileHandle) != NULL) (*StateCount)++;
	printf ("State Count: %d\n",(*StateCount));
	rewind(FileHandle);
	*X = malloc(StateLimit*(*StateCount)*sizeof(double));																
	if (*X == NULL) goto Error;
	i = 0;
	while (fscanf (FileHandle, "%lf", &(*X)[i]) == 1) {					//Keep scanning in floating points until there are no more
		i++;
		if (i == (StateLimit*(*StateCount))) {
			printf ("Aww Error: Your state file exceeds J=100, I don' want to deal with that\n");
			goto Error;
		}
	}
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont mater														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	*X = realloc(*X,i*sizeof(double));
	*StatePoints = (int) i/(*StateCount);
	*FileDelta = 2.0/(*StatePoints-1);
	free(TempString);
	return i;
Error:
	printf ("Error Loading file %s\n",FileName);
	return 0;
}

int Load_ETau_File2 (char *FileName, struct ETauStruct *StructToLoad, int *StateCount, int Verbose) 
{
//This loads the E_Tau file in by reading a single value from the file until there are no more values found in the file
//This unwraps the 2D file into a single long array for contiguousness, files should have a state's ET values in a single ROW, not column
int i,StateLimit;
FILE *FileHandle;
char *TempString;
	StateLimit = 10000;
	TempString = malloc(MAXLINESIZE*sizeof(char));
	if (TempString ==  NULL) {
		printf ("Error in Load_ETau_File2: Unable to allocate the read buffer. Consider lowering CHRLIMIT and recompiling\n");
		goto Error;
	}
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) {
		printf ("Error in Load_ETau_File2: Can't open file %s\n",FileName);
		goto Error;	
	}
	*StateCount = 0;
	while (fgets(TempString, MAXLINESIZE, FileHandle) != NULL) (*StateCount)++;	//Run through the file and keep going until we hit the end, track the number of lines/states
	rewind(FileHandle);	//Rewind to the start of the file
	
	//Currently we're just block allocating the max amount of space then deallocating what isnt actually used. 
	//At some point we should change this to work out the exact amount and only allocate what we need.
	(*StructToLoad).ETVals = malloc(StateLimit*(*StateCount)*sizeof(double));																
	if ((*StructToLoad).ETVals == NULL) {
		printf ("Error in Load_ETau_File2: Unable to allocate space for the data\n");
		goto Error;
	}
	i = 0;
	while (fscanf (FileHandle, "%lf", &((*StructToLoad).ETVals[i])) == 1) {					//Keep scanning in floating points until there are no more
		i++;
		if (i == (StateLimit*(*StateCount))) {
			printf ("Error: Your state file exceeds J=100, I don' want to deal with that\n");
			goto Error;
		}
	}
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont mater														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	(*StructToLoad).ETVals = realloc((*StructToLoad).ETVals,i*sizeof(double));
	(*StructToLoad).StatePoints = (int) i/(*StateCount);
	(*StructToLoad).Delta = 2.0/((*StructToLoad).StatePoints-1);
	if (i%(*StateCount) != 0) {
		printf ("Warning: There are an extra %d points in the ET file. This is likely an issue that needs to be resolved. Take results from this run with caution\n",i%(*StateCount));
	}
	if (Verbose) {
		printf ("=========Verbose Load_ETau_File2=========\n");
		printf("Loaded ET File %s with %d states, %d points per state or a delta kappa of %.2e\n",FileName,(*StateCount),(*StructToLoad).StatePoints,(*StructToLoad).Delta);
		printf ("=======================================\n");
	}
	free(TempString);
	return i;
Error:
	printf ("Error Loading file %s\n",FileName);
	return 0;
}

int Load_Base_Catalog (char *FileName, struct Transition **BaseCatalog,  int Verbose)
{
int i,FileLimit;
FILE *FileHandle;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;									//Confrim file opened without error			//Set line count to 0
	FileLimit = 200000;
	i = 0;
	
	*BaseCatalog = malloc(FileLimit*sizeof(struct Transition));
	while (fscanf (FileHandle, "%d %d %d", &(*BaseCatalog)[i].Upper, &(*BaseCatalog)[i].Lower, &(*BaseCatalog)[i].Type) == 3) {
		(*BaseCatalog)[i].Map = i;
		i++;
		if (i >= (FileLimit-1)) {
			printf ("Error: Base Catalog file exceeds %d transitions, I don't believe you\n", FileLimit);
			goto Error;
		}
		if (Verbose) printf ("Transition1: %d\tTransition2: %d\tType:%d\n", (*BaseCatalog)[i-1].Upper, (*BaseCatalog)[i-1].Lower,(*BaseCatalog)[i-1].Type);
	}
	if (Verbose) printf ("Loaded %d transitions from the Base Catalog File: %s\n",i,FileName);
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont matter														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	*BaseCatalog = realloc(*BaseCatalog,i*sizeof(struct Transition));
	return i;	//Return the number of states in the catalog
Error:
	printf ("Error loading Base Catalog file %s\n",FileName);
	return -1;
}

int Load_Base_Catalog_Dictionary (char *FileName, struct Level **DictIn,  int Verbose)
{
int i,FileLimit;
FILE *FileHandle;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;									//Confrim file opened without error	
	FileLimit = 200000;
	i = 0;																//Set line count to 0
	*DictIn = malloc(FileLimit*sizeof(struct Level));										
	while (fscanf (FileHandle, "%d %d %d %d", &(*DictIn)[i].Index, &(*DictIn)[i].J, &(*DictIn)[i].Ka, &(*DictIn)[i].Kc) == 4) {
		i++;
		if (i >= (FileLimit-1)) {
			printf ("Error: Base Catalog Dictionary file exceeds %d levels, I don't believe you\n", FileLimit);
			goto Error;
		}
		if (Verbose) printf ("Index: %d\tJ: %d\tKa: %d\tKc: %d\n", (*DictIn)[i-1].Index, (*DictIn)[i-1].J,(*DictIn)[i-1].Ka,(*DictIn)[i-1].Kc);
	}
	
	if (Verbose) printf ("Loaded %d levels from the Base Catalog Dictionary File: %s\n",i,FileName);
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont matter																					
	*DictIn = realloc(*DictIn,i*sizeof(struct Level));
	return i;
Error:
	printf ("Error loading Base Catalog Dictionary file %s\n",FileName);
	return -1;
}

int Load_Exp_File  (char *FileName, double **X, double **Y, int Verbose)
{
//Load the experimental data file and store it as an X-Y array
//We do not handle headers at all, any header on the file will ruin this function, probably
int i,FileLimit;
FILE *FileHandle;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;									//Confrim file opened without error															//Set line count to 0
	FileLimit = 200000;
	i = 0;
	*X = malloc(FileLimit*sizeof(double));									
	*Y = malloc(FileLimit*sizeof(double));						
	while (fscanf (FileHandle, "%lf %lf", &(*X)[i], &(*Y)[i]) == 2) {
		i++;
		if (i >= (FileLimit-1)) {
			if (i > 100000000) {
				printf ("Error: Experimental file exceeds 10 million points, dial that down a bit\n");
				goto Error;
			}
			FileLimit += 200000;
			*X = realloc (*X,FileLimit*sizeof(double));
			*Y = realloc (*Y,FileLimit*sizeof(double));
		}
	}
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont mater														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	*X = realloc(*X,i*sizeof(double));
	*Y = realloc(*Y,i*sizeof(double));
	if (Verbose) {
		int j;
		double YMax,YMin,Avg;
		Avg = 0.0;
		YMax = -1.0E+307;
		YMin = +1.0E+307;
		for (j=0;j<i;j++) {
			if ((*Y)[j] > YMax) YMax = 	(*Y)[j];
			if ((*Y)[j] < YMin) YMin = 	(*Y)[j];
			Avg += (*Y)[j];
		}	
		printf ("=========Verbose Load_Exp_File=========\n");
		printf ("Spectra X max/min are : %.4f/%.4f\n",(*X)[0],(*X)[i-1]);
		printf ("Spectra Y max/min are : %.4f/%.2e\n",YMax,YMin);
		printf ("Average amplitude: %.3f\n",Avg/i);
		printf ("=======================================\n");
	}
	return i;
Error:
	printf ("Error Loading file %s\n",FileName);
	return 0;
}

int Load_Str_File (char *FileName, double ***Data, int Verbose) 
{
//Function to load in a line strength/Sij file

int i,j,PointsPerState,StateCount;
double Test;
FILE *FileHandle;
int CharLimit = MAXLINESIZE;	
char *TempString;	
	TempString = malloc(CharLimit*sizeof(char));
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) {
		printf ("Error loading Sij data: Cannot open file %s\n",FileName);
		goto Error;
	}
	StateCount = 0;
	while (fgets(TempString, CharLimit, FileHandle) != NULL) StateCount++;	//Run through the file and keep going until we hit the end, track the number of lines/states
	rewind(FileHandle);	//Rewind to the start of the file
	i=0;
	while ((fscanf(FileHandle, "%lf", &Test) == 1)) i++;
	rewind(FileHandle);	//Rewind to the start of the file
	PointsPerState = i/StateCount;
	if (i%StateCount != 0) {
		printf ("Warning: There are an extra %d points in the Sij file. This is likely an issue that needs to be resolved. Take results from this run with caution\n",i%StateCount);
	}
	if (*Data == NULL) {
		*Data = malloc (StateCount*sizeof(double *));
		for (i=0;i<StateCount;i++) {	
			(*Data)[i] = malloc(PointsPerState*sizeof(double));
			if ((*Data)[i] == NULL) printf ("Null alloc\n");
		}
	} else {
		//Ive set this up to do the allocation here, so if it was done elsewhere theres no way to be sure it was set correctly. so we send it back if it was done ahead of time
		printf ("Please send a clean unallocated pointer to this function\n");
		goto Error;
	}
	//Loop over the file loading stuff in
	for (i=0;i<StateCount;i++) {
		for (j=0;j<PointsPerState;j++) {
			fscanf (FileHandle,"%lg",&(*Data)[i][j]);
		}
	}
	free(TempString);
	if (Verbose) {
		printf ("======Load_Str_File======");
		printf ("Loaded %d States with %d points per state\n",StateCount,PointsPerState);
		printf ("=========================");
	}
	return StateCount;
Error:
	return 0;
}

////////////////////////////////////
double Get_Kappa (double A, double B, double C) 
{
	//Computes Ray's asymmetry parameter
	return ((2.0*B-A-C)/(A-C));		
}

double Get_Delta (double A, double B, double C)
{
	//Function to get the inertial defect 
	return 505379.009*((1.0/C)-(1.0/B)-(1.0/A));
}

int Get_J (int TransitionIndex, struct Level *MyDictionary)
{
	return MyDictionary[TransitionIndex].J; 
}

double Partition_Function (double *Constants, double Temperature) 
{
//Gordy and Cook partition function approximation function. 
//Not currently used but left here in case it becomes useful
	return (5.34E+6*sqrt(pow(Temperature,3.0)/(Constants[0]*Constants[1]*Constants[2])));
}

double E_tau (int TransitionIndex, double Kappa, struct ETauStruct ETStruct) 
{
//A function to return the E_tau() value of a rigid rotor Hamiltonian
//Really this just fetches a value for the appropriate value of J/Ka/Kc
//Values are explicitly calculated when E_tau has an analytic form, for all other values we use a look up table calculated elsewhere
int Index;
	//Mapping trick to keep Kappa in line in bad fits
	//if we exceed 1 or -1 we can jump into the wrong state and get off track fast or just segfault
	//This remaps Kappa onto a normalized arctangent function so it can never exceed |1|, but asymptotically approaches it as the program ramps the real kappa
	//The 100 is to get it in a reasonable range for the function
	//This has the potential to return poor results on its own, but is better than the alternative
	if (Kappa > 1.0) Kappa = atan(100.0*Kappa)/1.570796326794896619231321691639;
	if (Kappa < -1.0) Kappa = atan(100.0*Kappa)/1.570796326794896619231321691639;
	Index = (int) ((Kappa+1.0)/ETStruct.Delta);
	switch (TransitionIndex) {	//This version just switches by a transition index, aka lookup table
		case 0:																																//000 (0)
			return 0.0;	//The easiest one
		//============ J = 1 ============
		case 1:
			return (Kappa-1.0);																												//101 (-1)
		case 2:
			return 0.0;																														//111 (0)			
		case 3:
			return (Kappa+1);																												//110 (+1)
		//============ J = 2 ============
		case 4:
			return (2.0*(Kappa-sqrt(Kappa*Kappa+3.0)));																						//202 (-2)
		case 5:
			return (Kappa-3.0);																												//212 (-1)
		case 6:
			return	(4.0*Kappa);																											//211 (0)
		case 7:
			return 	(Kappa+3.0);																											//221 (+1)
		case 8:
			return (2.0*(Kappa+sqrt(Kappa*Kappa+3.0)));																						//220 (+2)
		//============ J = 3 ============
		case 9:
			return (5.0*Kappa-3.0-2.0*sqrt(4.0*Kappa*Kappa+6.0*Kappa+6.0));																	//303 (-3)
		case 10:
			return 2.0*(Kappa-sqrt(Kappa*Kappa+15.0));																						//313 (-2)
		case 11:
			return (5.0*Kappa+3.0-2.0*sqrt(4.0*Kappa*Kappa-6.0*Kappa+6.0));																	//312 (-1)
		case 12:
			return  4.0*Kappa;																												//322 (0)
		case 13:
			return (5.0*Kappa-3.0+2.0*sqrt(4.0*Kappa*Kappa+6.0*Kappa+6.0));																	//321 (+1)
		case 14:
			return 2.0*(Kappa+sqrt(Kappa*Kappa+15.0));																						//331 (+2)
		case 15:
			return (5.0*Kappa+3.0+2.0*sqrt(4.0*Kappa*Kappa-6.0*Kappa+6.0));																	//330 (+3)
		//============ J = 4 ============	
		case 17: 
			return 5.0*Kappa-5.0-2.0*sqrt(4.0*Kappa*Kappa+10.0*Kappa+22.0);																	//414 (-3)	
		case 18:
			return 10.0*Kappa-2.0*sqrt(9.0*Kappa*Kappa+7.0);																				//413 (-2)	
		case 19:
			return 5.0*Kappa+5.0-2.0*sqrt(4.0*Kappa*Kappa-10.0*Kappa+22.0);																	//423 (-1)
		case 21:
			return 5.0*Kappa-5.0+2.0*sqrt(4.0*Kappa*Kappa+10.0*Kappa+22.0);																	//432 (+1)
		case 22:
			return 10.0*Kappa+2.0*sqrt(9.0*Kappa*Kappa+7.0);																				//431 (+2)
		case 23:
			return 5.0*Kappa+5.0+2.0*sqrt(4.0*Kappa*Kappa-10.0*Kappa+22.0);																	//441 (+3)
		//============ J = 5 ============	
		case 28:
			return 10.0*Kappa-6.0*sqrt(Kappa*Kappa+3.0);																					//524 (-2)	
		case 32:
			return 10.0*Kappa+6.0*sqrt(Kappa*Kappa+3.0);																					//542 (+2)														
		//============ J >= 6 ============
		//For any non hardcoded state we do the actual math
		default: 
			return ETStruct.ETVals[Index+ETStruct.StatePoints*TransitionIndex]+((ETStruct.ETVals[Index+ETStruct.StatePoints*TransitionIndex+1]-ETStruct.ETVals[Index+ETStruct.StatePoints*TransitionIndex])/ETStruct.Delta)*(Kappa-((Index*ETStruct.Delta)-1.0));
	}	
	printf ("Bad %f\n", Kappa);
	return -1.0;
}

double Rigid_Rotor (double A, double C, int J, int Index, double Kappa, struct ETauStruct ETStruct)
{
//Ease of use function to compute the energy of a single rigid rotor level
//The first term is basically just total angular momentum, the second is an asymmetric rotor/state-dependent correction
//Happily frequencies/energies are computed in whatever units the rotational constants are given in
	return 0.5*(A+C)*J*(J+1.0)+0.5*(A-C)*E_tau(Index,Kappa,ETStruct);
}

double Rigid_Rotor_Error (double A, double C, int J, int Index, double Kappa, struct ETauStruct ETStruct, double *Hkappa)
{
//Ease of use function to compute the energy of a single rigid rotor level
//The first term is basically just total angular momentum, the second is an asymmetric rotor/state-dependent correction
//Happily frequencies/energies are computed in hbar units the rotational constants are given in
	*Hkappa = E_tau(Index,Kappa,ETStruct);	
	return 0.5*(A+C)*J*(J+1.0)+0.5*(A-C)*(*Hkappa);
}

double Get_Frequency (int J_Up, int J_Low, int IndexUp, int IndexLow, double *Constants, struct ETauStruct ETStruct)
{
//Ease of use function to compute the frequency of an asymmetric rotor
double Kappa;
	Kappa = Get_Kappa (Constants[0],Constants[1],Constants[2]);
	return fabs(Rigid_Rotor(Constants[0],Constants[2],J_Up,IndexUp,Kappa,ETStruct)-Rigid_Rotor(Constants[0],Constants[2],J_Low,IndexLow,Kappa,ETStruct));
}

double Get_Frequency_Error (int J_Up, int J_Low, int IndexUp, int IndexLow, double *Constants, struct ETauStruct ETStruct, double *TransitionError, double *ConstantsError)
{
//Ease of use function to compute the frequency of an asymmetric rotor
double Kappa,H1,H2,E1,E2,AC;
	Kappa = Get_Kappa (Constants[0],Constants[1],Constants[2]);
	E1 = Rigid_Rotor_Error(Constants[0],Constants[2],J_Up,IndexUp,Kappa,ETStruct,&H1);
	E2 = Rigid_Rotor_Error(Constants[0],Constants[2],J_Low,IndexLow,Kappa,ETStruct,&H2);
	AC = sqrt(ConstantsError[0]*ConstantsError[0]+ConstantsError[2]*ConstantsError[2]);
	*TransitionError = sqrt((AC*0.5*(J_Up*(J_Up+1.0)+H1))*(AC*0.5*(J_Up*(J_Up+1.0)+H1))+(AC*0.5*(J_Low*(J_Low+1.0)+H2))*(AC*0.5*(J_Low*(J_Low+1.0)+H2)));
	return fabs(E1-E2);
}

int Get_Catalog (struct Transition *restrict CatalogtoFill, double *restrict Constants, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary)
{
//Utility function for calculating frequencies of a catalog
int i;	//Declaring i here because I like it, and apparently learned C pre C99
	for (i=0;i<CatLines;i++) {	
		(CatalogtoFill)[i].Frequency = Get_Frequency (	MyDictionary[(CatalogtoFill)[i].Upper].J,
														MyDictionary[(CatalogtoFill)[i].Lower].J,
														(CatalogtoFill)[i].Upper,
														(CatalogtoFill)[i].Lower,
														Constants,  
														ETStruct
		);
		if (Verbose) printf ("%d %d ",i,(CatalogtoFill)[i].Type);		
		if (Verbose) print_Transition ((CatalogtoFill)[i],MyDictionary);
	} 
	return 1;
}

int Get_Catalog_Error (struct Transition *restrict CatalogtoFill, double *restrict Constants, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary, double *ConstantsError)
{
//Utility function for calculating frequencies of a catalog
int i;	//Declaring i here because I like it, and apparently learned C pre C99
	for (i=0;i<CatLines;i++) {	
		(CatalogtoFill)[i].Frequency = Get_Frequency_Error (	MyDictionary[(CatalogtoFill)[i].Upper].J,
																MyDictionary[(CatalogtoFill)[i].Lower].J,
																(CatalogtoFill)[i].Upper,
																(CatalogtoFill)[i].Lower,
																Constants,  
																ETStruct,
																&((CatalogtoFill)[i].Error),
																ConstantsError
																
		);
		if (Verbose) printf ("%d ",(CatalogtoFill)[i].Type);
		if (Verbose) printf ("%.4f ",(CatalogtoFill)[i].Error);		
		if (Verbose) print_Transition ((CatalogtoFill)[i],MyDictionary);
	} 
	return 1;
}

double Get_Str (double Kappa, int Transition, int PointsPerSate, double **StrData)
{
// Function to pull out the correct(ish) Sij value 
// This uses integer math and no rounding to get the index of the array
// This will give a slightly inaccurate result, but should be close enough, at some point I could add an interpolation scheme
double RetStr;
int Index;		
	RetStr = 0.0;
	Index = (int) ((Kappa+1.0)*PointsPerSate*0.5);
	RetStr = StrData[Transition][Index];
	return RetStr;	
}

////////////////////////////////////
void print_Transition (struct Transition TransitionToPrint, struct Level *MyDictionary)
{
//Utility function for printing transitions
	printf ("%.3f %i %i %i -- %i %i %i\n",TransitionToPrint.Frequency, MyDictionary[TransitionToPrint.Upper].J,MyDictionary[TransitionToPrint.Upper].Ka,MyDictionary[TransitionToPrint.Upper].Kc,MyDictionary[TransitionToPrint.Lower].J,MyDictionary[TransitionToPrint.Lower].Ka,MyDictionary[TransitionToPrint.Lower].Kc);
}

void Sort_Catalog (struct Transition *CatalogtoSort, int TransitionCount, int SortMethod, int SortType) 
{
//Wrapper function for sorting a catalog
	int (*Catalog_Comparator)(const void *,const void *);
	switch (SortType) {
		case 0:
			Catalog_Comparator = &Catalog_Comparator_Frequency;
			break;
		case 1:
			Catalog_Comparator = &Catalog_Comparator_Intensity;	
			SortMethod = 1;
			break;
		case 2:
			Catalog_Comparator = &Catalog_Comparator_Index_Upper;
			SortMethod = 1;	
			break;
		case 3:
			Catalog_Comparator = &Catalog_Comparator_Index_Lower;
			SortMethod = 1;
			break;
		default:
			Catalog_Comparator = &Catalog_Comparator_Frequency;	
			break;
	}
	switch (SortMethod) {
		case 1:
			//Quick Sort - Fast for general sorting
			//This should be slower than insertion sort, but is here as a backup
			qsort(CatalogtoSort, TransitionCount, sizeof(struct Transition), Catalog_Comparator);
			break;
		default:
			//Insertion Sort - Use this one
			//Faster than quicksort when the target is already somewhat sorted
			//The program should in general produce fairly ordered catalogs already, so this ought to be the better choice
			insertionSort(CatalogtoSort, TransitionCount); 
			break;
	}
}

int Catalog_Comparator_Frequency (const void *a, const void *b) 
{
//Comparison function for qsort sorting of the catalog
	struct Transition A = *(struct Transition *) a;
	struct Transition B = *(struct Transition *) b;
	if (A.Frequency > B.Frequency) return 1;
	else if (A.Frequency < B.Frequency) return -1;
	else return 0;
}

int Catalog_Comparator_Intensity (const void *a, const void *b) 
{
//Comparison function for qsort sorting of the catalog
	struct Transition A = *(struct Transition *) a;
	struct Transition B = *(struct Transition *) b;
	if (A.Intensity < B.Intensity) return 1;
	else if (A.Intensity > B.Intensity) return -1;
	else return 0;
}

int Catalog_Comparator_Index_Upper (const void *a, const void *b) 
{
//Comparison function for qsort sorting of the catalog
	struct Transition A = *(struct Transition *) a;
	struct Transition B = *(struct Transition *) b;
	if (A.Upper > B.Upper) return 1;
	else if (A.Upper < B.Upper) return -1;
	else return 0;
}

int Catalog_Comparator_Index_Lower (const void *a, const void *b) 
{
//Comparison function for qsort sorting of the catalog
	struct Transition A = *(struct Transition *) a;
	struct Transition B = *(struct Transition *) b;
	if (A.Lower > B.Lower) return 1;
	else if (A.Lower < B.Lower) return -1;
	else return 0;
}

void insertionSort(struct Transition *CatalogtoSort, int TransitionCount) 
{ 
//Insertion sort for sorting the catalog
struct Transition Key;
int i, j; 
	for (i = 1; i < TransitionCount; i++) { 
		Key = CatalogtoSort[i]; 
        j = i - 1; 
        while (j >= 0 && CatalogtoSort[j].Frequency > Key.Frequency) { 
            CatalogtoSort[j + 1] = CatalogtoSort[j]; 
            j = j - 1; 
        } 
        CatalogtoSort[j + 1] = Key; 
    } 
} 

int Fill_Catalog_Restricted_Frequency (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, double FrequencyLow, double FrequencyHigh, int Verbose, struct Level *MyDictionary)
{
int i, CatLinesOut;     
	*CatalogtoFill = malloc(CatLines*sizeof(struct Transition));
	CatLinesOut = 0; 
	for (i=0;i<CatLines;i++) {
		if (((SourceCatalog)[i].Frequency > FrequencyLow) && ((SourceCatalog)[i].Frequency < FrequencyHigh)) {
			(*CatalogtoFill)[CatLinesOut] = (SourceCatalog)[i]; 	//Fairly certain this is a value copy not an address copy
			if (Verbose) printf ("%d ",(*CatalogtoFill)[CatLinesOut].Type);				
			if (Verbose) print_Transition ((*CatalogtoFill)[CatLinesOut],MyDictionary);
			CatLinesOut++;
		}
     } 
	*CatalogtoFill = realloc (*CatalogtoFill,CatLinesOut*sizeof(struct Transition));   
	return CatLinesOut;
}

int Fill_Catalog_Restricted_J (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, int JMin, int JMax, int Verbose, struct Level *MyDictionary)
{
int i, CatLinesOut;     
	*CatalogtoFill = malloc(CatLines*sizeof(struct Transition));
	CatLinesOut = 0; 
	for (i=0;i<CatLines;i++) {
		
		if ((MyDictionary[(SourceCatalog)[i].Upper].J >= JMin) && (MyDictionary[(SourceCatalog)[i].Upper].J <= JMax)) {
			(*CatalogtoFill)[CatLinesOut] = (SourceCatalog)[i];
			if (Verbose) printf ("%d ",(*CatalogtoFill)[CatLinesOut].Type);				
			if (Verbose) print_Transition ((*CatalogtoFill)[CatLinesOut],MyDictionary);
			CatLinesOut++;
		}
     } 
	*CatalogtoFill = realloc (*CatalogtoFill,CatLinesOut*sizeof(struct Transition));   
	return CatLinesOut;
}

int Fill_Catalog_Restricted_Intensity (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, double Percentile, int Verbose, struct Level *MyDictionary)
{
int i, CatLinesOut;     
	if (Percentile > 1.0) {
		printf ("Error: you need to specify a value between 0 and 1 for catalog intensity trimming");
		goto Error;
	}
	CatLinesOut = (int) CatLines*Percentile; 
	*CatalogtoFill = malloc(CatLinesOut*sizeof(struct Transition));
	Sort_Catalog (SourceCatalog, CatLines, 1, 1); 
	for (i=0;i<CatLinesOut;i++) {
		(*CatalogtoFill)[CatLinesOut] = (SourceCatalog)[i];
		if (Verbose) printf ("%d ",(*CatalogtoFill)[CatLinesOut].Type);				
		if (Verbose) print_Transition ((*CatalogtoFill)[CatLinesOut],MyDictionary);
     } 
	return CatLinesOut;
Error:
	return -1;
}

int Fill_Catalog_Restricted_Intensity_Count (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, int CatLinesOut, int Verbose, struct Level *MyDictionary)
{
//Count-based variant of the Fill Catalog by intensity function
//Returns the top-x most intense lines in a catalog
int i;     
	if (CatLinesOut > CatLines) CatLinesOut = CatLines;
	*CatalogtoFill = malloc(CatLinesOut*sizeof(struct Transition));
	if (*CatalogtoFill == NULL) goto Error;
	Sort_Catalog (SourceCatalog, CatLines, 1, 1); 
	for (i=0;i<CatLinesOut;i++) {
		(*CatalogtoFill)[i] = (SourceCatalog)[i];
		if (Verbose) printf ("%d ",i);
		if (Verbose) printf ("%f ",(*CatalogtoFill)[i].Frequency);
		if (Verbose) printf ("%d ",(*CatalogtoFill)[i].Type);				
		if (Verbose) print_Transition ((*CatalogtoFill)[i],MyDictionary);
     }
	Sort_Catalog (*CatalogtoFill, CatLinesOut, 0, 0); 
	return CatLinesOut;
Error:
	return -1;
}

int Fill_Catalog_Restricted_Ka (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, int KaMin, int KaMax, int Verbose, struct Level *MyDictionary)
{
int i, CatLinesOut;     
	*CatalogtoFill = malloc(CatLines*sizeof(struct Transition));
	CatLinesOut = 0; 
	for (i=0;i<CatLines;i++) {
		
		if ((MyDictionary[(SourceCatalog)[i].Upper].Ka >= KaMin) && (MyDictionary[(SourceCatalog)[i].Upper].Ka <= KaMax)) {
			(*CatalogtoFill)[CatLinesOut] = (SourceCatalog)[i];
			if (Verbose) printf ("%d ",(*CatalogtoFill)[CatLinesOut].Type);				
			if (Verbose) print_Transition ((*CatalogtoFill)[CatLinesOut],MyDictionary);
			CatLinesOut++;
		}
     } 
	*CatalogtoFill = realloc (*CatalogtoFill,CatLinesOut*sizeof(struct Transition));   
	return CatLinesOut;
}

void Calculate_State_Energies (struct Level *MyDictionary, struct Transition *SourceCatalog, int CatalogTransitions, int DictionaryLevels, int Verbose)
{
//Function to get the energies of levels in a catalog
int i;
	(MyDictionary[0]).Energy = 0.0;
	for (i=0;i<CatalogTransitions;i++) (MyDictionary[SourceCatalog[i].Upper]).Energy = (MyDictionary[SourceCatalog[i].Lower]).Energy+SourceCatalog[i].Frequency;	//Keeping the calculation in MHz because it's easier
	for (i=0;i<DictionaryLevels;i++) {
		(MyDictionary[i]).Energy *= 4.8E-5;
	}
	if (Verbose) for (i=0;i<CatalogTransitions;i++) printf ("%d State Energy:%f Upper Index: %d LowerIndex: %d\n", i, (MyDictionary[SourceCatalog[i].Upper]).Energy, SourceCatalog[i].Upper, SourceCatalog[i].Lower);
}

void Calculate_Intensities (struct Transition *SourceCatalog, int CatalogTransitions, struct Level *MyDictionary, double T, double *Dipoles, int Verbose)
{
int i;
	for (i=0;i<CatalogTransitions;i++) {
		(SourceCatalog[i]).Intensity = Dipoles[SourceCatalog[i].Type-1]*Dipoles[SourceCatalog[i].Type-1]*SourceCatalog[i].Frequency*fabs(	exp(-1.0*(MyDictionary[SourceCatalog[i].Lower]).Energy/T)-exp(-1.0*(MyDictionary[SourceCatalog[i].Upper].Energy)/T));
		if (Verbose) printf ("%e %f %f\n",(SourceCatalog[i]).Intensity,(MyDictionary[SourceCatalog[i].Lower]).Energy,(MyDictionary[SourceCatalog[i].Upper]).Energy);
	}
}

void Calculate_Intensities_Sij (struct Transition *SourceCatalog, int CatalogTransitions, struct Level *MyDictionary, double T, double *Dipoles, double **StrData, double Kappa, double *Constants, int Verbose)
{
int i;
double S,Q;
	Q = Partition_Function (Constants,T);
	for (i=0;i<CatalogTransitions;i++) {
		S = Get_Str (Kappa, SourceCatalog[i].Map, 4000, StrData);
		(SourceCatalog[i]).Intensity = S*(1.0/Q)*4.16231E-5*Dipoles[SourceCatalog[i].Type-1]*Dipoles[SourceCatalog[i].Type-1]*SourceCatalog[i].Frequency*fabs(exp(-1.0*(MyDictionary[SourceCatalog[i].Lower]).Energy/T)-exp(-1.0*(MyDictionary[SourceCatalog[i].Upper].Energy)/T));
		if (Verbose) printf ("%e %f %f\n",(SourceCatalog[i]).Intensity,(MyDictionary[SourceCatalog[i].Lower]).Energy,(MyDictionary[SourceCatalog[i].Upper]).Energy);
	}
}

int Standardize_Catalog (struct Transition *CatalogtoStandardize, int CatLines, int Verbose)
{
//Utility function for standardizing a catalog
int i;	//Declaring i here because I like it, and apparently learned C pre C99Constants
double Sum,Mean,StdDev;
	Sum = 0.0;
	StdDev = 0.0;
	for(i=0; i<CatLines; ++i) Sum += CatalogtoStandardize[i].Frequency;
    Mean = Sum/CatLines;
    for(i=0; i<CatLines; ++i) StdDev += pow(CatalogtoStandardize[i].Frequency - Mean, 2.0);
    StdDev = sqrt(StdDev/CatLines);
	for (i=0;i<CatLines;i++) CatalogtoStandardize[i].Frequency = (CatalogtoStandardize[i].Frequency-Mean)/StdDev;
	if (Verbose) printf ("Catalog Mean: %f\tCatalog Standard Deviation: %f\n",Mean,StdDev);
	return 1;
}

////////////////////////////////////
int Peak_Find (double **LineList, double Max, double Min, double *X, double *Y, int ArraySize, int Verbose)
{
//Pretty standard implementation of a peakfinding algorithm
//If a point is within the min/max amplitude and higher than the two points on either side we consider it a peak		
int i,PeakCount;	
	*LineList = malloc((int)(ArraySize/3)*sizeof(double));		//Allocate an array as large as one third the size of the experimental array, this is essentially the theoretical max # of peaks and we should never get close to that
	PeakCount = 0;
	for (i=1;i<ArraySize-1;i++) {
		if (((Y)[i] > Min) & ((Y)[i] < Max)) {		//Checking for min max first, not sure if this is faster or if this really matters at all, probably not
			if (((Y)[i-1] <(Y)[i]) & ((Y)[i+1] <(Y)[i])) {
				(*LineList)[PeakCount] = X[i];
				PeakCount++;
			}
		}
	}	
	if (Verbose) {
		printf ("=====Peak_Find=====\n");
		printf ("Found %d peaks\n",PeakCount);
		printf ("===================\n");
	}
	*LineList = realloc(*LineList,PeakCount*sizeof(double));	
	return PeakCount;
}

int Find_Triples (struct Triple *TripletoFit, double *LineFrequencies, double Window, int LineCount, int Verbose)
{
int i,Count;	
	TripletoFit->TriplesList = malloc(LineCount*3*sizeof(double));	//Allocate an array that is 3x the total number of lines we could fit, thats the max, assuming every line is within the search window of all three transitions
	TripletoFit->TriplesCount[0] = 0;
	TripletoFit->TriplesCount[1] = 0;
	TripletoFit->TriplesCount[2] = 0;
	i=0;
	Count = 0;
	while (LineFrequencies[i] < TripletoFit->TransitionList[0].Frequency-Window) {
		i++;
	}
	while ((LineFrequencies[i] < TripletoFit->TransitionList[0].Frequency+Window) & (i<LineCount)) {	//As long as the line is lower than the max frequency of the window and we aren't at the end of the array, add it to the list
		TripletoFit->TriplesList[Count] = LineFrequencies[i];
		Count++;
		i++;
	}	
	TripletoFit->TriplesCount[0] = Count;
	i=0;	//Reset the index above, since we have no way of knowing if these are in order of frequency or if they over lap. This is the safest, but slowest. Sorting could speed this up, but it is more work, something to add later if bigggg files start giving this issues
	Count = 0;
	while (LineFrequencies[i] < TripletoFit->TransitionList[1].Frequency-Window) {
		i++;
	}
	while ((LineFrequencies[i] < TripletoFit->TransitionList[1].Frequency+Window) & (i<LineCount)) {
		TripletoFit->TriplesList[Count+TripletoFit->TriplesCount[0]] = LineFrequencies[i];
		Count++;
		i++;
	}
	TripletoFit->TriplesCount[1] = Count;
	i=0;	//Reset the index above again
	Count = 0;
	while (LineFrequencies[i] < TripletoFit->TransitionList[2].Frequency-Window) {
		i++;
	}
	while ((LineFrequencies[i] < TripletoFit->TransitionList[2].Frequency+Window) & (i<LineCount)) {
		TripletoFit->TriplesList[Count+TripletoFit->TriplesCount[0]+TripletoFit->TriplesCount[1]] = LineFrequencies[i];
		Count++;
		i++;
	}
	TripletoFit->TriplesCount[2] = Count;
	//Check this alloc
	TripletoFit->TriplesList = realloc(TripletoFit->TriplesList,TripletoFit->TriplesCount[0]*TripletoFit->TriplesCount[1]*TripletoFit->TriplesCount[2]*sizeof(double));	
	printf ("=====Find_Triples=====\n");
	printf ("%d experimental lines sent to function\n",LineCount);
	printf ("Found %d lines for %f\n",TripletoFit->TriplesCount[0],TripletoFit->TransitionList[0].Frequency);
	printf ("Found %d lines for %f\n",TripletoFit->TriplesCount[1],TripletoFit->TransitionList[1].Frequency);
	printf ("Found %d lines for %f\n",TripletoFit->TriplesCount[2],TripletoFit->TransitionList[2].Frequency);
	printf ("======================\n");
	return 1;
}

double Fit_Triples (double *Guess, int Verbose, struct Opt_Bundle GSLOptBundle)
{
//Fit a set of three lines to A/B/C
//This function is semi-obsolete. Fitting a single set of transitions like this is hugely wasteful as you have to reallocate and free all the variables for each fit. Better to fit a full set of triples, or even to use a single workspace for all fits
	int info;
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  	gsl_multifit_nlinear_workspace *Workspace;
  	gsl_multifit_nlinear_fdf fdf;
  	gsl_vector_view x;
  	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  	const size_t n = 3;
  	const size_t p = 3;
  	gsl_vector *f;
  	const double xtol = 1e-8;
  	const double gtol = 1e-8;
  	const double ftol = 1e-1;
 	Workspace = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);
  	fdf.f = OptFunc_gsl;
  	fdf.df = NULL;   	
  	fdf.fvv = NULL;     
  	fdf.n = n;
  	fdf.p = p;
	fdf.params = &GSLOptBundle;
	fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  	x = gsl_vector_view_array (Guess, p);
   	gsl_multifit_nlinear_init (&x.vector, &fdf, Workspace);     	
	f = gsl_multifit_nlinear_residual(Workspace);
  	gsl_multifit_nlinear_driver(20, xtol, gtol, ftol, NULL, NULL, &info, Workspace);
  	gsl_vector *Final = gsl_multifit_nlinear_position(Workspace);
  	printf ("%10.4f %10.4f %10.4f\n",gsl_vector_get(Final, 0),gsl_vector_get(Final, 1),gsl_vector_get(Final, 2));
  	gsl_multifit_nlinear_free (Workspace);
	return 1.0;
}

void Initialize_Triples_Fitter (struct GSL_Bundle *FitBundle, struct Opt_Bundle *MyOpt_Bundle)
{
size_t p = 3;	//These are the sizes of the parameters and data points
size_t n = 3;	//Theyre hard coded because all triples fits are 3 parameters and 3 unknowns
	FitBundle->fdf_params = gsl_multifit_nlinear_default_parameters();
 	FitBundle->fdf.f = OptFunc_gsl;
   	FitBundle->fdf.df = NULL;   //Finite difference Jacobian because there is no general analytic version	
   	FitBundle->fdf.fvv = NULL;	//No geodesic acceleration, early tests showed no real improvement in using it
   	FitBundle->fdf.n = n;
  	FitBundle->fdf.p = p;
 	FitBundle->fdf_params.trs = gsl_multifit_nlinear_trs_lm;
 	FitBundle->T = gsl_multifit_nlinear_trust;
 	FitBundle->Workspace = gsl_multifit_nlinear_alloc (FitBundle->T, &(FitBundle->fdf_params), n, p);
 	FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);
	if (MyOpt_Bundle->TransitionsGSL == NULL) {
		printf ("Allocating transition memory space (%d transitions)\n",MyOpt_Bundle->TransitionCount);
		MyOpt_Bundle->TransitionsGSL = malloc(MyOpt_Bundle->TransitionCount*sizeof(struct Transition));
	}
 	else {
 		printf ("Reallocating transition memory space (%d transitions)\n",MyOpt_Bundle->TransitionCount);
 		MyOpt_Bundle->TransitionsGSL = realloc(MyOpt_Bundle->TransitionsGSL,MyOpt_Bundle->TransitionCount*sizeof(struct Transition));	
	}
}

void Initialize_Triples_Fitter_Alloc (struct GSL_Bundle *FitBundle, struct Opt_Bundle *MyOpt_Bundle)
{
size_t p = 3;	//These are the sizes of the parameters and data points
size_t n = 3;	//Theyre hard coded because all triples fits are 3 parameters and 3 unknowns
	FitBundle->fdf_params = gsl_multifit_nlinear_default_parameters();
 	FitBundle->fdf.f = OptFunc_gsl;
   	FitBundle->fdf.df = NULL;   //Finite difference Jacobian because there is no general analytic version	
   	FitBundle->fdf.fvv = NULL;	//No geodesic acceleration, early tests showed no real improvement in using it
   	FitBundle->fdf.n = n;
  	FitBundle->fdf.p = p;
 	FitBundle->fdf_params.trs = gsl_multifit_nlinear_trs_lm;
 	FitBundle->T = gsl_multifit_nlinear_trust;
 	FitBundle->Workspace = gsl_multifit_nlinear_alloc (FitBundle->T, &(FitBundle->fdf_params), n, p);
 	if (FitBundle->Workspace == NULL) printf ("Null alloc of FitBundle Workspace\n");
 	FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);
	printf ("Allocating transition memory space (%d transitions)\n",MyOpt_Bundle->TransitionCount);
	MyOpt_Bundle->TransitionsGSL = malloc(MyOpt_Bundle->TransitionCount*sizeof(struct Transition));
	if (MyOpt_Bundle->TransitionsGSL == NULL) printf ("Null alloc of transition space\n");

}

int OptFunc_gsl (const gsl_vector *x, void *params, gsl_vector *f)
{
double GSLConstants[3];	//Declare some doubles to hold our constants
	//I dont love doing this, but GSL really only wants one pointer to void for the function parameters. So this is a struct to hold all the things we need for calculating frequencies
	struct Opt_Bundle *p = (struct Opt_Bundle *) params;
	GSLConstants[0] = gsl_vector_get(x, 0);	//Pull these values from the vector
	GSLConstants[1] = gsl_vector_get(x, 1);
	GSLConstants[2] = gsl_vector_get(x, 2);
	gsl_vector_set (f, 0, Get_Frequency(p->MyDictionary[p->TransitionsGSL[0].Upper].J,
										p->MyDictionary[p->TransitionsGSL[0].Lower].J,
										p->TransitionsGSL[0].Upper,
										p->TransitionsGSL[0].Lower,
										GSLConstants,p->ETGSL) - p->TransitionsGSL[0].Frequency);
	
	gsl_vector_set (f, 1, Get_Frequency(p->MyDictionary[p->TransitionsGSL[1].Upper].J,
										p->MyDictionary[p->TransitionsGSL[1].Lower].J,
										p->TransitionsGSL[1].Upper,
										p->TransitionsGSL[1].Lower,
										GSLConstants,p->ETGSL) - p->TransitionsGSL[1].Frequency);
	gsl_vector_set (f, 2, Get_Frequency(p->MyDictionary[p->TransitionsGSL[2].Upper].J,
										p->MyDictionary[p->TransitionsGSL[2].Lower].J,
										p->TransitionsGSL[2].Upper,
										p->TransitionsGSL[2].Lower,
										GSLConstants,p->ETGSL) - p->TransitionsGSL[2].Frequency);
	return GSL_SUCCESS;
}

int Fit_Triples_Bundle (struct Triple TransitionstoFit, double *Guess, double **FitResults, struct Transition **MyFittingCatalog, int CatalogLines, struct GSL_Bundle *FitBundle, struct Opt_Bundle MyOpt_Bundle, ScoreFunction TriplesScoreFunction, void *ScoringParameters)
{
int i,j,k,info,Count,Iterations,Wins,Errors;
const double xtol = 1e-8;
const double gtol = 1e-8;
const double ftol = 1e-1;
const size_t p = 3;
struct Transition Transitions[3];
gsl_vector *Final;
gsl_vector_view x;
  	x = gsl_vector_view_array (Guess, p);							//Set the guess	
  	*FitResults = malloc(3*TransitionstoFit.TriplesCount[0]*TransitionstoFit.TriplesCount[1]*TransitionstoFit.TriplesCount[2]*sizeof(double));
  	//Extract the transitions we'll use for the fit
  	Transitions[0].Upper = TransitionstoFit.TransitionList[0].Upper;
  	Transitions[0].Lower = TransitionstoFit.TransitionList[0].Lower;
  	Transitions[1].Upper = TransitionstoFit.TransitionList[1].Upper;
  	Transitions[1].Lower = TransitionstoFit.TransitionList[1].Lower;
	Transitions[2].Upper = TransitionstoFit.TransitionList[2].Upper;
  	Transitions[2].Lower = TransitionstoFit.TransitionList[2].Lower;  
  	//printf ("%i %i %i %i %i %i\n",TransitionstoFit.TransitionList[0].Upper,TransitionstoFit.TransitionList[0].Lower,TransitionstoFit.TransitionList[1].Upper,TransitionstoFit.TransitionList[1].Lower,TransitionstoFit.TransitionList[2].Upper,TransitionstoFit.TransitionList[2].Lower);
  	Wins = 0;		//Track the total number of wins for the current scoring system
  	Count = 0;		//Track the total number of constants (A+B+C) in the fit results
  	Iterations = 0;	//Variable to track the total number of iterations throughout the fit, just a bookeeping thing for me to see how the fitter is operating
  	Errors = 0;		//A count of the number of unconverged fits, another metric for me to track the fitting
  	double MyConstants[3];
  	FitBundle->fdf.params = &MyOpt_Bundle;
  	for (i=0;i<TransitionstoFit.TriplesCount[0];i++) {
  		for (j=0;j<TransitionstoFit.TriplesCount[1];j++) {
  			for (k=0;k<TransitionstoFit.TriplesCount[2];k++) {
  				//This is why the transitions were extracted. We set the 3 transition's frequencies to those of the triple we're working on and hand it off to the fitter
				MyOpt_Bundle.TransitionsGSL[0].Frequency = TransitionstoFit.TriplesList[i];			
				MyOpt_Bundle.TransitionsGSL[1].Frequency = TransitionstoFit.TriplesList[j+TransitionstoFit.TriplesCount[0]];
				MyOpt_Bundle.TransitionsGSL[2].Frequency =  TransitionstoFit.TriplesList[k+TransitionstoFit.TriplesCount[0]+TransitionstoFit.TriplesCount[1]];								
				FitBundle->fdf.params = &MyOpt_Bundle;															//Set the parameters for this run					
   				gsl_multifit_nlinear_init (&x.vector, &(FitBundle->fdf), FitBundle->Workspace);	//reInitialize the workspace incase this isnt the first run of the loop  	
				FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);								//compute initial cost function
  				gsl_multifit_nlinear_driver(50, xtol, gtol, ftol, NULL, NULL, &info, FitBundle->Workspace);		//solve the system with a maximum of 20 iterations
  				Iterations += gsl_multifit_nlinear_niter (FitBundle->Workspace); 	//Track the iterations
  				Final = gsl_multifit_nlinear_position(FitBundle->Workspace);		//Snag the results
  				if (gsl_multifit_nlinear_niter (FitBundle->Workspace) == 50) {		//Check for an error, currently only considering non-convergence
  					//printf ("Error: Unconverged Fit: %.4f %.4f %.4f\n",Transitions[0].Frequency,Transitions[1].Frequency,Transitions[2].Frequency);
  					Errors++;
  				}
   				(*FitResults)[Count] = gsl_vector_get(Final, 0);	//Get the final A
  				MyConstants[0] = gsl_vector_get(Final, 0);				//Update the constants
  				Count++;											//Track the number of items in FitResults
  				(*FitResults)[Count] = gsl_vector_get(Final, 1);	
  				MyConstants[1] = gsl_vector_get(Final, 1);
  				Count++;
  				(*FitResults)[Count] = gsl_vector_get(Final, 2);
  				MyConstants[2] = gsl_vector_get(Final, 2);
  				Count++;
  				Get_Catalog (*MyFittingCatalog, MyConstants, CatalogLines,0,MyOpt_Bundle.ETGSL,MyOpt_Bundle.MyDictionary);	//Now we recompute the full catalog, this isnt necessary to complete the fit, but has to be done to score the fit
  				Sort_Catalog (*MyFittingCatalog,CatalogLines,0,0);					//Catalog is not necessarily sorted, so we sort it
				//TriplesScoreFunction (*MyFittingCatalog,ScoringParameters);
  			} 
  		} 
  	}
	printf ("%d\n",Count);
  	printf ("%e %e %e\n",MyConstants[0],MyConstants[1],MyConstants[2]);
  	printf ("%i average iterations, %i Errors\n",Iterations/(TransitionstoFit.TriplesCount[0]*TransitionstoFit.TriplesCount[1]*TransitionstoFit.TriplesCount[2]),Errors);
	free(*FitResults);
	return 1;
}

void callback (const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
//GSL function included in case it's needed for debugging
gsl_vector *f = gsl_multifit_nlinear_residual(w);
gsl_vector *x = gsl_multifit_nlinear_position(w);
double rcond;
	gsl_multifit_nlinear_rcond(&rcond, w);
	fprintf (stderr, "iter %2zu: A = %.4f, B = %.4f, C = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n", iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), 1.0 / rcond, gsl_blas_dnrm2(f));
}

int Initialize_SBFIT (struct GSL_Bundle *FitBundle, struct Opt_Bundle *MyOpt_Bundle)
{
size_t p = 3;	//These are the sizes of the parameters and data points
size_t n = MyOpt_Bundle->TransitionCount;	//The only major variation between the triples and non triples call
	FitBundle->fdf_params = gsl_multifit_nlinear_default_parameters();
 	FitBundle->fdf.f = SBFIT_OptFunc_gsl;
   	FitBundle->fdf.df = NULL;   //Finite difference Jacobian because there is no general analytic version	
   	FitBundle->fdf.fvv = NULL;	//No geodesic acceleration, early tests showed no real improvement in using it
   	FitBundle->fdf.n = n;
  	FitBundle->fdf.p = p;
 	FitBundle->fdf_params.trs = gsl_multifit_nlinear_trs_lm;
 	FitBundle->T = gsl_multifit_nlinear_trust;
 	FitBundle->Workspace = gsl_multifit_nlinear_alloc (FitBundle->T, &(FitBundle->fdf_params), n, p);
 	FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);
 	if (MyOpt_Bundle->TransitionsGSL == NULL) MyOpt_Bundle->TransitionsGSL = malloc(MyOpt_Bundle->TransitionCount*sizeof(struct Transition));
 	else MyOpt_Bundle->TransitionsGSL = realloc(MyOpt_Bundle->TransitionsGSL,MyOpt_Bundle->TransitionCount*sizeof(struct Transition));
 	return 1;
}

int Initialize_SBFIT_Alloc (struct GSL_Bundle *FitBundle, struct Opt_Bundle *MyOpt_Bundle)
{
size_t p = 3;	//These are the sizes of the parameters and data points
size_t n = MyOpt_Bundle->TransitionCount;	//The only major variation between the triples and non triples call
	printf ("Hello, %d\n",MyOpt_Bundle->TransitionCount);
	FitBundle->fdf_params = gsl_multifit_nlinear_default_parameters();
 	FitBundle->fdf.f = SBFIT_OptFunc_gsl;
   	FitBundle->fdf.df = NULL;   //Finite difference Jacobian because there is no general analytic version	
   	FitBundle->fdf.fvv = NULL;	//No geodesic acceleration, early tests showed no real improvement in using it
   	FitBundle->fdf.n = n;
  	FitBundle->fdf.p = p;
 	FitBundle->fdf_params.trs = gsl_multifit_nlinear_trs_lm;
 	FitBundle->T = gsl_multifit_nlinear_trust;
 	FitBundle->Workspace = gsl_multifit_nlinear_alloc (FitBundle->T, &(FitBundle->fdf_params), n, p);
 	FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);
 	printf ("Allocating transition memory space (%d transitions)\n",MyOpt_Bundle->TransitionCount);
 	MyOpt_Bundle->TransitionsGSL = malloc(MyOpt_Bundle->TransitionCount*sizeof(struct Transition));
 	if (MyOpt_Bundle->TransitionsGSL == NULL) printf ("Null alloc of transition space\n");
 	return 1;
}

int SBFIT (double *Guess, double *ChiSq, struct GSL_Bundle *FitBundle, struct Opt_Bundle MyOpt_Bundle, double *LineFrequencies, double FinalConstants[3])
{
int i,info,Success;
const double xtol = 1e-8;
const double gtol = 1e-8;
const double ftol = 1e-1;
const size_t p = 3;
gsl_vector *Final;
gsl_vector_view x;
  	//LineFrequencies are copied here and not in the initializer so that for cases where you want to fit the same number of lines over and over again you dont need to reinitialize the workspace 	
  	for (i=0;i<MyOpt_Bundle.TransitionCount;i++) MyOpt_Bundle.TransitionsGSL[i].Frequency = LineFrequencies[i];
  	x = gsl_vector_view_array (Guess, p);															//Set the guess	
	FitBundle->fdf.params = &MyOpt_Bundle;
	gsl_multifit_nlinear_init (&x.vector, &(FitBundle->fdf), FitBundle->Workspace);					//reInitialize the workspace incase this isnt the first run of the loop  	
	FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);								//compute initial cost function
	Success = 1;
	gsl_multifit_nlinear_driver(200, xtol, gtol, ftol, NULL, NULL, &info, FitBundle->Workspace);		//solve the system with a maximum of 50 iterations
	Final = gsl_multifit_nlinear_position(FitBundle->Workspace);									//Snag the results
	if (gsl_multifit_nlinear_niter (FitBundle->Workspace) == 200) {
		//printf ("Error: Unconverged Fit\n");
		Success = 0;
	}
	gsl_blas_ddot(FitBundle->f, FitBundle->f, ChiSq);
 	(FinalConstants)[0] = gsl_vector_get(Final,0);
 	(FinalConstants)[1] = gsl_vector_get(Final,1);
 	(FinalConstants)[2] = gsl_vector_get(Final,2);
	//printf ("A:%.4f B:%.4f C:%.4f\n",gsl_vector_get(Final, 0),gsl_vector_get(Final, 1),gsl_vector_get(Final, 2));
	return Success;
}

int Get_SBFIT_Error (struct GSL_Bundle *FitBundle, double Errors[3])
{
const size_t p = 3;
gsl_matrix *J;
gsl_matrix *covar = gsl_matrix_alloc (p, p);
	J = gsl_multifit_nlinear_jac(FitBundle->Workspace);
  	gsl_multifit_nlinear_covar (J, 0.0, covar);
	(Errors)[0] = sqrt(gsl_matrix_get(covar,0,0));
	(Errors)[1] = sqrt(gsl_matrix_get(covar,1,1));
	(Errors)[2] = sqrt(gsl_matrix_get(covar,2,2));
	gsl_matrix_free (covar);
	return 1;
}

int SBFIT_OptFunc_gsl (const gsl_vector *x, void *params, gsl_vector *f)
{
double GSLConstants[3];	//Declare some doubles to hold our constants
int i;
	//I dont love doing this, but GSL really only wants one pointer to void for the function parameters. So this is a struct to hold all the things we need for calculating frequencies
	struct Opt_Bundle *p = (struct Opt_Bundle *) params;
	GSLConstants[0] = fabs(gsl_vector_get(x, 0));	//Pull these values from the vector
	GSLConstants[1] = fabs(gsl_vector_get(x, 1));
	GSLConstants[2] = fabs(gsl_vector_get(x, 2));
	for (i=0;i<p->TransitionCount;i++) {
		gsl_vector_set (f, i, Get_Frequency(p->MyDictionary[p->TransitionsGSL[i].Upper].J,
											p->MyDictionary[p->TransitionsGSL[i].Lower].J,
											p->TransitionsGSL[i].Upper,
											p->TransitionsGSL[i].Lower,
											GSLConstants,p->ETGSL) - p->TransitionsGSL[i].Frequency);
	}
	return GSL_SUCCESS;
}

////////////////////////////////////
int Timing_Test(void)
{
/*
	Test Code - Speed test
	The code below runs the same catalog calculation over and over again as a speed test to see how fast it's running
	Current Benchmarks (Full J=25 Catalog):
		2013 Core i7 Macbook Pro - ~8 seconds for 100000 calculate+sorts
		2017 Core i7 4970 - Ubuntu ~6 seconds for 100000 calculate+sors
		Updated 7/22/19 
*/
double Constants[3], Dipoles[3];	
struct ETauStruct ETStruct;
struct Transition *BaseCatalog, *SortedCatalog;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
int CatalogTransitions,DictionaryLevels;		//Number of transitions in the catalogs
int i,j;

	
	if (!Initialize_Stuff(&(ETStruct.ETVals),&CatalogTransitions,&DictionaryLevels,&(ETStruct.Delta),&(ETStruct.StatePoints),&BaseDict,&BaseCatalog)) {
		goto Error;
	}
	Constants[0] = 3000.0;
	Constants[1] = 2000.0;
	Constants[2] = 1000.0;
	
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
			);
	Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);


	int Loops = 100000;		//Number of loops in a single run
	int TimingLoops = 10;	//Number if timing runs, used to capture variance in the run time
	double *Timing;
	Timing = malloc (TimingLoops*sizeof(double));

	
	//Asymmetric Top
	Constants[0] = 2000.0;
	Constants[1] = 1510.0;
	Constants[2] = 1000.0;
	Dipoles[0] = 1.0;
	Dipoles[1] = 1.0;
	Dipoles[2] = 1.0;
	
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);

	
	
	SortedCatalog = malloc (CatalogTransitions*sizeof(struct Transition));
	for (i=0;i<CatalogTransitions;i++) {
		SortedCatalog[i] = BaseCatalog[i];
	}
	
	Get_Catalog (	SortedCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Calculate_State_Energies (BaseDict, SortedCatalog, CatalogTransitions,DictionaryLevels,0);
	Calculate_Intensities (BaseCatalog, CatalogTransitions, BaseDict, 3.0, Dipoles,0);	

	
	printf ("Starting timing test run. %d loops per run, %d runs\n",Loops,TimingLoops);
	for (j=0;j<TimingLoops;j++) {
		clock_t begin = clock();
		for (i=0;i<Loops;i++) {
			Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
							Constants, 			//Rotational constants for the calculation
							CatalogTransitions,	//# of transitions in the catalog
							0,					//Verbose
							ETStruct,
							BaseDict
			);
			Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);
			Calculate_State_Energies (BaseDict, BaseCatalog, CatalogTransitions, DictionaryLevels,0);
		}
		clock_t end = clock();
		Timing[j] = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("%f\n", Timing[j]);
		
	}
	double Sum,Mean,StdDev;
	Sum = 0.0;
	StdDev = 0.0;
	for(i=0; i<TimingLoops; ++i)
    {
        Sum += Timing[i];
    }
    Mean = Sum/TimingLoops;
    for(i=0; i<TimingLoops; ++i)
        StdDev += pow(Timing[i] - Mean, 2.0);
    StdDev = sqrt(StdDev/TimingLoops);
	printf ("Average Time:%f, Standard Deviation:%f\n",Mean,StdDev);	
	free(Timing);
	return 1;
Error:
	printf ("Error: Couldn't initialize the timing test\n");	
	return 0;
}

int Accuracy_Test(void)
{
/*
	Test Code - Accuracy test
	The code below builds three different catalogs at the oblate, prolate, and asymmetric limits and saves them. 
	Please note that the filenames are not dynamic and you should change them depending on the resolution of the solver used.
	
*/
double Constants[3];	
struct ETauStruct ETStruct;
struct Transition *BaseCatalog;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
int i,CatalogTransitions,DictionaryLevels;		//Number of transitions in the catalogs
FILE *FileHandle;	
	
	if (!Initialize_Stuff(&(ETStruct.ETVals),&CatalogTransitions,&DictionaryLevels,&(ETStruct.Delta),&(ETStruct.StatePoints),&BaseDict,&BaseCatalog)) {
		goto Error;
	}
	//Oblate Top
	Constants[0] = 2000.0;
	Constants[1] = 1990.0;
	Constants[2] = 1000.0;
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);
	FileHandle = fopen("OblateCat_dK2.txt","w");
	for (i=0;i<CatalogTransitions;i++) {
		struct Transition TransitionToPrint = BaseCatalog[i];
		fprintf (FileHandle,"%.12f %i %i %i  %i %i %i\n",TransitionToPrint.Frequency, BaseDict[TransitionToPrint.Upper].J,BaseDict[TransitionToPrint.Upper].Ka,BaseDict[TransitionToPrint.Upper].Kc,BaseDict[TransitionToPrint.Lower].J,BaseDict[TransitionToPrint.Lower].Ka,BaseDict[TransitionToPrint.Lower].Kc);
	}
	fclose(FileHandle);
	
	//Prolate Top
	Constants[0] = 2000.0;
	Constants[1] = 1010.0;
	Constants[2] = 1000.0;
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);
	FileHandle = fopen("ProlateCat_dK2.txt","w");
	for (i=0;i<CatalogTransitions;i++) {
		struct Transition TransitionToPrint = BaseCatalog[i];
		fprintf (FileHandle,"%.12f %i %i %i  %i %i %i\n",TransitionToPrint.Frequency, BaseDict[TransitionToPrint.Upper].J,BaseDict[TransitionToPrint.Upper].Ka,BaseDict[TransitionToPrint.Upper].Kc,BaseDict[TransitionToPrint.Lower].J,BaseDict[TransitionToPrint.Lower].Ka,BaseDict[TransitionToPrint.Lower].Kc);
	}
	fclose(FileHandle);

	//Asymmetric Top
	Constants[0] = 2000.0;
	Constants[1] = 1510.0;
	Constants[2] = 1000.0;
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);
	FileHandle = fopen("AsymmetricCat_dK2.txt","w");
	for (i=0;i<CatalogTransitions;i++) {
		struct Transition TransitionToPrint = BaseCatalog[i];
		fprintf (FileHandle,"%.12f %i %i %i  %i %i %i\n",TransitionToPrint.Frequency, BaseDict[TransitionToPrint.Upper].J,BaseDict[TransitionToPrint.Upper].Ka,BaseDict[TransitionToPrint.Upper].Kc,BaseDict[TransitionToPrint.Lower].J,BaseDict[TransitionToPrint.Lower].Ka,BaseDict[TransitionToPrint.Lower].Kc);
	}
	fclose(FileHandle);
	return 1;
	
Error:
	printf ("Error: Couldn't initialize the accuracy test\n");	
	return 0;
}

void Test_SBFIT (void)
{
double GuessConstants[3], RealConstants[3], Dipoles[3],FitConstants[3], ChiSqr;
double *FittingFrequencies;
int i,CatalogTransitions,DictionaryLevels,JRestrictedLines,IntRestrictedLines,FrequencyCount;
struct ETauStruct ETStruct;
struct Transition *BaseCatalog, *JRestricted, *IntRestricted;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
struct GSL_Bundle TestGSLBundle;
struct Opt_Bundle TestOptBundle;
	RealConstants[0] = 3333;
	RealConstants[1] = 2222;
	RealConstants[2] = 1111;
	
	GuessConstants[0] = 3000;
	GuessConstants[1] = 2000;
	GuessConstants[2] = 1000;
	
	Dipoles[0] = 1.0;
	Dipoles[1] = 1.0;
	Dipoles[2] = 1.0;
	
	FrequencyCount = 10;
	Initialize_Stuff(&(ETStruct.ETVals),&CatalogTransitions,&DictionaryLevels,&(ETStruct.Delta),&(ETStruct.StatePoints),&BaseDict,&BaseCatalog);
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					RealConstants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Calculate_State_Energies (BaseDict, BaseCatalog, CatalogTransitions, DictionaryLevels,0);
	Calculate_Intensities (BaseCatalog, CatalogTransitions, BaseDict, 2.0, Dipoles,0);	
	JRestrictedLines = Fill_Catalog_Restricted_J (BaseCatalog, &JRestricted, RealConstants, CatalogTransitions, 0, 10, 0, BaseDict);
	IntRestrictedLines = Fill_Catalog_Restricted_Intensity_Count (JRestricted, &IntRestricted, RealConstants, JRestrictedLines, FrequencyCount, 0, BaseDict);
	FittingFrequencies = malloc(FrequencyCount*sizeof(double));
	TestOptBundle.ETGSL = ETStruct;
	TestOptBundle.MyDictionary = BaseDict;
	TestOptBundle.TransitionCount = FrequencyCount;
	Initialize_SBFIT (&TestGSLBundle, &TestOptBundle);
	for (i=0;i<FrequencyCount;i++) {
		FittingFrequencies[i] = IntRestricted[i].Frequency;
		TestOptBundle.TransitionsGSL[i] = IntRestricted[i];
	}
  	RealConstants[0] = 3343;
  	SBFIT (GuessConstants, &ChiSqr, &TestGSLBundle, TestOptBundle, FittingFrequencies, FitConstants);
  	printf ("Fit to constants %f %f %f with a chi squared of %f\n", FitConstants[0],FitConstants[1],FitConstants[2],ChiSqr);
	free(FittingFrequencies);
}

void Test_Triples (char *FileName, struct Transition *FittingCatalog, struct Level *FittingDictionary, int CatalogLines, struct ETauStruct FittingETStruct)
{
double *ExpX, *ExpY, *PeakList,GuessConstants[3],*Results;
int ExperimentalPoints,PeakCount;
struct Triple TestTriple;
struct GSL_Bundle TestGSLBundle;
struct Opt_Bundle TestOptBundle;
ScoreFunction TestFunction;
	
	ExperimentalPoints = Load_Exp_File  (FileName, &ExpX, &ExpY, 1);
	PeakCount = Peak_Find (&PeakList, 100.0, 0.003, ExpX, ExpY, ExperimentalPoints,0);
	printf ("Found %d experimental peaks\n", PeakCount);
	TestTriple.TransitionList[0] = FittingCatalog[400];
	TestTriple.TransitionList[1] = FittingCatalog[500];
	TestTriple.TransitionList[2] = FittingCatalog[450];
	Find_Triples (&TestTriple, PeakList, 100.0, PeakCount,0);
	printf ("Found %d %d %d experimental lines for the proposed triple\n",TestTriple.TriplesCount[0],TestTriple.TriplesCount[1],TestTriple.TriplesCount[2]);
	TestOptBundle.ETGSL = FittingETStruct;
	TestOptBundle.MyDictionary = FittingDictionary;
	TestOptBundle.TransitionsGSL[0] = TestTriple.TransitionList[0];
	TestOptBundle.TransitionsGSL[1] = TestTriple.TransitionList[1];
	TestOptBundle.TransitionsGSL[2] = TestTriple.TransitionList[2];
	
	
	Initialize_Triples_Fitter (&TestGSLBundle,&TestOptBundle);
	GuessConstants[0] = 5002.0;
	GuessConstants[1] = 3004.0;
	GuessConstants[2] = 1998.0;
	TestFunction = &DummyFunction;
	Fit_Triples_Bundle (	TestTriple, 
							GuessConstants, 
							&Results, 
							&FittingCatalog, 
							CatalogLines, 
							&TestGSLBundle, 
							TestOptBundle, 
							TestFunction, 
							0 
						);

}

double DummyFunction (struct Transition *MyCatalog, void *Data)
{
	
	return 1.0;
}

int Timing_Test_Triples (void)
{
/*
	Test Code - Speed test
	The code below runs the same triples bundle calculation over and over again as a speed test to see how fast it's running
	Current Benchmarks (Full J=25 Catalog):

		Updated 7/22/19 
*/
double Constants[3], Dipoles[3], GuessConstants[3],*ExpX, *ExpY, *PeakList, *Results;	
struct ETauStruct ETStruct;
struct Transition *BaseCatalog, *SortedCatalog;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
int CatalogTransitions,DictionaryLevels;		//Number of transitions in the catalogs
int i,j,ExperimentalPoints,PeakCount;
struct Opt_Bundle TestOptBundle;
struct GSL_Bundle TestGSLBundle;
struct Triple TestTriple;
ScoreFunction TestFunction;
	
	if (!Initialize_Stuff(&(ETStruct.ETVals),&CatalogTransitions,&DictionaryLevels,&(ETStruct.Delta),&(ETStruct.StatePoints),&BaseDict,&BaseCatalog)) {
		goto Error;
	}
	Constants[0] = 3000.0;
	Constants[1] = 2000.0;
	Constants[2] = 1000.0;

	
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
			);
	Sort_Catalog (BaseCatalog,CatalogTransitions,2,0);
	
	TestOptBundle.ETGSL = ETStruct;
	TestOptBundle.MyDictionary = BaseDict;
	TestOptBundle.TransitionCount = 3;

	int Loops = 3;		//Number of loops in a single run
	int TimingLoops = 10;	//Number if timing runs, used to capture variance in the run time
	double *Timing;
	Timing = malloc (TimingLoops*sizeof(double));

	
	//Asymmetric Top
	Constants[0] = 2000.0;
	Constants[1] = 1510.0;
	Constants[2] = 1000.0;
	Dipoles[0] = 1.0;
	Dipoles[1] = 1.0;
	Dipoles[2] = 1.0;
	
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	TestTriple.TransitionList[0] = BaseCatalog[400];
	TestTriple.TransitionList[1] = BaseCatalog[500];
	TestTriple.TransitionList[2] = BaseCatalog[450];
	
	
	SortedCatalog = malloc (CatalogTransitions*sizeof(struct Transition));
	for (i=0;i<CatalogTransitions;i++) {
		SortedCatalog[i] = BaseCatalog[i];
	}
	
	Get_Catalog (	SortedCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					ETStruct,
					BaseDict
	);
	Calculate_State_Energies (BaseDict, SortedCatalog, CatalogTransitions,DictionaryLevels,0);
	Calculate_Intensities (BaseCatalog, CatalogTransitions, BaseDict, 3.0, Dipoles,0);	

	ExperimentalPoints = Load_Exp_File  ("../tests/ft2494.txt", &ExpX, &ExpY, 1);
	PeakCount = Peak_Find (&PeakList, 100.0, 0.003, ExpX, ExpY, ExperimentalPoints,0);
	Find_Triples (&TestTriple, PeakList, 100.0, PeakCount,1);

	TestOptBundle.TransitionsGSL = NULL;
	Initialize_Triples_Fitter (&TestGSLBundle,&TestOptBundle);
	TestOptBundle.TransitionsGSL[0] = TestTriple.TransitionList[0];
	TestOptBundle.TransitionsGSL[1] = TestTriple.TransitionList[1];
	TestOptBundle.TransitionsGSL[2] = TestTriple.TransitionList[2];

	GuessConstants[0] = 5002.0;
	GuessConstants[1] = 3004.0;
	GuessConstants[2] = 1998.0;
	TestOptBundle.ETGSL = ETStruct;
	TestFunction = &DummyFunction;
	
	print_Transition (BaseCatalog[400], BaseDict);
	print_Transition (BaseCatalog[500], BaseDict);
	print_Transition (BaseCatalog[450], BaseDict);

	printf ("Starting timing test run. %d loops per run, %d runs\n",Loops,TimingLoops);
	for (j=0;j<TimingLoops;j++) {
		clock_t begin = clock();
		for (i=0;i<Loops;i++) {
			Fit_Triples_Bundle (	TestTriple, 
									GuessConstants, 
									&Results, 
									&BaseCatalog, 
									CatalogTransitions, 
									&TestGSLBundle, 
									TestOptBundle, 
									TestFunction, 
									0 
								);	
		}
		free(Results);
		clock_t end = clock();
		Timing[j] = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("%f\n", Timing[j]);
		
	}
	double Sum,Mean,StdDev;
	Sum = 0.0;
	StdDev = 0.0;
	for(i=0; i<TimingLoops; ++i)
    {
        Sum += Timing[i];
    }
    Mean = Sum/TimingLoops;
    for(i=0; i<TimingLoops; ++i)
        StdDev += pow(Timing[i] - Mean, 2.0);
    StdDev = sqrt(StdDev/TimingLoops);
	printf ("Average Time:%f, Standard Deviation:%f\n",Mean,StdDev);	
	free(Timing);
	return 1;
Error:
	printf ("Error: Couldn't initialize the timing test\n");	
	return 0;
}

////////////////////////////////////
int Search_DR_Hits (int DRPairs, double ConstStart, double ConstStop, double Step, double *DRFrequency, double Tolerance, int ExtraLineCount, double *ExtraLines, int **DRLinks, int LinkCount, struct Transition *CatalogtoFill, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary, char *FileName)
{
/*Function to find rotational constants through an arbitrary set of DR links
Inputs:
DRPairs - The number of DR transitions to be searched
ConstStart (MHz) - Value to start searching in A, B, and C
ConstStop (MHz)  - Value to stop searching in A, B, and C 
Step (MHz) - Step to take in A, B, and C
DRFrequency (MHz) - List of DR frequencies to match, should have a length equal to DRPairs
Tolerance (MHz) - How much error to allow between an experimental frequency and a predicted catalog line to consider it a match. Since it's absolute value match condition is +/-Tolerance
ExtraLineCount -  User can feed extra experimental lines as a constraint even if not part of the DR match set. If given we do additional scoring based on number of matches to this 
ExtraLines (MHz) - List of the extra experimental lines. Should be the same length as ExtraLineCount
DRLinks - 2D array, listing the linkages between states. Array shape should be DRLinks[LinkCount][2] so that DRLink[x][0] is the first transition and DRLink[x][1] is the second, 0 and 1 should be integers between 0 and DRPairs that tell the program which two transitions it's matching
LinkCount -  Number of DR linkages to deal with
CatalogtoFill - Catalog used for matching
CatLines - Number of lines in the catalog
Verbose - Flag for printing more verbose inforamtion from the function
ETStruct - Eigenvalue struct passed so we can do fitting in the function
MyDictionary - Catalog dictionary passed for fitting and printing


Todo:
-Test Sorted searches

*/
double CurrentA,CurrentB,CurrentC,Count,ChiSqr;
double Constants[3],FitConstants[3];
int *Match,**MatchArrays,i,j,k,DRMatch,AllLinks,Wins,LocalLink,MatchLimit,MatchCount,StartJ,StartK;
int ***MatchRecord; //Record all of our matches in one place || MatchRecord[Match][Link][Upper/Lower]
struct GSL_Bundle MyGSLBundle;
struct Opt_Bundle MyOptBundle;
FILE *FileHandle;

	MatchLimit = 100;
	MatchRecord = malloc (MatchLimit*sizeof(int **));
	for (i=0;i<MatchLimit;i++) {
		MatchRecord[i] = malloc (LinkCount*sizeof(int *));
		for (j=0;j<LinkCount;j++) MatchRecord[i][j] = malloc(2*sizeof(int));
	}

	//Error checking of inputs, shouldnt be an issue but still
	if (ConstStart < 0.0) {
		if (Verbose) printf ("Starting constants set too low, defaulting to 1GHz\n");
		ConstStart = 1000.0;
	}
	if (ConstStop < 0.0) {
		if (Verbose) printf ("Stopping constants set too low, defaulting to 10GHz\n");
		ConstStart = 10000.0;
	}
	if (Step < 0.0) {
		if (Verbose) printf ("Step set too low, defaulting to 10MHz\n");
		ConstStart = 10.0;
	}
	if (Tolerance < 0.0) {
		if (Verbose) printf ("Tolerance set too low, defaulting to 20MHz\n");
		ConstStart = 20.0;
	}
	FileHandle = fopen(FileName,"a");
	if (FileHandle == NULL) goto Error;

	//Initialize the GSL fitter
	if (DRPairs >= 3) {
		MyOptBundle.ETGSL = ETStruct;
		MyOptBundle.MyDictionary = MyDictionary;
		MyOptBundle.TransitionsGSL = NULL;
		MyOptBundle.TransitionCount = DRPairs;
		Initialize_SBFIT (&MyGSLBundle, &MyOptBundle);	
	}

	//Verbose startup 
	if (Verbose) {
		printf ("Grid Search from %.2f MHz to %.2f MHz in %.2f MHz steps, %.2e total catalogs\n", ConstStart,ConstStop,Step, (double) pow((ConstStop-ConstStart)/Step,3.0)/6.0);
		if ((DRPairs <=3) && (ExtraLineCount < 1)) {
			printf ("Insufficient information given, either give extra lines to score against or more DR links\n");
			return 0;
		} else {
			if (DRPairs <= 3) printf ("Not enough DR pairs for fitting, will match only against extra lines\n");
			if (ExtraLineCount < 1) printf ("No extra lines supplied, can only score by fitting\n");
		}
	}

	fprintf (FileHandle, "Constants Start :%f Constants Stop :%f Constants Step :%f Tolerance %f\n",ConstStart,ConstStop,Step,Tolerance);
	fprintf (FileHandle, "%d DR lines supplied, %d links between them, %d extra lines supplied",DRPairs,ExtraLineCount,LinkCount);
	
	
	Match = malloc(DRPairs*sizeof(int));	//Array of yes/no to track if each of the DR frequencies has a matching catalog transition or not
	if (Match == NULL) goto Error;	//Basic but overkill error checking
	MatchArrays = malloc (DRPairs*sizeof(int *));	//Array to store all the potential matches for each DR frequency
	for (i=0;i<DRPairs;i++) MatchArrays[i] = malloc(100*sizeof(Match)); //Each of these arrays stores the matches for each DR line, can currently match up to 100 catalog lines per DR frequency. Hard coded limit for now cause over 100 is a lot
	CurrentA = ConstStart;
	Count = 0.0;	//Tracking the number of counts we perform, using doubles to prevent int overflow
	clock_t start = clock();
	while (CurrentA < ConstStop) {
		CurrentB = ConstStart;
		while (CurrentB < ConstStop) {
			CurrentC = ConstStart;
			while (CurrentC < ConstStop) {
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {
					Count+= 1.0;
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					//Predict and check spectrum
					Get_Catalog (	CatalogtoFill, //Catalog to compute frequencies for
									Constants, //Rotational constants for the calculation
									CatLines,	//# of transitions in the catalog
									0,	//Verbose
									ETStruct,
									MyDictionary
					);
					for (i=0;i<DRPairs;i++) Match[i] = 0;	//Reset our matches to 0
					//Brute force check to see if our DR lines exist in the catalog, may be faster to check after a sort, need to try
					for (i=0;i<CatLines;i++) {
						for (j=0;j<DRPairs;j++) {
							if (fabs(CatalogtoFill[i].Frequency-DRFrequency[j]) < Tolerance) {
								MatchArrays[j][Match[j]] = i;
								Match[j]++;
							}
						}
					}
					DRMatch = 1;
					for (i=0;i<DRPairs;i++) if (Match[i] <= 0) DRMatch = 0;	//If any of our transitions werent matched we bail  
					for (i=0;i<MatchLimit;i++) {
						for (j=0;j<LinkCount;j++) {
							MatchRecord[i][j][0] = 0;
							MatchRecord[i][j][1] = 0;
						}
					}
					MatchCount = 0;
					
					//Code in test
					AllLinks = 0;
					if (DRMatch) { 
						i=0;
						//Assuming >2 links, need to add an error check for this before we start
						while (i<LinkCount) {	
							//Iterate through all matched levels for this linkage
							LocalLink = 0;	//Start by assuming there is no link between lvls
							if (MatchCount > 0) {
								StartJ = MatchRecord[MatchCount-1][i][0];
								StartK = MatchRecord[MatchCount-1][i][1];
							}else {
								StartJ = 0;
								StartK = 0;
							}
							for (j=StartJ;j<Match[DRLinks[i][0]];j++) {
								for (k=StartK;k<Match[DRLinks[i][1]];k++) {
									//If any link works we count it as a win
									if (Match_Levels(MatchArrays[DRLinks[i][0]][j], MatchArrays[DRLinks[i][1]][k], CatalogtoFill) ) {
										LocalLink = 1;
										MatchRecord[MatchCount][i][0] = MatchArrays[DRLinks[i][0]][j];
										MatchRecord[MatchCount][i][1] = MatchArrays[DRLinks[i][1]][k];
										break;
									} 
								}
							}
							//If a match was found within the set, we proceed, if not 
							if (LocalLink) {
								i++;
								if (i == LinkCount) {
									AllLinks = 1;
									MatchCount++;
								}
							} else break;
						}
					}
					if (AllLinks) {
						if (DRPairs >=3) {
							for (i=0;i<LinkCount;i++) {
								MyOptBundle.TransitionsGSL[DRLinks[i][0]] = CatalogtoFill[MatchRecord[0][i][0]];
								MyOptBundle.TransitionsGSL[DRLinks[i][1]] = CatalogtoFill[MatchRecord[0][i][1]];
							}
							SBFIT (Constants, &ChiSqr, &MyGSLBundle, MyOptBundle, DRFrequency, FitConstants);
							if ((ChiSqr/DRPairs) < 0.1) {
								Get_Catalog (	CatalogtoFill, //Catalog to compute frequencies for
												Constants, //Rotational constants for the calculation
												CatLines,	//# of transitions in the catalog
												0,	//Verbose
												ETStruct,
												MyDictionary
												);
								Wins = 0;
								for (i=0;i<ExtraLineCount;i++) {
									for (j=0;j<CatLines;j++) {
										if (fabs(CatalogtoFill[j].Frequency-ExtraLines[i]) < (Tolerance/100.0)) {
											Wins++;
										}
									}
								}
								if (ExtraLineCount) {
									if (Wins > 3) {
										fprintf (FileHandle, "Fitted constants A:%f B:%f C:%f ChiSqr:%f\n",FitConstants[0],FitConstants[1],FitConstants[2],ChiSqr);
									}
								} else {
									fprintf (FileHandle, "Fitted constants A:%f B:%f C:%f ChiSqr:%f\n",FitConstants[0],FitConstants[1],FitConstants[2],ChiSqr);
								}	
							}
						} else {
							//Working with only 2 DR links
					
						}
					}
				}
				CurrentC += Step;
			}
			CurrentB += Step;
		}
		CurrentA += Step;
	}
	clock_t end = clock();
	double Timing = (double)(end - start) / CLOCKS_PER_SEC;
	printf ("%e individual fits performed in %.2fs\n",Count,Timing);
	free(Match);
	for (i=-0;i<DRPairs;i++) free(MatchArrays[i]);
	free(MatchArrays);
	fclose(FileHandle);
	return 1;
Error:
	printf("Error running matching program");
	return 0;
}
 
int Match_Levels (int Match1, int Match2, struct Transition *MatchCatalog) 
{
	if (((MatchCatalog)[Match1].Upper == (MatchCatalog)[Match2].Upper) || ((MatchCatalog)[Match1].Lower == (MatchCatalog)[Match2].Lower) || ((MatchCatalog)[Match1].Upper == (MatchCatalog)[Match2].Lower) || ((MatchCatalog)[Match1].Lower == (MatchCatalog)[Match2].Upper)) {
		return 1;
	} else {
		return 0;
	}
} 

#endif /* __FITTER_H__ */