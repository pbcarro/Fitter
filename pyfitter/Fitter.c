/*
Built by Brandon Carroll 4/7/19 All Rights Reserved

Generic version of the fitting program

Build Command:
gcc -Wall -o Go Fitter.c -lm -lgsl -lgslcblas -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Fitter.so -shared -fPIC -O3 -funroll-loops Fitter.c -lm -lgsl -lgslcblas
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_linalg.h>

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
	struct Transition TransitionsGSL[3];
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
int Initialize_Stuff (double **/*ETArray*/, int */*CatTransitions*/, int */*DictTransitions*/, double */*FileDelta*/, int */*StatePoints*/, struct Level **/*DictionaryIn*/, struct Transition **/*CatalogIn*/);
int Load_ETau_File (char */*FileName*/, double **/*X*/, double */*FileDelta*/, int */*StatePoints*/, int */*StateCount*/);
int Load_ETau_File2 (char */*FileName*/, struct ETauStruct */*StructToLoad*/, int */*StateCount*/, int /*Verbose*/);
int Load_Base_Catalog (char */*FileName*/, struct Transition **/*BaseCatalog*/,  int /*Verbose*/);
int Load_Base_Catalog_Dictionary (char */*FileName*/, struct Level **/*DictIn*/,  int /*Verbose*/);
int Load_Exp_File  (char */*FileName*/, double **/*X*/, double **/*Y*/, int /*Verbose*/);

//Frequency predicting functions
double Get_Kappa (double /*A*/, double /*B*/, double /*C*/);  
int Get_J (int /*TransitionIndex*/, struct Level */*MyDictionary*/);
double Partition_Function (double */*Constants*/, double /*Temperature*/);
double E_tau (int /*TransitionIndex*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
double Rigid_Rotor (double /*A*/, double /*C*/, int /*J*/, int /*Index*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
double Get_Frequency (int /*J_Up*/, int /*J_Low*/, int /*IndexUp*/, int /*IndexLow*/, double */*Constants*/, struct ETauStruct /*ETStruct*/);
int Get_Catalog (struct Transition */*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level */*MyDictionary*/);

//General catalog functions
void print_Transition (struct Transition /*TransitionToPrint*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_Frequency (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*FrequencyLow*/, double /*FrequencyHigh*/, int /*Verbose*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_J (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*JMin*/, int /*JMax*/, int /*Verbose*/, struct Level */*MyDictionary*/);
int Fill_Catalog_Restricted_Intensity (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*Percentile*/, int /*Verbose*/, struct Level */*MyDictionary*/);
void Sort_Catalog (struct Transition */*Catalog*/, int /*TransitionCount*/, int /*SortMethod*/, int /*SortType*/);
int Catalog_Comparator_Frequency (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Intensity (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Upper (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Lower (const void */*a*/, const void */*b*/);
void insertionSort(struct Transition */*CatalogtoSort*/, int /*TransitionCount*/);
void Calculate_State_Energies (struct Level */*MyCatalog*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/);
void Calculate_Intensities (struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/, struct Level */*MyDictionary*/, double /*T*/, double */*Dipoles*/);

//Triples fitting functions
int Peak_Find (double **/*LineList*/, double /*Max*/, double /*Min*/, double */*X*/, double */*Y*/, int /*ArraySize*/, int /*Verbose*/);
int Find_Triples (struct Triple */*TripletoFit*/, double */*LineFrequencies*/, double /*Window*/, int /*LineCount*/);
double Fit_Lines (double */*Guess*/, int /*Verbose*/, struct Opt_Bundle /*GSLOptBundle*/);
void Initialize_Triples_Fitter (struct GSL_Bundle */*FitBundle*/);
int OptFunc_gsl (const gsl_vector */*x*/, void */*params*/, gsl_vector */*f*/);
int Fit_Triples (struct Triple /*TransitionstoFit*/, double */*Guess*/, double **/*FitResults*/, struct Transition **/*Catalog*/, int /*CatalogLines*/, double */*ExperimentalLines*/, int /*ExperimentalLineCount*/);
int Fit_Triples_Bundle (struct Triple /*TransitionstoFit*/, double */*Guess*/, double **/*FitResults*/, struct Transition **/*Catalog*/, int /*CatalogLines*/, struct GSL_Bundle */*FitBundle*/, struct Opt_Bundle /*MyOpt_Bundle*/, ScoreFunction /*TriplesScoreFunction*/, void */*ScoringParameters*/);
void callback (const size_t /*iter*/, void */*params*/, const gsl_multifit_nlinear_workspace */*w*/);

//New Functions
int Search_DR_Hits (int /*DRPairs*/, double /*ConstStart*/, double /*ConstStop*/, double /*Step*/, double */*DRFrequency*/, double /*Tolerance*/, int /*ExtraLineCount*/, double */*ExtraLines*/, int **/*DRLinks*/, int /*LinkCount*/, struct Transition */*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level */*MyDictionary*/);

void Test_Triples (char *, struct Transition *, struct Level *, int /*CatalogLines*/, struct ETauStruct /*FittingETStruct*/);
double DummyFunction (struct Transition */*MyCatalog*/, void */*Data*/);

int main (int argc, char *argv[])  
{
double Constants[3], Dipoles[3];	
struct ETauStruct ETStruct;
struct Transition *BaseCatalog, *SortedCatalog;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
int CatalogTransitions,DictionaryLevels;		//Number of transitions in the catalogs
int i,j;
double *EnergyLevels,*IntensityVals;

	
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
	
	Test_Triples ("../tests/ft2494.txt",BaseCatalog,BaseDict,CatalogTransitions,ETStruct);
	
	/*
	Test Code - Accuracy test
	The code below builds three different catalogs at the oblate, prolate, and asymmetric limits and saves them. 
	Please note that the filenames are not dynamic and you should change them depending on the resolution of the solver used.
	
	*/
	
	/*
	FILE *FileHandle;	
	
	//Oblate Top
	Constants[0] = 2000.0;
	Constants[1] = 1990.0;
	Constants[2] = 1000.0;
	Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
					Constants, 			//Rotational constants for the calculation
					CatalogTransitions,	//# of transitions in the catalog
					0,					//Verbose
					TestStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,1);
	FileHandle = fopen("OblateCat_dK2.txt","w");
	int i;
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
					TestStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,1);
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
					TestStruct,
					BaseDict
	);
	Sort_Catalog (BaseCatalog,CatalogTransitions,1);
	FileHandle = fopen("AsymmetricCat_dK2.txt","w");
	for (i=0;i<CatalogTransitions;i++) {
		struct Transition TransitionToPrint = BaseCatalog[i];
		fprintf (FileHandle,"%.12f %i %i %i  %i %i %i\n",TransitionToPrint.Frequency, BaseDict[TransitionToPrint.Upper].J,BaseDict[TransitionToPrint.Upper].Ka,BaseDict[TransitionToPrint.Upper].Kc,BaseDict[TransitionToPrint.Lower].J,BaseDict[TransitionToPrint.Lower].Ka,BaseDict[TransitionToPrint.Lower].Kc);
	}
	fclose(FileHandle);
	*/
		
	/*
	Test Code - Speed test
	The code below runs the same catalog calculation over and over again as a speed test to see how fast it's running
	Current Benchmarks (Full J=25 Catalog):
		2013 Core i7 Macbook Pro - ~8 seconds for 100000 calculate+sorts
		2017 Core i7 4970 - Ubuntu ~6 seconds for 100000 calculate+sors
		Updated 7/22/19 
	*/
	int Loops = 100000;		//Number of loops in a single run
	int TimingLoops = 10;	//Number if timing runs, used to capture variance in the run time
	double *Timing;
	Timing = malloc (TimingLoops*sizeof(double));
	
	EnergyLevels = malloc(DictionaryLevels*sizeof(double));
	IntensityVals = malloc(CatalogTransitions*sizeof(double));
	
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
	qsort(SortedCatalog, CatalogTransitions, sizeof(struct Transition), Catalog_Comparator_Index_Upper);
	qsort(SortedCatalog, CatalogTransitions, sizeof(struct Transition), Catalog_Comparator_Index_Lower);
	Calculate_State_Energies (BaseDict, SortedCatalog, CatalogTransitions);
	Calculate_Intensities (BaseCatalog, CatalogTransitions, BaseDict, 3.0, Dipoles);	

	
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
	printf ("Error Initializing program, exiting\n");
	return 1;
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
	Find_Triples (&TestTriple, PeakList, 100.0, PeakCount);
	printf ("Found %d %d %d experimental lines for the proposed triple\n",TestTriple.TriplesCount[0],TestTriple.TriplesCount[1],TestTriple.TriplesCount[2]);
	Initialize_Triples_Fitter (&TestGSLBundle);
	GuessConstants[0] = 3002.0;
	GuessConstants[1] = 2004.0;
	GuessConstants[2] = 998.0;
	TestOptBundle.ETGSL = FittingETStruct;
	TestOptBundle.MyDictionary = FittingDictionary;
	TestOptBundle.TransitionsGSL[0] = TestTriple.TransitionList[0];
	TestOptBundle.TransitionsGSL[1] = TestTriple.TransitionList[1];
	TestOptBundle.TransitionsGSL[2] = TestTriple.TransitionList[2];
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

int Initialize_Stuff (double **ETArray, int *CatTransitions, int *DictLevels, double *FileDelta, int *StatePoints, struct Level **DictionaryIn, struct Transition **CatalogIn) 
{
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
	if (Load_ETau_File ("J0_25_dk3.dat", ETArray,FileDelta,StatePoints,&TempPoints)) printf ("Tables Loaded\n");
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
//This loads the E_Tau file in by reading a single value from the file until there are no more values found in the file
//This unwraps the 2D file into a single long array for contiguousness, files should have a state's ET values in a single ROW, not column
//Figuring out where each state starts and stops is done elsewhere, so the global variables need to be set properly to deal with this
//ET is preallocated based on the global variables, there are no checks to see if the file matches until later
int i,StateLimit;
FILE *FileHandle;
char * line = NULL;
size_t len = 0;
	
	StateLimit = 10000;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;	
	*StateCount = 0;
	while (getline(&line, &len, FileHandle) > 0) (*StateCount)++;
	rewind(FileHandle);
	*X = malloc(StateLimit*(*StateCount)*sizeof(double));																
	if (*X == NULL) goto Error;
	i = 0;
	while (fscanf (FileHandle, "%lf", &(*X)[i]) == 1) {					//Keep scanning in floating points until there are no more
		i++;
		if (i == (StateLimit*(*StateCount))) {
			printf ("Error: Your state file exceeds J=100, I don' want to deal with that\n");
			goto Error;
		}
	}
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont mater														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	*X = realloc(*X,i*sizeof(double));
	*StatePoints = (int) i/(*StateCount);
	*FileDelta = 2.0/(*StatePoints-1);
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
char * line = NULL;
size_t len = 0;
	
	StateLimit = 10000;
	FileHandle = NULL;
	FileHandle = fopen (FileName, "r");									//Open file read only
	if (FileHandle == NULL) goto Error;	
	*StateCount = 0;
	while (getline(&line, &len, FileHandle) > 0) (*StateCount)++;	//Run through the file and keep going until we hit the end, track the number of lines/states
	rewind(FileHandle);	//Rewind to the start of the file
	(*StructToLoad).ETVals = malloc(StateLimit*(*StateCount)*sizeof(double));																
	if ((*StructToLoad).ETVals == NULL) goto Error;
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
	if (Verbose) {
		printf ("=========Verbose Load_ETau_File2=========\n");
		printf("Loaded ET File %s with %d states, %d points per state or a delta kappa of %.2e\n",FileName,(*StateCount),(*StructToLoad).StatePoints,(*StructToLoad).Delta);
		printf ("=======================================\n");
	}
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
	return i;
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

double Get_Kappa (double A, double B, double C) 
{
	//Computes Ray's asymmetry parameter
	return ((2.0*B-A-C)/(A-C));		
}

int Get_J (int TransitionIndex, struct Level *MyDictionary)
{
	return MyDictionary[TransitionIndex].J; 
}

double Partition_Function (double *Constants, double Temperature) 
{
	return (5.34E+6*sqrt(pow(Temperature,3.0)/(Constants[0]*Constants[1]*Constants[2])));
}

double E_tau (int TransitionIndex, double Kappa, struct ETauStruct ETStruct) 
{
//A function to return the E_tau() value of a rigid rotor Hamiltonian
//Really this just fetches a value for the appropriate value of J/Ka/Kc
//Values are explicitly calculated when E_tau has an analytic form, for all other values we use a look up table calculated elsewhere
int Index;
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

double Get_Frequency (int J_Up, int J_Low, int IndexUp, int IndexLow, double *Constants, struct ETauStruct ETStruct)
{
//Ease of use function to compute the frequency of an asymmetric rotor
double Kappa;
	Kappa = Get_Kappa (Constants[0],Constants[1],Constants[2]);
	return fabs(Rigid_Rotor(Constants[0],Constants[2],J_Up,IndexUp,Kappa,ETStruct)-Rigid_Rotor(Constants[0],Constants[2],J_Low,IndexLow,Kappa,ETStruct));
}

int Get_Catalog (struct Transition *CatalogtoFill, double *Constants, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary)
{
//Utility function for calculating frequencies of a catalog
int i;	//Declaring i here because I like it, and apparently learned C pre C99Constants
	for (i=0;i<CatLines;i++) {	
		(CatalogtoFill)[i].Frequency = Get_Frequency (	MyDictionary[(CatalogtoFill)[i].Upper].J,
														MyDictionary[(CatalogtoFill)[i].Lower].J,
														(CatalogtoFill)[i].Upper,
														(CatalogtoFill)[i].Lower,
														Constants,  
														ETStruct
		);
		if (Verbose) printf ("%d ",(CatalogtoFill)[i].Type);		
		if (Verbose) print_Transition ((CatalogtoFill)[i],MyDictionary);
	} 
	return 1;
}

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
		case 1:
			Catalog_Comparator = &Catalog_Comparator_Intensity;	
			SortMethod = 1;
		case 2:
			Catalog_Comparator = &Catalog_Comparator_Index_Upper;
			SortMethod = 1;	
		case 3:
			Catalog_Comparator = &Catalog_Comparator_Index_Lower;
			SortMethod = 1;
		default:
			Catalog_Comparator = &Catalog_Comparator_Frequency;	
	}
	switch (SortMethod) {
		case 1:
			//Quick Sort - Fast for general sorting
			//This should be slower than insertion sort, but is here as a backup
			qsort(CatalogtoSort, TransitionCount, sizeof(struct Transition), Catalog_Comparator);
		default:
			//Insertion Sort - Use this one
			//Faster than quicksort when the target is already somewhat sorted
			//The program should in general produce fairly ordered catalogs already, so this ought to be the better choice
			insertionSort(CatalogtoSort, TransitionCount); 
			return;
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
	if (A.Intensity > B.Intensity) return 1;
	else if (A.Intensity < B.Intensity) return -1;
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
		
		if ((MyDictionary[(SourceCatalog)[i].Upper].J > JMin) && (MyDictionary[(SourceCatalog)[i].Upper].J < JMax)) {
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

int Find_Triples (struct Triple *TripletoFit, double *LineFrequencies, double Window, int LineCount)
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
	TripletoFit->TriplesList = realloc(TripletoFit->TriplesList,TripletoFit->TriplesCount[0]*TripletoFit->TriplesCount[1]*TripletoFit->TriplesCount[2]*sizeof(double));	
	return 1;
}

double Fit_Lines (double *Guess, int Verbose, struct Opt_Bundle GSLOptBundle)
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

void Initialize_Triples_Fitter (struct GSL_Bundle *FitBundle)
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
  	
  	double Temp;
  	x = gsl_vector_view_array (Guess, p);							//Set the guess	
  	//Extract the transitions we'll use for the fit
  	Transitions[0].Upper = TransitionstoFit.TransitionList[0].Upper;
  	Transitions[0].Lower = TransitionstoFit.TransitionList[0].Lower;
  	Transitions[1].Upper = TransitionstoFit.TransitionList[1].Upper;
  	Transitions[1].Lower = TransitionstoFit.TransitionList[1].Lower;
	Transitions[2].Upper = TransitionstoFit.TransitionList[2].Upper;
  	Transitions[2].Lower = TransitionstoFit.TransitionList[2].Lower;  
  	printf ("%i %i %i %i %i %i\n",TransitionstoFit.TransitionList[0].Upper,TransitionstoFit.TransitionList[0].Lower,TransitionstoFit.TransitionList[1].Upper,TransitionstoFit.TransitionList[1].Lower,TransitionstoFit.TransitionList[2].Upper,TransitionstoFit.TransitionList[2].Lower);
  	Wins = 0;		//Track the total number of wins for the current scoring system
  	Count = 0;		//Track the total number of constants (A+B+C) in the fit results
  	Iterations = 0;	//Variable to track the total number of iterations throughout the fit, just a bookeeping thing for me to see how the fitter is operating
  	Errors = 0;		//A count of the number of unconverged fits, another metric for me to track the fitting
  	double MyConstants[3];
  	printf ("%i\n",CatalogLines);
  	FitBundle->fdf.params = &MyOpt_Bundle;
  	for (i=0;i<TransitionstoFit.TriplesCount[0];i++) {
  		for (j=0;j<TransitionstoFit.TriplesCount[1];j++) {
  			for (k=0;k<TransitionstoFit.TriplesCount[2];k++) {
  				//This is why the transitions were extracted. We set the 3 transition's frequencies to those of the triple we're working on and hand it off to the fitter
  				Transitions[0].Frequency = TransitionstoFit.TriplesList[i];				
  				Transitions[1].Frequency = TransitionstoFit.TriplesList[j+TransitionstoFit.TriplesCount[0]];
  				Transitions[2].Frequency = TransitionstoFit.TriplesList[k+TransitionstoFit.TriplesCount[0]+TransitionstoFit.TriplesCount[1]];
				FitBundle->fdf.params = &MyOpt_Bundle;															//Set the parameters for this run
   				gsl_multifit_nlinear_init (&x.vector, &(FitBundle->fdf), FitBundle->Workspace);	//reInitialize the workspace incase this isnt the first run of the loop  	
				FitBundle->f = gsl_multifit_nlinear_residual(FitBundle->Workspace);								//compute initial cost function
  				gsl_multifit_nlinear_driver(50, xtol, gtol, ftol, NULL, NULL, &info, FitBundle->Workspace);		//solve the system with a maximum of 20 iterations
  				Iterations += gsl_multifit_nlinear_niter (FitBundle->Workspace); 	//Track the iterations
  				Final = gsl_multifit_nlinear_position(FitBundle->Workspace);		//Snag the results
  				if (gsl_multifit_nlinear_niter (FitBundle->Workspace) == 50) {		//Check for an error, currently only considering non-convergence
  					printf ("Error: Unconverged Fit: %.4f %.4f %.4f\n",Transitions[0].Frequency,Transitions[1].Frequency,Transitions[2].Frequency);
  					Errors++;
  				}
   				//(*FitResults)[Count] = gsl_vector_get(Final, 0);	//Get the final A
  				Temp = gsl_vector_get(Final, 0);	//Get the final A
  				//MyConstants[0] = (*FitResults)[Count];				//Update the constants
  				MyConstants[0] = Temp;
  				Count++;											//Track the number of items in FitResults
  				//(*FitResults)[Count] = gsl_vector_get(Final, 1);	
  				Temp = gsl_vector_get(Final, 1);
  				//MyConstants[1] = (*FitResults)[Count];
  				MyConstants[1] = Temp;
  				Count++;
  				//(*FitResults)[Count] = gsl_vector_get(Final, 2);
  				//MyConstants[2] = (*FitResults)[Count];
  				Temp = gsl_vector_get(Final, 2);
  				MyConstants[2] = Temp;
  				Count++;
  				Get_Catalog (*MyFittingCatalog, MyConstants, CatalogLines,0,MyOpt_Bundle.ETGSL,MyOpt_Bundle.MyDictionary);	//Now we recompute the full catalog, this isnt necessary to complete the fit, but has to be done to score the fit
  				Sort_Catalog (*MyFittingCatalog,CatalogLines,0,0);					//Catalog is not necessarily sorted, so we sort it
				TriplesScoreFunction (*MyFittingCatalog,ScoringParameters);
  			} 
  		} 
  	}
	printf ("%d\n",Wins);
  	printf ("%e %e %e\n",MyConstants[0],MyConstants[1],MyConstants[2]);
  	printf ("%i average iterations, %i Errors\n",Iterations/(TransitionstoFit.TriplesCount[0]*TransitionstoFit.TriplesCount[1]*TransitionstoFit.TriplesCount[2]),Errors);
	return 1;
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
		printf ("Found %d peaks\n",PeakCount);
	}
	*LineList = realloc(*LineList,PeakCount*sizeof(double));	
	return PeakCount;
}

void callback (const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
gsl_vector *f = gsl_multifit_nlinear_residual(w);
gsl_vector *x = gsl_multifit_nlinear_position(w);
double rcond;
	gsl_multifit_nlinear_rcond(&rcond, w);
	fprintf (stderr, "iter %2zu: A = %.4f, B = %.4f, C = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n", iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), 1.0 / rcond, gsl_blas_dnrm2(f));
}

void Calculate_State_Energies (struct Level *MyCatalog, struct Transition *SourceCatalog, int CatalogTransitions)
{
//Function to get the energies of levels in a catalog
int i;
	(MyCatalog[0]).Energy = 0.0;
	for (i=0;i<CatalogTransitions;i++) {
		(MyCatalog[SourceCatalog[i].Upper]).Energy = (MyCatalog[SourceCatalog[i].Lower]).Energy+SourceCatalog[i].Frequency;	//Keeping the calculation in MHz because it's easier
		(MyCatalog[SourceCatalog[i].Upper]).Energy *= 4.8E-5;
		//printf ("%d State Energy:%f Upper Index: %d LowerIndex: %d\n", i, (*Energies)[SourceCatalog[i].Upper], SourceCatalog[i].Upper, SourceCatalog[i].Lower);
	}
}

void Calculate_Intensities (struct Transition *SourceCatalog, int CatalogTransitions, struct Level *MyDictionary, double T, double *Dipoles)
{
int i;
	for (i=0;i<CatalogTransitions;i++) {
		(SourceCatalog[i]).Intensity = Dipoles[SourceCatalog[i].Type-1]*SourceCatalog[i].Frequency*fabs(exp(-1.0*(MyDictionary[SourceCatalog[i].Lower]).Energy/T)-exp(-1.0*(MyDictionary[SourceCatalog[i].Upper].Energy)/T));
		//printf ("%e %f %f\n",(*Intensity)[i],Energies[SourceCatalog[i].Lower],Energies[SourceCatalog[i].Upper]);
	}
}

int Search_DR_Hits (int DRPairs, double ConstStart, double ConstStop, double Step, double *DRFrequency, double Tolerance, int ExtraLineCount, double *ExtraLines, int **DRLinks, int LinkCount, struct Transition *CatalogtoFill, double *Constants, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary)
{
double CurrentA,CurrentB,CurrentC,Count;
int *Match,i,j,DRMatch,AllLinks,Wins;	

	Match = malloc(DRPairs*sizeof(int));
	if (Match == NULL) goto Error;
	CurrentA = ConstStart;
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
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
					Get_Catalog (	CatalogtoFill, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatLines,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									MyDictionary
					);
					//printf ("%f %f %f\n",CurrentA,CurrentB,CurrentC);
					Wins = 0;
					for (i=0;i<DRPairs;i++) Match[i] = -1;		//Reset our matches to 0
					for (i=0;i<CatLines;i++) {
						for (j=0;j<DRPairs;j++) {
							if (fabs(CatalogtoFill[i].Frequency-DRFrequency[j]) < Tolerance) {
								Match[j] = i;
								Wins++;
							}
						}
						for (j=0;j<ExtraLineCount;j++) if (fabs(CatalogtoFill[i].Frequency-ExtraLines[j]) < Tolerance) Wins++;	
					}
					for (i=0;i<DRPairs;i++) if (Match[i] < 0) DRMatch = 0;	//If any of our transitions werent matched we bail 
					DRMatch = 1;
					AllLinks = 1;
					if (DRMatch) { 
						for (i=0;i<LinkCount;i++) {
							if (!((CatalogtoFill)[Match[DRLinks[0][i]]].Upper == (CatalogtoFill)[Match[DRLinks[1][i]]].Upper) || ((CatalogtoFill)[Match[DRLinks[0][i]]].Lower == (CatalogtoFill)[Match[DRLinks[1][i]]].Lower) || ((CatalogtoFill)[Match[DRLinks[0][i]]].Upper == (CatalogtoFill)[Match[DRLinks[1][i]]].Lower) || ((CatalogtoFill)[Match[DRLinks[0][i]]].Lower == (CatalogtoFill)[Match[DRLinks[1][i]]].Upper)) {	
									AllLinks = 0;	
							}	
						}
						if (AllLinks) {
								//Do something cause at this point we have correctly matched the DR conditions
						}
					}					
				}
				CurrentC += Step;
			}
			CurrentB += Step;
		}
		CurrentA += Step;
	}
	printf ("%e individual fits performed",Count);
	free(Match);
	return 1;
Error:
	printf("Error running matching program");
	return 0;
}

double DummyFunction (struct Transition *MyCatalog, void *Data)
{
	
	return 1.0;
}





