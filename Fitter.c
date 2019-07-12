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


//=============Structures==============
struct Level
{
	unsigned int Index;
	unsigned int J;
	unsigned int Ka;
	unsigned int Kc;
	//double Energy;
};

struct Transition 
{
	double Frequency;
	unsigned int Upper;
	unsigned int Lower;
	unsigned int Type;
	//double Intensity;
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


//=============Function Prototypes==============

//Program setup functions
int Initialize_Stuff (double **/*ETArray*/, int */*CatTransitions*/, int */*DictTransitions*/, double */*FileDelta*/, int */*StatePoints*/, struct Level **/*DictionaryIn*/, struct Transition **/*CatalogIn*/);
int Load_ETau_File (char */*FileName*/, double **/*X*/, double */*FileDelta*/, int */*StatePoints*/, int */*StateCount*/);
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
int Fill_Catalog_Restricted (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*FrequencyLow*/, double /*FrequencyHigh*/, int */*Dipoles*/, int /*Verbose*/, struct Level */*MyDictionary*/);
void Sort_Catalog (struct Transition */*Catalog*/, int /*TransitionCount*/, int /*SortMethod*/);
int Catalog_Comparator (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Upper (const void */*a*/, const void */*b*/);
int Catalog_Comparator_Index_Lower (const void */*a*/, const void */*b*/);
void insertionSort(struct Transition */*CatalogtoSort*/, int /*TransitionCount*/);
void Calculate_State_Energies (double **/*UpperStateEnergies*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/);
void Calculate_Intensities (double **/*Intensity*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/, double */*Energies*/, double /*T*/, double */*Dipoles*/);


int main (int argc, char *argv[])  
{
double Constants[3], Dipoles[3];	
struct ETauStruct ETStruct;
struct Transition *BaseCatalog, *SortedCatalog;	//Model catalog used to save time and simplify	
struct Level	*BaseDict;						//The base catalog dictionary, translates from an index to J/Ka/Kc 
int CatalogTransitions,DictionaryLevels;		//Number of transitions in the catalogs
double *EnergyLevels,*IntensityVals;

	int Loops = 100000;
	int i,j;
	double Timing[100];	
	if (!Initialize_Stuff(&(ETStruct.ETVals),&CatalogTransitions,&DictionaryLevels,&(ETStruct.Delta),&(ETStruct.StatePoints),&BaseDict,&BaseCatalog)) {
		goto Error;
	}
	
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
	Calculate_State_Energies (&EnergyLevels, SortedCatalog, CatalogTransitions);
	Calculate_Intensities (&IntensityVals, BaseCatalog, CatalogTransitions, EnergyLevels, 3.0, Dipoles);	

	
	
	for (j=0;j<100;j++) {
		clock_t begin = clock();
		for (i=0;i<Loops;i++) {
			Get_Catalog (	BaseCatalog, 		//Catalog to compute frequencies for
							Constants, 			//Rotational constants for the calculation
							CatalogTransitions,	//# of transitions in the catalog
							0,					//Verbose
							ETStruct,
							BaseDict
			);
			Sort_Catalog (BaseCatalog,CatalogTransitions,2);
		}
		clock_t end = clock();
		Timing[j] = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("%f\n", Timing[j]);
		
	}
	double Sum,Mean,StdDev;
	Sum = 0.0;
	StdDev = 0.0;
	for(i=0; i<100; ++i)
    {
        Sum += Timing[i];
    }
    Mean = Sum/100.0;
    for(i=0; i<100; ++i)
        StdDev += pow(Timing[i] - Mean, 2.0);
    StdDev = sqrt(StdDev/100);
	printf ("Average Time:%f, Standard Deviation:%f\n",Mean,StdDev);	
	return 1;

Error:
	printf ("Error Initializing program, exiting\n");
	return 1;
}

int Initialize_Stuff (double **ETArray, int *CatTransitions, int *DictLevels, double *FileDelta, int *StatePoints, struct Level **DictionaryIn, struct Transition **CatalogIn) 
{
int TempPoints;
	TempPoints = 0;
	srand(time(NULL));								//Initializing the time calls elsewhere	
	*DictLevels = Load_Base_Catalog_Dictionary (	"Base Catalog Dictionary1.txt", 
													DictionaryIn,//&BaseDict, 
													0 //Verbose
	);	
	printf ("Dictionary Loaded\n");
	*CatTransitions = Load_Base_Catalog (	"Base Catalog.txt",
											CatalogIn, 		//List of transitions in the lower state
											0				//Verbose
	);	
	printf ("Catalog Loaded\n");

	//Load all of our ET values into the array and check that it worked, using static external files cause I'm lazy
	if (Load_ETau_File ("J0_25_dk2.dat", ETArray,FileDelta,StatePoints,&TempPoints)) printf ("Tables Loaded\n");
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
int i;	//Declaring i here because I like it, and apparently learned C pre C99
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

void Sort_Catalog (struct Transition *CatalogtoSort, int TransitionCount, int SortMethod) 
{
//Wrapper function for sorting a catalog
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

int Catalog_Comparator (const void *a, const void *b) 
{
//Comparison function for qsort sorting of the catalog
	struct Transition A = *(struct Transition *) a;
	struct Transition B = *(struct Transition *) b;
	if (A.Frequency > B.Frequency) return 1;
	else if (A.Frequency < B.Frequency) return -1;
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

int Fill_Catalog_Restricted (struct Transition *SourceCatalog, struct Transition **CatalogtoFill, double *Constants, int CatLines, double FrequencyLow, double FrequencyHigh, int *Dipoles, int Verbose, struct Level *MyDictionary)
{
int i, CatLinesOut;     
	*CatalogtoFill = malloc(CatLines*sizeof(struct Transition));
	CatLinesOut = 0; 
	for (i=0;i<CatLines;i++) {
		if (((SourceCatalog)[i].Frequency > FrequencyLow) && ((SourceCatalog)[i].Frequency < FrequencyHigh) &&  ((Dipoles[(SourceCatalog)[i].Type-1] == 1))) {
			(*CatalogtoFill)[CatLinesOut] = (SourceCatalog)[i];
			if (Verbose) printf ("%d ",(*CatalogtoFill)[CatLinesOut].Type);				
			if (Verbose) print_Transition ((*CatalogtoFill)[CatLinesOut],MyDictionary);
			CatLinesOut++;
		}
     } 

	*CatalogtoFill = realloc (*CatalogtoFill,CatLinesOut*sizeof(struct Transition));   
	return CatLinesOut;
}

int Peak_Find (double **LineList, double Max, double Min, double *X, double *Y, int ArraySize)
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
	*LineList = realloc(*LineList,i*sizeof(double));	
	return PeakCount;
}

void Calculate_State_Energies (double **Energies, struct Transition *SourceCatalog, int CatalogTransitions)
{
//Function to get the energies of levels in a catalog
//
int i;
	(*Energies)[0] = 0.0;
	for (i=0;i<CatalogTransitions;i++) {
		(*Energies)[SourceCatalog[i].Upper] = (*Energies)[SourceCatalog[i].Lower]+SourceCatalog[i].Frequency;	//Keeping the calculation in MHz because it's easier
		(*Energies)[SourceCatalog[i].Upper] *= 4.8E-5;
		//printf ("%d State Energy:%f Upper Index: %d LowerIndex: %d\n", i, (*Energies)[SourceCatalog[i].Upper], SourceCatalog[i].Upper, SourceCatalog[i].Lower);
	}
}

void Calculate_Intensities (double **Intensity, struct Transition *SourceCatalog, int CatalogTransitions, double *Energies, double T, double *Dipoles)
{
int i;
	for (i=0;i<CatalogTransitions;i++) {
		(*Intensity)[i] = Dipoles[SourceCatalog[i].Type-1]*SourceCatalog[i].Frequency*fabs(exp(-1.0*Energies[SourceCatalog[i].Lower]/T)-exp(-1.0*Energies[SourceCatalog[i].Upper]/T));
		//printf ("%e %f %f\n",(*Intensity)[i],Energies[SourceCatalog[i].Lower],Energies[SourceCatalog[i].Upper]);
	}
}









