/*
Built by Brandon Carroll 7/31/19 All Rights Reserved

C Brute force search extension to the Code

Build Command:
gcc -Wall -o Brute Brute\ Force\ Extension.c -lm -lgsl -lgslcblas -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Brute.so -shared -fPIC -O3 -funroll-loops Brute\ Force\ Extension.c -lm -lgsl -lgslcblas
*/

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
#include "Fitter.c"


//=============Structures==============
struct MultiSave 
{
	double Score;
	double A;
	double B;
	double C;
};

struct Axis
{
	double *Array;
	unsigned int Length;
};

struct Cube
{
	struct Axis AAxis;
	struct Axis BAxis;
	struct Axis CAxis;
};

typedef double (*WinCounter)(double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);	//Generic function pointer for the scoring function used in the triples fitter


//=============Function Prototypes==============

//Brute Force Functions
double Brute_Force (double /*CostantsStart*/, double /*CosntantsStop*/, double /*ConstantsStep*/, double * /*ExperimentalLines*/, int ExperimentalLineCount, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/);
double Brute_Force_ConstantsArray (double * /*ConstantsArray*/, int /*ConstantsSize*/, double * /*ExperimentalLines*/, int ExperimentalLineCount, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/);
double Brute_Force_Top_Results (double /*ConstantsStart*/, double /*ConstantsStop*/, double /*ConstantsStep*/, double * /*ExperimentalLines*/, int /*ExperimentalLineCount*/, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/, int /*SaveCount*/, struct MultiSave * /*Saves*/, int /*Verbose*/);

double Brute_Force_Cube (struct Cube /*SearchCube*/, double * /*ExperimentalLines*/, int /*ExperimentalLineCount*/, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, WinCounter /*WinFunction*/, int /*SaveCount*/, struct MultiSave * /*Saves*/, char * /*FileName*/, int /*Verbose*/);

//Scoring Functions
int CountWins (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_No_Double (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_No_Double_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);


void insertionSort_Saves(struct MultiSave * /*SavestoSort*/, int /*SaveCount*/);
int Load_Exp_Lines  (char * /*FileName*/, double ** /*X*/, int /*Verbose*/);
int Allocate_MultiSave (int /*Size*/, struct MultiSave ** /*SavestoAllocate*/);
int Save_MultiSave (char * /*FileName*/, int /*Size*/, struct MultiSave * /*SavestoSave*/);


//=============Functions========================
double Brute_Force (double CostantsStart, double CosntantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod)
{
double CurrentA,CurrentB,CurrentC,BestWins,BestA,BestB,BestC,Wins,Count;
double Constants[3];
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	BestWins = -1;
	CurrentA = CostantsStart;
	while (CurrentA < CosntantsStop) {
		CurrentB = CostantsStart;
		while (CurrentB < CosntantsStop) {
			CurrentC = CostantsStart;
			while (CurrentC < CosntantsStop) {
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					Get_Catalog (	SearchingCatalog, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatalogTransitions,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									SearchingDictionary
								);
					switch (ScoreMethod) {
						case 1:
							Wins = CountWins (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
					}
					if (Wins > BestWins) {
						BestA = CurrentA;
						BestB = CurrentB;
						BestC = CurrentC;
						BestWins = Wins;
						printf ("New Leader: A:%f B:%f C:%f Wins:%f Kappa: %f\n",BestA,BestB,BestC,BestWins,Get_Kappa(BestA,BestB,BestC));
					}
					Count+=1.0;
				}
				CurrentC += ConstantsStep;
			}
			CurrentB += ConstantsStep;
		}
		CurrentA += ConstantsStep;
		printf ("A:%f\n",CurrentA);
	}
	return BestWins;
}

double Brute_Force_ConstantsArray (double *ConstantsArray, int ConstantsSize, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double CurrentA,CurrentB,CurrentC,BestWins,BestA,BestB,BestC,Wins,Count;
double Constants[3];
int i,j,k;
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	BestWins = -1;
	for (i=0;i<ConstantsSize;i++) {
		CurrentA = ConstantsArray[i];
		for (j=0;j<ConstantsSize;j++) {
			CurrentB = ConstantsArray[i];
			for (k=0;k<ConstantsSize;k++) {
				CurrentC = ConstantsArray[i];
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					Get_Catalog (	SearchingCatalog, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatalogTransitions,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									SearchingDictionary
								);
					switch (ScoreMethod) {
						case 1:
							Wins = CountWins (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
					}
					if (Wins > BestWins) {
						BestA = CurrentA;
						BestB = CurrentB;
						BestC = CurrentC;
						BestWins = Wins;
						printf ("New Leader: A:%f B:%f C:%f Wins:%f Kappa: %f\n",BestA,BestB,BestC,BestWins,Get_Kappa(BestA,BestB,BestC));
					}
					Count+=1.0;
				}
			}
		}
		printf ("A:%f\n",CurrentA);
	}
	return BestWins;
}

double Brute_Force_Top_Results (double ConstantsStart, double ConstantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod, int SaveCount, struct MultiSave *Saves, int Verbose)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double CurrentA,CurrentB,CurrentC,Wins,Count,Timing;
double Constants[3];
int i;
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	CurrentA = ConstantsStart;
	clock_t begin = clock();
	while (CurrentA < ConstantsStop) {
		CurrentB = ConstantsStart;
		while (CurrentB < ConstantsStop) {
			CurrentC = ConstantsStart;
			while (CurrentC < ConstantsStop) {
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {
							
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					Get_Catalog (	SearchingCatalog, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatalogTransitions,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									SearchingDictionary
								);
					switch (ScoreMethod) {
						case 1:
							Wins = CountWins (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
						case 3:
							Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
						case 4:
							Wins = CountWins_No_Double_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);

					}
					if (Wins > Saves[0].Score) {
						Saves[0].Score = Wins;
						Saves[0].A = CurrentA;
						Saves[0].B = CurrentB;
						Saves[0].C = CurrentC;
						insertionSort_Saves(Saves, SaveCount); 
						if (Verbose > 1) printf ("New Good One -- %.2f %.2f %.2f %.2f Kappa:%f\n",Wins,CurrentA,CurrentB,CurrentC,Get_Kappa(CurrentA,CurrentB,CurrentC));
					}
					Count+=1.0;
				}
				CurrentC += ConstantsStep;
			}
			CurrentB += ConstantsStep;
		}
		CurrentA += ConstantsStep;
		if (Verbose) printf ("A:%f\n",CurrentA);
	}
	clock_t end = clock();
	Timing = (double)(end - begin) / CLOCKS_PER_SEC;
	if (Verbose) printf ("%.1f Fits in %.2f sec\n", Count,Timing);
	if (Verbose > 1) for (i=0;i<SaveCount;i++) printf ("%d: Score:%f %f %f %f\n",i,Saves[i].Score,Saves[i].A,Saves[i].B,Saves[i].C);
	return 1;
}

double Brute_Force_Pointer_Scoring (double ConstantsStart, double ConstantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, WinCounter WinFunction, int SaveCount, struct MultiSave *Saves, int Verbose)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double CurrentA,CurrentB,CurrentC,Wins,Count,Timing;
double Constants[3];
int i;
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	CurrentA = ConstantsStart;
	clock_t begin = clock();
	while (CurrentA < ConstantsStop) {
		CurrentB = ConstantsStart;
		while (CurrentB < ConstantsStop) {
			CurrentC = ConstantsStart;
			while (CurrentC < ConstantsStop) {
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {					
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					Get_Catalog (	SearchingCatalog, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatalogTransitions,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									SearchingDictionary
								);
					Wins = WinFunction (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
					if (Wins > Saves[0].Score) {
						Saves[0].Score = Wins;
						Saves[0].A = CurrentA;
						Saves[0].B = CurrentB;
						Saves[0].C = CurrentC;
						insertionSort_Saves(Saves, SaveCount); 
						if (Verbose > 1) printf ("New Good One -- %.2f %.2f %.2f %.2f Kappa:%f\n",Wins,CurrentA,CurrentB,CurrentC,Get_Kappa(CurrentA,CurrentB,CurrentC));
					}
					Count+=1.0;
				}
				CurrentC += ConstantsStep;
			}
			CurrentB += ConstantsStep;
		}
		CurrentA += ConstantsStep;
		if (Verbose) printf ("A:%f\n",CurrentA);
	}
	clock_t end = clock();
	Timing = (double)(end - begin) / CLOCKS_PER_SEC;
	if (Verbose) printf ("%.1f Fits in %.2f sec\n", Count,Timing);
	if (Verbose > 1) for (i=0;i<SaveCount;i++) printf ("%d: Score:%f %f %f %f\n",i,Saves[i].Score,Saves[i].A,Saves[i].B,Saves[i].C);
	return 1;
}

double Brute_Force_Cube (struct Cube SearchCube, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, WinCounter WinFunction, int SaveCount, struct MultiSave *Saves, char *FileName, int Verbose)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double Wins,Count,Timing;
double Constants[3];
int i,j,k;
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0.0;
	clock_t begin = clock();
	for (i=0;i<SearchCube.AAxis.Length;i++) {
		Constants[0] = SearchCube.AAxis.Array[i];
		for (j=0;k<SearchCube.BAxis.Length;j++) {
			Constants[1] = SearchCube.BAxis.Array[j];
			for (k=0;k<SearchCube.CAxis.Length;k++) {
				Constants[2] = SearchCube.CAxis.Array[k];
				if ((Constants[0] > Constants[1]) && (Constants[1] > Constants[2])) {					
					Get_Catalog (	SearchingCatalog, 		//Catalog to compute frequencies for
									Constants, 			//Rotational constants for the calculation
									CatalogTransitions,	//# of transitions in the catalog
									0,					//Verbose
									ETStruct,
									SearchingDictionary
								);
					Wins = WinFunction (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
					if (Wins > Saves[0].Score) {
						Saves[0].Score = Wins;
						Saves[0].A = Constants[0];
						Saves[0].B = Constants[1];
						Saves[0].C = Constants[2];
						insertionSort_Saves(Saves, SaveCount); 
						if (Verbose > 1) printf ("New Good One -- %.2f %.2f %.2f %.2f Kappa:%f\n",Wins,Constants[0],Constants[1],Constants[2],Get_Kappa(Constants[0],Constants[1],Constants[2]));
					}
					Count+=1.0;
				}
			}
		}
		if (Verbose) printf ("A:%f\n",Constants[0]);
	}
	clock_t end = clock();
	Timing = (double)(end - begin) / CLOCKS_PER_SEC;
	if (FileName != NULL) Save_MultiSave (FileName, SaveCount, Saves);
	if (Verbose) printf ("%.1e Fits in %.2f sec\n", Count,Timing);
	if (Verbose > 1) for (i=0;i<SaveCount;i++) printf ("%d: Score:%f %f %f %f\n",i,Saves[i].Score,Saves[i].A,Saves[i].B,Saves[i].C);
	return 1;
}

//Scoring functions
int CountWins (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance) 
{
int i,j, Wins;
	Wins = 0;
	for (i=0;i<CatalogTransitions;i++) {
		for (j=0;j<ExperimentalLines;j++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				Wins++;
				break;
			}
		}
	}
	return Wins;
}

double CountWins_Exp (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance) 
{
int i,j;
double Wins;
	Wins = 0.0;
	for (i=0;i<CatalogTransitions;i++) {
		for (j=0;j<ExperimentalLines;j++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				Wins += exp(-fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency));
				break;
			}
		}
	}
	return Wins;
}

double CountWins_No_Double (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance) 
{
double Wins,LastWin,LastExp;
int i,j;
	Wins = 0.0;
	LastWin = 0.0;
	LastExp = 0.0;
	for (j=0;j<ExperimentalLines;j++) {
		for (i=0;i<CatalogTransitions;i++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				if ((fabs(LastWin-SourceCatalog[i].Frequency) > Tolerance) && (fabs(LastExp-ExperimentalFrequencies[j]) > Tolerance)) {
					Wins++;
					LastWin = SourceCatalog[i].Frequency;
					LastExp = ExperimentalFrequencies[j];
					break;
				}
			}
		}
	}
	return  Wins;
}

double CountWins_No_Double_Exp (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance) 
{
double Wins,LastWin,LastExp;
int i,j;
	Wins = 0.0;
	LastWin = 0.0;
	LastExp = 0.0;
	for (j=0;j<ExperimentalLines;j++) {
		for (i=0;i<CatalogTransitions;i++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				if ((fabs(LastWin-SourceCatalog[i].Frequency) > Tolerance) && (fabs(LastExp-ExperimentalFrequencies[j]) > Tolerance)) {
					Wins += exp(-fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency));
					LastWin = SourceCatalog[i].Frequency;
					LastExp = ExperimentalFrequencies[j];
					break;
				}
			}
		}
	}
	return  Wins;
}

//General Functions
void insertionSort_Saves(struct MultiSave *SavestoSort, int SaveCount) 
{ 
//Insertion sort for sorting the saves
struct MultiSave Key;
int i, j; 
	SaveCount++;
	for (i=1;i<SaveCount;i++) { 
		Key = SavestoSort[i]; 
        j = i - 1; 
        while (j >= 0 && SavestoSort[j].Score > Key.Score) { 
            SavestoSort[j + 1] = SavestoSort[j]; 
            j = j - 1; 
        } 
        SavestoSort[j + 1] = Key; 
    } 
} 

int Load_Exp_Lines  (char *FileName, double **X, int Verbose)
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
	while (fscanf (FileHandle, "%lf", &(*X)[i]) == 1) {
		i++;
		if (i >= (FileLimit-1)) {
			if (i > 100000000) {
				printf ("Error: Experimental file exceeds 10 million points, dial that down a bit\n");
				goto Error;
			}
			FileLimit += 200000;
			*X = realloc (*X,FileLimit*sizeof(double));
		}
	}
	fclose (FileHandle);												//Not currently checking to see if fclose works, shouldnt affect file load and theres not much to be done if there is an error, the file stream should close when the program ends so hopefully this wont mater														//i is incremented at the end when the EOF occurs which doesnt represent real data so it gets fixed 								
	*X = realloc(*X,i*sizeof(double));
	if (Verbose) printf ("Loaded %d lines\n",i);
	return i;
Error:
	printf ("Error Loading file %s\n",FileName);
	return 0;
}

int Allocate_MultiSave (int Size, struct MultiSave **SavestoAllocate)
{
	*SavestoAllocate = malloc(Size*sizeof(struct MultiSave));
	if (*SavestoAllocate == NULL) goto Error;
	return 1;
Error:
	return 0;	
}

int Save_MultiSave (char *FileName, int Size, struct MultiSave *SavestoSave)
{
int i;
FILE *FileHandle;
	
	FileHandle = NULL;
	FileHandle = fopen (FileName, "w");									//Open file read only
	if (FileHandle == NULL) goto Error;
	for (i=0;i<Size;i++) {
		fprintf (FileHandle,"%.3f\t%.2f\t%.2f\t%.2f\n",SavestoSave[i].Score,SavestoSave[i].A,SavestoSave[i].B,SavestoSave[i].C);
	}
	fclose(FileHandle);
	return 1;
Error:
	printf ("Error: Unable to open save file");
	return 0;	
}

