/*
Built by Brandon Carroll 7/31/19 All Rights Reserved

C Brute force search extension to the Code

Build Command:
gcc -Wall -o Brute Brute\ Force\ Extension.c -lm -lgsl -lgslcblas -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Brute.so -shared -fPIC -O3 -funroll-loops Brute\ Force\ Extension.c -lm -lgsl -lgslcblas
*/

#ifndef __BRUTE_FORCE_H__
#define __BRUTE_FORCE_H__

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
#include "Fitter.h"


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
double Brute_Force_Pointer_Scoring (double /*ConstantsStart*/, double /*ConstantsStop*/, double /*ConstantsStep*/, double * /*ExperimentalLines*/, int /*ExperimentalLineCount*/, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, WinCounter /*WinFunction*/, int /*SaveCount*/, struct MultiSave * /*Saves*/, int /*Verbose*/);
double Brute_Force_Cube (struct Cube /*SearchCube*/, double * /*ExperimentalLines*/, int /*ExperimentalLineCount*/, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, WinCounter /*WinFunction*/, int /*SaveCount*/, struct MultiSave * /*Saves*/, char * /*FileName*/, int /*Verbose*/);
double Brute_Force_Fit_Four (double /*AStart*/, double /*AStop*/, double /*ConstantsStart*/, double /*ConstantsStop*/, double /*ConstantsStep*/, double * /*ExperimentalLines*/, int /*ExperimentalLineCount*/, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/, int /*SaveCount*/, struct MultiSave * /*Saves*/, int /*Verbose*/);



//Scoring Functions
int CountWins (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_No_Double (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_No_Double_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_No_Double_Nearest (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
int Find_Wins (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct Transition ** /*MatchedTransitions*/);
int Find_Wins_No_Double_Nearest (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct Transition ** /*FrequencyList*/, int /*Count*/, double **);

void insertionSort_Saves (struct MultiSave * /*SavestoSort*/, int /*SaveCount*/);
void insertionSort_Saves_Descending (struct MultiSave * /*SavestoSort*/, int /*SaveCount*/);
int Load_Exp_Lines  (char * /*FileName*/, double ** /*X*/, int /*Verbose*/);
int Allocate_MultiSave (int /*Size*/, struct MultiSave ** /*SavestoAllocate*/);
int Save_MultiSave (char * /*FileName*/, int /*Size*/, struct MultiSave * /*SavestoSave*/);
int Build_Axis_Linear (struct Axis * /*TargetAxis*/, double /*AxisStart*/, double /*AxisStepSize*/, unsigned int /*AxisSteps*/);
int Build_Cube_Linear (struct Cube * /*TargetCube*/, double * /*AxisStart*/, double * /*AxisStepSize*/, unsigned int * /*AxisSteps*/);



double Factorial (int /*Input*/);
int Pick_Four (int /*InputSize*/, int *** /*Return*/, int * /*ReturnSize*/);

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
							break;
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
							break;
						case 3:
							Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;
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
							break;
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
							break;
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
							break;
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
							break;
						case 3:
							Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;
						case 4:
							Wins = CountWins_No_Double_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;

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

double Brute_Force_Fit (double ConstantsStart, double ConstantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod, int SaveCount, struct MultiSave *Saves, int FrequencyCount, int Verbose)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double CurrentA,CurrentB,CurrentC,Wins,Count,Timing,ChiSqr;
double Constants[3];
double *FittingFrequencies,*FittedConstants;
int i,BadFits;
struct GSL_Bundle MyGSLBundle;
struct Opt_Bundle MyOptBundle;

	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	BadFits = 0;
	CurrentA = ConstantsStart;
	clock_t begin = clock();
	FittingFrequencies = malloc(sizeof(double));
	FittedConstants = malloc(3*sizeof(double));
	
	MyOptBundle.ETGSL = ETStruct;
	MyOptBundle.MyDictionary = SearchingDictionary;
	MyOptBundle.TransitionsGSL = NULL;
	
	for (i=0;i<SaveCount;i++) Saves[i].Score = 10000.0;
	
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
							break;
						case 2:
							Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
							break;
						case 3:
							Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;
						case 4:
							Wins = CountWins_No_Double_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;
						default:
							Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
							break;
					}
					if (Wins > 3) {
						MyOptBundle.TransitionCount = Wins;
						Initialize_SBFIT (&MyGSLBundle, &MyOptBundle);	
						FittingFrequencies = realloc(FittingFrequencies,Wins*sizeof(double));
						if (FittingFrequencies == NULL) goto Error;
						Find_Wins_No_Double_Nearest (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance, &(MyOptBundle.TransitionsGSL), Wins, &FittingFrequencies);
						if (!SBFIT (Constants, &ChiSqr, &MyGSLBundle, MyOptBundle, FittingFrequencies, &FittedConstants)) BadFits++;
						Wins = ChiSqr;
						if (Wins < Saves[0].Score) {
							Saves[0].Score = Wins;
							Saves[0].A = CurrentA;
							Saves[0].B = CurrentB;
							Saves[0].C = CurrentC;
							insertionSort_Saves_Descending(Saves, SaveCount); 
							if (Verbose > 1) printf ("New Good One -- %.2f %.2f %.2f %.2f Kappa:%f\n",Wins,CurrentA,CurrentB,CurrentC,Get_Kappa(CurrentA,CurrentB,CurrentC));
						}
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
	if (Verbose > 1) for (i=0;i<SaveCount;i++) printf ("%d: Score:%f %f %f %f\n",i,Saves[i].Score,Saves[i].A,Saves[i].B,Saves[i].C);
	if (Verbose) printf ("%.1f Fits in %.2f sec\n", Count,Timing);
	if (Verbose) printf ("%d Bad Fits\n",BadFits);
	free (FittingFrequencies);
	free (FittedConstants);
	return 1;
Error:
	printf ("Memory Error\n");
	return 0;	
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

double Brute_Force_Fit_Four (double AStart, double AStop, double ConstantsStart, double ConstantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod, int SaveCount, struct MultiSave *Saves, int Verbose)
{
/*
	-Newest variant of the brute force approach to fitting spectra.
	-The premise of this is that fitting three lines to three constants is essentially useless, but four lines will often have an unacceptable RMS, eliminating a lot of spurious fits
	-For each catalog a list of matching lines within the tolerance is produced, then an array of all possible combinations is made and fit, anything that is non ridiculous is kept

AStart (MHx) - Individual start for the A constant for subdividing this task into chunks for multithreading
AStop (MHz) - Individual stop for the A constant for subdividing this task into chunks for multithreading
ConstantsStart (MHz) - B/C Start
ConstantsStop (MHz) - B/C Stop
ConstantsStep (MHz) - Step for all constants
ExperimentalLines - List of experimental lines to match against
ExperimentalLineCount - # of lines in the experimental line list
SearchingCatalog - Catalog used for matching against experimental lines. For most cases, it's recommended this be fairly small catalog. Catalog should be built prior to passing it to the function
CatalogTransitions - # of transitions in the catalog, again unless you have a lot of processing power keep this low
Tolerance (MHz) - Absolute value Max amount a predicted line can differ from an experimental line and still be a match. Generally should be ~1-2x the stepsize
ETStruct - ETStruct used for fitting etc...
SearchingDictionary - Dictionary that goes with the catalog/ET
ScoreMethod - Method used for deciding how many matches we have, see the individual score methods for details, default is 1 which is a simple tolerance match, and generally the best to use since most filtering is done later
SaveCount - Max number of good saves to keep, currently not in use
Saves - Saves of good fits, currently not in use
Verbose - The standard verbosity flag, higher numbers produce higher levels of detail

*/
double CurrentA,CurrentB,CurrentC,Count,Timing,ChiSqr,Kappa,Delta,MaxKappa,MaxDelta,MinDelta,MaxChiSqr,BadFits;
double Constants[3];
double *FittingFrequencies,*FitConstants;
int i,FittableLines,BinomialSize,Wins;
int **PickFourArray;
struct GSL_Bundle MyGSLBundle;
struct Opt_Bundle MyOptBundle;
struct Transition *FoundLines;
FILE *FileHandle;

	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	BadFits = 0;
	MaxKappa = 0.95;
	MaxDelta = 50.0;
	MinDelta = -250.0;
	BadFits = 0.0;
	MaxChiSqr = 0.2;
	
	clock_t begin = clock();
	FoundLines = NULL;
	PickFourArray = malloc (50000000*sizeof(int *));
	for (i=0;i<50000000;i++) PickFourArray[i] = malloc(4*sizeof(int));
	FitConstants = malloc(3*sizeof(double));
	MyOptBundle.ETGSL = ETStruct;
	MyOptBundle.MyDictionary = SearchingDictionary;
	MyOptBundle.TransitionsGSL = NULL;
	MyOptBundle.TransitionCount = 4;
	Initialize_SBFIT (&MyGSLBundle, &MyOptBundle);	
	FittingFrequencies = malloc (4*sizeof(double));
	for (i=0;i<SaveCount;i++) Saves[i].Score = 10000.0;
	FoundLines = malloc (CatalogTransitions*sizeof(struct Transition)); //Array for the lines found to possibly match a set of constants, max number of lines we could match is the number of catalog transitions, realistically far fewer
	CurrentA = AStart;
	CurrentB = ConstantsStart;
	CurrentC = ConstantsStart;
	while (CurrentA < AStop) {
		CurrentB = ConstantsStart;
		while (CurrentB < ConstantsStop) {
			CurrentC = ConstantsStart;
			while (CurrentC < ConstantsStop) {
				if ((CurrentA > CurrentB) && (CurrentB > CurrentC)) {
							
					Constants[0] = CurrentA;
					Constants[1] = CurrentB;
					Constants[2] = CurrentC;
					Kappa = Get_Kappa(CurrentA,CurrentB,CurrentC);	//First we do the easy math to see if it's a structurally reasonable molecule
					Delta = Get_Delta(CurrentA,CurrentB,CurrentC);	//This is somewhat arbitrary, but grounded in experience, feel free to adjust as needed
					if ((Kappa < MaxKappa) && (Kappa > -1.0*MaxKappa) && (Delta < MaxDelta) && (Delta > MinDelta)) {	//Only move forward if we have a sane molecule
						Get_Catalog (	SearchingCatalog, 	//Catalog to compute frequencies for
										Constants, 			//Rotational constants for the calculation
										CatalogTransitions,	//# of transitions in the catalog
										0,					//Verbose
										ETStruct,			
										SearchingDictionary
									);		
						switch (ScoreMethod) {	//Score the catalog, really for this only the first method is decent
							case 1:
								Wins = CountWins (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
								break;
							case 2:
								Wins = CountWins_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);	
								break;
							case 3:
								Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
								break;
							case 4:
								Wins = CountWins_No_Double_Exp (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
								break;
							default:
								Wins = CountWins_No_Double (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance);
								break;
						}
						if (Wins > 3) {
							FittableLines = Find_Wins (ExperimentalLines, ExperimentalLineCount, SearchingCatalog, CatalogTransitions, Tolerance, &FoundLines);	//Separate function for actually pulling out the match transitions rather than just counting, could possibly combine with the score function above
							if ((FittingFrequencies == NULL) || (FittableLines < 4)) goto Error;	//If it didnt work for whatever reason we bail
							Pick_Four (FittableLines, &PickFourArray, &BinomialSize);	//Generate the array of possible options, easier to just remake it each time
							for (i=0;i<BinomialSize;i++) {	//Now we iterate through all possible sets of lines
								MyOptBundle.TransitionsGSL[0] = FoundLines[PickFourArray[i][0]];
								MyOptBundle.TransitionsGSL[1] = FoundLines[PickFourArray[i][1]];
								MyOptBundle.TransitionsGSL[2] = FoundLines[PickFourArray[i][2]];
								MyOptBundle.TransitionsGSL[3] = FoundLines[PickFourArray[i][3]];
								FittingFrequencies[0] = FoundLines[PickFourArray[i][0]].Frequency; //Assign the lines to the fitting setup
								FittingFrequencies[1] = FoundLines[PickFourArray[i][1]].Frequency;
								FittingFrequencies[2] = FoundLines[PickFourArray[i][2]].Frequency;
								FittingFrequencies[3] = FoundLines[PickFourArray[i][3]].Frequency;
								if (!SBFIT (Constants, &ChiSqr, &MyGSLBundle, MyOptBundle, FittingFrequencies, &FitConstants)) {	//Fit the four lines
									BadFits+=1.0;	//Count the bad fits for later
								} else {
									Kappa = Get_Kappa(FitConstants[0],FitConstants[1],FitConstants[2]);	//Recheck the structure now that weve fit, semi redundant but still can cut some junk out
									Delta = Get_Delta(FitConstants[0],FitConstants[1],FitConstants[2]);
									if ((ChiSqr < MaxChiSqr) && (Kappa < MaxKappa) && (Kappa > -1.0*MaxKappa) && (Delta < MaxDelta) && (Delta > MinDelta)) {	//Also recheck against chisqr, we need a converged fit with a sane chi sqr for a four line fit or theres no point in continuing
										Saves[0].Score = ChiSqr; 	//Save the best fit, but honestly, its pretty useless, this needs to be put through a much more rigorous 
										Saves[0].A = FitConstants[0];
										Saves[0].B = FitConstants[1];
										Saves[0].C = FitConstants[2];
										insertionSort_Saves_Descending(Saves, SaveCount); 
										FileHandle = fopen("TestLog.txt","a");
										fprintf (FileHandle,"%.3f %.3f %.3f\n",FitConstants[0],FitConstants[1],FitConstants[2]);
										fclose(FileHandle);
										if (Verbose > 1) {
											Count+=1.0;
											printf ("New Good One %.2e -- ChiSqr:%.2f A:%.2f B:%.2f C:%.2f Kappa:%f Delta:%f\n",Count,ChiSqr,FitConstants[0],FitConstants[1],FitConstants[2],Kappa,Delta);
										}
									} else {
										BadFits+=1.0;	//Count the bad fits for later
									}
								}
							}
						}
					}
				}
				CurrentC += ConstantsStep;
			}
			CurrentB += ConstantsStep;
		}
		CurrentA += ConstantsStep;
		printf ("A:%f\n",CurrentA);
	}
	clock_t end = clock();
	Timing = (double)(end - begin) / CLOCKS_PER_SEC;
	//if (Verbose > 1) for (i=0;i<SaveCount;i++) printf ("%d: Score:%f %f %f %f\n",i,Saves[i].Score,Saves[i].A,Saves[i].B,Saves[i].C);
	if (Verbose) printf ("%.1f Fits in %.2f sec\n", Count,Timing);
	if (Verbose) printf ("%f Bad Fits Ratio of bad to good fits: %.2f\n",BadFits,(BadFits/Count));
	free (FittingFrequencies);
	free (FitConstants);
	return 1;
Error:
	printf ("Memory Error\n");
	return 0;	
}


//Meta Functions


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
	return Wins;
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

int Find_Wins (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance, struct Transition **MatchedTransitions) 
{
int i,j, Wins,Max;
	Max = 100;
	Wins = 0;
	for (i=0;i<CatalogTransitions;i++) {
		for (j=0;j<ExperimentalLines;j++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				(*MatchedTransitions)[Wins] = SourceCatalog[i];
				(*MatchedTransitions)[Wins].Frequency = ExperimentalFrequencies[j];
				Wins++;
				if (Wins>Max-1) {
					return Max;
				}
			}
		}
	}
	return Wins;
}

int Find_Wins_No_Double_Nearest (double *ExperimentalFrequencies, int ExperimentalLines, struct Transition *SourceCatalog, int CatalogTransitions, double Tolerance, struct Transition **FrequencyList, int Count, double **FittingFrequencies) 
{
double LastWin,LastExp,Test;
int i,j,Wins;
	Wins = 0;
	LastWin = 0.0;
	LastExp = 0.0;
	for (i=0;i<CatalogTransitions;i++) {
		for (j=0;j<ExperimentalLines;j++) {
			if (fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency) < Tolerance) {
				if ((fabs(LastWin-SourceCatalog[i].Frequency) > Tolerance) && (fabs(LastExp-ExperimentalFrequencies[j]) > Tolerance)) {
					Test = fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency); 
					if (j<ExperimentalLines-1) j++;
					else break;
					while ((fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency)<Test) && (j<ExperimentalLines-1)) {
						Test = fabs(ExperimentalFrequencies[j]-SourceCatalog[i].Frequency);
						j++;
					}
					j--;
					(*FrequencyList)[Wins] = SourceCatalog[i];
					(*FittingFrequencies)[Wins] = ExperimentalFrequencies[j];
					LastWin = SourceCatalog[i].Frequency;
					LastExp = ExperimentalFrequencies[j];
					Wins++;
					break;
				}
			}
		}
	}
	//printf ("==\n");
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

void insertionSort_Saves_Descending(struct MultiSave *SavestoSort, int SaveCount) 
{ 
//Insertion sort for sorting the saves
struct MultiSave Key;
int i, j; 
	SaveCount++;
	for (i=1;i<SaveCount;i++) { 
		Key = SavestoSort[i]; 
        j = i - 1; 
        while (j >= 0 && SavestoSort[j].Score < Key.Score) { 
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

int Build_Axis_Linear (struct Axis *TargetAxis, double AxisStart, double AxisStepSize, unsigned int AxisSteps)
{
int i;
	TargetAxis->Length = AxisSteps;
	TargetAxis->Array = malloc(TargetAxis->Length*sizeof(double));
	if (TargetAxis->Array == NULL) goto Error;
	for (i=0;i<TargetAxis->Length;i++) TargetAxis->Array[i] = AxisStart+i*AxisStepSize;
	return 1;
Error:
	printf ("Error allocating axis\n");
	return 0;	
}

int Build_Cube_Linear (struct Cube *TargetCube, double *AxisStart, double *AxisStepSize, unsigned int *AxisSteps)
{
	Build_Axis_Linear (&(TargetCube->AAxis),AxisStart[0],AxisStepSize[0],AxisSteps[0]);
	Build_Axis_Linear (&(TargetCube->AAxis),AxisStart[1],AxisStepSize[1],AxisSteps[1]);
	Build_Axis_Linear (&(TargetCube->AAxis),AxisStart[2],AxisStepSize[2],AxisSteps[2]);
	return 1;	
}

double Factorial (int Input) 
{
//Basic function to compute the factorial of an integer input
double n;
int i;
	n = 1;
	for (i=1;i<=Input;i++) n= (double) n*i;
	return n;
}

int Pick_Four (int InputSize, int ***Return, int *ReturnSize)
{
/*
  
 Builds an array of all possible combinations of 4 input lines
  
*/
int i,j,k,l,Binomial,Count;
	Binomial = (int) (Factorial(InputSize)/(Factorial(InputSize-4)*24));
	if (Binomial > 50000000) {
		Binomial = 50000000;
		printf ("Warning binomial size exceeds 1GB, restricting to 50000000 sets\n");
	}
	Count = 0;
	for (i=0;i<InputSize-3;i++) {
		for (j=i+1;j<InputSize-2;j++) {
			for (k=j+1;k<InputSize-1;k++) {
				for (l=k+1;l<InputSize;l++) {
					if ((i==j) || (i==k) || (i==l) || (j==k) || (j==l) || (k==l)) break;
					else {
						//printf ("%d: %d %d %d %d\n",Count,i,j,k,l);
						(*Return)[Count][0] = i;
						(*Return)[Count][1] = j;
						(*Return)[Count][2] = k;
						(*Return)[Count][3] = l;
						Count++;
					}
				}
			}
		}
	}
	(*ReturnSize) = Binomial;
	return 1;
}

#endif /* __BRUTE_FORCE_H__ */