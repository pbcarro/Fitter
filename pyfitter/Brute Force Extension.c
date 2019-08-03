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

//=============Function Prototypes==============

//Brute Force Functions
double Brute_Force (double /*CostantsStart*/, double /*CosntantsStop*/, double /*ConstantsStep*/, double * /*ExperimentalLines*/, int ExperimentalLineCount, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/);
double Brute_Force_ConstantsArray (double * /*ConstantsArray*/, int /*ConstantsSize*/, double * /*ExperimentalLines*/, int ExperimentalLineCount, struct Transition * /*SearchingCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/, struct ETauStruct /*ETStruct*/, struct Level * /*SearchingDictionary*/, int /*ScoreMethod*/);


//Scoring Functions
int CountWins (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);
double CountWins_Exp (double * /*ExperimentalFrequencies*/, int /*ExperimentalLines*/, struct Transition * /*SourceCatalog*/, int /*CatalogTransitions*/, double /*Tolerance*/);


void insertionSort_Saves(struct MultiSave * /*SavestoSort*/, int /*SaveCount*/);

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

double Brute_Force_Top_Results (double ConstantsStart, double ConstantsStop, double ConstantsStep, double *ExperimentalLines, int ExperimentalLineCount, struct Transition *SearchingCatalog, int CatalogTransitions, double Tolerance, struct ETauStruct ETStruct, struct Level *SearchingDictionary, int ScoreMethod, int SaveCount, struct MultiSave *Saves)
{
//Variant of the Brute Force search that takes an array of values for the constants so you can use nonlinear steps for more effective searches
double CurrentA,CurrentB,CurrentC,Wins,Count;
double Constants[3];
	Count = 0.0;		//Tracking the number of counts we perform, using doubles to prevent int overflow
	Wins = 0;
	CurrentA = ConstantsStart;
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
					}
					if (Wins > Saves[SaveCount-1].Score) {
						Saves[SaveCount-1].Score = Wins;
						Saves[SaveCount-1].A = CurrentA;
						Saves[SaveCount-1].B = CurrentB;
						Saves[SaveCount-1].C = CurrentC;
						insertionSort_Saves(Saves, SaveCount); 
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




