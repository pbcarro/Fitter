/*
Built by Brandon Carroll 4/7/19 All Rights Reserved

C Parser Extension to the Code

Build Command:
gcc -Wall -o Parse Parser\ Extension.c -lm -lgsl -lgslcblas -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Parser.so -shared -fPIC -O3 -funroll-loops Parser\ Extension.c -lm -lgsl -lgslcblas
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
struct FullState
{
	double Frequency;
	unsigned int UpperState;
	unsigned int LowerState;
	unsigned int Type;
	unsigned int J_Upper;
	unsigned int Ka_Upper;
	unsigned int Kc_Upper;
	unsigned int J_Lower;
	unsigned int Ka_Lower;
	unsigned int Kc_Lower;
};

//=============Function Prototypes==============
int Allocate_FullStates (int /*StateCount*/, struct FullState ** /*StateArray*/);
int Free_FullStates (struct FullState ** /*StateArray*/);
int Parse_Catalog (int StateCount, struct FullState * /*StateArray*/, struct Transition * /*SourceCatalog*/, struct Level * /*SourceDictionary*/, int /*Verbose*/);

//Functions
int Allocate_FullStates (int StateCount, struct FullState **StateArray)
{

	*StateArray = malloc (StateCount*sizeof(struct FullState));
	if (*StateArray == NULL) goto Error;
return 1;

Error:
	printf ("Error: Unable to allocate states");
	return 0;
}

int Free_FullStates (struct FullState **StateArray)
{
	free(*StateArray);
return 1;
}

int Parse_Catalog (int StateCount, struct FullState *StateArray, struct Transition *SourceCatalog, struct Level *SourceDictionary, int Verbose)
{
int i;
	for (i=0;i<StateCount;i++) {
		StateArray[i].Frequency = SourceCatalog[i].Frequency;
		StateArray[i].UpperState = SourceCatalog[i].Upper;
		StateArray[i].LowerState = SourceCatalog[i].Lower;
		StateArray[i].Type = SourceCatalog[i].Type;
		StateArray[i].J_Upper = SourceDictionary[SourceCatalog[i].Upper].J;
		StateArray[i].Ka_Upper = SourceDictionary[SourceCatalog[i].Upper].Ka;
		StateArray[i].Kc_Upper = SourceDictionary[SourceCatalog[i].Upper].Kc;
		StateArray[i].J_Lower = SourceDictionary[SourceCatalog[i].Lower].J;
		StateArray[i].Ka_Lower = SourceDictionary[SourceCatalog[i].Lower].Ka;
		StateArray[i].Kc_Lower = SourceDictionary[SourceCatalog[i].Lower].Kc;
		
		
		if (Verbose) {
			printf ("Index:%d\tFrequency:%f\tUpper State Index:%d\tLower State Index:%d\tType:%d\tJ Upper:%d\tKa Upper:%d Kc Upper:%d\tJ Lower:%d\tKa Lower:%d\tKc Lower:%d\n",i,	StateArray[i].Frequency,
																										StateArray[i].UpperState,
																										StateArray[i].LowerState,
																										StateArray[i].Type,
																										StateArray[i].J_Upper,
																										StateArray[i].Ka_Upper,
																										StateArray[i].Kc_Upper,
																										StateArray[i].J_Lower,
																										StateArray[i].Ka_Lower,
																										StateArray[i].Kc_Lower
																										);
		}
	}			
	return 1;
}







