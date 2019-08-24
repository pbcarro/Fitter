/*
Built by Brandon Carroll 7/31/19 All Rights Reserved

C Brute force search extension to the Code

Build Command:
gcc -Wall -o Brute Brute\ Force\ Extension.c -lm -lgsl -lgslcblas -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Brute.so -shared -fPIC -O3 -funroll-loops Brute\ Force\ Extension.c -lm -lgsl -lgslcblas
*/

#ifndef __DISTORTION_H__
#define __DISTORTION_H__

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

double Rigid_Rotor_Distort (double /*A*/, double /*C*/, int /*J*/, double /*DJ*/ int /*Index*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
double Get_Frequency_Distort (int /*J_Up*/, int /*J_Low*/, int /*IndexUp*/, int /*IndexLow*/, double * /*Constants*/, double /*DJ*/, struct ETauStruct /*ETStruct*/);
int Get_Catalog_Distort (struct Transition * /*CatalogtoFill*/, double * /*Constants*/, double /*DJ*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level * /*MyDictionary*/);

double Rigid_Rotor_Distort (double A, double C, int J, double DJ, int Index, double Kappa, struct ETauStruct ETStruct)
{
//Ease of use function to compute the energy of a single rigid rotor level
//The first term is basically just total angular momentum, the second is an asymmetric rotor/state-dependent correction
//Happily frequencies/energies are computed in whatever units the rotational constants are given in
	return 0.5*(A+C)*J*(J+1.0)+0.5*(A-C)*E_tau(Index,Kappa,ETStruct)-J*(J+1.0)*DJ;
}

double Get_Frequency_Distort (int J_Up, int J_Low, int IndexUp, int IndexLow, double *Constants, double DJ, struct ETauStruct ETStruct)
{
//Ease of use function to compute the frequency of an asymmetric rotor
double Kappa;
	Kappa = Get_Kappa (Constants[0],Constants[1],Constants[2]);
	return fabs(Rigid_Rotor_Distort(Constants[0],Constants[2],J_Up,DJ,IndexUp,Kappa,ETStruct)-Rigid_Rotor_Distort(Constants[0],Constants[2],J_Low,DJ,IndexLow,Kappa,ETStruct));
}

int Get_Catalog_Distort (struct Transition *CatalogtoFill, double *Constants, double DJ, int CatLines, int Verbose, struct ETauStruct ETStruct, struct Level *MyDictionary)
{
//Utility function for calculating frequencies of a catalog
int i;	//Declaring i here because I like it, and apparently learned C pre C99Constants
	for (i=0;i<CatLines;i++) {	
		(CatalogtoFill)[i].Frequency = Get_Frequency (	MyDictionary[(CatalogtoFill)[i].Upper].J,
														MyDictionary[(CatalogtoFill)[i].Lower].J,
														(CatalogtoFill)[i].Upper,
														(CatalogtoFill)[i].Lower,
														Constants,  
														DJ,
														ETStruct
		);
		if (Verbose) printf ("%d ",(CatalogtoFill)[i].Type);		
		if (Verbose) print_Transition ((CatalogtoFill)[i],MyDictionary);
	} 
	return 1;
}

#endif /* __BRUTE_FORCE_H__ */