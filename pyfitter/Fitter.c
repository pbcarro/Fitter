/*
Built by Brandon Carroll 4/7/19 All Rights Reserved

Generic version of the fitting program

Build Command:
gcc -Wall -o Go Fitter.c -lm -lgsl -lgslcblas -O3 -funroll-loops
gcc -Wall -o Go Fitter.c -lm -O3 -funroll-loops

Ctypes Shared Lib:
gcc -Wall -o Fitter.so -shared -fPIC -O3 -funroll-loops Fitter.c -lm -lgsl -lgslcblas
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
// #include "GSL_Include/vector/gsl_vector.h"
// #include "GSL_Include/blas/gsl_blas.h"
// #include "GSL_Include/matrix/gsl_matrix.h"
// #include "GSL_Include/multifit_nlinear/gsl_multifit_nlinear.h"
// #include "GSL_Include/linalg/gsl_linalg.h"
#include "Fitter.h"


//Functions
int main (int argc, char *argv[])
{
	//Test_SBFIT();
	//Timing_Test();
	Accuracy_Test();
	//Timing_Test_Triples();
		
	return 1;
}

