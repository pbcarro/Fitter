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
#include "Fitter.h"
#include "Brute Force Extension.h"


//=============Functions========================
int main (int argc, char *argv[])
{

	return 1;
}



