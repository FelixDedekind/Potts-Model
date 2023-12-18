#ifndef variables_defined
#define variables_defined

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "macros.h"

#define PI 3.141592653589793238

extern site *sitelist; // Declare the array as an external variable
extern double J;           //this is the constant of energy in the hamiltonian
extern double kB;           //boltzmann constant
extern double T;            //temperature

void malloc_sitelist();
void free_sitelist();

#endif
