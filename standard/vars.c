#include "vars.h"


site *sitelist; // Define the array in this source file
double *acc_rates;
double J = 2;
double kB = 1;
double T = 1.9;

void malloc_sitelist() {
    sitelist = malloc(N * sizeof(site)); // Allocate memory for the sitelist
    acc_rates = malloc((nei_num+1) * sizeof(double));

    if (sitelist == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
    if (acc_rates == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
}


void free_sitelist() {
    free(sitelist); // Free the allocated memory when done
    free(acc_rates);
}