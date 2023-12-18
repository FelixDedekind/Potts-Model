#include "vars.h"


site *sitelist; // Define the array in this source file
double J = 1;
double kB = 1;
double T = 1;

void malloc_sitelist() {
    sitelist = malloc(N * sizeof(site)); // Allocate memory for the sitelist

    if (sitelist == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
}


void free_sitelist() {
    free(sitelist); // Free the allocated memory when done
}