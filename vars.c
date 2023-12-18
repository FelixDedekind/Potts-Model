#include "vars.h"


site *sitelist; // Define the array in this source file

void mallocSitelist() {
    sitelist = malloc(N * sizeof(site)); // Allocate memory for the sitelist

    if (sitelist == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
}


void freeSitelist() {
    free(sitelist); // Free the allocated memory when done
}