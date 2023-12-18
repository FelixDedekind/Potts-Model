#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm

int main() {
    mallocSitelist();   //mallocs sitelist
    initiateSites();    //initializes with spin up and writes neighbour list

    

    freeSitelist();     //frees sitelist
    return 0;
}
