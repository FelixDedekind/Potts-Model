#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm

int main() {
    srand(time(NULL));      //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list

    print_config();


    int mc_timesteps = 100;

    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        print_config();
    }

    free_sitelist();     //frees sitelist
    return 0;
}
