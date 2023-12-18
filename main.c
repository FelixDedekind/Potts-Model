#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm

int main() {
    srand(time(NULL));      //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list
    initiate_energy_table();
    print_config();


    FILE* out = fopen("energy_over_time.txt", "w");

    int mc_timesteps = 100;

    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        fprintf(out, "%d %f \n", ii, calc_energy());
        print_config();
    }


    fclose(out);
    free_sitelist();     //frees sitelist
    return 0;
}
