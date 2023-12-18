#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm

double test_temp(double Temp) {
    FILE* out = fopen("energy_over_temp.txt", "a");
    T = Temp;
    int mc_timesteps = 5000;
    int ii;
    int naverage = 1000;
    double Eavg = 0;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        if(ii > mc_timesteps - naverage) {
            Eavg += calc_energy()/(double)naverage;
        }
    }
    printf("%f ",calc_energy());
    fprintf(out, "%f %f \n",Temp, Eavg);
    fclose(out);
}

int main() {
    srand(time(NULL));   //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_energy_table();
    initiate_sites();    //initializes with spin up and writes neighbour list
    
    FILE* out = fopen("energy_over_temp.txt", "w");
    fclose(out);


    double T_init = 0.6;
    double T_final = 5;
    int T_steps = 100;
    int tt;
    for(tt = 0; tt < T_steps; tt++) {
        initiate_sites();
        test_temp(T_init+(T_final-T_init)*((double)tt/(double)T_steps));
    }


    /*FILE* out = fopen("energy_over_time.txt", "w");

    int mc_timesteps = 5000;

    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        fprintf(out, "%d %f \n", ii, calc_energy());
        //print_config();
    }


    print_config();

    fclose(out); */
    free_sitelist();     //frees sitelist
    return 0;
}
