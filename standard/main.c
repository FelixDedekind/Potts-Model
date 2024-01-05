#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm

double test_temp(double Temp) {
    FILE* out_en = fopen("energy_over_temp.txt", "a");
    FILE* out_mag = fopen("mag_over_temp.txt", "a");
    //printf("opened output files \n");
    T = Temp;
    //printf("set temperatures \n");
    int mc_timesteps = 20000;
    int ii;
    int naverage = 2000;
    double mag_avg = 0;
    double E_avg = 0;
    //printf("defined variables in test_temp \n");
    for(ii = 0; ii < mc_timesteps; ii++) {
        //printf("trying timestep: %d/%d \n",ii,mc_timesteps);
        mc_timestep();
        //printf("timestep done \n");
        if(ii > mc_timesteps - naverage) {
            E_avg += calc_energy()/(double)naverage;
            mag_avg += calc_magnetization()/(double)naverage;
        }
    }
    fprintf(out_en, "%f %f \n",Temp, E_avg);
    fprintf(out_mag, "%f %f \n",Temp, mag_avg);
    //printf("printed to files \n");
    fclose(out_en);
    fclose(out_mag);
    //printf("closed files \n");
}

int main() {
    srand(time(NULL));   //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list
    
    
    FILE* out_en = fopen("energy_over_temp.txt", "w");
    fclose(out_en);

    FILE* out_mag = fopen("mag_over_temp.txt", "w");
    fclose(out_mag);


    double T_init = 1.9;
    double T_final = 2.1;
    int T_steps = 100;
    int tt;
    for(tt = 0; tt < T_steps; tt++) {
        printf("Progress: %d/%d \n", tt, T_steps);
        initiate_sites();
        //printf("initiated sites \n");
        test_temp(T_init+(T_final-T_init)*((double)tt/(double)T_steps));
        //printf("tested temperature \n");
        //print_config();
    }

    /*
    FILE* out = fopen("energy_over_time.txt", "w");

    int mc_timesteps = 20000;

    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        fprintf(out, "%d %f \n", ii, calc_energy());
        //print_config();
    } 


    print_config();

    fclose(out);  */
    free_sitelist();     //frees sitelist
    return 0;
}
