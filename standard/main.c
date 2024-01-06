#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm -O2

double test_temp(double Temp) {
    FILE* out_en = fopen("energy_over_temp.txt", "a");
    FILE* out_mag = fopen("mag_over_temp.txt", "a");
    T = Temp;
    int mc_timesteps = 15000;
    int ii;
    int naverage = 2000;
    double mag_avg = 0;
    double E_avg = 0;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        if(ii > mc_timesteps - naverage) {
            E_avg += calc_energy()/(double)naverage;
            mag_avg += calc_magnetization()/(double)naverage;
        }
    }
    fprintf(out_en, "%f %f \n",Temp, E_avg);
    fprintf(out_mag, "%f %f \n",Temp, mag_avg);
    fclose(out_en);
    fclose(out_mag);
}

int main() {
    srand(time(NULL));   //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list
    print_config();

    
    FILE* out_en = fopen("energy_over_temp.txt", "w");
    fclose(out_en);

    FILE* out_mag = fopen("mag_over_temp.txt", "w");
    fclose(out_mag);


    double T_init = 1.4;
    double T_final = 2.5;
    int T_steps = 50;
    int runs = 5;
    int tt;
    for(int cc = 0; cc < runs; cc++) {
        for(tt = 0; tt < T_steps; tt++) {
            printf("Progress: %d/%d: %d/%d \n", cc, runs, tt, T_steps);
            initiate_sites();
            test_temp(T_init+(T_final-T_init)*((double)tt/(double)T_steps));
            print_config();
            printf("mag = %f \n",calc_magnetization());
        }
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
