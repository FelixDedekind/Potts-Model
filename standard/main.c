#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm -O2

void test_temp(double Temp, double* mag, double* E, double* lcs) {
    //setting global temperature
    T = Temp;
    //number of timesteps and number of timesteps to average over
    int mc_timesteps = 15000;
    int naverage = 2000;
    //perform timesteps
    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        if(ii > mc_timesteps - naverage) {
            *E += calc_energy()/(double)naverage;
            *mag += calc_magnetization()/(double)naverage;
            *lcs += calc_largest_cluster()/(double)naverage;
        }
    }
}

int main() {
    srand(time(NULL));   //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list
    init_acc_rates();    //initialize array of possible acceptance rates

    double T_init = 1.;
    double T_final = 2.5;
    int T_steps = 50;
    int runs = 5;

    double mags[runs];
    double Es[runs];
    double lcs[runs];

    FILE* out_en = fopen("./energy_over_temp.txt", "w");
    FILE* out_mag = fopen("./mag_over_temp.txt", "w");
    FILE* out_lcs = fopen("./lcs_over_temp.txt", "w");

    int tt,cc,jj;
    for(tt = 0; tt < T_steps; tt++) {
        for(cc = 0; cc < runs; cc++) {
            printf("Progress: %d/%d: %d/%d \n", tt, T_steps, cc, runs);
            randomize_phis();
            init_acc_rates();
            mags[cc] = 0;
            Es[cc] = 0;
            lcs[cc] = 0;
            test_temp(T_init+(T_final-T_init)*((double)tt/(double)T_steps),&mags[cc],&Es[cc],&lcs[cc]);
            print_config();
            printf("mag = %f \n",mags[cc]);
        }
        double mean_mag = 0;
        double std_mag = 0;
        for(jj = 0; jj < runs; jj++) mean_mag += mags[jj]/runs;
        for(jj = 0; jj < runs; jj++) std_mag += pow(mags[jj]-mean_mag,2)/N;
        std_mag = sqrt(std_mag);
        fprintf(out_mag,"%f %f %f \n", T, mean_mag, std_mag);


        double mean_E = 0;
        double std_E = 0;
        for(jj = 0; jj < runs; jj++) mean_E += Es[jj]/runs;
        for(jj = 0; jj < runs; jj++) std_E += pow(Es[jj]-mean_E,2)/N;
        std_E = sqrt(std_E);
        fprintf(out_en,"%f %f %f \n", T, mean_E, std_E);

        double mean_lcs = 0;
        double std_lcs = 0;
        for(jj = 0; jj < runs; jj++) mean_lcs += lcs[jj]/runs;
        for(jj = 0; jj < runs; jj++) std_lcs += pow(lcs[jj]-mean_lcs,2)/N;
        std_lcs = sqrt(std_lcs);
        fprintf(out_lcs,"%f %f %f \n", T, mean_lcs, std_lcs);
    }
     
    fclose(out_en);
    fclose(out_mag);
    fclose(out_lcs);

    /*
    FILE* out = fopen("energy_over_time.txt", "w");

    int mc_timesteps = 10000;

    int ii;
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
        fprintf(out, "%d %f \n", ii, calc_energy());
    } 
    print_config();
    update_labels();
    print_labels();



    fclose(out); */
    free_sitelist();     //frees sitelist
    return 0;
}
