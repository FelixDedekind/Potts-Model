#include "vars.h"
#include "func.h"

//compile with gcc -o pottsmodel main.c vars.c func.c -lm -O2

int correlation_function_and_time(double* correlations, int test_timesteps, int samples) {
    //get average and std first
    randomize_phis();
    int equi_time = 15000;
    int period = 2500;
    int measurements = 10;
    double mean = 0;
    double* explicits = malloc((measurements+1)*sizeof(double));
    int ii, tt;
    for(ii = 0; ii < equi_time; ii++) {
        mc_timestep();
    }
    explicits[0] = calc_magnetization();
    mean += explicits[0];
    for(ii = 0; ii < measurements; ii++) {
        for(tt = 0; tt < period; tt++) {
            mc_timestep();
        }
        explicits[ii] = calc_magnetization();
        mean += explicits[ii];
    }
    mean /= (measurements+1);
    randomize_phis();

    //read out fluctuations at equidistant times
    int sample_period = test_timesteps/(samples-1);
    double * temp_deltam = malloc(samples*sizeof(double));
    int jj = 0;
    for(ii = 0; ii < test_timesteps; ii++) {
        mc_timestep();
        if(ii == jj*sample_period) {
            temp_deltam[jj] = (calc_magnetization() - mean);
            jj++;
        }
    }

    //calculate correlation function
    double temporary_correlation;
    for(tt = 0; tt < samples/2; tt++) {
        temporary_correlation = 0.;
        for(ii = 0; ii < samples-tt; ii++) {
            temporary_correlation += temp_deltam[ii]*temp_deltam[ii+tt];
        }
        correlations[tt] = (double)temporary_correlation/(samples-tt+1);
        correlations[tt] /= correlations[0];
    }

    //calculate correlation time
    int closest = samples/2;
    double max = correlations[0];
    for(tt = 0; tt < samples/2; tt++) {
        if(correlations[tt]-0.5*max < 0) {
            closest = tt;
            break;
        }
    }

    //print time-dependent correlation function into file
    FILE* out_corr = fopen("./correlation.txt", "w");
    for(tt = 0; tt - samples/2; tt++) fprintf(out_corr, "%i %f \n",tt*sample_period, correlations[tt]);
    fclose(out_corr);
    free(temp_deltam);
    free(explicits);
    return closest*sample_period;
}

int get_correlation_time() {
    printf("probing correlation time at temperature: T=%f \n", T);
    int test_timesteps = 4000; 
    int tests = 5;
    int samples = 201;
    double * correlations = malloc(samples/2*sizeof(double));
    double tau = 0;
    int ii;
    for(ii = 0; ii < tests; ii++) {
        tau += correlation_function_and_time(correlations,test_timesteps,samples);
    }
    tau /= tests;
    printf("average correlation time: %f \n", tau);
    free(correlations);

    return tau;
}

void test_temp(int tau, double* mag, double* E, double* lcs) {
    printf("probing temperature: T=%f \n", T);
    randomize_phis();
    int mc_timesteps = 5*tau;
    int measurements = 10;
    int ii,tt;
    //perform timesteps
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
    }
    for(ii = 0; ii < measurements; ii++) {
        for(tt = 0; tt < tau; tt++) {
            mc_timestep();
        }
        *E += calc_energy()/(double)measurements;
        *mag += calc_magnetization()/(double)measurements;
        *lcs += calc_largest_cluster()/(double)measurements;
    }
}

void test_temp_percolation(int tau) {
    printf("probing temperature: T=%f \n", T);
    randomize_phis();
    FILE* out_perc_prob = fopen("./perc_prob_over_temp.txt", "a");
    int mc_timesteps = 5*tau;
    int measurements = 1000;
    int ii,tt;
    //perform timesteps
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
    }
    for(ii = 0; ii < measurements; ii++) {
        for(tt = 0; tt < tau; tt++) {
            mc_timestep();
        }
        if(calc_largest_cluster()>=cluster_threshold) {
            fprintf(out_perc_prob,"1");
        } else {
            fprintf(out_perc_prob,"0");
        }
    }
    fprintf(out_perc_prob,"\n");
    fclose(out_perc_prob);
}

int main() {
    srand(time(NULL));   //seed rng
    malloc_sitelist();   //mallocs sitelist
    initiate_sites();    //initializes with spin up and writes neighbour list
    init_acc_rates();    //initialize array of possible acceptance rates

    double phase_transition = 1/(log(1+sqrt(q)));
    printf("q = %i \n", q);
    printf("phase transition predicted at %f \n", phase_transition);
    double T_init = phase_transition-0.3-0.01;
    double T_final = phase_transition+0.3+0.01;
    int T_steps = 20;
    int runs = 1;
    int tau = 0;

    double mags[runs];
    double Es[runs];
    double lcs[runs];

    FILE* out_en = fopen("./energy_over_temp.txt", "w");
    FILE* out_mag = fopen("./mag_over_temp.txt", "w");
    FILE* out_lcs = fopen("./lcs_over_temp.txt", "w");
    FILE* out_perc_prob = fopen("./perc_prob_over_temp.txt", "w");
    fclose(out_perc_prob);


    int tt,cc,jj;
    for(tt = 0; tt < T_steps; tt++) {
        T = T_init+(T_final-T_init)*((double)tt/(double)T_steps);
        init_acc_rates();
        tau = get_correlation_time();
        if(tau < 500.0) tau = 500.0;
        for(cc = 0; cc < runs; cc++) {
            printf("Progress: %d/%d: %d/%d \n", tt, T_steps, cc, runs);
            randomize_phis();
            mags[cc] = 0;
            Es[cc] = 0;
            lcs[cc] = 0;
            //test_temp(tau, &mags[cc],&Es[cc],&lcs[cc]);
            test_temp_percolation(tau);
            //print_config();
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
    T = 1.1;
    int test_timesteps = 7000; 
    int samples = 201;
    double * correlations = malloc(samples/2*sizeof(double));
    correlation_function_and_time(correlations,test_timesteps,samples);
    */

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
