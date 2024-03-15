#include "vars.h"

void site_to_coords(int* coords, int sitenum);
int coords_to_site(int* coords);
void initiate_sites();
//void initiate_energy_table();
double calc_energy();
void randomize_phis();
void try_change_spin(int index);
void print_config();
void mc_timestep();
void init_acc_rates();
double calc_magnetization();
int calc_aligned_neighbours_difference(int index, int newphi);
void print_labels();
void update_labels();
double calc_largest_cluster();
void check_labeling();
int correlation_function_and_time(double* correlations, int test_timesteps, int samples);
int get_correlation_time();
void test_temp(int tau, double* mag, double* E, double* lcs);
void test_temp_percolation(int tau);