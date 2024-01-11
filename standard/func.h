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
int calc_energy_difference(int index, int newphi);
