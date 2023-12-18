#include "vars.h"

void site_to_coords(int* coords, int sitenum);
int coords_to_site(int* coords);
void initiate_sites();
double calc_energy();
void try_change_spin(int index);
void print_config();
void mc_timestep();