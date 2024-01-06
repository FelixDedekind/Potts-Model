#include "func.h"
#include "vars.h"

int real_mod(int a, int b) {            //mod operation for positive and negative numbers
    return (a%b+b)%b;
}

void site_to_coords(int* coords, int sitenum) {     //writes coordinates of site-number in n^dim grid into coords -> important for finding neighbours
    int rest = sitenum;
    int dd;


    for(dd = 0; dd < dim-1; dd++) {
        coords[dd] = rest / (int)pow(n,dim-dd-1);
        rest %= (int)pow(n,dim-dd-1);
    }
    coords[dim-1] = rest;
}

int coords_to_site(int* coords) {     //converts coordinates in n^dim grid into site-number
    int sitenum = 0;
    for(int dd = 0; dd < dim; dd++) {
        sitenum += coords[dd]*pow(n,dim-dd-1);
    }
    return sitenum;
}


void initiate_sites() {
    if(sitelist == NULL) {
        printf("error: sitelist not defined! \n");
        exit(-1);
    } else {
        int cc,dd;
        int coords[dim];
        int neicoords[dim];
        int nei1,nei2;
        int rest = 0;
        int rand_spin;
        for(cc = 0; cc < N; cc++) {
            int rand_spin = rand()%q;
            sitelist[cc].phi = rand_spin;       //all spins randomized

            site_to_coords(coords, cc);

            for(dd = 0; dd < dim; dd++) {           //first copy the coordinates
                neicoords[dd] = coords[dd];
            }

            for(dd = 0; dd < dim; dd++) {           //then calculate coordinates of neighbours in each direction
                neicoords[dd] = real_mod(coords[dd]+1,n);
                nei1 = coords_to_site(neicoords);
                neicoords[dd] = real_mod(coords[dd]-1,n);
                nei2 = coords_to_site(neicoords);
                printf("i: %d, nei1: %d, nei2: %d \n",cc,nei1,nei2);
                sitelist[cc].neis[2*dd] = nei1;        //and write coordinates to respective neighbour list
                sitelist[cc].neis[2*dd+1] = nei2;
            }
          
        }
    }
}


//void randomize_phis() <-- to be added


double calc_energy() {
    double energy = 0;
    int cc,dd;
    for(cc = 0; cc < N; cc++) {
        for(dd = 0; dd < nei_num; dd++) 
        {
            if(sitelist[cc].phi-sitelist[sitelist[cc].neis[dd]].phi==0)
            {
                energy -=J;
            }
        }
    }
    return energy/2;    //avoid double-counting
}

double calc_energy_difference(int index, int newphi) {
    double e0 = 0,e1 = 0;
    int dd;
    for(dd = 0; dd < nei_num; dd++) {
        if(sitelist[index].phi==sitelist[sitelist[index].neis[dd]].phi)
        {
            e0-=J;
        }
    }
    int oldphi = sitelist[index].phi;
    sitelist[index].phi = newphi;
    for(dd = 0; dd < nei_num; dd++) {
        if(sitelist[index].phi==sitelist[sitelist[index].neis[dd]].phi)
        {
            e1-=J;
        }
    }
    sitelist[index].phi = oldphi;
    return e1-e0;
}


void try_change_spin(int index) {
    int rand_phi = rand()%q;
    double de = calc_energy_difference(index, rand_phi);
    int flipable = 1;
    if(de>0.)
    {
        if(exp(-de/(kB*T)) < (double)rand()/RAND_MAX) flipable = 0;           // make shorter
    }
    if(flipable==1) sitelist[index].phi = rand_phi;
}


void print_config() {           //2d cross section of configuration printed to terminal
    int cc,dd;
    printf("printing 2d cross section of config to terminal: \n");
    for(cc = 0; cc < n; cc++) {
        for(dd = 0; dd < n; dd++) {
            printf("%d ", sitelist[cc*n+dd].phi);
        }
        printf("\n");
    }
}



void mc_timestep() {
    int ii;
    int rand_site;
    for(ii = 0; ii < N; ii++) {
        rand_site = rand()%N;
        if(rand_site<0 || rand_site>=N) printf("BUNGER BUNGER BUNGER BUNGER BUNGER!!!! rand_site: %d \n \n", rand_site);
        try_change_spin(rand_site);
    }
}


double calc_magnetization() {           // this doesnt make a lot of sense
    double mag = 0;
    int cc;
    for(cc = 0; cc < N; cc++) {
        if(sitelist[cc].phi == 0) 
        {
            mag += 1;
        } else {
            mag -= 1/(q-1);
        }
    }
    return mag;
}