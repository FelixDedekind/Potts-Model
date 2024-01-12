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
            sitelist[cc].label = cc;
            sitelist[cc].above = -1;
            sitelist[cc].below = -1;

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
                sitelist[cc].neis[2*dd] = nei1;        //and write coordinates to respective neighbour list
                sitelist[cc].neis[2*dd+1] = nei2;
            }
          
        }
    }
}


void randomize_phis() {
    int cc;
    int rand_spin;
    for(cc = 0; cc < N; cc++) {
            rand_spin = rand()%q;
            sitelist[cc].phi = rand_spin;
    }
}

void init_acc_rates()
{
    int ii;
    for(ii = 0; ii < nei_num+1; ii++)
    {
        acc_rates[ii]=exp(-(ii*J)/(kB*T));
    }
}

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

int calc_energy_difference(int index, int newphi) {
    int n0 = 0,n1 = 0;
    int dd;
    int oldphi = sitelist[index].phi;
    for(dd = 0; dd < nei_num; dd++) {
        if(oldphi==sitelist[sitelist[index].neis[dd]].phi)
        {
            n0+=1;
        }
        if(newphi==sitelist[sitelist[index].neis[dd]].phi)
        {
            n1+=1;
        }
    }
    return n1-n0;
}


void try_change_spin(int index) {
    int rand_phi = rand()%q;
    int de = -calc_energy_difference(index, rand_phi);
    int flipable = (int) ((de<=0) || (acc_rates[de] >= (double)rand()/RAND_MAX));
    sitelist[index].phi = flipable*rand_phi + (1-flipable) * sitelist[index].phi;
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
    double maxmag = 0;
    int cc,dd;
    for(dd=0; dd < q; dd++) {
        mag = 0;
        for(cc = 0; cc < N; cc++) {       
            if(sitelist[cc].phi == dd) 
            {
                mag += (double)1;
            } else {
                mag -= (double)1/(q-1);
            }
        }
        if(mag>maxmag) maxmag=mag;
    }
    return maxmag/N;
}

void print_labels() {
    int cc,dd;
    printf("printing 2d cross section of labels to terminal: \n");
    for(cc = 0; cc < n; cc++) {
        for(dd = 0; dd < n; dd++) {
            printf("%d ", sitelist[cc*n+dd].label);
        }
        printf("\n");
    }
}

void update_labels() {
    int searchlabel;
    int replacelabel;
    int cc, dd;

    //reset labels
    for(cc = 0; cc < N; cc++) {
            sitelist[cc].label = cc;
            sitelist[cc].above = -1;
            sitelist[cc].below = -1;
    }
    

    //check for connections
    for(cc = 0; cc < N; cc++) {
        for(dd = 0; dd < nei_num; dd++) {
            if(sitelist[cc].label==sitelist[sitelist[cc].neis[dd]].label) {
                continue;
            }
            if(sitelist[cc].phi==sitelist[sitelist[cc].neis[dd]].phi)  {
                //define label to replace other labels with
                replacelabel = sitelist[cc].label;

                //get to lowest site in upper array
                int bottom_of_up = cc;
                while(sitelist[bottom_of_up].below!=-1) {
                    bottom_of_up = sitelist[bottom_of_up].below;
                } 

                //get to top site in lower array
                int top_of_bot = sitelist[cc].neis[dd];
                sitelist[top_of_bot].label = replacelabel;
                while(sitelist[top_of_bot].above!=-1) {
                    top_of_bot = sitelist[top_of_bot].above;
                }

                //combine clusters by setting new branch neighbours
                sitelist[bottom_of_up].below = top_of_bot;
                sitelist[top_of_bot].above = bottom_of_up;

                //update labels
                int counterlabel = bottom_of_up;
                while(sitelist[counterlabel].below!=-1) {
                    counterlabel = sitelist[counterlabel].below;
                    sitelist[counterlabel].label = replacelabel;
                }
            }
        }
    }
}

void calc_percolation() {
    update_labels();
}