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


void initiateSites() {
    if(sitelist == NULL) {
        printf("error: sitelist not defined! \n");
        exit(-1);
    } else {
        int cc,dd;
        int coords[dim];
        int neicoords[dim];
        int nei1,nei2;
        int rest = 0;
        for(cc = 0; cc < N; cc++) {
            sitelist[cc].phi = 0;       //all spins pointing up

            site_to_coords(coords, cc);

            for(dd = 0; dd < dim; dd++) {           //first copy the coordinates
                neicoords[dd] = coords[dd];
            }

            for(dd = 0; dd < dim; dd++) {
                neicoords[dd] = real_mod(coords[dd]+1,n);
                nei1 = coords_to_site(neicoords);
                neicoords[dd] = real_mod(coords[dd]-1,n);
                nei2 = coords_to_site(neicoords);
                neicoords[dd] = coords[dd];
                sitelist[cc].neis[2*dd] = nei1;
                sitelist[cc].neis[2*dd+1] = nei2;
            }
          
        }
    }
}

