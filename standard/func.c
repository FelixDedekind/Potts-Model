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
    //error if sitelist not defined yet
    if(sitelist == NULL) {
        printf("error: sitelist not defined! \n");
        exit(-1);
    } else {
        int cc,dd;

        //array to save coordinates of points on grid
        int coords[dim];
        //array to save coordinates of NEIGHBOURING points on grid
        int neicoords[dim];

        int nei1,nei2;

        int rand_spin;

        for(cc = 0; cc < N; cc++) {
            //initiate label to number of site
            sitelist[cc].label = cc;
            //initiate without any bonds to neighbours
            sitelist[cc].above = -1;
            sitelist[cc].below = -1;        

            //all spins randomized
            int rand_spin = rand()%q;
            sitelist[cc].phi = rand_spin;       

            //calculate coordinates of site in n^dim space to easily find neighbours
            site_to_coords(coords, cc);         

            for(dd = 0; dd < dim; dd++) {
                //copy the coordinates of the current site to neicoords           
                neicoords[dd] = coords[dd];
            }

            for(dd = 0; dd < dim; dd++) {                      //calculate coordinates of neighbours in each direction of each dimension
                //go positive step in dimension dd
                neicoords[dd] = real_mod(coords[dd]+1,n);
                //save site number of this neighnour
                nei1 = coords_to_site(neicoords);

                //go negative step in dimension dd
                neicoords[dd] = real_mod(coords[dd]-1,n);
                //save site number of this neighbour
                nei2 = coords_to_site(neicoords);
                
                //write coordinates to respective neighbour list
                sitelist[cc].neis[2*dd] = nei1;        
                sitelist[cc].neis[2*dd+1] = nei2;

                //reset neicoords array for calculation in next dimension
                neicoords[dd] = coords[dd];             
            }
          
        }
    }
}


void randomize_phis() {             //re-randomizes spins without re-calculating neighbours
    int cc;
    int rand_spin;
    for(cc = 0; cc < N; cc++) {
            rand_spin = rand()%q;
            sitelist[cc].phi = rand_spin;
    }
}

void init_acc_rates()               //initiates list of all possible acceptance rates
{                                   //this way, we do not have to call exp in the monte carlo steps
    int ii;
    for(ii = 0; ii < nei_num+1; ii++)
    {
        acc_rates[ii]=exp(-(ii*J)/(kB*T));
    }
}

double calc_energy() {              //calculates the energy of the entire system for the standard potts model
    double energy = 0;
    int cc,dd;
    for(cc = 0; cc < N; cc++) {
        for(dd = 0; dd < nei_num; dd++) 
        {
            if(sitelist[cc].phi-sitelist[sitelist[cc].neis[dd]].phi==0)     //<--- the Hamiltonian uses a Kronecker delta
            {
                energy -=J;
            }
        }
    }
    //avoid double-counting by taking half
    return energy/2;    
}

int calc_aligned_neighbours_difference(int index, int newphi) {             //calculates difference of how many neighbours are aligned before and after a spin flip
    int n0 = 0,n1 = 0;
    int dd;
    int oldphi = sitelist[index].phi;
    for(dd = 0; dd < nei_num; dd++) {
        //aligned before flipping spin
        if(oldphi==sitelist[sitelist[index].neis[dd]].phi)  n0+=1;

        //aligned after flipping spin
        if(newphi==sitelist[sitelist[index].neis[dd]].phi)  n1+=1;
    }
    return n1-n0;
}


void try_change_spin(int index) {               //try to flip a single spin
    int rand_phi = rand()%q;

    //if more neighbours are aligned after the spin flip, dn is positive
    int dn = calc_aligned_neighbours_difference(index, rand_phi);   

    //flipable = 0 means dont flip, and the other way around
    int flipable = (int) ((dn>=0) || (acc_rates[abs(dn)] >= (double)rand()/RAND_MAX));

    //flip if flipable is 1, dont flip if it is 1
    sitelist[index].phi = flipable*rand_phi + (1-flipable) * sitelist[index].phi;       
}


void print_config() {           //2d cross section of configuration printed to terminal - works well for n<150 and q<=10 or so, any more than that and it gets janky
    //can only print 2d crosssection if dim is at least 2
    if(dim>=2) {            
        int cc,dd;
        printf("printing 2d cross section of config to terminal: \n");
        for(cc = 0; cc < n; cc++) {
            for(dd = 0; dd < n; dd++) {
                printf("%d ", sitelist[cc*n+dd].phi);
            }
            printf("\n");
        }    
    }
    
}



void mc_timestep() {                //performs a mc timestep, i.e. tries N spin flips
    int ii;
    int rand_site;
    for(ii = 0; ii < N; ii++) {
        rand_site = rand()%N;
        try_change_spin(rand_site);
    }
}


double calc_magnetization() {           // MAKE AN EXCEPTION FOR Q = 1
    //if q > 1, calculate magnetization in direction of each basis vector
    if(q>1) {               
        double mag = 0;
        //we will calculate the magnetization in direction of each basis vector and take the maximum magnetization in the end
        double maxmag = 0;      
        int cc,dd;
        for(dd=0; dd < q; dd++) {
            mag = 0;
            for(cc = 0; cc < N; cc++) {       
                if(sitelist[cc].phi == dd) 
                {
                    mag += 1.;
                } else {
                    //if all spins are random, the total magnetization is 0, thus 1./(q-1)
                    mag -= 1./(q-1);            
                }
            }
            if(mag>maxmag) maxmag=mag;
        }
        return maxmag/N;    
    } else {                // if only q=1 exists (1 direction), the magnetization per site will always be 1
        return 1;
    }
    
}

void print_labels() {           //prints a cross section of the labels to the terminal
    int cc,dd,ee;
    printf("printing 2d cross section of labels to terminal: \n");
    for(cc = 0; cc < n; cc++) {
        for(dd = 0; dd < n; dd++) {
            printf("%d ", sitelist[cc*n+dd].label);
        }
        printf("\n");
    }
}

void update_labels() {      //function to label clusters
    int replacelabel;
    int cc, dd;

    //reset labels and cluster neighbours
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
                int bot_of_top = cc;
                while(sitelist[bot_of_top].below!=-1) bot_of_top = sitelist[bot_of_top].below;

                //get to top site in lower array
                int top_of_bot = sitelist[cc].neis[dd];
                while(sitelist[top_of_bot].above!=-1) top_of_bot = sitelist[top_of_bot].above;

                //combine clusters by setting new branch neighbours
                sitelist[bot_of_top].below = top_of_bot;
                sitelist[top_of_bot].above = bot_of_top;

                //update labels
                int counterlabel = bot_of_top;
                while(sitelist[counterlabel].below!=-1) {
                    counterlabel = sitelist[counterlabel].below;
                    sitelist[counterlabel].label = replacelabel;
                }
            }
        }
    }
    //check_labeling();
}

void calc_percolation() {       //function to calculate the percolation
    //first update labels
    update_labels();        
}

void check_labeling()
{
    printf("Checking labelling\n");
    int cc, dd;
    for(cc = 0; cc < N; cc++) 
    {
        for(dd = 0; dd < nei_num; dd++) 
        {
            if(sitelist[cc].phi==sitelist[sitelist[cc].neis[dd]].phi)
            {
                if(sitelist[cc].label!=sitelist[sitelist[cc].neis[dd]].label)
                {
                    printf("Error1: %d\n",cc);
                }
            }
            if(sitelist[cc].phi!=sitelist[sitelist[cc].neis[dd]].phi)
            {
                if(sitelist[cc].label==sitelist[sitelist[cc].neis[dd]].label)
                {
                    printf("Error2: %d\n",cc);
                }
            }
        }    
    
    }
}
