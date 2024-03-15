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


void initiate_sites() {             // assign neighbours, initialize spins
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


void init_acc_rates()  {             //initiates list of all possible acceptance rates
    int ii;                         //this way, we do not have to call exp in the monte carlo steps
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
    //can only print 2d crosssection if dim is at least 2!
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


double calc_magnetization() {       // calculates magnetization
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


void print_labels() {           //prints a cross section of the labels to the terminal (i.e. clustering)
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
}


double calc_largest_cluster() {       //function to calculate the percolation
    //first update labels
    update_labels();   
    //variables
    int cc;
    int maxsize = 0;
    int * sizes;
    //list of all possible labels
    sizes = malloc(N*sizeof(int));
    for(cc = 0; cc < N; cc++) sizes[cc]=0;
    //count how often labels appear
    for(cc = 0; cc < N; cc++) sizes[sitelist[cc].label]+=1;
    //find maximum size of clusters
    for(cc = 0; cc < N; cc++) {
        if(maxsize < sizes[cc]) maxsize = sizes[cc];
    }
    free(sizes);
    //return maximum size
    return (double)maxsize/N;     
}


void check_labeling() {      //prints if any neighbouring sites have same spin but are labelled differently (as a cross-check)
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


int correlation_function_and_time(double* correlations, int test_timesteps, int samples) {      // calculates and prints correlation function, calculates correlation time
    //get average first
    randomize_phis();
    //large amount of timesteps to ensuredly get to equilibrium
    int equi_time = 30000;
    // do 10 measurements with 2500 mc timesteps inbetween
    int period = 2500;
    int measurements = 10;

    double mean = 0;
    int ii, tt;
    // relax first
    for(ii = 0; ii < equi_time; ii++) {
        mc_timestep();
    }

    // save explicit calculated magnetizations at different times in an array
    double* explicits = malloc((measurements+1)*sizeof(double)); 

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

    // re-randomize
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
        //normalize correlations to get normalized correlation function
        correlations[tt] /= correlations[1];
    }

    //calculate correlation time
    // do this using first point where C(t) < C(0)/2. This works since correlation function behaves nicely!
    int first = samples/2;
    double max = correlations[0];
    for(tt = 0; tt < samples/2; tt++) {
        if(correlations[tt]-0.5*max < 0) {
            first = tt;
            break;
        }
    }

    //print time-dependent correlation function into file
    FILE* out_corr = fopen("./correlation.txt", "w");
    for(tt = 0; tt - samples/2; tt++) fprintf(out_corr, "%i %f \n",tt*sample_period, correlations[tt]);
    fclose(out_corr);

    // free temporary arrays
    free(temp_deltam);
    free(explicits);
    return first*sample_period;
}


int get_correlation_time() {        // averages correlation times and returns mean
    printf("probing correlation time at temperature: T=%f \n", T);
    // tested time interval
    int test_timesteps = 4000; 
    int tests = 5;
    // samples taken in interval
    int samples = 401;
    double * correlations = malloc(samples/2*sizeof(double));
    double tau = 0;
    int ii;
    // average correlation time
    for(ii = 0; ii < tests; ii++) {
        tau += correlation_function_and_time(correlations,test_timesteps,samples);
    }
    tau /= tests;
    printf("average correlation time: %f \n", tau);
    free(correlations);

    return tau;
}


void test_temp(int tau, double* mag, double* E, double* lcs) {      // test temperatures for energy, magnetization, largest cluster size
    printf("probing temperature: T=%f \n", T);
    // randomize angles
    randomize_phis();
    // measurement properties
    int mc_timesteps = 5*tau;
    int measurements = 10;
    int ii,tt;
    // perform timesteps
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
    }
    // perform measurements every \tau timesteps
    for(ii = 0; ii < measurements; ii++) {
        for(tt = 0; tt < tau; tt++) {
            mc_timestep();
        }
        *E += calc_energy()/(double)measurements;
        *mag += calc_magnetization()/(double)measurements;
        *lcs += calc_largest_cluster()/(double)measurements;
    }
}


void test_temp_percolation(int tau) {       // this function was implemented but NEVER USED! It sets an arbitrary threshold to 0.7
    printf("probing temperature: T=%f \n", T);      // tests 500 configurations, if more than 70 % of sites are in one cluster, prints "1" to file ("0" else)
    randomize_phis();
    FILE* out_perc_prob = fopen("./perc_prob_over_temp.txt", "a");
    // relaxation timesteps and measurements
    int mc_timesteps = 5*tau;
    int measurements = 500;
    int ii,tt;
    //perform timesteps
    for(ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep();
    }
    // measurements and printing
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