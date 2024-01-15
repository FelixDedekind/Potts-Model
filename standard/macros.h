#ifndef macros_defined
#define macros_defined

#ifndef dim             /*dim will be the dimension of calculations. If no calculation is defined, we set to 2*/
#define dim 2           
#endif 

#define nei_num 2*dim



#ifndef q               
#define q 10
#endif

#ifndef n
#define n 50
#endif
#define N (int)pow(n,dim)

typedef struct site {
    int phi;                //is int in {0,1,...,q-1}
    int neis[nei_num];      //list of neighbours
    int label;              //label for percolation
    int above, below;       //bond pointers for percolation
}site;

#endif