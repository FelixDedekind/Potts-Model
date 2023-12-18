#ifndef macros_defined
#define macros_defined

#ifndef dim             /*dim will be the dimension of calculations. If no calculation is defined, we set to 2*/
#define dim 2
#endif

#define nei_num 2*dim



#ifndef q               //define number of possible spin configurations
#define q 4
#endif

#ifndef n
#define n 10
#endif
#define N pow(n,dim)

typedef struct site {
    int phi;              //spin *configuration* given by 0,1,2,...,q. The angle is given by phi/q*2*PI. 
    int neis[nei_num];      //list of neighbours
}site;

#endif