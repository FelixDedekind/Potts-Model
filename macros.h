#ifndef macros_defined
#define macros_defined

#ifndef dim             /*dim will be the dimension of calculations. If no calculation is defined, we set to 2*/
#define dim 2
#endif

#define size 2*dim



#ifndef q
#define q 4
#endif

#ifndef n
#define n 10
#endif
#define N pow(n,dim)

typedef struct site {
    float phi;
    int neis[size];
}site;

#endif