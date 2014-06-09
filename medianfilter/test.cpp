#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define TAU 9
#define SZ_DATA 25 //1001*15000

//forward function definitions
double* medianfilter( double* data, const int Npoints, const int szWindow );


int main()
{
    int i;
    double data[SZ_DATA];
    
    //generate "random" data for input
    srand(42);
    
    for(i=0; i<SZ_DATA; ++i)
        data[i] = rand() % 100 + 1;
        //data[i] = SZ_DATA-i;
    
    //median filter the data
    medianfilter(data, SZ_DATA, TAU);
    
    for(i=0; i<SZ_DATA; ++i)
        printf("%.0f ",data[i]);
    
    return 0;
}












