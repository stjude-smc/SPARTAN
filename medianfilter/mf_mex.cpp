/*
MEDIANFILTER

DESCRIPTION  Linear median filter (operates on rows)
AUTHOR       Daniel Terry
COPYRIGHT    2008, Scott Blanchard
*/

using namespace std;

//#include <algorithm>
//#include <vector>
#include <iostream>
#include <string.h> //memcpy
//#include <math.h> //ceil

#include "mex.h"



//forward function declarations (see mfilt.cpp)
double* medianfilter( double* data, const int Npoints );


//Matlab entry point
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    //Check input arguments
    if (nlhs!=1 || nrhs!=1)
        mexErrMsgTxt( "Invalid arguments" );
    
    
    //Load input matrix
    const mxArray* input = prhs[0];  //MxN matrix of doubles

    mwSize Npoints = mxGetM(input);  //num rows
    mwSize Ntraces = mxGetN(input);  //num columns
    
    //if (Ntraces != 1)
    //    mexErrMsgTxt( "Only one trace at a time!!");
    
    
    //Create output matrix
    mxArray* output = mxDuplicateArray(input);
    double*  p_output = mxGetPr( output );
    
    for( int i=0; i<Ntraces; ++i)
        medianfilter( p_output+(i*Npoints), Npoints );
    
    plhs[0] = output;
    return;
}








