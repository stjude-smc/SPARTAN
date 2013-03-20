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
double* medianfilter( double* data, const int Npoints, const int TAU );


//Matlab entry point
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    //Check input arguments
    if (nlhs!=1 || nrhs!=2)
        mexErrMsgTxt( "Invalid arguments" );

	//Transpose column-major (MATLAB-style) input tino row-major for easier
    //access within C.
	mxArray* RHS[1] = { (mxArray*)prhs[0] };
	mxArray* LHS[1];
	mexCallMATLAB( 1,LHS, 1,RHS, "transpose" );
    
    //Load input matrix
	const int TAU = mxGetScalar(prhs[1]);
    const mxArray* input = LHS[0];  //MxN matrix of doubles

    mwSize Npoints = mxGetM(input);  //num rows
    mwSize Ntraces = mxGetN(input);  //num columns

	if (Npoints<TAU*2)
		mexErrMsgTxt( "Traces are too small to filter!" );
    
    //Create output matrix
    mxArray* output = LHS[0];//mxDuplicateArray(input);
    double*  p_output = mxGetPr( output );
    
    for( int i=0; i<Ntraces; ++i)
        medianfilter( p_output+(i*Npoints), Npoints, TAU );
    
	//Transpose to convert output back to column-major.
	RHS[0] = output;
	mexCallMATLAB( 1,LHS, 1,RHS, "transpose" );

    plhs[0] = LHS[0];
    return;
}
