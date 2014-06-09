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


#define	INPUT	prhs[0]
#define OUTPUT  plhs[0]


//forward function declarations (see mfilt.cpp)
double* medianfilter( double* data, const int nPoints, const int TAU );


//Matlab entry point
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    //Check input arguments
    if (nlhs!=1 || nrhs!=2)
        mexErrMsgTxt( "Invalid arguments" );

    mwSize nTraces = mxGetM(INPUT);  //number of rows
    mwSize nPoints = mxGetN(INPUT);  //number of columns
	const int TAU = mxGetScalar(prhs[1]);
    
    //mexPrintf("nTraces=%d, nPoints=%d, TAU=%d\n",nTraces,nPoints,TAU);
    
	if (nPoints<TAU*2)
		mexErrMsgTxt( "Traces are too small to filter!" );
    
    // Allocate output array memory. and copy the input data there.
    // Since MATLAB is column major and C is row-major, we need to
    // transpose for direct array access to be efficient.
    // This automatically allocates memory in OUTPUT.
    //OUTPUT = mxCreateDoubleMatrix(nTraces, nPoints, mxREAL);
    mxArray *intermediate[1];
    mexCallMATLAB( 1,intermediate, 1,(mxArray**)&INPUT, "transpose" );
    
    if( mxGetClassID(prhs[0]) != mxDOUBLE_CLASS )
        mexErrMsgTxt( "mex medianfilter only works on doubles." );
    
    //Create output matrix
    double* p_output = mxGetPr( intermediate[0] ); //get pointer to actual data.
    
    for( int i=0; i<nTraces; ++i)
       medianfilter( p_output+(i*nPoints), nPoints, TAU );
    
    // Transpose again to return data as column-major to MATLAB.
    mexCallMATLAB( 1,&OUTPUT, 1,intermediate, "transpose" );
    
    mxDestroyArray(intermediate[0]);
    return;
}






