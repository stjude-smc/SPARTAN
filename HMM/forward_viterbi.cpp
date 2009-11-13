/*

 *WARNING: this method may fail if multiple states share the same conductance!

 * TODO: raise an error OR SOMETHING if the tree cannot be saved...
 */



#include "mex.h"
#include <string.h>
#include <math.h>


//
int forwardViterbi(
        double* start_p, double* trans_p, double* emit_p, \
        const int nStates, const int nObs, \
        int*& vPath, double& vLL );  //output parameters


//
inline mxArray* callMatlab( const mxArray* input, const char* command );


//Matlab entry point
//FORMAT: [vPath, vLL] = forward_viterbi(start_p, trans_p, emit_p)
//
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
    if( nrhs!=3 )
        mexErrMsgTxt( "Not enough input arguments" );
    
    if( nlhs>2 )
        mexErrMsgTxt( "Too many output arguments" );
        
    //Verify input argument dimensions
    int M,N, nStates, nObs;
    
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    mxAssert( M>1 && N==1, "Argument 1 (p0) has wrong size" );
    nStates = M;
    
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
    mxAssert( M==nStates && N==nStates, \
                "Argument 1 (p0) and 2 (A) have incompatible sizes" );
    
    M = mxGetM(prhs[2]);
    N = mxGetN(prhs[2]);
    mxAssert( M==nStates, \
                "Argument 2 (A) and 3 (B) have incompatible sizes" );
    mxAssert( N>1, "B matrix is empty?" );
    nObs = N;
    
    //mexPrintf("Viterbi: %d states and %d datapoints\n",nStates,nObs);
    
    //Prepare matrices for use in C code
    //NOTE: we have to tell matlab to free these somehow...
    mxArray* start_p = mxDuplicateArray( prhs[0] );
    mxArray* trans_p = callMatlab(prhs[1],"transpose");
    mxArray* emit_p = mxDuplicateArray( prhs[2] );
    
    mxAssert( start_p>0 && trans_p>0 && emit_p>0, "alloc error" );
    
    //Run Viterbi
    int* vPath = 0;
    double  vLL = 987654;
    
    forwardViterbi( mxGetPr(start_p), mxGetPr(trans_p), mxGetPr(emit_p), \
                   nStates, nObs, vPath, vLL );
    vLL = vLL/nObs;
    
    //Pass the results as output arguments to MATLAB
    if( nlhs > 0 )
    {
        plhs[0] = mxCreateNumericMatrix(nObs,1,mxINT32_CLASS,mxREAL);
        mxAssert( sizeof(int)==mxGetElementSize(plhs[0]), \
                   "size assignment mismatch..." );
        if( vPath!=0 )
           memcpy( mxGetData(plhs[0]), vPath, nObs*sizeof(int) );
    }
    
    if( nlhs > 1 )
        plhs[1] = mxCreateDoubleScalar(vLL);
    
    
    //Cleanup
    mxDestroyArray(start_p);
    mxDestroyArray(trans_p);
    mxDestroyArray(emit_p);
    if( vPath!=0 )
       mxFree(vPath);
}

//
inline mxArray* callMatlab( const mxArray* input, const char* command )
{
    mxArray *rhs[1] = {(mxArray*)input};
    mxArray *output[1];
    
    mexCallMATLAB(1, output, 1, rhs, command);
    
    return output[0];
}

//
inline double valArgMax( const double* pArray, const int N, int* index )
{
    const double* start = pArray;
    const double* end   = start+N;
    double val = *pArray;
    *index = 1;
    
    while( ++pArray < end )
    {
        if( *pArray>val )
        {
            *index = pArray-start+1;
            val = *pArray;
        }
    }
    
    return val;
}

/*
//
inline void vecAddTo( double* a, double* b, const int N )
{
    const double* end = a+N;
    
    do
    {
        *a+=*b;
        ++b;
    } while( ++a < end );
        
}

inline void vecAddTo( double* a, const double b, const int N )
{
    const double* end = a+N;
    
    do
    {
        *a+=b;
    } while( ++a < end );
        
}

//
inline void vecCopy( double* a, double* b, const int N )
{
    const double* end = a+N;
    
    do
    {
       *a=*b;
       ++b;
    } while( ++a < end );
}
*/

//allocate vPath to nObs before starting...
int forwardViterbi(
        double* lsp, double* ltp, double* lep, \
        const int nStates, const int nObs, \
        int*& vPath, double& vLL )
{
    // Maximal probability (and best path) that ends in state=i at time=t
    // "If I am here, by what route is it most likely I arrived?"
    
    
    double* delta; //[nObs][nStates]; //partial probabilities
    int   * psi;   //[nObs][nStates]; //back pointers (most likely previous state)
    double* pCurr;
    int argmax, i,j,t;
        
    delta = (double*)mxCalloc(nObs*nStates,sizeof(double));
    pCurr = (double*)mxCalloc(nStates,sizeof(double));
    psi   = (int*)mxCalloc(nObs*nStates,sizeof(int));
    vPath = (int*)mxCalloc(nObs,sizeof(int));
    
    mxAssert(delta&&pCurr&&psi&&vPath,"alloc error");
    
    memset( delta, 0, sizeof(double)*nObs*nStates );
    memset( pCurr, 0, sizeof(double)*nStates );
    memset( psi,   0, sizeof(int)*nObs*nStates );
    memset( vPath, 0, sizeof(int)*nObs ); 
        
    //log transform input arguments
    double* ptr;
    
    mxAssert(lsp&&ltp&&lep,"bad args");
    
    ptr = lsp+nStates;
    while( (--ptr)>=lsp )
        *ptr = (*ptr>0 ? log(*ptr) : -1.0e10);
    
    ptr = ltp+(nStates*nStates);
    while( (--ptr)>=ltp )
        *ptr = (*ptr>0 ? log(*ptr) : -1.0e10);
    
    ptr = lep+(nObs*nStates);
    while( (--ptr)>=lep )
        *ptr = (*ptr>0 ? log(*ptr) : -1.0e10);
        
        
    // Initialization
    //delta(:,1) = lsp + lep(:,1);
//     vecCopy( (double*)delta, lsp, nStates );
//     vecAddTo( (double*)delta, lep, nStates );
    for( i=0; i<nStates; ++i )
       delta[i] = lsp[i] + lep[i];
    
    // Induction: calcualte probability at each timepoint
    for( t=1; t<nObs; ++t )
    {
        //delta_t(j) = MAX_i[ delta_t-1(i) * A(i,j) ]   * B(j)

        for( j=0; j<nStates; ++j )  //current state
        {
            //pCurr = delta(:,t-1) + ltp(:,j) + lep(j,t);
            //[delta(j,t),psi(j,t)] = max(pCurr);
            
            for( i=0; i<nStates; ++i )
                pCurr[i] = delta[(t-1)*nStates+i] + \
                           ltp[i*nStates+j] + lep[t*nStates+j];
            
            //partial prob. of state=j at time=t (delta);
            //most likely state (psi)
            delta[t*nStates+j] = valArgMax( pCurr, nStates, &psi[t*nStates+j] );
        }
    }
    //mexPrintf("done. ");
    
     
    // Termination
    vLL = valArgMax( delta+nStates*(nObs-1), nStates, &argmax );
    //tLL = sum(delta(:));
    
    //mexPrintf("after valmax. ");

    // Backtrace to determine most likely path
    vPath[nObs-1] = argmax;

    for( t=nObs-2; t>=0; --t )
        vPath[t] = psi[ (t+1)*nStates + vPath[t+1]-1 ];
    
    //cleanup
    mxFree( delta );
    mxFree( pCurr );
    mxFree( psi );
    //mexPrintf("really done.\n");
    return 1;
}
