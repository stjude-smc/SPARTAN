/*

 *WARNING: this method may fail if multiple states share the same conductance!

 * TODO: raise an error OR SOMETHING if the tree cannot be saved...
 */



#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>


// Macros for 2-dimensional array access. GETO for observations (time
// across the rows) and GETS for states. A is the matrix variable name.
#define GETS(A,I,J)  A[(I) + (J)*nStates]


// Forward function definitions.
double forwardViterbi(  double* lsp, double* ltp, double* lep,          \
                        const int nStates, const int nObs, int* vPath  );


//Matlab entry point
//FORMAT: [vPath, vLL] = forward_viterbi(start_p, trans_p, emit_p)
//
//   start_p is a  Sx1 vector of state initial probabilities.
//   trans_p is an SxS matrix of transition probabilities.
//   emit_p  is a  NxS matrix of emission probabilities for each state
//                     at each time (N observations).
//
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
    mwSize nStates, nObs;
    mxArray* vPath; //optimal state sequence (idealization)
    double vLL;     //optimal path log-likelihood score
    
    //Verify input argument dimensions
    if( nrhs!=3 )
        mexErrMsgTxt( "Not enough input arguments" );
    
    if( nlhs>2 )
        mexErrMsgTxt( "Too many output arguments" );
    
    nStates = mxGetM(prhs[0]);
    if( nStates<=1 || mxGetN(prhs[0])!=1 )
        mexErrMsgTxt( "Argument 1 (p0) has wrong size" );
    if( !mxIsDouble(prhs[0]) )
        mexErrMsgTxt( "Argument 1 (p0) must be of type double" );
    
    if( mxGetM(prhs[1])!=nStates || mxGetN(prhs[1])!=nStates )
        mexErrMsgTxt( "Argument 1 (p0) and 2 (A) have incompatible sizes" );
    if( !mxIsDouble(prhs[1]) )
        mexErrMsgTxt( "Argument 2 (A) must be of type double" );
    
    nObs = mxGetN(prhs[2]);
    if( mxGetM(prhs[2])!=nStates)
    	mexErrMsgTxt( "Argument 2 (A) and 3 (B) have incompatible sizes" );
    if( nObs<=1 )
        mexErrMsgTxt( "B matrix is empty?" );
    if( !mxIsDouble(prhs[2]) )
        mexErrMsgTxt( "Argument 3 (B) must be of type double!" );
    
    //mexPrintf("Viterbi: %d states and %d datapoints\n",nStates,nObs);
    
    
    //Allocate output memory for viterbi path (idealization).
    vPath = mxCreateNumericMatrix(1,nObs,mxINT32_CLASS,mxREAL);
    
    //Run Viterbi
    vLL = forwardViterbi( mxGetPr(prhs[0]), mxGetPr(prhs[1]), \
                mxGetPr(prhs[2]), nStates, nObs, (int*)mxGetData(vPath) );
    vLL = vLL/nObs;
    
    
    //Pass the results as output arguments to MATLAB
    if( nlhs > 0 )
        plhs[0] = vPath;
    
    if( nlhs > 1 )
        plhs[1] = mxCreateDoubleScalar(vLL);
}



//
inline double valArgMax( const double* pArray, const int N, int* index )
{
    //This function returns a 1-based index into the array that has the
    //maximum value.
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


//allocate vPath to nObs before starting...
// Maximal probability (and best path) that ends in state=i at time=t
// "If I am here, by what route is it most likely I arrived?"
double forwardViterbi(  double* lsp, double* ltp, double* lep,          \
                        const int nStates, const int nObs, int* vPath  )
{
    int endState, i,j,t;  //iterators. i is "from", j is "to", t is "time".
    double vLL;  //optimal path's log likelihood
    
    //[nStates]: //probabilities of each state at current time
    double* pCurr = (double*) mxCalloc( nStates, sizeof(double) );
    //[nObs][nStates]; //partial probabilities
    double* delta = (double*) mxCalloc( nObs*nStates, sizeof(double) );
    //[nObs][nStates]; //back pointers (most likely previous state)
    int* psi = (int*) mxCalloc( nObs*nStates, sizeof(int) );
    
    
    if( lsp==0 || ltp==0 || lep==0 )
        mexErrMsgTxt("bad args");
    if( delta==0 || psi==0 || vPath==0 )
        mexErrMsgTxt("alloc error");
    
    
    //mexPrintf("zeroing temporary memory. ");
    memset( delta, 0, sizeof(double)*nObs*nStates );
    memset( pCurr, 0, sizeof(double)*nStates );
    memset( psi,   0, sizeof(int)*nObs*nStates );
    memset( vPath, 0, sizeof(int)*nObs ); 
    
    
    // Initialization
    //delta(:,1) = lsp + lep(:,1);
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
            
            // How likely am I to traverse all the way to step t-1,
            // then transition to j, and see the current observation?
            for( i=0; i<nStates; ++i )  //next state
            {
                pCurr[i] = GETS(delta,i,t-1) + GETS(ltp,i,j) + GETS(lep,j,t);
                
                if(pCurr[i]>0)
                {
                    mexPrintf("Delta=%f, ltp=%f, lep=%f\n",GETS(delta,i,t-1), GETS(ltp,i,j), GETS(lep,j,t) );
                    mexPrintf("State: t=%d,j=%d,i=%d, LL=%f\n",t,j,i,pCurr[i]);
                    mexErrMsgTxt( "Invalid pCurr" );
                }
            }
            
            // Which of the possible previous (source) states is most likely
            // given this information?
            // partial prob. of state=j at time=t (delta);
            // most likely state (psi)
            GETS(delta,j,t) = valArgMax( pCurr, nStates, &GETS(psi,j,t) );
        }
    }
    //mexPrintf("done. ");
    
     
    // Termination: Find most likely endpoint.
    vLL = valArgMax( &GETS(delta,0,nObs-1), nStates, &endState );    
    
    if( vLL>0 )
        mexErrMsgTxt( "Invalid vLL" );
    

    // Backtrace to determine most likely path
    vPath[nObs-1] = endState;

    for( t=nObs-2; t>=0; --t )
    {
        vPath[t] = GETS(psi,vPath[t+1]-1,t+1);
        if( vPath[t+1]-1 < 0 || vPath[t+1]-1 > nStates )
            mexErrMsgTxt( "Backtrace bounds check failed." );
    }
    
    //cleanup
    mxFree( pCurr );
    mxFree( delta );
    mxFree( psi );
    //mexPrintf("really done.\n");
    return vLL;
}
