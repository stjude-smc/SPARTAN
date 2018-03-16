/* Gillespie stochastic simulation algorithm.
 * 
 *   
 * 
 */

#include "mex.h"
#include "matrix.h"
#include <string.h>  //memcpy
#include <math.h>    //log
#include <vector>
#include "twister.h"  //Mersenne Twister PRNG implementation

using namespace std;

// Macros for 2-dimensional array access.
#define GETS(A,I,J)  A[(I) + (J)*nStates]

// Uniform random number (Mersenne Twister) over [0,1]-real-interval
#define RAND randomMT()*(1.0/4294967295.0)



// Gillespie algorithm (direct method) for stochastic simulation.
void gillespie(  const int nStates, const int startState, const double endTime, \
                   const double* Qtau, const double* Qcumsum, \
                   vector<short>& states, vector<double>& times  ) //outputs
{
    double cumTime = 0;
    short nextState, curState=startState-1;  //zero-based
    double r;
    
    while( cumTime<endTime )
    {
        // Time until next transition, drawn from exponential distribution
        // with time constants of 1000/sum(Q(s,:)).
        times.push_back(  Qtau[curState] * log(RAND)  );
        
        // Randomly sample final state with probabilities calculated as the
        // fraction of all possible rate constants exiting current state.
        r = RAND;
        nextState = 0;
        while(  r > GETS(Qcumsum,curState,nextState)  )
            ++nextState;
        
        //mxAssert( r>=0.0 & r<=1.0, "PRNG out of range!" );
        //mxAssert( nextState>=1 & nextState<=nStates, "Invalid state number" );
        
        // Update iterators
        states.push_back( curState+1 );
        curState = nextState;
        cumTime += times.back();
    }
    
    return;
}



/* Matlab entry point
 * FORMAT: [states, dwells] = gillespie(Qtau, Qcumsum, startState, endTime)
 *
 * startState: Initial state of the system before iteration.
 * endTime:    Maximum trace length in milliseconds.
 * Qtau:       Expected lifetime of each state: sum(Q(s,:))
 * Qcumsum:    cumsum(s,:)/Qtau(s) for randomly choosing the next state in sequence.
 *
 * states:     int16 array of state numbers (1..nStates).
 * times:      double array of dwell times (in ms), same size as states.
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Check number of input and output parameters
    if( nrhs!=4 )
        mexErrMsgTxt( "Not enough input arguments" );
    
    if( nlhs>2 )
        mexErrMsgTxt( "Too many output arguments" );
    
    // Veryify input parameter sizes
    int nStates = mxGetNumberOfElements(prhs[2]);
    
    if( nStates<=1 )
        mexErrMsgTxt("Qtau cannot be empty");
    
    if( mxGetM(prhs[3])!=nStates || mxGetN(prhs[3])!=nStates )
        mexErrMsgTxt("Qcumsum must be a symmetric matrix in the number of states");
    
    
    // Run Gillespie algorithm
    vector<short> states;
    vector<double> times;  //outputs from gillespie(), passed by reference
    
    gillespie( nStates, mxGetScalar(prhs[0]), mxGetScalar(prhs[1]), \
               mxGetPr(prhs[2]), mxGetPr(prhs[3]), states, times);
    
    
    // Return results to back to MATLAB
    if( nlhs >= 1 ) {
        plhs[0] = mxCreateNumericMatrix(1, states.size(), mxINT16_CLASS, mxREAL);
        memcpy( mxGetPr(plhs[0]), &states[0], states.size()*sizeof(short) );
    }
    if( nlhs >= 2 ) {
        plhs[1] = mxCreateDoubleMatrix(1, times.size(), mxREAL);
        memcpy( mxGetPr(plhs[1]), &times[0], times.size()*sizeof(double) );
    }
    
//     mexPrintf("DONE\n");
    return;
}


