/* forwardBackward  Forward-backward algorithm for hidden Markov modeling
 * 
 * [LL,alpha,beta,gamma,E] = forwardBackward(p0,A,B)
 * Returns the log likelihood (LL) of one FRET trace given model parameters.
 * p0 are initial state probabilities, A are transition probabilities between
 * states, and B are Guassian observation probabilities (frames in rows and
 * states across columns). For parameter optimization in Baum-Welch, the
 * forward and backward partial probabilities (alpha and beta), state
 * probabilities at each time (gamma), and average transition probabilities
 * (E) are calculated if requested.
 *
 * NOTE: this function requires Eigen3 template library.
 * http://eigen.tuxfamily.org/dox/
 *
 * See also: forwardBackward.m, bwOptimize, mplOptimize, batchKinetics.
 *
 *
 * Copyright 2008-2018 Cornell University All Rights Reserved.
 */


/* Notes on Eigen idioms and idiosyncrasies:
 * 
 * 1) Eigen variables can be arrays, which only support element-wise
 *    operations, or matrices, which support linear algebra operations.
 *    The two types can be interconverted with .array() and .matrix().
 *
 * 2) An set of operations on Eigen variables builds an expression, which
 *    is only evaluated when necessary. This enables many optimizations.
 * 
 * 3) Map<> variables create a Eigen wrapper around existing data allocated
 *    for example on the heap or Matlab's address space to avoid copying.
 *
 * 2) VectorsXd is a column vector by default. Use RowVectorXd for rows.
 *
 */


#include <stdio.h>
#include <csignal>
#include <iostream>
#include <cstdlib>

#include <string.h>
#include <math.h>

#ifdef MEX_FCN
#include "mex.h"
#include "matrix.h"
#define NUMEL(X) mxGetNumberOfElements(X)
#define ISVECTOR(X) ( mxGetN(X)!=1 && mxGetN(X)!=1 )
#endif

// Raise an error when Eigen encounters an invalid expression, rather than
// an abort interrupt, which is the default behavior. This allows us to
// catch the exception rather than crashing.
#ifndef NDEBUG
#include <stdexcept>
#define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);
#endif

#include <Eigen/Dense>
using namespace Eigen;


//const double PI=3.141592653589793238463;



// Calculate forward probabilities
//double* pAlpha, double* pBeta, double* pGamma, double* pE )  //outputs
//(data, mu, sigma, p0, A);
double forwardBackward( const double* pp0, const double* pA, const double* pB, \
                        const int nStates,   const int nObs, \
                        double* pAlpha, double* pBeta, double* pGamma, double* pE)
                         
{
    double LL=0;
    
    // Create Eigen matrix wrappers around input pointers to Matlab data.
    // Note that both Eigen and Matlab are column major.
    Map<const RowVectorXd> p0(pp0, nStates);
    Map<const MatrixXd> A(pA, nStates, nStates);
    Map<const MatrixXd> B(pB, nObs, nStates);
    
    
    // Calculate emission probabilities (B matrix).
    //Map<const ArrayXd> data(pData, nObs);
    //MatrixXd B(nObs, nStates);
    //for( int i=0; i<nStates; ++i )
    //    B.col(i) = exp(-0.5 * square((data-pMu[i])/pSigma[i])) / (sqrt(2*PI) * pSigma[i])  + 1e-15;
    
    
    // Calculate forward probability for each timepoint in the series.
    // alpha(t,i) = P( observations 1..t & state(t)=i | model )
    // alpha(t,j) = SUM_i[ alpha(t-1,i) * A(i,j) * B(t,j) ]
    if(pAlpha==NULL) return LL;
    Map<MatrixXd> alpha( pAlpha, nObs,nStates );
    ArrayXd nrm(nObs);  //scaling coefficients
    
    alpha.row(0) = p0.cwiseProduct( B.row(0) );
    nrm(0) = 1/alpha.row(0).sum();
    alpha.row(0) *= nrm(0);

    for( int t=1; t<nObs; ++t )
    {
        // Matlab: alpha(t,:) = alpha(t-1,:) * A * diag(Bx(t,:));
        alpha.row(t) = alpha.row(t-1) * A * B.row(t).asDiagonal();
        
        // Normalize, keeping coefficient for scaling of backward probabilities
        nrm(t) = 1/alpha.row(t).sum();
        alpha.row(t) *= nrm(t);
    }
    LL = -log(nrm).sum();
    
    
    // Calculate backward probabilities
    // beta(t,i) = P( observations t+1..end | state(t)=i & model )
    // beta(t,i) = SUM_j[  A(i,j) * B(t+1,j) * beta(t+1,j)  ]
    // Matlab: beta(t,:) = A * beta(t+1,:) * diag(B(t+1,:)) * nrm(t);
    if(pBeta==NULL) return LL;
    Map<MatrixXd> beta(pBeta, nObs,nStates);
    beta.row(nObs-1).array() = 1;
    
    for( int t=nObs-2; t>=0; --t )
        beta.row(t) = A * B.row(t+1).asDiagonal() * beta.row(t+1).transpose() * nrm(t);
    
    
    // Calculate probability of being in each state at each time:
    // gamma(t,i) = P( state(t)=i | all obseravtions & model )
    // SUM_t(gamma) is the expected number of times each state is occupied.
    // NOTE: converting these statements to the equivalent ArrayXXd version
    // produces slightly different results for unknown reasons.
    if(pGamma==NULL) return LL;
    
    Map<ArrayXXd> gamma(pGamma, nObs,nStates);
    gamma = alpha.cwiseProduct( beta );
    gamma.colwise() /= gamma.rowwise().sum();
    
    
    // Calculate instantaneous transition probabilities.
    // Used by Baum Welch for estimating the transition probability matrix A.
    // NOTE: the inner loop below can be done in a single line, but the
    // correct syntax for eigen with implicit expansion isn't clear.
    if(pE==NULL) return LL;
    Map<MatrixXd> Etot(pE, nStates,nStates);
    MatrixXd es(nStates, nStates);
    
    for( int t=0; t<nObs-1; ++t )
    {
        //E(t,i,j) = alpha(t,i) * A(i,j) * B(t+1,j) * beta(t+1,j) / norm(E(i,j))
        for( int i=0; i<nStates; ++i )
            es.row(i) = alpha(t,i) * A.row(i).array() * B.row(t+1).array() * beta.row(t+1).array();
        
        Etot += es / es.sum();   //normalize
    }
    
    
    return LL;
    
}  //function forwardBackward



//Matlab entry point
//FORMAT: [LL,alpha,beta,gamma,E] = forwardBackward( p0, A, B )
//
// FIXME: does not support degenerate states!
// FIXME: for now we require alpha output. but may only want LL.
//
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
#ifdef MEX_FCN
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double LL = 0;
    
    //Verify number of input/output arguments
    if( nrhs!=3 )
        mexErrMsgTxt( "Incorrect number of input arguments" );
    
    if( nlhs>5 || nlhs<2 )
        mexErrMsgTxt( "Incorrect number of output arguments" );
    
    mwSize nObs    = mxGetM(prhs[2]);  //number of frames in input data
    mwSize nStates = mxGetN(prhs[2]);
    
    
    //Verify all arguments are valid and have matching dimensions
    for( int i=0; i<nrhs; ++i )
        if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) )
            mexErrMsgTxt( "All arguments must be real doubles." );
    
    if( !ISVECTOR(prhs[0]) )
        mexErrMsgTxt( "First argument (p0) must be a vector." );
    
    if( nStates<1 || nObs<1 )
        mexErrMsgTxt( "Data or model empty?" );
    
    if( mxGetM(prhs[1])!=nStates || mxGetN(prhs[1])!=nStates  )
        mexErrMsgTxt( "Parameter size mismatch: number of states not equal. (Degenerate states not supported)" );

    
    //Allocate output arrays
    if( nlhs>1 )  plhs[1] = mxCreateDoubleMatrix(nObs,nStates,mxREAL);    //alpha
    if( nlhs>2 )  plhs[2] = mxCreateDoubleMatrix(nObs,nStates,mxREAL);    //beta
    if( nlhs>3 )  plhs[3] = mxCreateDoubleMatrix(nObs,nStates,mxREAL);    //gamma
    if( nlhs>4 )  plhs[4] = mxCreateDoubleMatrix(nStates,nStates,mxREAL); //E
    
    double* pAlpha = (nlhs>1) ? mxGetPr(plhs[1]) : NULL;
    double* pBeta  = (nlhs>2) ? mxGetPr(plhs[2]) : NULL;
    double* pGamma = (nlhs>3) ? mxGetPr(plhs[3]) : NULL;
    double* pE     = (nlhs>4) ? mxGetPr(plhs[4]) : NULL;
    
    
    //Run the forward backward algorithm, save results into output pointers.
    try {
        //data, mu, sigma, p0, A, nStates, nObs, alpha, beta, gamma, E.
        LL = forwardBackward( mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), \
                              nStates, nObs, pAlpha, pBeta, pGamma, pE );
    } catch (const std::exception& e) {
        // Gracefully catch errors in Eigen (presumably only in debug mode).
        // Requires the eigen_assert declaration above.
        std::ostringstream fmt;
        fmt << "forwardBackward.mex internal failure: " << e.what();
        mexErrMsgTxt( fmt.str().c_str() );
    }
    
    plhs[0] = mxCreateDoubleScalar(LL);
    return;
}

#endif



