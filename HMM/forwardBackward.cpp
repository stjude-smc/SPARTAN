/* forwardBackward  Forward-backward algorithm probabilities
 * 
 * [alpha,beta,gamma,E] = forwardBackward(p0,A,B)
 * Calculates the forward (alpha) and backward (beta) probabilities.
 * p0 is a vector of initial probabilities for each state.
 * A is the transition probability matrix (rows sum to 1).
 * B is the probability of observing the data at time t (rows) if the hidden
 *   state is i (columns).
 *
 * WARNING: assumes no degenerate states!
 *
 * See also: forwardBackward.m, bwOptimize, mplOptimize.
 *
 * NOTE: this function requires Eigen3 template library to be installed.
 * 
 *
 */

/* A few notes on Eigen idioms and idiosyncrasies:
 * 
 * 1) Eigen tends to raise abort interrupts instead of exceptions, for example
 *    when you try to multiply mismatched matrices. I added the eigen_assert()
 *    below to make these regular exceptions.
 *
 * 2) Eigen has no equivalent operators for element-wise multiplication, etc.
 *    Instead, one must use .array() to inform eigen that we don't want to do
 *    any linear algebra.
 *
 * 3) VectorsXd is a column vector by default. Use RowVectorXd for rows.
 *
 * 4) 
 *
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
#endif

#ifndef NDEBUG
#include <stdexcept>
#define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);
#endif

#include <Eigen/Dense>
using namespace Eigen;
//prevent Eigen from allocating anything for debugging.
//Eigen::set_is_malloc_allowed(false);




// Calculate forward probabilities
//double* pAlpha, double* pBeta, double* pGamma, double* pE )  //outputs
double forwardBackward(  const double* pp0, const double* pA, const double* pB, \
                         const int nStates, const int nObs, \
                         double* pAlpha, double* pBeta, double* pGamma, double* pE)
                         
{
    double LL=0;
    
    // Create Eigen matrix wrappers around input pointers to Matlab data.
    // Note that both Eigen and Matlab are column major.
    Map<const RowVectorXd> p0(pp0, nStates);
    Map<const MatrixXd> A(pA, nStates, nStates);
    Map<const MatrixXd> B(pB, nObs, nStates);
    
    
    // Calculate forward probability for each timepoint in the series.
    // alpha(t,i) = P( observations 1..t & state(t)=i | model )
    // alpha(t,j) = SUM_i[ alpha(t-1,i) * A(i,j) * B(t,j) ]
    if(pAlpha==NULL) return LL;
    Map<MatrixXd> alpha( pAlpha, nObs,nStates );
    VectorXd nrm(nObs);  //scaling coefficients
    
    alpha.row(0) = p0.cwiseProduct( B.row(0) );
    nrm(0) = 1/alpha.row(0).sum();
    alpha.row(0) = alpha.row(0) * nrm(0);

    for( int t=1; t<nObs; ++t )
    {
        // Matlab: alpha(t,:) = alpha(t-1,:) * A * diag(Bx(t,:));
        alpha.row(t) = alpha.row(t-1) * A * B.row(t).asDiagonal();
        
        // Normalize, keeping coefficient for scaling of backward probabilities
        nrm(t) = 1/alpha.row(t).sum();
        alpha.row(t) = alpha.row(t) * nrm(t);
    }
    LL = -log(nrm.array()).sum();  //log likelihood = -sum( log(nrm) )
    
    
    // Calculate backward probabilities
    // beta(t,i) = P( observations t+1..end | state(t)=i & model )
    // beta(t,i) = SUM_j[  A(i,j) * B(t+1,j) * beta(t+1,j)  ]
    // Matlab: beta(t,:) = A * beta(t+1,:) * diag(B(t+1,:)) * nrm(t);
    if(pBeta==NULL) return LL;
    Map<MatrixXd> beta(pBeta, nObs,nStates);
    beta.row(nObs-1).array() = 1;
    
    for( int t=nObs-2; t>=0; --t )
    {
        VectorXd temp = A * B.row(t+1).asDiagonal() * beta.row(t+1).transpose();
        beta.row(t) = temp.transpose() * nrm(t);
    }
    
    
    // Calculate probability of being in each state at each time:
    // gamma(t,i) = P( state(t)=i | all obseravtions & model )
    // SUM_t(gamma) is the expected number of times each state is occupied.
    if(pGamma==NULL) return LL;
    Map<MatrixXd> gamma(pGamma, nObs,nStates);
    
    gamma = alpha.cwiseProduct( beta );
    VectorXd gnorm = gamma.array().rowwise().sum();
    gamma = gamma.cwiseQuotient( gnorm.replicate(1,nStates) );  //slow??
    
    
    // Calculate instantaneous transition probabilities.
    // Used by Baum Welch for estimating the transition probability matrix A.
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
    if( nStates<1 || nObs<1 )
        mexErrMsgTxt( "Data or model empty?" );
    
    if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || \
        mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) || mxIsComplex(prhs[2]) )
        mexErrMsgTxt( "All argument must be real doubles" );
    
    if( mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=nStates )
        mexErrMsgTxt( "Argument 1 (p0) size is not correct. Should be a 1 x nStates row vector" );
    
    if( mxGetM(prhs[1])!=nStates || mxGetN(prhs[1])!=nStates )
        mexErrMsgTxt( "Argument 2 (transition probabilities) has an invalid size. Should be nStates x nStates" );

    
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
        LL = forwardBackward( mxGetPr(prhs[0]), mxGetPr(prhs[1]), \
                              mxGetPr(prhs[2]), nStates, nObs, \
                              pAlpha, pBeta, pGamma, pE );
    } catch (const std::exception& e) {
        // Gracefully catch errors in Eigen (presumably only in debug mode).
        // Requires the eigen_assert declaration above.
        std::ostringstream fmt;
        fmt << "forwardBackward internal failure: " << e.what();
        mexErrMsgTxt( fmt.str().c_str() );
    }
    
    plhs[0] = mxCreateDoubleScalar(LL);
    return;
}
#endif



