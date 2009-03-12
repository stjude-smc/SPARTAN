/*

 *WARNING: this method may fail if multiple states share the same conductance!
*/


#include "qubmatlab.h"



//Matlab entry point
//FORMAT: tree (filename) -> struct representing tree
//
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
    //Check input arguments
    if (nlhs!=1 || nrhs!=1)
        mexErrMsgTxt( "Invalid arguments" );
    
    //Parse input parameters
    char* modelFilename = mxArrayToString(prhs[0]);
    
    //Load QuB model as a QUB_Tree
    QUB_Tree model( QUB_Tree::Open(modelFilename, true) );
    if ( model.isNull() )
        mexErrMsgTxt("Not a valid model file.");
    mxFree( modelFilename );
    
    plhs[0] = treeToStruct( model );
    
    return;
}

    



