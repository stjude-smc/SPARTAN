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
    if (nlhs!=0 || nrhs!=2)
        mexErrMsgTxt( "Invalid arguments" );
    
    //if( ~mxIsStruct(prhs[0]) )
    //    mexErrMsgTxt("Argument 1 must be a structure!");
    
    //Construct a QUB_Tree object from the given structure
    mxArray* modelStruct  = mxDuplicateArray( prhs[0] );  //not necessary?
    QUB_Tree outputTree = structToTree( modelStruct, "ModelFile" );
    
    //Save the tree to file
    char* modelFilename = mxArrayToString(prhs[1]);
    outputTree.saveAs(modelFilename);
    outputTree.close();
    
    //Cleanup
    mxFree( modelFilename );
    mxDestroyArray( modelStruct );
    return;
}
