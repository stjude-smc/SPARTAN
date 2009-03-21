/*

 *WARNING: this method may fail if multiple states share the same conductance!

 * TODO: raise an error OR SOMETHING if the tree cannot be saved...
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
    if (nlhs!=0 || nrhs<2)
        mexErrMsgTxt( "Invalid arguments" );
    
    //if( ~mxIsStruct(prhs[0]) )
    //    mexErrMsgTxt("Argument 1 must be a structure!");
    
    string rootName;
    if( nrhs>2 )
    {
        char* pname = mxArrayToString(prhs[2]);
        rootName = pname;
        mxFree( pname );
    }
    else
        rootName = "Root";
    
    //Construct a QUB_Tree object from the given structure
    mxArray* modelStruct  = mxDuplicateArray( prhs[0] );  //not necessary?
    QUB_Tree outputTree = structToTree( modelStruct, rootName );
    
    //Save the tree to file
    char* modelFilename = mxArrayToString(prhs[1]);
    bool result = outputTree.saveAs(modelFilename);
    outputTree.close();
    
    if( !result )
        mexErrMsgTxt( "Can't save model file" );
    
    //Cleanup
    mxFree( modelFilename );
    mxDestroyArray( modelStruct );
    return;
}
