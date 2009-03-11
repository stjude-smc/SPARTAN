/*

 *WARNING: this method may fail if multiple states share the same conductance!
*/


#include "mex.h"
#include "qubtree/QUB_Tree.h"

#include "string.h"

#include <vector>
#include <string>

using namespace std;


#define numel mxGetNumberOfElements
#define isfield(A,B) mxGetFieldNumber(A,B)!=-1


int countChildren( QUB_Tree tree, string name );
vector<string> fieldNames( mxArray* structure );
// mxArray* traverseNode( QUB_Tree node, int depth=0 );
QUB_Tree structToTree( mxArray* structure, const char* rootName );
vector<string> listNames( QUB_Tree tree );

//Shortcut for accessing elements in a 2D array.
//A is the matrix, m is a mwSize object, i is row, j is column (as in Matlab).
#define getArrayElement(A,M,i,j) A[(i)+(j)*M]


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
    mxArray* modelTree  = mxDuplicateArray( prhs[0] );  //not necessary?
    QUB_Tree outputTree = structToTree( modelTree, "ModelFile" );
    
    //Save the tree to file
    char* modelFilename = mxArrayToString(prhs[1]);
    outputTree.saveAs(modelFilename);
    outputTree.close();
    
    //Cleanup
    mxFree( modelFilename );
    mxDestroyArray( modelTree );
    return;
}


//Now...how do we deal with struct arrays? -- for now just use first element.
QUB_Tree structToTree( mxArray* structure, const char* rootName )
{
    QUB_Tree output = QUB_Tree::Create(rootName);
    int nFields = mxGetNumberOfFields(structure);
    int fieldID,j;
    
    //for each field, ...
    for( fieldID=0; fieldID<nFields; ++fieldID )
    {
        const char* fieldName = mxGetFieldNameByNumber(structure,fieldID);
        
        //for( j=0; j<numel(structure); ++j )
        //{
            j = 0;
            mxArray* field = mxGetFieldByNumber(structure,j,fieldID);
            
            if( strcmp(fieldName,"dataType")==0 )
                continue;
            
            //If data element found, save as node data.
            if( strcmp(fieldName,"data")==0 && ~mxIsEmpty(field) )
            {
                int M = mxGetM(field);
                int N = mxGetN(field);
                
                mexPrintf("%s: %d x %d ",rootName,M,N);
                
                if( mxIsDouble(field) )
                {
                    double* fieldData = mxGetPr(field);
                    output.setNumData( QTR_TYPE_DOUBLE, M,N, fieldData );
                    mexPrintf("double\n");
                }
                else if( mxIsInt32(field) )
                {
                    int* fieldData = (int*)mxGetData(field);
                    output.setNumData( QTR_TYPE_INT, M,N, fieldData );
                    mexPrintf("int32\n");
                }
                else if( mxIsChar(field) )
                {
                    char* str = mxArrayToString(field);
                    output.setData( str );
                    mxFree(str);
                    mexPrintf("string\n");
                }
                else
                    mexPrintf("unknown type\n");
            }
            
            //Otherwise, add it as a child node...
            else
                output.appendChild( structToTree(field,fieldName) );
        //}
    }
    
    return output;
}
