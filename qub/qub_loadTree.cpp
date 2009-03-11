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
mxArray* traverseNode( QUB_Tree node, int depth=0 );

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
    if (nlhs!=1 || nrhs!=1)
        mexErrMsgTxt( "Invalid arguments" );
    
    //Parse input parameters
    char* modelFilename = mxArrayToString(prhs[0]);
    
    //Load QuB model as a QUB_Tree
    QUB_Tree model( QUB_Tree::Open(modelFilename, true) );
    if ( model.isNull() )
        mexErrMsgTxt("Not a valid model file.");
    mxFree( modelFilename );
    
    plhs[0] = traverseNode( model );
    
    return;
}

    
    
    
mxArray* traverseNode( QUB_Tree node, int depth )
{
    int i;
    
    
    //Create structure for adding fields
    mxArray* structure = mxCreateStructMatrix(1,1, 0,NULL);
    
    /*
    mxClassID mxTypeLookup[32];
    mxTypeLookup[ QTR_TYPE_CHAR   ] = mxINT8_CLASS;
    mxTypeLookup[ QTR_TYPE_UCHAR  ] = mxUINT8_CLASS;
    mxTypeLookup[ QTR_TYPE_SHORT  ] = mxINT16_CLASS;
    mxTypeLookup[ QTR_TYPE_USHORT ] = mxUINT16_CLASS;
    mxTypeLookup[ QTR_TYPE_INT    ] = mxINT32_CLASS;
    mxTypeLookup[ QTR_TYPE_UINT   ] = mxUINT32_CLASS;
    mxTypeLookup[ QTR_TYPE_LONG   ] = mxINT32_CLASS;
    mxTypeLookup[ QTR_TYPE_ULONG  ] = mxUINT32_CLASS;
    mxTypeLookup[ QTR_TYPE_FLOAT  ] = mxSINGLE_CLASS;
    mxTypeLookup[ QTR_TYPE_DOUBLE ] = mxDOUBLE_CLASS;
    */
    
    //Parse data entry...
    mxArray* data = 0;
    mwSize M,N;
    mxClassID mxtype;
    double* pddata;
    int* pidata;
    unsigned char* pcdata;
    
    
    switch( node.dataType() )
    {
    case QTR_TYPE_EMPTY:
        break; //no data
        
    case QTR_TYPE_POINTER:
        mexWarnMsgTxt("Pointer fields not supported.");
        break;
        
    case QTR_TYPE_STRING:
        //string
        data = mxCreateString( node.dataAsString(true).c_str() );
        break;
    
    case QTR_TYPE_UCHAR:
    case QTR_TYPE_CHAR:
    case QTR_TYPE_USHORT:
    case QTR_TYPE_SHORT:
    case QTR_TYPE_UINT:
    case QTR_TYPE_INT:
    case QTR_TYPE_ULONG:
    case QTR_TYPE_LONG:
        M = node.dataCols();
        N = node.dataRows();
        
        data = mxCreateNumericMatrix( M,N,mxINT32_CLASS,mxREAL );
        pidata = (int*)mxGetData( data );
        node.getDataAsInts(pidata,0,N-1);
        break;
        
    case QTR_TYPE_FLOAT:
    case QTR_TYPE_DOUBLE:
        M = node.dataCols();
        N = node.dataRows();
        
        data = mxCreateDoubleMatrix( M,N,mxREAL );
        pddata = mxGetPr( data );
        node.getDataAsDoubles(pddata,0,N-1);
        break;
        
    default:
        mexPrintf("%d",node.dataType());
        mexWarnMsgTxt("Unsupported field type");
    }
    
    if( data )
    {
        int fid = mxAddField(structure,"data");
        mxSetFieldByNumber(structure,0,fid,data);
        
        fid = mxAddField(structure,"dataType");
        data = mxCreateNumericMatrix( 1,1,mxUINT8_CLASS,mxREAL );
        pcdata = (unsigned char*)mxGetData(data);
        *pcdata = node.dataType();
        mxSetFieldByNumber(structure,0,fid,data);
    }
    //otherwise, it will just be an empty matrix
    
    
    
    //Get a list of all unique
    vector<string> childNames = listNames(node);
    
    //For each unique name, create a struct
    QUB_TreeMonoIter tci;
    
    for( i=0; i<childNames.size(); ++i ) //for each child...
    {
        string childName = childNames[i];
        int nTwins = countChildren( node, childName );
        
        //mexPrintf("%d* %s (%d x)\n", depth, childName.c_str(), nTwins);
        
        mxArray* twins = mxCreateStructMatrix(nTwins,1, 0, NULL);
        
        if( twins<=0 )
            mexErrMsgTxt("bad twins..."); 
        
        //for each twin node, compile all fields into a single structure
        int tci_i = 0;
        for (tci=node.find(childName); !tci->isNull(); tci.nextSameName())
        {
            if( tci_i>=nTwins )
                mexErrMsgTxt("too far"); 
            
            mxArray* twin = traverseNode(*tci,depth+1);
            
            //for each field in the current twin node
            for( int fid=0; fid<mxGetNumberOfFields(twin); ++fid ) //for each field
            {
                int newFid = mxAddField( twins, mxGetFieldNameByNumber(twin,fid) );
                if( newFid<0 )
                    mexErrMsgTxt("invalid field...");
                
                mxArray* field = mxDuplicateArray( mxGetFieldByNumber(twin,0,newFid) );
                
                if( field<=0 )
                    mexErrMsgTxt("invalid field..."); 
                
                //mexPrintf(" %d*     test %d %d %d\n",depth,fid,tci_i,mxGetNumberOfFields(field));
                mxSetFieldByNumber( twins, tci_i, newFid, field );
            }
            
            mxDestroyArray(twin);
            
            ++tci_i;
        }
        
        //Add field from above structure to output structure        
        int fieldID = mxAddField(structure,childName.c_str());
        mxSetFieldByNumber(structure,0,fieldID,twins);
    }
    
    return structure;
}



vector<string> fieldNames( mxArray* structure )
{
    vector<string> output;
    int nFields = mxGetNumberOfFields(structure);
    
    for( int i=0; i<nFields; ++i )
        output.push_back( mxGetFieldNameByNumber(structure,i) );
    
    return output;
}

int countChildren( QUB_Tree tree, string name )
{
    int nChildren = 0;
    QUB_TreeMonoIter tci;
    
    for (tci=tree.find(name); !tci->isNull(); tci.nextSameName())
        ++nChildren;
    
    return nChildren;
}

vector<string> listNames( QUB_Tree tree )
{
    vector<string> output;
    QUB_TreeMonoIter tci;
    
    for (tci=tree.children(); !tci->isNull(); ++tci)
        output.push_back( tci->name() );
    
    //remove duplicates
    vector<string>::iterator new_end = unique(output.begin(), output.end());
   // delete all elements past new_end 
   output.erase(new_end, output.end());

    
    return output;
}


