

#include "qubmatlab.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cstring>

using namespace std;



//Now...how do we deal with struct arrays? -- for now just use first element.
void structToTree( mxArray* structure, string rootName, QUB_Tree& parent )
{
    int nFields = mxGetNumberOfFields(structure);
    int fieldID,i;
    
        
    for( i=0; i<numel(structure); ++i )
    {
        QUB_Tree newChild = QUB_Tree::Create(rootName);
        
        //for each field, ...
        for( fieldID=0; fieldID<nFields; ++fieldID )
        {
            string fieldName = mxGetFieldNameByNumber(structure,fieldID);
            replace( fieldName.begin(),fieldName.end(),'_',' ' );

            mxArray* field = mxGetFieldByNumber(structure,i,fieldID);

            if( fieldName=="dataType" )
                continue;

            //If data element found, save as node data.
            if( fieldName=="data" && !mxIsEmpty(field) )
            {
                int M = mxGetM(field);
                int N = mxGetN(field);

                //mexPrintf("%s: %d x %d ",rootName.c_str(),M,N);

                if( mxIsDouble(field) )
                {
                    double* fieldData = mxGetPr(field);
                    newChild.setNumData( QTR_TYPE_DOUBLE, M,N, fieldData );
                    //mexPrintf("double\n");
                }
                else if( mxIsInt32(field) )
                {
                    int* fieldData = (int*)mxGetData(field);
                    newChild.setNumData( QTR_TYPE_INT, M,N, fieldData );
                    //mexPrintf("int32\n");
                }
                else if( mxIsChar(field) )
                {
                    char* str = mxArrayToString(field);
                    newChild.setData( str );
                    mxFree(str);
                    //mexPrintf("string\n");
                }
                //else
                //  mexPrintf("unknown type\n");
             
            }

            //Otherwise, add it as a child node...
            else
                structToTree(field,fieldName,newChild);
        }
    
        //Add new child node to parent
        parent.appendChild( newChild );
    }
    return;
}




//
QUB_Tree structToTree( mxArray* structure, string rootName )
{
    QUB_Tree outputTree = QUB_Tree::Create(rootName);
    structToTree( structure, rootName, outputTree );
    return outputTree[rootName];
}





    
//
mxArray* treeToStruct( QUB_Tree node, int depth )
{
    //char pattern[100]; //for sprintf
    //memset( pattern, ' ', depth*2 );
    //pattern[depth*2] = 0;

    //Create structure for adding fields
    mxArray* structure = mxCreateStructMatrix(1,1, 0,NULL);
    
    if( !structure )
        mexErrMsgTxt("can't alloc struct"); 
    
    
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
    
    
    //Parse data entry...
    int i;
    mxArray* data = 0;
    mwSize M,N;
    mxClassID mxtype;
    void* pdata;
    unsigned char* pcdata;
    
    
    M = node.dataCols();
    N = node.dataRows();
    
    switch( node.dataType() )
    {
    case QTR_TYPE_EMPTY:
        break; //no data
        
    case QTR_TYPE_POINTER:
        mexWarnMsgIdAndTxt("qubtree:PointerFieldsNotSupported","Pointer fields not supported.");
        break;
        
    case QTR_TYPE_STRING:
        //string
        //mexPrintf("%s - %s\n", pattern, node.dataAsString(true).c_str());
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
    case QTR_TYPE_FLOAT:
    case QTR_TYPE_DOUBLE:
        //mexPrintf("%s - %d x %d matrix\n", pattern, M,N);

        if( M<1 || N<1 )
        {
            mexWarnMsgTxt("invalid data size");
            break;
        }
        if( !(M==1 || N==1) )
            mexWarnMsgIdAndTxt("qubtree:MatrixFieldsNotSupported","Matrix data not supported!");
        
        mxtype = mxTypeLookup[node.dataType()];
        data = mxCreateNumericMatrix( M,N,mxtype,mxREAL );
        pdata = mxGetData( data );
        memcpy( pdata, node.data(), M*N*mxGetElementSize(data) );
        break;
        
    default:
        //mexPrintf("%s - %d (Unknown type)",pattern, node.dataType());
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

        //mexPrintf("%s%d* %s (%d)\n", pattern, depth, childName.c_str(), nTwins);
        
        mxArray* twins = mxCreateStructMatrix(nTwins,1, 0, NULL);
        if( twins<=(void*)0 )
            mexErrMsgTxt("can't alloc twins"); 
        
        //for each twin node, compile all fields into a single structure
        int tci_i = 0;
        for (tci=node.find(childName); !tci->isNull(); tci.nextSameName())
        {
            if( tci_i>=nTwins )
                mexErrMsgTxt("too far"); 
            
            mxArray* twin = treeToStruct(*tci, depth+1);
            
            //for each field in the current twin node
            for( int fid=0; fid<mxGetNumberOfFields(twin); ++fid ) //for each field
            {
                int newFid = mxAddField( twins, mxGetFieldNameByNumber(twin,fid) );
                if( newFid<0 )
                    mexErrMsgTxt("invalid field...");
                
                mxArray* field = mxDuplicateArray( mxGetFieldByNumber(twin,0,newFid) );
                
                if( field<=(void*)0 )
                    mexErrMsgTxt("invalid field..."); 
                
                //mexPrintf(" %d*     test %d %d %d\n",depth,fid,tci_i,mxGetNumberOfFields(field));
                mxSetFieldByNumber( twins, tci_i, newFid, field );
            }
            
            mxDestroyArray(twin);
            
            ++tci_i;
        }
        
        //Add field from above structure to output structure
        replace( childName.begin(),childName.end(),' ','_' );     
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




