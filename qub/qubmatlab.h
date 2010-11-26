/*

 *WARNING: this method may fail if multiple states share the same conductance!
*/


#include "mex.h"
#include "qubtree/QUB_Tree.h"


//useful shortcuts to make MEX/C++ look like Matlab
#define numel mxGetNumberOfElements
#define isfield(A,B) mxGetFieldNumber(A,B)!=-1

//Shortcut for accessing elements in a 2D array.
//A is the matrix, m is a mwSize object, i is row, j is column (as in Matlab).
#define getArrayElement(A,M,i,j) A[(i)+(j)*M]

#define MIN(a,b)  (((a) > (b)) ? (a) : (b))
#define MAX(a,b)  (((a) < (b)) ? (a) : (b))


//Generate a structure array representing a QUB_Tree object
mxArray* treeToStruct( QUB_Tree tree, int depth=0 );

//Convert a structure array created by treeToStruct
//into a QUB_Tree object with root node named <rootName>
QUB_Tree structToTree( mxArray* structure, string rootName );

//QuB's Maximum-Interval-Likelihood algorithm for direct
//optimization a kinetic model based on dwell-time information.
QUB_Tree milOptimize( string dataFilename, string modelFilename );
QUB_Tree milOptimize( QUB_Tree data, QUB_Tree model );




//accessory functions (see treestruct.cpp)
int countChildren( QUB_Tree tree, string name );
vector<string> fieldNames( mxArray* structure );
vector<string> listNames( QUB_Tree tree );
