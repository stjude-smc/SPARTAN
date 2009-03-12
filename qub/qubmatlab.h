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
mxArray* treeToStruct( QUB_Tree tree );

//Convert a structure array created by treeToStruct
//into a QUB_Tree object with root node named <rootName>
QUB_Tree structToTree( mxArray* structure, const char* rootName );
