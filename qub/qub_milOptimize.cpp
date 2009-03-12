/*
qubMIL_fromFile

DESCRIPTION  Optimizes a kinetic model given dwell-times data
             using QuB's MIL algorithm
AUTHOR       Daniel Terry
COPYRIGHT    2009, Scott Blanchard
*/


#include "qubmatlab.h"


// QUB header files (qub/qubopt)
#ifdef QUB_DLL
#include "qubopt.h"
#include "mdlrep.h"
#include "QUB_IdlSeg.h"	 //for loading idealization (.dwt) files
#include "QUB_Tree.h"
#endif

//QTR_Impl* statusCallback( QTR_Impl* input );


#include <string>
#include <iostream>

using namespace std;




#ifdef MAIN_FCN

#include <cstdlib>

#define ERROR(M) cerr << M; exit(-1);
#define WARN(M) cout << "Warning: " << M;

int main(int argc, char* argv[])
{
	int dummy;

	string base("Z:\\qub\\testdata\\");
	string dataPath     = base+ "snr8.qub.dwt";
	string modelPath    = base+ "model.qmf";
	string outputModel  = base+ "outputModel.qmf";
	string outputResult = base+ "outputResult.qtr";

	cout << "Running..." << endl;

	//Run MIL
	QUB_Tree result = milOptimize(dataPath,modelPath);

	//Save the optimized model
	QUB_Tree model = result["ModelFile"];
	model.saveAs( outputModel );
	model.close();

	//Save the result data
	result.saveAs( outputResult );
	result.close();

	cout << result["LL"].dataAsDouble(-7.0) << endl;
	
	cout << "Finished.";

	cin >> dummy;
	return 0;
}
#endif



#ifdef MEX_FCN
#include "mex.h"

#define ERROR(M) mexErrMsgTxt(M); return
#define WARN(M) mexWarnMsgTxt(M)


#define numel mxGetNumberOfElements
#define isfield(A,B) mxGetFieldNumber(A,B)!=-1

//Shortcut for accessing elements in a 2D array.
//A is the matrix, m is a mwSize object, i is row, j is column (as in Matlab).
#define getArrayElement(A,M,i,j) A[(i)+(j)*M]


//Matlab entry point
//FORMAT: [dwtFilename,modelFilename] -> resultsTree
//
//<plhs> contains left-hand-side (<nlhs> of them), for return value
//<prhs> contains right-hand-side (<nrhs> of them), for parameters
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
    //Check input arguments
    if (nlhs!=1 || nrhs!=2)
        mexErrMsgTxt( "Invalid arguments" );
    
    //if( ~mxIsStruct(prhs[0]) )
    //    mexErrMsgTxt("Argument 1 must be a structure!");
    
    //Parse input parameters
    char* dataPath = mxArrayToString( prhs[0] );
    char* modelPath = mxArrayToString( prhs[1] );
    
    //Run MIL (via miltreeiface.exe)
    QUB_Tree result = milOptimize(dataPath,modelPath);

    //Convert the result to a struct for processing in Matlab
    plhs[0] = treeToStruct( result );
    
    //Cleanup
    mxFree(dataPath);
    mxFree(modelPath);
    return;
}

#endif




//Use the idealized DATA to optimize the starting MODEL for best LL.
//Dead-time correction is hard-coded at 0.5*framerate.  FIXME
QUB_Tree milOptimize(string dwtFilename, string modelFilename)
{
    QUB_Tree result;
    int retval = -9;
	int max_iter = 100;
    
    //Load starting model
	QUB_Tree model( QUB_Tree::Open(modelFilename, true) );
	if ( model.isNull() )
		mexErrMsgTxt("No model");

	//Setup MIL configuration from givenp arameters
    model.setData("Initial");
    
	QUB_Tree config( QUB_Tree::Create("Properties") );
	config.children().insert( model.clone() );
	config["GroupViterbi"].setData(QTR_TYPE_INT, 0);
	config["join segments"].setData(QTR_TYPE_INT, 0);
	config["ChannelIndex"].setData(QTR_TYPE_INT, 0);
	config["ThreadCount"].setData(QTR_TYPE_INT, 1);
	config["Mode"].setData((max_iter > 0) ? "optimize" : "evaluate");
	config["use segments"].setData("together");
	config["SearchLimit"].setData(QTR_TYPE_DOUBLE, 10.0);
	config["DFP"]["MaxStep"].setData(QTR_TYPE_DOUBLE, 0.1);
	config["DFP"]["ConvLL"].setData(QTR_TYPE_DOUBLE, 0.0001);
	config["DFP"]["ConvGrad"].setData(QTR_TYPE_DOUBLE, 0.01);
	config["DFP"]["MaxIterations"].setData(QTR_TYPE_INT, max_iter);
	config["DFP"]["MaxRestarts"].setData(QTR_TYPE_INT, 0);
	
    //WIN32 MIL version - directly call qubopt.dll (if available)
    #ifdef QUB_DLL
        //Load dwell-times file into a tree
        QUB_Tree data( QUB_Idl_ReadDWT(dataFilename.c_str()) );

        //Set deadtime
        double sampling = data["sampling"].dataAsDouble(0);
        double deadtime = sampling/2;
        data["ExpCond"]["DeadTime"].setData(QTR_TYPE_DOUBLE, deadtime);

        //Create data list
        QTR_Impl* adata[2];
        adata[0] = data.getImpl();
        adata[1] = NULL;

        //Run MIL via qubopt.dll
        result = miltreeiface(config.getImpl(), adata, NULL, NULL, NULL);
        QTR_DECREF( result.getImpl() );

    //Generic MIL version - use a bridge interface...
    //Under Linux, use Wine to run the program.
    #else
        //Save MIL input parameter tree
        config.saveAs( ".milconfig.qtr" );
        config.close();

        //Run MIL via command line
        string copyCmd = "\"" + dwtFilename + "\" .mildata.dwt";
        string milCmd = "miltreeiface.exe .milconfig.qtr .mildata.dwt .milresult.qtr";
        #ifdef WIN32
        copyCmd = "copy /Y " + copyCmd; 
        #else
        copyCmd = "cp -f " + copyCmd; 
        milCmd = "wine " + milCmd;
        #endif
        //milCmd = "/home/dsterry/cornell/code/cascade_git/qub/"+milCmd;
        system( copyCmd.c_str() );
        retval = system( milCmd.c_str() );
    
    #endif
    
    switch( retval )
    {
    case -1:
        mexErrMsgTxt("miltreeiface: incorrect args");
        break;
    case -2:
        mexErrMsgTxt("miltreeiface: invalid or missing DWT file");
        break;
    case 0:
        result = QUB_Tree::Open(".milresult.qtr").clone(true);
        break;
    default:
        mexPrintf("(%d) ",retval);
        mexErrMsgTxt("miltreeiface: unknown error");
    }
    
	return result;
}
