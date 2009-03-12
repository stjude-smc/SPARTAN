/*
qubMIL_fromFile

DESCRIPTION  Optimizes a kinetic model given dwell-times data
             using QuB's MIL algorithm
AUTHOR       Daniel Terry
COPYRIGHT    2009, Scott Blanchard
 *
 *
 * NOTE: Idealization generally leaves a very long dwell in the dark
 *   state at the end of each trace.  This greatly confuses MIL -- 
 *   it will take forever to converge and its results are inaccurate.
 *   TODO: take steps to remove this dwell here or insure that all
 *   idealizers remove the final dwell in dark state.
*/

#include "qubmatlab.h"


// QUB header files (qub/qubopt)
#include "qubopt.h"
#include "mdlrep.h"
#include "QUB_IdlSeg.h"	 //for loading idealization (.dwt) files

//QTR_Impl* statusCallback( QTR_Impl* input );


#include <string>
#include <iostream>

using namespace std;



//#define MAIN_FCN 1
#define MEX_FCN 1


#ifdef MAIN_FCN
/*
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
*/
#endif


#ifdef MEX_FCN
#include "mex.h"

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
    
    //Run MIL (from qubopt.dll)
    QUB_Tree result = milOptimize(dataPath,modelPath);

    //Convert the result to a struct for processing in Matlab
    plhs[0] = treeToStruct( result );
    
    //Cleanup
    mxFree(dataPath);
    mxFree(modelPath);
    return;
}

#endif



//This appears to be a simple C/command line interface to MIL.
//
//It will:
//   1. Run MIL-Together on the given data and model,
//   2. Save an execution report as BASE_report.txt
//   3. Save resulting model as BASE_model.mdl
//
//PARAMETERS:
//   char*  data_path   : path/filename of idealization data (.dwt file)
//   char*  model_path  : path/filename of model file (.mdl) containing
//                        initial rates, constraints, etc.
//   int    max_iter    : max. number of iterations before terminating
//   double tdead_ms    : dead-time setting (in ms) = 0.5*dt
//   char*  output_path : path/filename of MIL results as a model (.mdl) file.
//
//ERROR CODES (double):
//  -11:   data file could not be loaded
//  -13:   model file could not be loaded
//  -17:   MIL execution failed
//   otherwise, the return value is log-likelihood of the data given optimized model
//
QUB_Tree milOptimize(string dataFilename, string modelFilename  )
{
	QUB_Tree result;

    //Load dwell-times data from file
	QUB_Tree data( QUB_Idl_ReadDWT(dataFilename.c_str()) );
	if ( data.isNull() )
		cerr << "No data\n";

    //Remove final dwell in dark state (assuming class=0)
    //(or else MIL may get confused...)
    /*
    sampling = data["sampling"].dataAsDouble();
    
    QUB_TreeMonoIter tci;
    for( tci=node.find("Segment"); !tci->isNull(); tci.nextSameName() )
    {
        if( (*tci)["Classes"].dataAsInt() )
    }
     **/
    
	//Load starting model
	QUB_Tree model( QUB_Tree::Open(modelFilename, true) );
	if ( model.isNull() )
		cerr << "No model\n";

	//Run MIL
	if( !data.isNull() && !model.isNull() )
		result = milOptimize(data,model);

	return result;
}


//Use the idealized DATA to optimize the starting MODEL for best LL.
//Dead-time correction
QUB_Tree milOptimize(QUB_Tree data, QUB_Tree model)
{
	QUB_Tree result;
	int max_iter = 100;
	double sampling = data["sampling"].dataAsDouble(0);
	double deadtime = sampling/2;

	//data.saveAs("data.qtr");
	//data.close();

	//Setup MIL configuration from givenp arameters
	data["ExpCond"]["DeadTime"].setData(QTR_TYPE_DOUBLE, deadtime);

	QUB_Tree config( QUB_Tree::Create("Properties") );
	config.children().insert( model.clone() );
	config["MaxIterations"].setData(QTR_TYPE_INT, max_iter); 
	config["Mode"].setData((max_iter > 0) ? "optimize" : "evaluate");
	config["use segments"].setData("together");
	config["SearchLimit"].setData(QTR_TYPE_DOUBLE, 10.0);
	config["DFP"]["MaxStep"].setData(QTR_TYPE_DOUBLE, 0.1);
	config["DFP"]["ConvLL"].setData(QTR_TYPE_DOUBLE, 0.0001);
	config["DFP"]["ConvGrad"].setData(QTR_TYPE_DOUBLE, 0.01);

	//Create data list
	QTR_Impl* adata[2];
	adata[0] = data.getImpl();
	adata[1] = NULL;

	//Run MIL
	result = miltreeiface(config.getImpl(), adata, NULL, NULL, NULL);
	QTR_DECREF( result.getImpl() );

	return result;
}

