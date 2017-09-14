


// QUB header files (qub/qubopt)
#include "qubopt.h"
#include "mdlrep.h"
#include "QUB_Tree.h"
#include "mil_eval_tree.h"

#include <string>
#include <iostream>

using namespace std;


int main(int argc, char* argv[])
{
	//Parse input parameter list
	if( argc!= 4 )
	{
		cerr << "Incorrect number of arguments" << endl;
		return -1;
	}

	string configFilename(argv[1]);
	string dataFilename(argv[2]);
	string resultFilename(argv[3]);

	//Load input parameters
	QUB_Tree data( QUB_Idl_ReadDWT(dataFilename.c_str()) );
	if ( data.isNull() )
	{
		cerr << "Could not load DWT data" << endl;
		return -2;
	}

	QUB_Tree config = QUB_Tree::Open(configFilename);
	QUB_Tree result;


	cout << "Running";

	//Run MIL
	
	//Set deadtime
	double sampling = data["sampling"].dataAsDouble(0);
	double deadtime = sampling/2;
	data["ExpCond"]["DeadTime"].setData(QTR_TYPE_DOUBLE, deadtime);

	//Create data list
	QTR_Impl* adata[2];
	adata[0] = data.getImpl();
	adata[1] = NULL;
	//Run MIL
	result = miltreeiface(config.getImpl(), adata, NULL, NULL, NULL);
	

	//finished.
	cout << "Saving results..." << endl << fflush;
	

	//Save the result data.
    //Saved first to a temporary filename, then renamed. This prevents MATLAB
    //from openning the file before it is completely written.
    string tempName = resultFilename + ".tmp";
    
	result.saveAs( tempName );
	result.close();
	QTR_DECREF( result.getImpl() );

    rename( tempName.c_str(), resultFilename.c_str() );


	cout << result["LL"].dataAsDouble(-7.0) << endl;
	cout << "Finished." << fflush;

	return 0;
}
