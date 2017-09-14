/* Copyright 1998-2012 Research Foundation State University of New York */

/* This file is part of QUB Express.                                     */

/* QUB Express is free software; you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation, either version 3 of the License, or     */
/* (at your option) any later version.                                   */

/* QUB Express is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */

/* You should have received a copy of the GNU General Public License,    */
/* named LICENSE.txt, in the QUB Express program directory.  If not, see */
/* <http://www.gnu.org/licenses/>.                                       */

#ifndef MIL_EVAL_TREE_H
#define MIL_EVAL_TREE_H

#include "qubopt.h"
#include "CountedPtr.h"
#include "dfp_optimize.h"
#include "modelTree.h"

extern "C" QUBOPT_API
QTR_Impl * QUB_Idl_ReadDWT(const char *filename);

	class mil_model{
		public:
			mil_model( QUB_Tree modelNode, double searchLimit ); // changes a copy of modelNode

			void reset(); // re-copy original

			void updateRates(); // from mi.x into node
			void updateRates( fq::vector<double>& xsd );
			
			QUB_Tree  original;
			QUB_Tree  node; // the one used by tree
			ModelTree tree;
			double    searchLimit;
			mdlinf    mi;
			int       nconstraint;
			int       nclass;
			int       nchannel;

		private:
			void normalizePr();
		};

	class mil_metamodel{
		public:
			mil_metamodel( mil_model& mdl, int nchan );

			mil_model&       model;
			int              nchannel;

			// generated on construction:
			fq::vector<int>    idrug, ivolt; // idrug[i] = mi.path[i][2], ivolt[i] = mi.path[i][3]

			int              nmetastate, nmetaclass, nmetapath;
			fq::vector<int>    metaclass, nmetagroup;
			fq::vector<float>  metaamp;
			fq::vector<double> mprx;
			fq::matrix<int>      metastate, metapath, subpath, metaindex;
			int              maxnmetastate, maxnmetapath;
		};

	void UpdateResult( QUB_Tree result, double ll, fq::vector<double>& gz );
	void BlessResult( QUB_Tree result );

	class mil_workspace { 
			// can evaluate segments and make histograms for one file (dataset with same conditions)
		public:
			mil_workspace( QUB_Tree dataset, mil_model& mdl, bool usePequil, bool smoothHistoBins, bool sqrtHisto, int *StopFlag, std::ostream& errstream );

			void reset(); // back to original rates

			int segCount();
			int classCount(); // in metamodel
			int classCountInData();
			int componentsInClass( int cls );

			// evaluate(...) is where the action is.
			// evaluate multi workspaces = sum( eachWorkspace.evaluate() )
			// call updateZ() first if z has changed
			// returns number of failed segments
			int evaluate( fq::vector<double>& gz, double& ll, int flagGrad, int seg=-1 ); // default: all segs

			// provide the nodes that hold the "Histogram" and "TimeConstants" elements.
			// you can use them repeatedly.
			void getHistograms( QUB_Tree histGroup, int seg=-1 ); // bins n' bars
			// you provide # of bins in histGroup["BinCount"]
			void getPDFs( QUB_Tree histGroup, QUB_Tree timeConstContainer ); // pdf and TCs for given bins

			QUB_Tree getModelStats(double tMax, int tBins);

			bool resetZ( bool alsoResetModel=true ); // if false, you get the modified rates
			bool updateZ( fq::vector<double>& newz );
			void updateModel(); // update tree from most recent z
			void updateSD( fq::matrix<double>& hessian ); // into the tree

			fq::vector<double> z;      // constrained parameters
			bool             zValid; // false if the initial guess didn't meet constraints

			QUB_Tree         data; // with dead time applied and subtracted
			std::vector<QUB_Tree> segments;

		protected:
			void copyMetaData( QUB_Tree orig );
			void copyData( QUB_Tree orig );
			void copyDataJoined( QUB_Tree orig );

			void setupArrays();

			void setTimeConstHint( QUB_Tree hist, QUB_Tree tc );

			void calHist( int& nbin, float *bins, float *hist, int cls, int seg=-1 ); // can decrease nbin
			bool calPDFandTC( int nbin, float *bins, float **pdf, double *tau, double *amp, double& meanTau, int cls );

			void peqm(double **q, int ic, double *pe);
			bool pdfcof(double **q, int ic, double *pe, double *c, double *wr);
			void pdfbin(int ic, double *w, double *c, double td, int nbin, float* bin, float **pdf);

			void getModelStats(double tMax, int tBins, double **q, string suffix, QUB_Tree result);

			void setupDQDQ();          // path            -> dQdq
			void calcStartProbs();     // q, dQdq         -> pr, dpr_q
			void calcMultiProbs();     // pr, metastate   -> mpr, dmpr_pr
			void calc_dlogLL_q( int cls, fq::vector<double>& dLL_q );

			int              maxndwt;
			std::vector<int>      ndwt;
			std::vector<int*>     idwt;
			std::vector<float*>   tdwt;
			int              classesInData;

			mil_metamodel    model;
			bool             usePeq; // for start probs

			int              *Stopped;
			std::ostream&    milerr;

			fq::vector<double> C; // , V;  C now has concs and voltages intermingled
			double           td;

			// constraints: x = mtx * z + vct
			fq::matrix<double>   mtx;
			fq::vector<double> vct;

			int              posieigencount;
			fq::vector<double> df_q, df_x, df_z, df_qq;
			fq::matrix<double>   q, qq, qqe, dx_z, dq_x, dqq_q, df_qqe;
			fq::tensor<double>   dqqe_qq;
			fq::matrix<double*>  dqqe_qqAsPtr;
			fq::matrix<double>   dmpr_pr, dpr_q, dmpr_q;
			fq::vector<double> v;
			fq::matrix<double>   w, a, c, c1;
			fq::tensor<double>   p, p1;
			fq::tensor<double>   dQdq; // as in MIP, for dpr_q

			fq::matrix<double>   alpha, beta;
			fq::vector<double> scale;

			bool             smoothHistBins, sqrtHist;
			double           histBinR;
		};

	typedef CountedPtr<mil_workspace> mil_worksptr;


	class mil_eval_tree : public dfp_optimizable{
		public:
			mil_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
							QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB );

			QUB_Tree execute();

			virtual bool getStartingZ( fq::vector<double> &z );
			virtual bool evaluate( fq::vector<double> &z, fq::vector<double> &gz,
									double &ll, int &nf, int &ndf, bool flagGrad=true );
			virtual bool checkContinue( int iter, int nf, int ndf, double ll,
										fq::vector<double> &z, fq::vector<double> &gz );

		protected:
			void initResults();
			void addTCSelectVars( QUB_Tree seg );
			void addRateSelectVars( QUB_Tree seg );
			void finishResults();

			void updatePDFs();
			void sendResult( QUB_Tree resNode, int iter, bool final=false, int errCode = 0 );

			bool canRun();
			bool evaluateAll( fq::vector<double> &z, fq::vector<double> &gz,
							  double &ll, int &nf, int &ndf, bool flagGrad=true );

			dfp_result runDFP();

			mil_model model;
			std::vector<mil_worksptr> workspaces;
			QUB_Tree cfg;
			QTR_Callback reportCB, resultCB, percentCB;

			mil_reportstream milerr;

			bool opt, together, pdfsEachIter;
			int maxIter, maxRestarts;
			double maxStep, convLL, convGrad;
			
			int *Stopped;

			mil_worksptr currSpace; // for together == false
			int currSeg, iSeg, nSeg;
			int iterations;
			
			QUB_Tree result;
			QUB_Tree pctNode;
			int *pctPtr;
		};

#ifdef __cplusplus
extern "C" {
#endif

	QUBOPT_API QTR_Impl * miltreeiface( QTR_Impl *cfg, QTR_Impl **data,
				  QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB );

	typedef struct {
		int              nclass;
		int              *classOfState; // aka mi->clazz

		int              nmetastate;
		int              nmetaclass;
		int              nmetapath;
		int              maxnmetastate;
		int              maxnmetapath;

		int              *metaclass;
		int              *nmetagroup;
		double           *metapr;
		int              **metastate;
		int              **metapath;
		int              **subpath;
		int              **metaindex;

		mil_model        *model;
		mil_metamodel    *metamodel;
		} mil_export_model;

	QUBOPT_API mil_export_model * new_mil_model( QTR_Impl *model, int nchannel );
	QUBOPT_API void del_mil_model( mil_export_model *mdl );

#ifdef __cplusplus
}
#endif



#endif
