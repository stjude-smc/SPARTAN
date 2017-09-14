/* Copyright 1998-2014 Research Foundation State University of New York */

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

#include "mil_eval_tree.h"

#include "histpdf.h"
#include "mdlrep.h"
#include "qmatrixutil.h"
#include "qublib.h"

#include <ctime>
#include <float.h>
#include <stdexcept>
#include <string.h>
#include <iostream>

using namespace fq;

#ifndef _WIN32
  #define _isnan isnan
  #define _finite finite
#endif

//#include <boost/math/special_functions/next.hpp> // compiler trouble; replacement in milutil
#define TDEAD_DELTA_FLOATS 5

//----- local functions 
void mil_expQ(int ic, int n, double** w, tensor<double>& p, tensor<double>& p1,  double t, double** a);
void ApplySearchLimit( mdlinf& mi, double searchLimit );
void CleanResultFloats( QUB_Tree node );

void mil_xtoq( int nstate, int npath, int **path, int *idrug, double *C, 
			   int *ivolt, double *V, double *x, double **q, double **dq_x ){
	int	i,j,k,l;
	double nQ;
	
	for (i=0; i<nstate; i++)			// Fill q with 0's
		for (j=0; j<nstate; j++) 
			q[i][j]=0.0;

	for (k=0; k<npath; k++) {			// Set q[state][state] from path array 
	    i=path[k][0];
	    j=path[k][1];
		nQ = C[idrug[k]] * exp( x[2*k] + V[ivolt[k]] * x[2*k+1] );

	    q[i][i] -= (q[i][j]=nQ) ;		// diagonal = -sum(row)
		for (l=0; l<npath; l++)			// Set dq_x to 0's 
			dq_x[k][2*l] = dq_x[k][2*l+1] = 0;
		dq_x[k][2*k] = nQ ;				// Set Derivatives for diagonal 
		dq_x[k][2*k+1] = V[ivolt[k]] *	nQ; 
		}


	}

static QTR_Impl * QTR_NullCallback(QTR_Impl *) {
  return NULL;
}

static QUB_Tree QTR_DoCallback(QTR_Callback cb, QUB_Tree arg) {
  QTR_Impl *impl = NULL;
  QUB_Tree result;
  if ( cb ) {
    impl = cb(arg.getImpl());
    if ( impl ) {
      result = impl;
      QTR_DECREF(impl);
    }
  }
  return result;
}

// *************************************************************************

mil_model::mil_model( QUB_Tree modelNode, double search_limit )
	: original( modelNode ), node( modelNode.clone() ), tree( node ),
	  searchLimit( search_limit )
{
	node.setData( "Final" );
	tree.condenseClasses();
	nconstraint = tree.toMdlinf( mi );
	nclass = tree.classesPresent();
	nchannel = node.find("ChannelCount")->dataAsInt(1);
	normalizePr();
	ApplySearchLimit( mi, searchLimit );
	tree.restoreClasses();
}

void mil_model::reset()
{
	node = original.clone();
	node.setData( "Final" );
	tree = ModelTree( node );
	tree.condenseClasses();
	nconstraint = tree.toMdlinf( mi );
	normalizePr();
	ApplySearchLimit( mi, searchLimit );
	tree.restoreClasses();
}

void mil_model::updateRates()
{
	tree.updateRates( mi );
}

void mil_model::updateRates( fq::vector<double>& xsd )
{
	tree.updateRates( mi, xsd );
}

// If sum<=0 set all = 1/n otherwaise scale all to sum=1
void mil_model::normalizePr() {
	double sumprob = 0.0;
	int i;

	for ( i=0; i<mi.pr.size(); ++i )
		sumprob += fabs(mi.pr[i]);

	if ( sumprob <= 0.0 ) {
		sumprob = mi.pr.size();
		for ( i=0; i<mi.pr.size(); ++i )
			mi.pr[i] = 1.0/sumprob;
		}
	else 
		for( i=0; i<mi.pr.size(); ++i )
			mi.pr[i] = fabs(mi.pr[i]) / sumprob;
	}

// *************************************************************************

mil_metamodel::mil_metamodel( mil_model& mdl, int nchan )
		: model( mdl ), 
		  nchannel( nchan ), 
		  idrug( mdl.mi.npath ), 
		  ivolt( mdl.mi.npath ) {

	mdlinf& mi = mdl.mi;
	int i, j, k;
	
	for ( i=0; i<mi.npath; ++i ) {
		idrug[i] = mi.path[i][2];
		ivolt[i] = mi.path[i][3];
		}

	int nstate = mi.nstate;

	// generate multi-channel model
	double tmp = 0.0;
	int nn = nchannel + 1;
	for ( ; nn <= nstate + nchannel - 1; ++nn )
		tmp += log( (double) nn );
	for ( nn=1; nn<=nstate-1; ++nn )
		tmp -= log( (double) nn );
	nn = (int) exp( tmp ) + 10;

	maxnmetastate = 0;
	maxnmetapath = 0;

	if ( nstate > 0 ) {
		// iterate through combinations of states / channels 
		matrix<int> cc( nn, nstate );
		CChannelCombine oMeta(nchannel, nstate);
		do {
			for ( j=0; j<nstate; ++j ){
				cc[maxnmetastate][j] = oMeta.State0(j);
				}
			++maxnmetastate;
			} while( oMeta.GetNextCombine() );

		int s;
		for ( i=0; i<maxnmetastate; ++i ) {
			for ( j=0; j<maxnmetastate; ++j ) {
				for ( s=k=0; k<nstate; ++k )
					s += abs( cc[i][k] - cc[j][k] );
				if ( s == 2 )
					++maxnmetapath;
			}
		}
	}

	datinf di;
	mdl.tree.toDatinf(di);
	fq::vector<float> ratio(mdl.nclass);
	for ( i=0; i<mdl.nclass; ++i )
		ratio[i] = float(di.i[i] - di.i[0]);

	metastate.resize( maxnmetastate, nstate );
	metaclass.resize( maxnmetastate );
	metapath.resize( maxnmetapath, 2 );
	subpath.resize( maxnmetapath, 2 );
	nmetagroup.resize( maxnmetastate );
	metaindex.resize( maxnmetastate, 2*maxnmetastate );
	mprx.resize( maxnmetastate );
	metaamp.resize( maxnmetastate );

	metaMarkov( nchannel, nstate, mi.npath, mi.path, mi.clazz, ratio, &nmetastate,
			    &nmetapath, &nmetaclass, metaamp, metastate, metapath, metaclass, subpath );
	aggregate( nmetastate, nmetaclass, metaclass, nmetagroup, metaindex );

	for ( tmp=0.0, i=0; i<nstate; ++i )
		tmp += mi.pr[i];
	if ( tmp < 1.0e-5 ) {
		for ( i=0; i<nstate; ++i )
			mi.pr[i] = 1.0/(double)nstate;
	}
	CalcMultiPr( nchannel, nstate, nmetastate, metastate, mi.pr, mprx );
}

void ApplySearchLimit( mdlinf& mi, double searchLimit ){
	// on the importance of the search limit:
	// the dfp optimizer gets a gradient, then goes lnsrch in that direction.
	// lnsrch: try a few stabs in that direction until you find something good enough
	// i've heard it zeroes in by assuming the likelihood space is parabolic.
	// when it tries rates that are so far off that the likelihood space is way not parabolic,
	// lnsrch gets confused and gets stuck in the wrong direction.
	// thus, rates outside the search limit are rejected out of hand, in mil_ztox.
	double logSL = log( searchLimit );

	int n = mi.x.size();
	for ( int i=0; i<n; i += 2 ) {
		mi.xlimit[0][i] = mi.x[i] - logSL;
		mi.xlimit[1][i] = mi.x[i] + logSL;
		mi.xlimit[0][i+1] = mi.x[i+1] / searchLimit; // should this be "- logSL" instead?
		mi.xlimit[1][i+1] = mi.x[i+1] * searchLimit;
		}
	}


// *************************************************************************

void mil_expQ(int ic, int n, double** w, tensor<double>& p, tensor<double>& p1, 
				 double t, double** a) {
	int	i,j,k;
	double	tmp;
	
	for (i=0; i<n; i++) 
		for (j=0; j<n; j++) 
			a[i][j]=0.0;
	
	for (k=0; k<n; k++) {
		tmp=exp(w[ic][k]*t);
	    for (i=0; i<n; i++) 
		    for (j=0; j<n; j++) 
				a[i][j] += tmp*p[ic][i][k]*p1[ic][k][j];
		}
	}

// similar to above, extract full A matrix for model report, not submatrix as above
void mil_expQ_full(int n, double* w, double** p, double** p1, double t, double** a){
	int i, j, k;
	double tmp;

	for ( i=0; i<n; ++i )
		for ( j=0; j<n; ++j )
			a[i][j] = 0.0;

	for ( k=0; k<n; ++k ) {
		tmp = exp(w[k]*t);
		for ( i=0; i<n; ++i )
			for ( j=0; j<n; ++j )
				a[i][j] += tmp * p[i][k] * p1[k][j];
		}
	}
	    

double mil_gama(int ic, int k, int l, double t, double** w) {
   /* the calculation here assumes the two eigenvalues are not equal,
      but they may actually be. We add a second check "w1==w2" in addition
      to k=l, this needs to be proved. */
	double y;
	
	double w1 = w[ic][k];
	double w2 = w[ic][l];
	if (k==l) 
		y=t*exp(w1*t);
	else 
		y=(exp(w1*t)-exp(w2*t))/(w1-w2);
	
	return y;
	}

void mil_fbward( int ndwell, int *idwell, float *tdwell, double** q,
				double** w, tensor<double>& p, tensor<double>& p1,
				double** alpha, double** beta,
				double* scale, int flagBack, mil_metamodel& model, int *Stopped ){
	int i,j,k,m,n,ic,jc,k2;
	double tmp;
	
	fq::vector<double> xV(model.nmetastate);
	fq::vector<double> yV(model.nmetastate);
	matrix<double> aM(model.nmetastate,model.nmetastate);
	fq::vector<double> pstartV(model.nmetastate);
	
	double *x = xV;
	double *y = yV;
	double **a = aM;
	double *pstart = pstartV;
	
	/* scaled alpha[i][0:2*ndwell] and scale[0:2*ndwell] */    
	
	ic = idwell[0];
	jc = 0;
	m  = model.nmetagroup[ic];
	
	/* starting pr of the states in the first dwell */ 
	for (tmp=0, i=0; i<m; i++)
		tmp += model.mprx[model.metaindex[ic][i]];
	for (i=0; i<m; i++)
		pstart[i] = model.mprx[model.metaindex[ic][i]] / tmp;
	
	for (scale[0]=0, i=0; i<m; i++) { 
		alpha[i][0] = pstart[i];
		scale[0]   += alpha[i][0];
		}
	scale[0]=1.0/scale[0];
	for (i=0; i<m; i++) 
		alpha[i][0] *= scale[0];
	
	for (k=0; k<ndwell; k++) {
		if (Stopped && *Stopped)
			break;
		
		k2 = 2*k;
		mil_expQ(ic,m,w,p,p1,tdwell[k],a);
		for (j=0; j<m; j++) 
			x[j]=alpha[j][k2];
		vxm(m,m,a,x,y);
		tmp = dsumv(y,m);
		if (tmp<1.0e-10)
			tmp = 1.0e-10;
		scale[k2+1]=1.0/tmp;
		for (j=0; j<m; j++) 
			alpha[j][k2+1]=y[j]*scale[k2+1];
		
		if (k<ndwell-1) {
			jc = idwell[k+1];
			n  = model.nmetagroup[jc];
			for (i=0; i<m; i++)
				for (j=0; j<n; j++) 
					a[i][j]=q[model.metaindex[ic][i]][model.metaindex[jc][j]];
			}
		else {
			n = model.nmetastate - model.nmetagroup[ic];
			for (i=0; i<m; i++)
				for (j=0; j<n; j++) 
					a[i][j]=q[model.metaindex[ic][i]][model.metaindex[ic][model.nmetastate+j]];
			}
		for (j=0; j<m; j++) 
			x[j]=alpha[j][k2+1];
		vxm(m,n,a,x,y);
		tmp = dsumv(y,n);
		if (tmp<1.0e-10)
			tmp = 1.0e-10;
		scale[k2+2]=1.0/tmp;
		for (j=0; j<n; j++) 
			alpha[j][k2+2]=y[j]*scale[k2+2];
		
		ic = jc;
		m  = n;
		}
	
	
	/* scaled beta[i][t]  */
	if (flagBack==1) {
		jc = idwell[ndwell-1];
		n  = model.nmetastate - model.nmetagroup[jc];
		for (j=0; j<n; j++) beta[j][2*ndwell]=scale[2*ndwell];
		
		for (k=ndwell-1; k>=0; k--) {
			if (Stopped && *Stopped)
				break;
			
			k2 = 2*k;
			for (j=0; j<n; j++) x[j]=beta[j][k2+2];
			ic = idwell[k];
			m  = model.nmetagroup[ic];
			if (k<ndwell-1) {
				for (i=0; i<m; i++)
					for (j=0; j<n; j++)
						a[i][j]=q[model.metaindex[ic][i]][model.metaindex[jc][j]];
				}
			else {
				for (i=0; i<m; i++)
					for (j=0; j<n; j++)
						a[i][j]=q[model.metaindex[ic][i]][model.metaindex[ic][model.nmetastate+j]];
				}
			mxv(m,n,a,x,y);
			for (j=0; j<m; j++) { 
				if ( (y[j] > 1.0e120) && (scale[k2+1] > 1.0) ) // Chris attempt to avoid overflow 11-19-02
					beta[j][k2+1] = 1.0e120;
				else if ( (y[j] < 1.0e-120) && (scale[k2+1] < 1.0) )
					beta[j][k2+1] = 1.0e-120;						// Chris attempt to avoid overflow 5.21.02 rev 11-19-02
				else
					beta[j][k2+1]=y[j]*scale[k2+1];               // float overflow
				}
			
			for (j=0; j<m; j++) x[j]=beta[j][k2+1];
			mil_expQ(ic,m,w,p,p1,tdwell[k],a);
			mxv(m,m,a,x,y);
			for (j=0; j<m; j++) {
				if ( (y[j] > 1.0e120) && (scale[k2] > 1.0) ) // Chris attempt to avoid overflow 11-19-02
					beta[j][k2] = 1.0e120;
				else if ( (y[j] < 1.0e-120) && (scale[k2] < 1.0) )
					beta[j][k2] =1.0e-120;						// Chris attempt to avoid overflow 5.21.02 rev 11-19-02
				else
					beta[j][k2]=y[j]*scale[k2];                   // float overflow
				}
			
			jc = ic;
			n  = m;
			}
		} 
	}

// working theory 11-19-02:
// the overflows come from wildly asymmetrical rates:
//   e.g. k0 1->2 is big, k0 2->1 is small
// when they're similar, the same scale factors apply forward and backward,
// but discrepancies accumulate as beta is generated.
// multiply in 1.0e3 at each event, and before you know it you've passed the expo limit of IEEE double precision
void mil_dlogLL_q( double** df_qe, tensor<double>& dqe_q, double* df_q, mil_metamodel& model ){
	int i,j,k;
	double s;
	
	for (k=0; k<model.nmetapath; k++) {
		s = 0.0;
		for (i=0; i<model.nmetastate; i++) 
			for (j=0; j<model.nmetastate; j++)
				s += df_qe[i][j]*dqe_q[i][j][k];
			df_q[k]=s;
		}
	}

void mil_dlogLL_qe( double** w, tensor<double>& p, tensor<double>& p1,
			  int ndwell, int *idwell, float *tdwell,
			  double** alpha, double** beta, double** df_q,
			  mil_metamodel& model, int *Stopped ){
	int	 i,j,k,l,m,n,t,ic,jc;
	double tmp;

	// care-free stack-based memory management
	matrix<double> aM (model.nmetastate, model.nmetastate);
	matrix<double> bM (model.nmetastate, model.nmetastate);
	matrix<double> cM (model.nmetastate, model.nmetastate);
	matrix<double> skM(model.nmetastate, model.nmetastate);
	matrix<double> slM(model.nmetastate, model.nmetastate);

	// fast access
	double **a = aM;
	double **b = bM;
	double **c = cM;
	double **sk = skM;
	double **sl = slM;
  
	for (ic=0; ic<model.nmetaclass; ic++) {
		m=model.nmetagroup[ic];
		for (dzerom(m,m,sk),k=0; k<m; k++) {
			for (dzerom(m,m,sl),l=0; l<m; l++) {
				if (Stopped && *Stopped)
					break;
				for (dzerom(m,m,a),t=0; t<ndwell; t++) 
					if (idwell[t]==ic) {
						tmp=mil_gama(ic,k,l,tdwell[t],w);
						for (i=0; i<m; i++) 
							for (j=0; j<m; j++) 
								a[i][j] += alpha[i][2*t]*beta[j][2*t+1]*tmp;
						}
					for (i=0; i<m; i++) 
						for (j=0; j<m; j++) 
							b[i][j]=p[ic][j][l]*p1[ic][l][i];           
						mxm(m,m,m,a,b,c);
						for (i=0; i<m; i++) 
							for (j=0; j<m; j++) 
								sl[i][j] += c[i][j];
				}
			for (i=0; i<m; i++) 
				for (j=0; j<m; j++) 
					b[i][j]=p[ic][j][k]*p1[ic][k][i]; 
				mxm(m,m,m,b,sl,c);
				for (i=0; i<m; i++) 
					for (j=0; j<m; j++) 
						sk[i][j] += c[i][j];
			}
		for (i=0; i<m; i++) 
			for (j=0; j<m; j++) {
				k=model.metaindex[ic][i];
				l=model.metaindex[ic][j];
				df_q[k][l] = sk[i][j];
				}
			
			for (jc=0; jc<model.nmetaclass; jc++) {
				if (jc==ic) continue;
				n=model.nmetagroup[jc];
				if (Stopped && *Stopped)
					break;
				for (dzerom(m,n,a),t=0; t<ndwell-1; t++) {
					if (idwell[t]==ic && idwell[t+1]==jc) {
						for (i=0; i<m; i++) 
							for (j=0; j<n; j++) 
								a[i][j]+=alpha[i][2*t+1]*beta[j][2*t+2];
						}
					}
				if (idwell[ndwell-1]==ic) {
					for (i=0; i<m; i++) 
						for (j=0; j<n; j++) 
							a[i][j]+=alpha[i][2*t+1]*beta[j][2*t+2];
					}
				for (i=0; i<m; i++) 
					for (j=0; j<n; j++) {
						k=model.metaindex[ic][i];
						l=model.metaindex[jc][j];
						df_q[k][l] = a[i][j];
						}
				}
		}
	}

// *************************************************************************
mil_workspace::mil_workspace( QUB_Tree dataset, mil_model& mdl, bool usePequil, bool smoothHistoBins, bool sqrtHisto, int *StopFlag, std::ostream& errstream )
	: data( QUB_Tree::Create("DataSet") ),
	  classesInData( -1 ), 
	  model( mdl, (int) dataset["ExpCond"].find("ChannelCount")->dataAsDouble( mdl.nchannel ) ),
	  usePeq( usePequil ),
	  Stopped( StopFlag ), 
	  milerr( errstream ),
	  mtx( 2 * mdl.mi.npath + 1, 2 * mdl.mi.npath + 1 ), 
	  vct( 2 * mdl.mi.npath + 1 ),
	  posieigencount( 0 ),
      smoothHistBins( smoothHistoBins ),
	  sqrtHist( sqrtHisto ) {

	// V    = dataset["ExpCond"]["Voltage"].dataAsDouble( 0.0 );
	td   = dataset["ExpCond"]["DeadTime"].dataAsDouble( 1.0e-9 );

	//C[0] = 1.0;
	//C[1] = dataset["Concentration"].dataAsDouble( 1.0 );

	C.resize( model.model.tree.nparam + 2 );
	C[0] = 1.0;
	C[1] = 0.0;
	for ( QUB_TreeMonoIter ligi = dataset["ExpCond"].children(); ! (*ligi).isNull(); ++ligi ) {
		int ix = model.model.tree.indexOfParam[ ligi->name() ];
		if ( ix )
			C[ix] = ligi->dataAsDouble( 1.0 );
		}

	if ( (* dataset.find("Join Segments") ).isNull() )
		copyData( dataset );
	else
		copyDataJoined( dataset );

	setupArrays();
	}

void mil_workspace::reset(){
	resetZ();
	}

int mil_workspace::segCount(){
	return (int) ndwt.size();
	}

int mil_workspace::classCount(){ // in metamodel
	return model.nmetaclass;
	}

int mil_workspace::classCountInData(){
	if ( classesInData < 0 ) {
		for ( int i=0; i<int(ndwt.size()); ++i ) {
			int *idwell = idwt[i];
			for ( int j=ndwt[i]-1; j>=0; --j ) {
				if ( idwell[j] > classesInData )
					classesInData = idwell[j];
				}
			}
		++classesInData; // max class present -> # of classes
		}
	return classesInData;
	}

int mil_workspace::componentsInClass( int cls ){
	return model.nmetagroup[ cls ];
	}

#pragma warning (disable: 4101)	// Disable unreferenced variables message 
int mil_workspace::evaluate( fq::vector<double>& gz, double& ll, int flagGrad, int seg ){
	int rtnVal = 0;
	
	std::vector<int> segs;
	if ( seg < 0 )
		for ( int i=0; i<int(ndwt.size()); ++i )
			segs.push_back( i );
	else if ( seg < int(ndwt.size()) )
		segs.push_back( seg );

	ll = 0.0;
	if( flagGrad )
		dzerov( gz.size(), &(gz[0]) );

	int i, j;
	//std::clock_t before = std::clock();
	double sum;
	int nx = model.model.mi.x.size();
	int nz = gz.size();

	for ( std::vector<int>::iterator segi = segs.begin(); segi != segs.end(); ++segi ) {
		try {
			int ndwell = ndwt[ *segi ];
			if ( ndwell == 0 ) continue;

			int *idwell = idwt[ *segi ];
			float *tdwell = tdwt[ *segi ];

			/* check the total starting pr of the states int the first dwell */
			for (sum=0, i=0; i<model.nmetagroup[idwell[0]]; i++)
				sum += model.mprx[model.metaindex[idwell[0]][i]];
			if (sum < 1.0e-5) {  /* discard this segment */
				milerr << "Segment " << *segi << " starts in a class with zero start probability...discarding it." << endl;
				continue;
			}

			mil_fbward(ndwell,idwell,tdwell,qqe,w,p,p1,alpha,beta,scale,flagGrad,model,Stopped);

			for (j=0; j<=2*ndwell; j++)
				ll -= -log(scale[j]); // -=: to maximize using the minimizer // -log: log(1/x)

			if ( flagGrad ) {
				mil_dlogLL_qe(w,p,p1,ndwell,idwell,tdwell,alpha,beta,df_qqe,model,Stopped);
				mil_dlogLL_q(df_qqe,dqqe_qq,df_qq,model);
				vxm(model.nmetapath,model.model.mi.npath,dqq_q,df_qq,df_q);


				//fq::vector<double> df_q_fromPeq;
				//calc_dlogLL_q( idwell[0], df_q_fromPeq );
				//segments[ *segi ]["df_q_peq"].setNumData(QTR_TYPE_DOUBLE, model.model.mi.npath, 1, &(df_q_fromPeq[0]));
				//for ( j=model.model.mi.npath-1; j>=0; --j )
				//	df_q[j] += df_q_fromPeq[j];


				vxm(model.model.mi.npath,nx,dq_x,df_q,df_x);
				vxm(nx,nz,dx_z,df_x,df_z);
				for (j=0; j<nz; j++)
					gz[j] -= df_z[j];
			}
		}
		catch (const std::overflow_error &oe) {
#ifdef _WIN32
			unsigned int fstat = _clearfp(); // == 0x80007; should be 0
			_fpreset();
#endif

			milerr << "Numerical Overflow: segment " << (*segi+1) << endl;
			++rtnVal;
		}
		catch (const std::underflow_error &oe) {
#ifdef _WIN32
			unsigned int fstat = _clearfp();
			_fpreset();
#endif

			milerr << "Numerical Underflow: segment " << (*segi+1) << endl;
			++rtnVal;
		}
		catch (const std::exception &ee) {
#ifdef _WIN32
			unsigned int fstat = _clearfp(); // just in case
			_fpreset();
#endif
			
			milerr << "Exception: " << ee.what() << ": segment " << (*segi+1) << endl;
			++rtnVal;
		}
		catch (...) {
#ifdef _WIN32
			unsigned int fstat = _clearfp(); // just in case
			_fpreset();
			
			if ( fstat & _SW_OVERFLOW )
				milerr << "Numerical Overflow: segment " << (*segi+1) << endl;
			else if ( fstat & _SW_UNDERFLOW )
				milerr << "Numerical Underflow: segment " << (*segi+1) << endl;
			else if ( fstat & _SW_ZERODIVIDE )
				milerr << "Division by Zero: segment " << (*segi+1) << endl;
			else
#endif
				milerr << "Exception in MIL Core: segment " << (*segi+1) << endl;
			
			++rtnVal;
		}
	}
/*	
	if ( model.nchannel > 1 || seg < 0 ) // guard against 1channel only test code below
		return rtnVal;

	milerr << "old mil: " << ( ( std::clock() - before ) / (double)CLOCKS_PER_SEC ) << endl;
	before = std::clock();

	int Ns = model.nmetastate;
	fq::matrix<double> K0(Ns, Ns), K1(Ns, Ns);
	fq::matrix<int> Ligand(Ns, Ns), Voltage(Ns, Ns);
	for (int k=0; k<model.model.mi.npath; k++) {			// Set q[state][state] from path array 
	    int st_i = model.model.mi.path[k][0];
	    int st_j = model.model.mi.path[k][1];
		K0[st_i][st_j] = exp(model.model.mi.x[2*k]); // warning: can't handle channelcount > 1
		K1[st_i][st_j] = model.model.mi.x[2*k+1];
		Ligand[st_i][st_j] = model.model.mi.path[k][2];
		Voltage[st_i][st_j] = model.model.mi.path[k][3];
	}

	double new_ll = -666;
	int new_rtn = inter_ll(Ns, model.metaclass, model.mprx,
		K0, K1, Ligand, Voltage, C,
		td*1e-3, 1, &(ndwt[seg]), &(idwt[seg]), &(tdwt[seg]),
		&new_ll);

	milerr << "new mil: " << ( ( std::clock() - before ) / (double)CLOCKS_PER_SEC ) << '[' << new_ll << ']' << endl;
*/
	return rtnVal;
}
#pragma warning ( default : 4101)

void mil_workspace::getHistograms( QUB_Tree histGroup, int seg )
{
	QUB_TreeIter histi;
	int i;

	// clear existing
	while ( ! (* (histi = histGroup.find("Histogram")) ).isNull() )
		histi.remove();

	int nbin = histGroup["BinCount"].dataAsInt( 32 );

	for ( i = classCount() - 1; i >= 0; --i ) {
		QUB_Tree hist = QUB_Tree::Create("Histogram");
		histi.insert( hist ); // insert _before_ the prev. one
		
		ostringstream ost;
		ost << "Class " << (i+1);
		hist["Title"].setData( ost.str() );

		hist["XLabel"].setData( "duration [log10 ms]" );
		if ( sqrtHist )
			hist["YLabel"].setData( "sqrt(count / total)" );
		else
			hist["YLabel"].setData( "count / total" );

		fq::vector<float> tmpBins( nbin ), tmpBars( nbin );
		int actualBins = nbin;
		calHist( actualBins, tmpBins, tmpBars, i, seg );

		hist["BinCount"].setData( QTR_TYPE_INT, actualBins );
		hist["ActualBins"].setNumData( QTR_TYPE_FLOAT, actualBins, 1, &(tmpBins[0]) );
		hist["Bars"].setNumData( QTR_TYPE_FLOAT, actualBins, 1, &(tmpBars[0]) );

		for ( int j=0; j<actualBins; ++j )
			tmpBins[j] = float(log10( tmpBins[j] * 1.0e3 ));
		hist["Bins"].setNumData( QTR_TYPE_FLOAT, actualBins, 1, &(tmpBins[0]) );
	}
}

class VecForSortingByFirstElem
{
public:
	std::vector<double> vec;
	/*
	bool operator< (const VecForSortingByFirstElem& vec2) {
		return vec[0] < vec2.vec[0];
	}  */
};

bool operator<(const VecForSortingByFirstElem& vec1, const VecForSortingByFirstElem& vec2) {
     return vec1.vec[0] < vec2.vec[0];
}

void mil_workspace::getPDFs( QUB_Tree histGroup, QUB_Tree timeConstContainer )
{
	QUB_TreeMonoIter histi = histGroup.find("Histogram");
	QUB_TreeIter tci = timeConstContainer.find("TimeConstants");
	int i, j;

	for ( i=0; (i<classCount()) && ! (*histi).isNull(); ++i, histi.nextSameName(), tci.nextSameName() ) {
		QUB_Tree hist = *histi;
		QUB_Tree tc = *tci;

		int ncomp = componentsInClass( i );
		if ( tc.isNull() ) {
			tc = QUB_Tree::Create("TimeConstants");
			tc["Class"].setData( QTR_TYPE_INT, i );
			tc["Tau"].setNumData( QTR_TYPE_DOUBLE, ncomp, 1 );
			tc["Amp"].setNumData( QTR_TYPE_DOUBLE, ncomp, 1 );
			tci.insert( tc );
		}

		int nbins = hist["BinCount"].dataAsInt( 32 );
		if ( nbins ) {
			float *bins = (float *) hist["ActualBins"].data();
			matrix<float> pdf( ncomp, nbins );
			std::vector<float> sumpdf( nbins, 0.0 );

			/*
			QUB_Tree line = hist["Line"];
			if ( line.dataCount() != nbins ) {
				line.setNumData( QTR_TYPE_FLOAT, nbins, 1 );
				line["Color"].setData( QTR_TYPE_INT, i );
			}
			float *pdf = (float *) line.data();
			*/

			double *tau = (double *) tc["Tau"].data();
			double *amp = (double *) tc["Amp"].data();
			double meanTau;

			try {
				calPDFandTC( nbins, bins, pdf, tau, amp, meanTau, i );
				tc["MeanTau"].setData(QTR_TYPE_DOUBLE, meanTau);

				QUB_TreeIter lines;

				if ( ncomp > 1 ) {
					lines = hist.find("Component");
					for ( j=0; j<ncomp; ++j ) {
						if ( lines->isNull() ) {
							lines.insert("Component");
						}
						lines->setNumData( QTR_TYPE_FLOAT, nbins, 1, pdf[j] );
						lines.next("Component");
					}
				}

				if ( sqrtHist ) {
					for ( int b=0; b<nbins; ++b ) {
						for ( j=0; j<ncomp; ++j )
							sumpdf[b] += pdf[j][b] * pdf[j][b];
						sumpdf[b] = (float) sqrt(sumpdf[b]);
					}
				} else {
					for ( int b=0; b<nbins; ++b )
						for ( j=0; j<ncomp; ++j )
							sumpdf[b] += pdf[j][b];
				}
				lines = hist.find("Line");
				if ( lines->isNull() ) {
					lines.insert("Line");
					(*lines)["Color"].setData( QTR_TYPE_INT, i );
				}
				lines->setNumData( QTR_TYPE_FLOAT, nbins, 1, &(sumpdf[0]) );

				// sort time constants ascending
				std::vector<VecForSortingByFirstElem> tcvec( ncomp );
				for ( j=0; j<ncomp; ++j ) {
					tcvec[j].vec.push_back( tau[j] );
					tcvec[j].vec.push_back( amp[j] );
				}
				sort( tcvec.begin(), tcvec.end() );
				for ( j=0; j<ncomp; ++j ) {
					tau[j] = tcvec[j].vec[0];
					amp[j] = tcvec[j].vec[1];
				}

				setTimeConstHint( hist, tc );

				// add a VLine per time constant
				/*
				for ( j=0; j<ncomp; ++j ) {
					lines.next("VLine");
					if ( lines->isNull() ) {
						lines.insert("VLine");
						lines->setData( QTR_TYPE_DOUBLE, log10(tau[j]) );
					}
					else {
						lines->dataAs(0, (double)0.0) = log10(tau[j]);
					}
				}
				*/
			}
			catch (...) {
#ifdef _WIN32
				unsigned int fstat = _clearfp(); // just in case
				_fpreset();
#endif
				
				for ( j=0; j<nbins; ++j ) {
					for ( int k=0; k<ncomp; ++k )
						pdf[k][j] = 0.0f;
					sumpdf[j] = 0.0f;
				}
				for ( j=0; j<ncomp; ++j ) {
					tau[j] = amp[j] = 0.0;
				}
			}
		}
	}
}

// sort time constants decreasing
typedef pair<double, double> TimeConst;
bool TimeConstByDecAmp(const TimeConst& a, const TimeConst& b)
{
	return a.second > b.second;
}

// build a qubtree with model report info 
void mil_workspace::getModelStats(double tMax, int tBins, double **q, string suffix, QUB_Tree result)
{
	// big hypothesis:  Q* (qqe) can be manipulated just like an ordinary Q matrix

	double sampling = data["sampling"].dataAsDouble() * 1.0e-3;
	int n = model.nmetastate;
	int ic, i;//, j, k;

	result["Q"+suffix].setNumData(QTR_TYPE_DOUBLE, n, n, (void**)(double**)q);

	matrix<double> a(n+1, n+1);
	if ( mdlrep_AMatrix(q, n, sampling, a, milerr) )
		result["A"+suffix].setNumData(QTR_TYPE_DOUBLE, n, n, (void**)(double**)a);

	fq::vector<double> evals(n);
	matrix<double> evects(n, n);
	mdlrep_Eigen(q, n, evals, evects);
	result["Q_eigen"+suffix].setNumData(QTR_TYPE_DOUBLE, n, 1, evals);

	tensor<double> Ak(n, n, n);
	mdlrep_Spectrum(evals, evects, n, Ak, milerr);

	fq::vector<double> currB(n-1), currT(n-1), stateAmp(n);
	for ( i=0; i<n; ++i )
		stateAmp[i] = model.metaamp[ model.metaclass[i] ];
	mdlrep_Current(evals, Ak, model.mprx, stateAmp, n, currB, currT);
	result["I_B"+suffix].setNumData(QTR_TYPE_DOUBLE, n-1, 1, currB);
	result["I_T"+suffix].setNumData(QTR_TYPE_DOUBLE, n-1, 1, currT);
	result["I_T"+suffix].setLineComment("ms");

	result["P0"].setNumData(QTR_TYPE_DOUBLE, n, 1, model.mprx);

	fq::vector<double> Peq(n);
	mdlrep_Peq(q, n, Peq, milerr);
	result["Peq"+suffix].setNumData(QTR_TYPE_DOUBLE, n, 1, Peq);

	//matrix<double> FL(n, n);
	//MeanFirstPassage(a, Peq, sampling*1e3, n, FL);
	//result["MeanFirstPassage"+suffix].setNumData(QTR_TYPE_DOUBLE, n, n, (void**) (double**) FL);

	result["Ieq"+suffix].setData(QTR_TYPE_DOUBLE, mdlrep_Ieq(Peq, stateAmp, n));


	int nc = model.nmetagroup[0];
	fq::vector<double> flBins(tBins), flVals(tBins), flAmp(nc), flTau(nc);
	mdlrep_FirstLatency(q, model.mprx, model.metaclass, n, tMax, tBins, flBins, flVals, flAmp, flTau, milerr);
	result["FL_Amp"+suffix].setNumData(QTR_TYPE_DOUBLE, nc, 1, flAmp);
	result["FL_Tau"+suffix].setNumData(QTR_TYPE_DOUBLE, nc, 1, flTau);
	result["FL_Tau"+suffix].setLineComment("sec");
	result["FirstLatency"+suffix].setNumData(QTR_TYPE_DOUBLE, tBins, 1, flVals);
	QUB_Tree hist = result.appendChild("HistogramGroup")["Histogram"];
	hist["Title"].setData("First Latency"+suffix);
	hist["XLabel"].setData("t");
	hist["YLabel"].setData(" "); // ?
	hist["BinCount"].setData(QTR_TYPE_INT, tBins);
	hist["Bins"].setNumData(QTR_TYPE_DOUBLE, tBins, 1, flBins);
	hist["Bars"].setNumData(QTR_TYPE_DOUBLE, tBins, 1, fq::vector<double>(tBins, 0.0));
	hist["Line"].setNumData(QTR_TYPE_DOUBLE, tBins, 1, flVals);



	result["MeanTimeConst"+suffix].setNumData(QTR_TYPE_DOUBLE, model.nmetaclass, 1, (void*)0);
	double* meanTC = (double *) result["MeanTimeConst"+suffix].data();

	fq::vector<double> pe(n), w(n), c(n);

	for ( ic=0; ic<model.nmetaclass; ++ic ) {
		int ncomp = componentsInClass(ic);

		peqm(q,ic,pe); 
		pdfcof(q,ic,pe,c,&(w[0]));

		meanTC[ic] = 0.0;
		std::vector<TimeConst> tcs;
		for ( i=0; i<ncomp; ++i ) {
			TimeConst tc;
			tc.first = (w[i] == 0.0) ? 0.0 : (-1000.0 / w[i]);
			tc.second = (w[i] == 0.0) ? 0.0 : (-c[i] / w[i]);
			tcs.push_back(tc);
			meanTC[ic] += tc.first * tc.second;
		}
		sort(tcs.begin(), tcs.end(), TimeConstByDecAmp);

		QUB_Tree TCNode = result.appendChild("TimeConst"+suffix);
		TCNode.setNumData(QTR_TYPE_DOUBLE, ncomp, 2, (void*)0);
		double *TC = (double *) TCNode.data();
		for ( i=0; i<ncomp; ++i ) {
			TC[2*i+0] = tcs[i].first;
			TC[2*i+1] = tcs[i].second;
		}
	}
}

QUB_Tree mil_workspace::getModelStats(double tMax, int tBins)
{
	// big hypothesis:  Q* (qqe) can be manipulated just like an ordinary Q matrix

	QUB_Tree result = QUB_Tree::Create("Model Stats");
	
	QUB_Tree dataset = result["DataSet"];
	dataset.appendClone(data["FileName"]);
	dataset.appendClone(data["sampling"]);
	dataset.appendClone(data["ExpCond"]);

	getModelStats(tMax, tBins, qq, "", result);  // qq is multi channel, no dead time 
	getModelStats(tMax, tBins, qqe, "*", result);  // qqe * is with dead time 
	
	return result;
}

bool mil_workspace::resetZ( bool alsoResetModel )
{
	if ( alsoResetModel )
		model.model.reset();

	int i, j;
	int nconstraint = model.model.nconstraint;
	mdlinf& mi = model.model.mi;

	// generate initial z, change mtx and vct from [mtx * x + vct = 0] to [x = mtx * z + vct]
	for ( i=0; i<nconstraint; ++i ) {
		for ( j=0; j<2*mi.npath; ++j )
			mtx[i][j] = mi.mtx[i][j];
		vct[i] = mi.vct[i];
	}

	int nz = 2 * mi.npath;
	z.resize( nz );
	zValid = (0!=freePar( nconstraint, 2 * mi.npath, mi.x, mtx, vct, &nz, z, milerr ));
	z.resize( nz );

	return zValid;
}

bool mil_workspace::updateZ( fq::vector<double>& newz )
{
	try {
		mdlinf& mi = model.model.mi;

		z = newz;
		mil_ztox( z.size(), mi.x.size(), mtx, vct, z, mi.x, dx_z, mi.xlimit[0], mi.xlimit[1] );
		mil_xtoq( mi.nstate, mi.npath, mi.path, model.idrug, C, model.ivolt, C, mi.x, q, dq_x );
		qtometaq( mi.nstate, mi.npath, mi.path, model.nmetastate, model.nmetapath, model.metapath,
				  model.metastate, model.subpath, q, qq, dqq_q );
		qtoqe( model.nmetastate, model.nmetapath, model.metapath, model.metaclass, model.nmetaclass, td * 1.0e-3,
			   qq, qqe, dqqe_qqAsPtr );
		//calcStartProbs();
		//calcMultiProbs();

		int i, j, k, n;
		for (k=0; k<model.nmetaclass; k++) {
			n = model.nmetagroup[k]; 

			for (i=0; i<n; i++) {
				for (j=0; j<n; j++) {
					a[i][j]=qqe[model.metaindex[k][i]][model.metaindex[k][j]];
					
					if ( _isnan(a[i][j]) || a[i][j] == DBL_MAX || a[i][j] > 10e20) {
						throw a[i][j];
					}
				}
			}

			if ( qspctrm(n,a,v,c,c1,posieigencount,milerr) != 1 ) {
				return false;
			}
			for (i=0; i<n; i++) {
				w[k][i]=v[i];
				for (j=0; j<n; j++) {
					p[k][i][j] =c[i][j];
					p1[k][i][j]=c1[i][j];
				}
			}
		}
	}
	catch (double derr) {
		//unsigned int fstat = _clearfp(); // just in case
		//_fpreset();
        milerr << " Invalid A matrix (" << derr << ").  Try changing the dead time." << endl;
		return false;
	}
	catch (...) {
#ifdef _WIN32
		unsigned int fstat = _clearfp(); // just in case
		_fpreset();
		
		if ( fstat & _SW_OVERFLOW )
			milerr << "Numerical Overflow while applying rates" << endl;
		else if ( fstat & _SW_UNDERFLOW )
			milerr << "Numerical Underflow while applying rates" << endl;
		else if ( fstat & _SW_ZERODIVIDE )
			milerr << "Division by Zero while applying rates" << endl;
		else
#endif
			milerr << "Exception in MIL Core while applying rates" << endl;

		return false;
	}


	return true;
}

void mil_workspace::updateModel()
{
	model.model.updateRates();
}

void mil_workspace::updateSD( matrix<double>& hessian )
{ // after updateZ(); writes to node
	mdlinf& mi = model.model.mi;

	int i, j;
	double da, db;
	
	// we already have x (from updateZ()), and mtx == dx_z
	// mil_ztox( z.size(), _x.size(), _Mtx, _vct,
	// 	  z, _x, dx_z, _xmin, _xmax );

	matrix<double>& dx_z = mtx;
	fq::vector<double>& x = mi.x;
	fq::vector<double> sd( x.size(), 0.0 );

	int npath = mi.path.nr;
	int nz = z.size();

	for ( i=0; i<npath; i++ ) {
		for ( da=db=0., j=0; j<nz; j++ ) {
			da += DSQR( dx_z[2*i  ][j] ) * hessian[j][j];
			db += DSQR( dx_z[2*i+1][j] ) * hessian[j][j];
		}
		if ( da < 0.0 ) da = 0.0; // this line and next 9/25/01 Chris
		if ( db < 0.0 ) db = 0.0; // to prevent sqrt( - whatever )
		da = sqrt( da ) * exp( x[2*i] );
		db = sqrt( db );
		sd[2*i  ] = da;
		sd[2*i+1] = db;
	}
	model.model.updateRates( sd );
}

int CopyIdlWithTDead( QUB_Tree seg, int *idwt, float *tdwt, double tdead, bool continuing=false )
// tdead in millisec; output in sec
// continuing (as in joining segments): idwt[-1] and tdwt[-1] are valid
// returns # of dwells after applying tdead
{
	float td = float(tdead * 1.0e-3);
	
	int ndwells = seg["DwellCount"].dataAsInt();
	int *idwt_in = (int *) seg["Classes"].data();
	float *tdwt_in = (float *) seg["Durations"].data();

	int i; // outgoing index
	int j; // incoming index
	int cls;
	float tm;

	// first, drop initial too-short events
	for ( j=0; j<ndwells; j++ ) {
		tm = tdwt_in[j] * 1.0e-3f;
		if ( (float_distance(tm, td) < TDEAD_DELTA_FLOATS) || (tm > td) )
			break;
	}
	// then add the time-values of other too-shorts to the prev. dwell
	for ( i=0; j<ndwells; j++ ) {
		cls = idwt_in[j];
		tm = tdwt_in[j] * 1.0e-3f;

		if ( (continuing || (i > 0)) && (idwt[i-1] == cls) )
			tdwt[i-1] = (float)(tdwt[i-1] + tm);
		else if ( (continuing || (i > 0)) && ((float_distance(tm, td) >= TDEAD_DELTA_FLOATS) && (tm < td)) )
			tdwt[i-1] = (float)(tdwt[i-1] + tm);
		else {
			idwt[i] = cls;
			tdwt[i] = (float) (tm - td);
			i++;
		}
	}
	ndwells = i;

	return ndwells;
}

void mil_workspace::copyMetaData( QUB_Tree orig )
{
	QUB_TreeIter metaout = data.end();
	
	for ( QUB_TreeMonoIter metai = orig.children(); ! (*metai).isNull(); ++metai ) {
		QUB_Tree meta = *metai;
		string name = meta.name();
		if ( name != "Segment" ) {
			metaout.insert( meta.clone() );
			++metaout;
		}
	}
}

void mil_workspace::copyData( QUB_Tree orig )
{
	int ndwell, *idwell;
	float *tdwell;

	maxndwt = 0;

	// data = orig.clone();

	copyMetaData( orig );

	QUB_TreeIter outsegs = data.end();
	int segix = 0;

	for ( QUB_TreeMonoIter segi = orig.find("Segment");
	      ! (*segi).isNull(); segi.nextSameName() ) {

		QUB_Tree seg( (*segi).clone() );
		seg["Index"].setData( QTR_TYPE_INT, segix++ );
		segments.push_back( seg );
		outsegs.insert( seg );
		++outsegs;

		QUB_Tree durs( seg["Durations"] );
		// durs.enforceFloat( durs.dataCount() );

		idwell = (int *) seg["Classes"].data();
		tdwell = (float *) durs.data();
		ndwell = CopyIdlWithTDead( seg, idwell, tdwell, td, false );

		seg["DwellCount"].setData( QTR_TYPE_INT, ndwell );
		// *( (int *) seg["DwellCount"].data() ) = ndwell;

		ndwt.push_back( ndwell );
		idwt.push_back( idwell );
		tdwt.push_back( tdwell );
		if ( ndwell > maxndwt )
			maxndwt = ndwell;
	}
}

void mil_workspace::copyDataJoined( QUB_Tree orig )
{
	QUB_TreeMonoIter segi;

	copyMetaData( orig );
	
	for ( segi = orig.find("Segment"), maxndwt = 0;
	      ! (*segi).isNull(); segi.nextSameName() )
		maxndwt += (*segi)["Classes"].dataCount();

	QUB_Tree seg = data["Segment"];
	seg["Index"].setData( QTR_TYPE_INT, 0 );
	segments.push_back( seg );

	seg["DwellCount"].setData( QTR_TYPE_INT, maxndwt );
	seg["Classes"].setNumData( QTR_TYPE_INT, maxndwt, 1 );
	seg["Durations"].setNumData( QTR_TYPE_FLOAT, maxndwt, 1 );

	int *idwell = (int *) seg["Classes"].data();
	float *tdwell = (float *) seg["Durations"].data();

	int i;

	int bounds[2];

	for ( segi = orig.find("Segment"), maxndwt = 0, i = 0;
	      ! (*segi).isNull(); segi.nextSameName(), ++i ) {
		QUB_Tree seg_in( *segi );
		QUB_Tree durs( seg_in["Durations"] );
		// durs.enforceFloat( durs.dataCount() );

		maxndwt += CopyIdlWithTDead( seg_in, idwell+maxndwt, tdwell+maxndwt, td, (i > 0) );

		if ( i == 0 ) {
			bounds[0] = seg_in.dataAs( 0, (int) 0 );
			seg["start"].setData( QTR_TYPE_DOUBLE, seg_in["start"].dataAsDouble() );
		}
		bounds[1] = seg_in.dataAs( 1, (int) 0 );
	}
	seg.setNumData( QTR_TYPE_INT, 2, 1, bounds );

	*( (int *) seg["DwellCount"].data() ) = maxndwt;

	ndwt.push_back( maxndwt );
	idwt.push_back( idwell );
	tdwt.push_back( tdwell );
}

void mil_workspace::setupArrays()
{
	mdlinf& mi = model.model.mi;
	zValid = resetZ( false );

	int nx = mi.x.size();
	int nz = z.size();

	q.resize( mi.nstate, mi.nstate );
	qq.resize( model.nmetastate, model.nmetastate );
	qqe.resize( model.nmetastate, model.nmetastate );
	dx_z.resize( nx, nz );
	dq_x.resize( mi.npath, nx );
	dqq_q.resize( model.nmetapath, mi.npath );
	dqqe_qq.resize( model.nmetastate, model.nmetastate, model.nmetapath );
	df_qqe.resize( model.nmetastate, model.nmetastate);
	df_qq.resize( model.nmetapath );
	df_q.resize( mi.npath );
	df_x.resize( nx );
	df_z.resize( nz );
	w.resize( model.nmetaclass, model.nmetastate );
	p.resize( model.nmetaclass, model.nmetastate, model.nmetastate );
	p1.resize( model.nmetaclass, model.nmetastate, model.nmetastate );
	a.resize( model.nmetastate + 1, model.nmetastate + 1 );
	c.resize( model.nmetastate + 1, model.nmetastate + 1 );
	c1.resize( model.nmetastate + 1, model.nmetastate + 1 );
	v.resize( model.nmetastate + 1 );

	alpha.resize( model.nmetastate, 2 * maxndwt + 1 );
	beta.resize( model.nmetastate, 2 * maxndwt + 1 );
	scale.resize( 2 * maxndwt + 1 );

	dqqe_qqAsPtr.resize( model.nmetastate, model.nmetastate );
	for ( int i=0; i<model.nmetastate; ++i )
		for ( int j=0; j<model.nmetastate; ++j )
			dqqe_qqAsPtr[i][j] = &(dqqe_qq[i][j][0]);

	//setupDQDQ();
}

void mil_workspace::setTimeConstHint( QUB_Tree hist, QUB_Tree tc )
{
	int ntc = tc["Tau"].dataCount();
	double *tau = (double *) tc["Tau"].data();
	double *amp = (double *) tc["Amp"].data();

	ostringstream ost;
	ost << setw(12) << "Tau" << setw(12) << "Amp";
	for ( int i=0; i<ntc; ++i )
		ost << "\r\n" << setw(12) << setprecision(4) << tau[i]
					  << setw(12) << setprecision(4) << amp[i];

	hist["Hint"].setData( ost.str() );
}

void mil_workspace::calHist( int& nbin, float *bins, float *hist, int cls, int seg )
{
	int nseg, segoff;
	if ( seg < 0 ) {
		nseg = (int) ndwt.size();
		segoff = 0;
	}
	else if ( seg < int(ndwt.size()) ) {
		nseg = 1;
		segoff = seg;
	}
	else {
		nseg = 0;
		segoff = 0;
	}

	double sampling = data["sampling"].dataAsDouble() * 1.0e-3;

	binning( 0.0, nseg, &(ndwt[segoff]), &(tdwt[segoff]), &nbin, bins, (float)(td * 1.0e-3), histBinR );
	if ( smoothHistBins )
		calhist_smooth( cls, nseg, &(ndwt[segoff]), &(idwt[segoff]), &(tdwt[segoff]), nbin, bins, hist, (float)(td * 1.0e-3), (float)sampling, (float)histBinR );
	else
		calhist( 0.0, cls, nseg, &(ndwt[segoff]), &(idwt[segoff]), &(tdwt[segoff]), nbin, bins, hist, (float)(td * 1.0e-3) );

	int ndwell_tot = 0; // prepare to normalize histogram by number of events in this class
	for ( int i=0; i<nseg; ++i ) {
		int *idwell = idwt[segoff+i];
		for ( int j=ndwt[segoff+i]-1; j>=0; --j ) {
			if ( idwell[j] == cls )
				++ndwell_tot;
		}
	}

	if ( sqrtHist ) {
		for ( int i=0; i<nbin; ++i )
		  hist[i] = ndwell_tot ? (float) (sqrt( hist[i] / (float) ndwell_tot )) : 0.0f;
	}
	else {
		for ( int i=0; i<nbin; i++)
		  hist[i] = ndwell_tot ? (float) (hist[i] / (float) ndwell_tot) : 0.0f;
	}
}

bool mil_workspace::calPDFandTC( int nbin, float *bins, float **pdf, double *tau, double *amp, double& meanTau, int cls )
{
	fq::vector<double> pe( model.nmetastate ),
					 c(  model.nmetastate ),
					 w(  model.nmetastate );
	int ncomp = componentsInClass( cls );

	peqm(qqe,cls,pe);  
	if ( pdfcof(qqe,cls,pe,c,w) ) {
		pdfbin(cls,w,c,td,nbin,bins,pdf);
	}
	else {
		for ( int j=0; j<ncomp; ++j )
			for ( int k=0; k<nbin; ++k )
				pdf[j][k] = 0.0;
	}
   
	meanTau = 0.0;
	for ( int i=0; i<ncomp; ++i ) {
		tau[i] = (w[i] == 0.0) ? 0.0 : (-1000.0 / w[i]);
		amp[i] = (w[i] == 0.0) ? 0.0 : (-c[i] / w[i]);
		meanTau += tau[i] * amp[i];
	}

	return true;
}

void mil_workspace::peqm(double **q, int ic, double *pe){
	int i,k,l,m,n; 
	
	matrix<double> qmn(model.nmetastate+1, model.nmetastate+1);
	matrix<double> qnn(model.nmetastate+1, model.nmetastate+1);
	matrix<double> qnm(model.nmetastate+1, model.nmetastate+1);
	matrix<double> qmm(model.nmetastate+1, model.nmetastate+1);
	matrix<double> a(model.nmetastate+1, model.nmetastate+1);
	fq::vector<double> wr(model.nmetastate+1);
	fq::vector<double> wi(model.nmetastate+1);
	
	matrix<double> amm(model.nmetastate+1, model.nmetastate+1); 
	matrix<double> ann(model.nmetastate+1, model.nmetastate+1); 
	
	m = model.nmetagroup[ic];
	n = model.nmetastate - model.nmetagroup[ic];
	
	for (k=0; k<m; k++) 
		for (l=0; l<m; l++)
			qmm[k][l]= q[model.metaindex[ic][k]]
						[model.metaindex[ic][l]];
       
	for (k=0; k<m; k++) 
		for (l=0; l<n; l++)
			qmn[k][l]= q[model.metaindex[ic][k]]
						[model.metaindex[ic][model.nmetastate+l]];

	for (k=0; k<n; k++) 
		for (l=0; l<n; l++)
			qnn[k][l]= q[model.metaindex[ic][model.nmetastate+k]]
						[model.metaindex[ic][model.nmetastate+l]];
	
	for (k=0; k<n; k++) 
		for (l=0; l<m; l++)
			qnm[k][l]=q[model.metaindex[ic][model.nmetastate+k]][model.metaindex[ic][l]];
	   
	// Should A be initialized to I before each gaussj0?

	gaussj0(qmm,m,a,1,milerr);
		
	gaussj0(qnn,n,a,1,milerr);

	mxm(m,m,n,qmm,qmn,a);
	mxm(m,n,n,a,qnn,qmn);

	mxm(m,n,m,qmn,qnm,a);

	// qmm= transpose(a)
	for( k=0; k<m; k++ )
		for( l=0; l<m; l++ ) 
			qmm[k][l]=a[l][k];

	eigen(qmm,m,wr+1,wi+1,a);

	double dev = 1e8, tmp;
	k=-1;
	for( i=1; i<=m; i++) {
		tmp = fabs(wr[i]-1);
		if (tmp<dev) {
			k = i-1;
			dev = tmp;
		}
		// changed by Chris 10-16-5
		// most upside-down pdfs go away when you pick the eigenvector whose eigenvalue is _closest_ to 1
		//if (fabs(wr[i]-1)<1.0e-1 && fabs(wi[i])<1.0e-2) {
		// k=i-1; 
		//}
	}

	// ASSERT(k>=0), otherwise, Err in finding pe 
	// if (k < 0) milerr << "peqm: error in finding pe" << endl;
	tmp=0.0;
	if( k >= 0 )
		for( i=0; i<m; i++)
			tmp += a[i][k];

	if ( tmp != 0.0 )
		for (i=0; i<m; i++)
			pe[i] = a[i][k]/tmp;
	}


bool mil_workspace::pdfcof(double **q, int ic, double *pe, double *c, double *wr){
	int     i,j,k,l,m,n; 
	double  tmp;
	
	matrix<double> qmm(model.nmetastate+1, model.nmetastate+1); 
	matrix<double> qmn(model.nmetastate+1, model.nmetastate+1);
	matrix<double> a(model.nmetastate+1, model.nmetastate+1); 
	matrix<double> s(model.nmetastate+1, model.nmetastate+1);
	matrix<double> s1(model.nmetastate+1, model.nmetastate+1);

	m = model.nmetagroup[ic];
	n = model.nmetastate - model.nmetagroup[ic];
	
	for (k=0; k<m; k++) 
		for (l=0; l<m; l++) {
			i = model.metaindex[ic][k];
			j = model.metaindex[ic][l];
			qmm[k][l] = q[i][j];
			if ( _isnan(a[i][j]) || a[i][j] == DBL_MAX || qmm[k][l] > 10e20 )
				return false;
			}
		
	if ( ! qspctrm(m,qmm,wr,s,s1,posieigencount,milerr) )
		return false;
	
	for (k=0; k<m; k++) 
		for (l=0; l<n; l++) {
			i = model.metaindex[ic][k];
			j = model.metaindex[ic][model.nmetastate+l];
			qmn[k][l] = q[i][j];
			} 

	for (k=0; k<m; k++) {
		for (i=0; i<m; i++) 
			for (j=0; j<m; j++)
				qmm[i][j]=s[i][k]*s1[k][j];
		mxm(m,m,n,qmm,qmn,a);
		for (tmp=0.0,i=0; i<m; i++) 
			for (j=0; j<n; j++) 
				tmp += pe[i]*a[i][j];
		c[k]=tmp;
		}

	return true;
	}  


void mil_workspace::pdfbin(int ic, double *w, double *c, double td, int nbin, float* bin, float **pdf)
   // td in ms, bin[] in sec
{
	int     i,j;
	float   a,b,eps=1.0e-5f;
	
	//a=float(td);
	a = float(bin[0] * 1.0e3 / histBinR);
	for (i=0; i<nbin; i++) 
   {
	   b = bin[i] * float(1.0e3);
	   //pdf[i] = 0.0;
	   for (j=0; j<model.nmetagroup[ic]; j++) 
      {
		  if (w[j] > 10)
			 pdf[j][i] = 0.0;
	      else if (fabs(w[j])<eps) 
	         pdf[j][i]=float(c[j]*(b-a)*1.0e-3);
	      else 
	         pdf[j][i]=float(c[j]/w[j]*(exp(w[j]*(b-td)*1.0e-3)-exp(w[j]*(a-td)*1.0e-3)));

		  if ( sqrtHist )
			 pdf[j][i] = (float) sqrt( pdf[j][i] );
	   }
	   //if ( pdf[i] < 0.0 ) {
		//   milerr << "skipping class " << (ic+1) << " pdf: found negative probability" << endl;
		//   throw 1;
	   //}
	   a = b;
	}
}   



void mil_workspace::setupDQDQ()          // path            -> dQdq
{
	int i, j, qq;
	mdlinf& mi = model.model.mi;
	dQdq.resize( mi.npath, mi.nstate, mi.nstate );

	for ( qq=0; qq<mi.npath; ++qq ) {
		for ( i=0; i<mi.nstate; ++i ) {
			for ( j=0; j<mi.nstate; ++j ) {
				dQdq[qq][i][j] = 0;
			}
		}
		dQdq[ qq ][ mi.path[qq][0] ][ mi.path[qq][0] ] = -1;
		dQdq[ qq ][ mi.path[qq][0] ][ mi.path[qq][1] ] = 1;
	}
}

void mil_workspace::calcStartProbs()     // q, dQdq         -> pr, dpr_q
{
	mdlinf& mi = model.model.mi;
	int i, j, k, qq;

	dpr_q.resize(mi.npath, mi.nstate);

	if ( ! usePeq ) {
		for ( qq=0; qq<mi.npath; ++qq )
			for ( i=0; i<mi.nstate; ++i )
				dpr_q[qq][i] = 0.0;
		return;
	}
			
	matrix<double> R(mi.nstate, mi.nstate+1);
	matrix<double> RT(mi.nstate+1, mi.nstate);
	matrix<double> RRT(mi.nstate, mi.nstate);
	matrix<double> dRRTdq(mi.nstate, mi.nstate);
	matrix<double> Mat1(mi.nstate, mi.nstate);
	fq::vector<double> U(mi.nstate);

	//Peq = 1.(R.RT)-1

    //U
	for ( i=0; i<mi.nstate; ++i )
		U[i] = 1;

    //R, RT
	for ( i=0; i<mi.nstate; ++i ) {
		for ( j=0; j<mi.nstate; ++j ) {
			R[i][j] = q[i][j];
			RT[j][i] = q[i][j];
		}
	}
	for ( k=0; k<mi.nstate; ++k ) {
		R[k][mi.nstate] = 1;
		RT[mi.nstate][k] = 1;
	}

    //R.RT
	for ( i=0; i<mi.nstate; ++i ) {
		for ( j=0; j<mi.nstate; ++j ) {
			RRT[i][j] = 0.0;
			for ( k=0; k<=mi.nstate; ++k ) {
				RRT[i][j] += R[i][k] * RT[k][j];
			}
		}
	}

    //RRTI = (R.RT)-1
	gaussj_invert(RRT, mi.nstate, milerr);
	matrix<double>& RRTI = RRT; // since that's what it's become

    //MulVecMatVec(U, Mat2, Peq, StateCount);
	vxm(mi.nstate, mi.nstate, RRTI, U, mi.pr);

	// dpr_q:

	for ( i=0; i<mi.nstate; ++i )
		U[i] = -1;

	for ( qq=0; qq<mi.npath; ++qq ) {
		//d(R.RT)/dq = dRdq.RT + R.dRTdq
		for ( i=0; i<mi.nstate; ++i ) {
			for ( j=0; j<mi.nstate; ++j ) {
				dRRTdq[i][j] = 0.0;
				for ( k=0; k<=mi.nstate; ++k ) {
					//dRdq.RT
					if ( k < mi.nstate )
						dRRTdq[i][j] += dQdq[qq][i][k] * RT[k][j];
					//R.dRTdq
					if ( k < mi.nstate )
						dRRTdq[i][j] += R[i][k] * dQdq[qq][j][k];
				}
			}
		}

		//Mat1 = RRTI.dRRTdq
		// MulMatMat(RRTI, dRRTdq, Mat1, StateCount);
		mxm(mi.nstate, mi.nstate, mi.nstate, RRTI, dRRTdq, Mat1);

		//dRRTdq = Mat1.RRTI
		//MulMatMat(Mat1, RRTI, dRRTdq, StateCount);
		mxm(mi.nstate, mi.nstate, mi.nstate, Mat1, RRTI, dRRTdq);

		//MulVecMatVec(U, dRRTdq, TPMDouble(dPeqdq[0])[qq], StateCount);
		vxm(mi.nstate, mi.nstate, dRRTdq, U, dpr_q[qq]);
	}
}

void mil_workspace::calcMultiProbs()     // pr, metastate   -> mpr, dmpr_pr
{
   mdlinf& mi = model.model.mi;
   int i, j;
   double tmp, dum;

   dmpr_pr.resize(model.nmetastate, mi.nstate);
   dmpr_q.resize(model.nmetastate, mi.npath);

   for (i=0; i<model.nmetastate; i++)
   {
      tmp = lnfactrl(model.nchannel);
      for (j=0; j<mi.nstate; j++)
      {
         dum = (mi.pr[j] < 1.0e-5) ? 1.0e-10 : mi.pr[j];
         tmp += model.metastate[i][j]*log(dum) - lnfactrl(model.metastate[i][j]);
      }
      model.mprx[i] = exp(tmp);
   }

   for ( i=0; i<model.nmetastate; ++i ) {
	   for ( j=0; j<mi.nstate; ++j ) {
		   dmpr_pr[i][j] = model.mprx[i] * model.metastate[i][j] / mi.pr[j];
	   }
   }

   mxm(model.nmetastate, mi.nstate, mi.npath, dmpr_pr, dpr_q, dmpr_q);
}

void mil_workspace::calc_dlogLL_q( int cls, fq::vector<double>& dLL_q )
{
	mdlinf& mi = model.model.mi;
	int i, qq;
	fq::vector<double> beta0(model.nmetastate);

	dLL_q.resize(mi.npath);

	// extract beta[0] as full nmetastate vector
	for ( i=0; i<model.nmetagroup[cls]; ++i )
		beta0[ model.metaindex[cls][i] ] = beta[i][0];

	vxm(model.nmetastate, mi.npath, dmpr_q, beta0, dLL_q);
	for ( qq=0; qq<mi.npath; ++qq )
		dLL_q[qq] *= -1;
	//for ( qq=0; qq<mi.npath; ++qq )
	//	dLL_q[qq] = - (beta0 * dmpr_q[qq]);
}


// *************************************************************************

mil_eval_tree::mil_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
							  QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB )
  : model( config["ModelFile"], config["SearchLimit"].dataAsDouble(1000.0) ),
    cfg( config ), reportCB( rptCB ), resultCB( rsltCB ), percentCB( pctCB ),
	milerr( rptCB )
{
	QUB_Tree stopFlag = * (config.find("StopFlag"));
	if ( stopFlag.isNull() )
		Stopped = 0;
	else
		Stopped = (int *) stopFlag.data();

	bool smoothHistBins = config["Histograms"]["SmoothBins"].dataAsInt();
	bool sqrtHist = config["Histograms"]["Sqrt"].dataAsInt();

	bool usePeq = config["StartProb"].dataAsString() == "equilibrium";
	for ( std::vector<QUB_Tree>::iterator dsi = dataSets.begin(); dsi != dataSets.end(); ++dsi )
		workspaces.push_back( mil_worksptr( new mil_workspace( *dsi, model, usePeq, smoothHistBins, sqrtHist, Stopped, milerr ) ) );

	opt = config["Mode"].dataAsString() == "optimize";
	together = config["use segments"].dataAsString() == "together";
	pdfsEachIter = (0!=config["Histograms"]["PDFs Each Iteration"].dataAsInt( 0 ));
	
	QUB_Tree dfp = cfg["DFP"];
	maxIter = dfp["MaxIterations"].dataAsInt( 100 );
	maxRestarts = dfp["MaxRestarts"].dataAsInt( 0 );
	maxStep = dfp["MaxStep"].dataAsDouble( 1.0 );
	convLL = dfp["ConvLL"].dataAsDouble( 0.0001 );
	convGrad = dfp["ConvGrad"].dataAsDouble( 0.00005 );

	pctNode = QUB_Tree::Create("");
	pctNode.setData( QTR_TYPE_INT, 0 );
	pctPtr = (int *) pctNode.data();

	nSeg = 0;
	for ( int i=0; i<int(workspaces.size()); ++i )
		nSeg += workspaces[i]->segCount();
}

void mil_eval_tree::initResults()
{
	result = QUB_Tree::Create("MIL Result");
	result.appendClone( cfg );

	QUB_Tree bincount = cfg["Histograms"]["BinCount"];

	if ( together ) {
		result["LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
		result["Gradient"].setData( QTR_TYPE_DOUBLE, 0.0 );
		result["Iterations"].setData( QTR_TYPE_INT, 0 );
	}

	for ( int i=0; i<int(workspaces.size()); ++i ) {
		mil_workspace& wspace = * workspaces[i];
		result.appendChild( wspace.data );

		string filename = wspace.data["FileName"].dataAsString();
		size_t pos; // of final '\'
		if ( (pos = filename.rfind("\\")) != string::npos )
			filename = filename.substr( pos+1 );

		if ( together ) {
			QUB_Tree histGroup = wspace.data["HistogramGroup"];
			histGroup["Title"].setData( filename );

			histGroup.appendClone( bincount );
			wspace.getHistograms( histGroup );
		}
		for ( int j=0; j<wspace.segCount(); ++j ) {
			QUB_Tree seg = wspace.segments[j];
			seg["LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
			seg["Gradient"].setData( QTR_TYPE_DOUBLE, 0.0 );
			seg["Iterations"].setData( QTR_TYPE_INT, 0 );

			if ( ! together ) {
				QUB_Tree histGroup = seg["HistogramGroup"];
				ostringstream ost;
				ost << filename << " Segment " << (j+1);
				histGroup["Title"].setData( ost.str() );

				histGroup.appendClone( bincount );
				wspace.getHistograms( histGroup, j );
			}
		}
	}
}

void mil_eval_tree::addTCSelectVars( QUB_Tree seg )
{
	QUB_TreeIter seginsert = seg.end();

	int tci = 1;
	QUB_TreeMonoIter tcni = seg["HistogramGroup"].find("TimeConstants");
	for ( ; ! (*tcni).isNull(); ++tci, tcni.nextSameName() ) {
		ostringstream tauName;
		tauName << "Tau" << tci;
		QUB_Tree tau = (*tcni)["Tau"].clone();
		tau.setName( tauName.str() );
		seginsert.insert( tau );
		++seginsert;

		ostringstream ampName;
		ampName << "TauAmp" << tci;
		QUB_Tree amp = (*tcni)["Amp"].clone();
		amp.setName( ampName.str() );
		seginsert.insert( amp );
		++seginsert;
	}
}

QUB_Tree MakeRateSelectVar( int from, int to, double rate, bool isK1=false )
{
	ostringstream ost;
	ost << (from+1) << "->" << (to+1);
	if ( isK1 )
		ost << "b";
	
	QUB_Tree rateNode = QUB_Tree::Create( ost.str() );

	rateNode.setData( QTR_TYPE_DOUBLE, rate );
	return rateNode;
}

void mil_eval_tree::addRateSelectVars( QUB_Tree seg )
{
	/*
	QUB_TreeIter seginsert = seg.end();

	QUB_TreeMonoIter ri = seg["ModelFile"]["Rates"].find("Rate");
	for ( ; ! (*ri).isNull(); ri.next("Rate") ) {
		QUB_Tree rate = *ri;
		int from = ( (int *) rate["States"].data() )[0];
		int to   = ( (int *) rate["States"].data() )[1];
		double k0a = ( (double *) rate["k0"].data() )[0];
		double k1a = ( (double *) rate["k1"].data() )[0];
		double k0b = ( (double *) rate["k0"].data() )[1];
		double k1b = ( (double *) rate["k1"].data() )[1];

		seginsert.insert( MakeRateSelectVar( from, to, k0a ) );
		++seginsert;
		seginsert.insert( MakeRateSelectVar( to, from, k0b ) );
		++seginsert;
		if ( abs(k1a) > 1.0e-3 ) {
			seginsert.insert( MakeRateSelectVar( from, to, k1a, true ) );
			++seginsert;
		}
		if ( abs(k1b) > 1.0e-3 ) {
			seginsert.insert( MakeRateSelectVar( to, from, k1b, true ) );
			++seginsert;
		}
	}
	*/
}

void mil_eval_tree::finishResults()
{
	double sumLL = 0.0;
	int sumNevent = 0;

	for ( QUB_TreeMonoIter dsi = result.find("DataSet"); ! (*dsi).isNull(); dsi.nextSameName() ) {
		for ( QUB_TreeMonoIter segi = (*dsi).find("Segment"); ! (*segi).isNull(); segi.nextSameName() ) {
			QUB_Tree seg = *segi;
			seg.find("Classes").remove();
			seg.find("Durations").remove();
			seg.find("Firsts").remove();
			seg.find("Lasts").remove();
			seg.find("TempLL").remove();
			seg.find("TempGradient").remove();
			QUB_Tree ll = *(seg.find("LL"));
			if ( ! ll.isNull() ) {
				sumLL += ll.dataAsDouble();
				int dwellCount = seg["DwellCount"].dataAsInt();
				sumNevent += dwellCount;
				seg["LL per event"].setData( QTR_TYPE_DOUBLE, ll.dataAsDouble() / (dwellCount ? dwellCount : 1) );
			}

			if ( ! together ) {
				try {
					addTCSelectVars( seg );
					addRateSelectVars( seg );
				}
				catch (...) {
#ifdef _WIN32
					_clearfp();
					_fpreset();
#endif
				}
			}
		}
	}
	result.find("TempLL").remove();
	result.find("TempLL per event").remove();
	result.find("TempGradient").remove();

	if ( sumNevent > 0 )
		result["LL per event"].setData( QTR_TYPE_DOUBLE, sumLL / sumNevent );

	CleanResultFloats( result );
}

void mil_eval_tree::updatePDFs()
{
	if ( together ) {
		for ( int i=0; i<int(workspaces.size()); ++i ) {
			QUB_Tree histGroup = workspaces[i]->data["HistogramGroup"];
			workspaces[i]->getPDFs( histGroup, histGroup );
		}
	}
	else {
		QUB_Tree histGroup = currSpace->segments[currSeg]["HistogramGroup"];
		currSpace->getPDFs( histGroup, histGroup );
	}
}

void mil_eval_tree::sendResult( QUB_Tree resNode, int iter, bool final, int errCode )
{
	CleanResultFloats(resNode);
	if ( final ) {
		resNode["Final"];
		resNode["ErrorCode"].setData( QTR_TYPE_INT, errCode );
	}
	else {
		resNode["Iterations"].setData( QTR_TYPE_INT, iter );
	}
	QTR_DoCallback(resultCB, resNode);
}

bool mil_eval_tree::canRun()
{
	bool can_run = true;

	std::vector<double> concs;
	int maxnchannel = 1;
	
	// all of the following must be true:
	//   [# cls in data] <= [# cls in metamodel]
	for ( int i=0; i<int(workspaces.size()); ++i ) {
		// some groundwork for later checks
		concs.push_back( workspaces[i]->data["ExpCond"]["Ligand"].dataAsDouble( 1.0 ) );
		
		int nchannel = (int) workspaces[i]->data["ExpCond"].find("ChannelCount")->dataAsDouble( model.nchannel );
		if ( nchannel > maxnchannel )
			maxnchannel = nchannel;

		if ( workspaces[i]->classCountInData() > workspaces[i]->classCount() ) {
			milerr << "Not enough classes in the model (file " << (i+1) << ")" << endl;
			can_run = false;
		}
	}

	//   if there's more than one channel, there can't be substates
	/*
	set<int> classesPresent; // in single-channel model
	for ( int abc=0; abc<model.mi.clazz.size(); ++abc )
		classesPresent.insert( model.mi.clazz[abc] );

	if ( (maxnchannel > 1) && (classesPresent.size() > 2) ) {
		milerr << "Can't optimize multiple channels with substates" << endl;
		can_run = false;
	}
	*/

	// make sure loops in detailed balance are stimulus-friendly
	QUB_Tree loopMsg = DoCheckLoopStimuli(model.node);
	if ( loopMsg.dataCount() ) {
	  milerr << loopMsg.dataAsString() << endl;
	  can_run = false;
	}

	return can_run;
}

bool mil_eval_tree::evaluateAll( fq::vector<double> &z, fq::vector<double> &gz,
								double &ll, int &nf, int &ndf, bool flagGrad )
{
	fq::vector<double> seggz( gz.size() );
	double segll;

	int i;
	// int nerr = 0, nseg = 0;

	for ( i=0; i<gz.size(); ++i )
		gz[i] = 0.0;
	ll = 0.0;

	for ( i=0; i<int(workspaces.size()); ++i ) {
		mil_workspace& wspace = * workspaces[i];

		if ( ! wspace.updateZ( z ) ) {
			return false;
			//nseg += wspace.segCount();
			//nerr += wspace.segCount();
			//continue;
		}

		for ( int j=0; j<wspace.segCount(); ++j/*, ++nseg*/ ) {
			/*nerr +=*/ if ( wspace.evaluate( seggz, segll, flagGrad, j ) != 0 )
				return false;

			UpdateResult( wspace.segments[j], segll, seggz );

			ll += segll;
			for ( int k=0; k<gz.size(); ++k )
				gz[k] += seggz[k];

			if ( Stopped && *Stopped )
				return false;
		}

		if ( Stopped && *Stopped )
			return false;
	}

	UpdateResult( result, ll, gz ); // only on check?

	++nf;
	if ( flagGrad )
		++ndf;

	return true; // ( nerr < nseg );
}

// i tried allowing some segs to have errors, but it reached a worse LL in more iterations, probably slower too

dfp_result mil_eval_tree::runDFP()
{
	iterations = 0;
	++maxRestarts;
	bool doRestart = true;

	dfp_result res;

	while ( doRestart && maxRestarts ) {
		doRestart = false;
		--maxRestarts;

		res = dfp_optimize( this, maxIter, convLL, convGrad, maxStep );

		switch ( res.err ) {
		case 0:
		case -1: // killed in evaluate
		case -2: // killed in check
			break;
		case -3: // max iterations
			milerr << "Exceeded Maximum Iterations" << endl;
			doRestart = true;
			break;
		case -666: // bad constraints
			milerr << "Impossible to satisfy all constraints" << endl;
			break;
		default:
			milerr << "Unknown Error in MIL Core" << endl;
		}
	}

	return res;
}

QUB_Tree mil_eval_tree::execute()
{
	result = QUB_Tree::Create("");
	string mode = cfg["Mode"].dataAsString();
	if ( mode == "model report" ) {
		if ( workspaces.size() && workspaces[0]->updateZ( workspaces[0]->z ) ) {
			result = workspaces[0]->getModelStats( cfg["FirstLatMax"].dataAsDouble(), cfg["FirstLatBins"].dataAsInt() );
			result.appendClone( cfg );
		}
	}
	else {
		initResults();
		int err = 123; // not even run

		if ( mode == "evaluate" ) {
			fq::vector<double>& z = workspaces[0]->z;
			fq::vector<double> gz( z.size() );
			double ll;
			int nf, ndf;

			if ( canRun() )
				evaluateAll( z, gz, ll, nf, ndf, true );
			BlessResult( result );

			if ( together ) {
				updatePDFs();
			}
			else {
				for ( int i=0; i<int(workspaces.size()); ++i ) {
					currSpace = workspaces[i];
					for ( currSeg=0; currSeg<currSpace->segCount(); ++currSeg )
						updatePDFs();
				}
			}
		}
		else if ( mode == "optimize" ) {
			if ( together ) {
				result.appendChild( model.node );
				if ( canRun() ) {
					dfp_result dfpres = runDFP();
					workspaces[0]->updateSD( dfpres.hessian );
					err = dfpres.err;
					//milerr << "Covariance:\n" << dfpres.hessian << endl;
				}
				if ( ! pdfsEachIter && ! err )
					updatePDFs();
				sendResult( result, -1, true, err );
			}
			else {
				bool can_run = canRun();
				iSeg = 0;
				for ( int i=0; i<int(workspaces.size()); ++i ) {
					currSpace = workspaces[i];
					for ( currSeg=0; currSeg<currSpace->segCount(); ++currSeg, ++iSeg ) {
						currSpace->reset();
						currSpace->segments[currSeg].appendChild( model.node );
						
						ostringstream mdlname;
						mdlname << "Final " << (currSpace->segments[currSeg]["Index"].dataAsInt() + 1);
						model.node.setData( mdlname.str() );

						if ( can_run ) {
							dfp_result dfpres = runDFP();
							currSpace->updateSD( dfpres.hessian );
							err = dfpres.err;
							//milerr << "Covariance:\n" << dfpres.hessian << endl;
						}
						if ( ! pdfsEachIter && ! err )
							updatePDFs();
						sendResult( currSpace->segments[currSeg], -1, true, err );

						if ( Stopped && *Stopped )
							break;
					}

					if ( Stopped && *Stopped )
						break;
				}
			}
		}
		else { // hist
			err = 0;
			for ( int i=0; i<int(workspaces.size()); ++i )
				err = err | (int)workspaces[i]->updateZ( workspaces[i]->z ); // generate metaq etc.

			if ( together && ! err ) {
				updatePDFs();
			}
			else if ( ! err ) {
				for ( int i=0; i<int(workspaces.size()); ++i ) {
					currSpace = workspaces[i];
					for ( currSeg=0; currSeg<currSpace->segCount(); ++currSeg ) {
						updatePDFs();
					}
				}
			}
		}

		finishResults();
	}
	result.clone().saveAs("mil_result.qtr");
	return result;
}

bool mil_eval_tree::getStartingZ( fq::vector<double> &z )
{
	z = workspaces[0]->z;
	return workspaces[0]->zValid;
}

bool mil_eval_tree::evaluate( fq::vector<double> &z, fq::vector<double> &gz,
								double &ll, int &nf, int &ndf, bool flagGrad ) // iter, pdf ??
{
	bool rtnVal = false;
	if ( together ) {
		rtnVal = evaluateAll( z, gz, ll, nf, ndf, flagGrad );
	}
	else {
		if ( currSpace->updateZ( z ) )
			rtnVal = (currSpace->evaluate( gz, ll, flagGrad, currSeg ) == 0);
		else
			rtnVal = false;

		UpdateResult( currSpace->segments[currSeg], ll, gz );
	}

	if ( rtnVal )
		rtnVal = ! (Stopped && *Stopped);

	return rtnVal;
}

bool mil_eval_tree::checkContinue( int iter, QUBOPT_VAR_NOT_USED int nf, QUBOPT_VAR_NOT_USED int ndf, QUBOPT_VAR_NOT_USED double ll,
				   QUBOPT_VAR_NOT_USED fq::vector<double> & z, QUBOPT_VAR_NOT_USED fq::vector<double> & gz )
{
	++iterations;

	workspaces[0]->updateModel();

	if ( together ) {
		if ( iter == 1 || pdfsEachIter )
			updatePDFs();
		BlessResult( result );
		sendResult( result, iterations );
	}
	else {
		if ( iter == 1 || pdfsEachIter )
			updatePDFs();
		BlessResult( currSpace->segments[currSeg] );
		sendResult( currSpace->segments[currSeg], iterations );
	}

	*pctPtr = (100 * iter) / maxIter;
	if ( ! together )
		*pctPtr = (*pctPtr + 100 * iSeg) / nSeg;
	if ( *pctPtr > 100 )
		*pctPtr = 100;

	cout << "." << flush;
	QTR_DoCallback(percentCB, pctNode);

	return ! (Stopped && *Stopped);
}

void UpdateResult( QUB_Tree result, double ll, fq::vector<double>& gz )
{
	ll *= -1;

	result["TempLL"].setData( QTR_TYPE_DOUBLE, ll );

	QUB_Tree grad = result["TempGradient"];
	if ( int(grad.dataCount()) != gz.size() )
		grad.setNumData( QTR_TYPE_DOUBLE, gz.size(), 1, &(gz[0]) );
	else
		memcpy( grad.data(), &(gz[0]), gz.size() * sizeof(double) );
}

void CleanResultFloats( QUB_Tree node )
{
	int n = node.dataCount();
	if ( node.dataType() == QTR_TYPE_FLOAT ) {
		float *data = (float *) node.data();
		for ( int i=0; i<n; ++i )
			if ( ! _finite( data[i] ) )
				data[i] = 0.0;
	}
	else if ( node.dataType() == QTR_TYPE_DOUBLE ) {
		double *data = (double *) node.data();
		for ( int i=0; i<n; ++i )
			if ( ! _finite( data[i] ) )
				data[i] = 0.0;
	}
	for (QUB_TreeMonoIter ci = node.children(); ! ci->isNull(); ++ci)
		CleanResultFloats( *ci );
}

void BlessResult( QUB_Tree result )
{
	QUB_Tree tmpGrad = *(result.find("TempGradient"));
	if ( ! tmpGrad.isNull() ) {
		double ll = result["TempLL"].dataAsDouble();
		QUB_Tree blessGrad = result["Gradient"];

		if ( blessGrad.dataCount() != tmpGrad.dataCount() ) {
			blessGrad.setNumData( QTR_TYPE_DOUBLE, tmpGrad.dataCount(), 1, (double *) tmpGrad.data() );
			
			if ( _isnan(ll) || ll == DBL_MAX || ll > 10e20)
				result["LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
			else
				result["LL"].setData( QTR_TYPE_DOUBLE, ll );
			if ( _isnan(ll) || ll == DBL_MAX || ll > 10e20)
				result["Initial LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
			else
				result["Initial LL"].setData( QTR_TYPE_DOUBLE, ll );
		}
		else {
			memcpy( blessGrad.data(), tmpGrad.data(), tmpGrad.dataCount() * sizeof(double) );
			if ( _isnan(ll) || ll == DBL_MAX || ll > 10e20)
				* ((double *) result["LL"].data() ) = 0.0;
			else
				memcpy( result["LL"].data(), &ll, sizeof(double) );
		}
		
		double *gz = (double *) blessGrad.data();
		for ( int i=0; i<int(tmpGrad.dataCount()); ++i ) {
			if ( _isnan(gz[i]) || gz[i] == DBL_MAX || gz[i] > 10e20)
				gz[i] = 0.0;
		}
	}
	for (QUB_TreeMonoIter dsi = result.find("DataSet"); ! (*dsi).isNull(); dsi.nextSameName() )
		for ( QUB_TreeMonoIter segi = (*dsi).find("Segment"); ! (*segi).isNull(); segi.nextSameName() )
			BlessResult( *segi );
}

// *************************************************************************

string ExtractFileName( string fname )
{
	size_t slash = fname.rfind("\\");
	if ( slash == string::npos )
		slash = -1;
	size_t dot = fname.rfind(".");
	if ( dot == string::npos || dot < slash )
		dot = fname.size();
	return fname.substr( slash+1, dot - slash - 1 );
}

string SpaceCat( string pre, string text )
{
	ostringstream ost;
	ost << pre << ' ' << text;
	string out = ost.str();
	return out;
}

void mil_eval_renameOptVars( QUB_Tree node, string& fname )
{
	node.find("LL")->setName( SpaceCat("LL", fname) );
	node.find("LL per event")->setName( SpaceCat("LL per event", fname) );
	node.find("Initial LL")->setName( SpaceCat("Initial LL", fname) );
	node.find("Gradient")->setName( SpaceCat("Gradient", fname) );
	node.find("Iterations")->setName( SpaceCat("Iterations", fname) );
	node.find("ErrorCode")->setName( SpaceCat("ErrorCode", fname) );
	
	QUB_Tree x = * node.find("HistogramGroup")->find("Title");
	x.setData( SpaceCat( x.dataAsString(), fname ) );

	x = * node.find("ModelFile");
	x.setData( SpaceCat( fname, x.dataAsString() ) );
}

void mil_eval_renameVars( QUB_Tree results, QUB_Tree model )
{
	string fname = ExtractFileName( model["FileName"].dataAsString() );

	QUB_Tree node;
	
	node = results["Properties"]["ModelFile"];
	node.setData( SpaceCat( fname, node.dataAsString() ) );

	mil_eval_renameOptVars( results, fname );
	for ( QUB_TreeMonoIter dsi = results.find("DataSet"); ! dsi->isNull(); dsi.nextSameName() ) {
		mil_eval_renameOptVars( *dsi, fname );
		for ( QUB_TreeMonoIter segi = dsi->find("Segment"); ! segi->isNull(); segi.nextSameName() ) {
			mil_eval_renameOptVars( *segi, fname );
		}
	}
}

void mil_eval_renameAndReparentOptVars( QUB_Tree node, string& fname, QUBOPT_VAR_NOT_USED QUB_Tree sumNode )
{
	QUB_TreeIter inserter = node.end();
	QUB_Tree var;

	var = * node.find("LL");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("LL", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("LL per event");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("LL per event", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("Initial LL");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("Initial LL", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("Gradient");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("Gradient", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("Iterations");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("Iterations", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("ErrorCode");
	if ( ! var.isNull() ) {
		var.setName( SpaceCat("ErrorCode", fname) );
		inserter.insert( var );
		++inserter;
	}
		
	var = * node.find("HistogramGroup");
	if ( ! var.isNull() ) {
		QUB_Tree x = * var.find("Title");
		x.setData( SpaceCat( x.dataAsString(), fname ) );
		inserter.insert( var );
		++inserter;
	}

	var = * node.find("ModelFile");
	if ( ! var.isNull() ) {
		var.setData( SpaceCat( fname, var.dataAsString() ) );
		inserter.insert( var );
		++inserter;
	}
}

void mil_eval_renameAndReparentVars( QUB_Tree results, QUB_Tree model, QUB_Tree sumResults )
{
	string fname = ExtractFileName( model["FileName"].dataAsString() );

	QUB_Tree node;
	
	node = results["Properties"]["ModelFile"];
	node.setData( SpaceCat( fname, node.dataAsString() ) );
	sumResults["Properties"].appendChild( node );

	mil_eval_renameAndReparentOptVars( results, fname, sumResults );
	QUB_TreeMonoIter dsi, sumdsi, segi, sumsegi;

	for ( dsi = results.find("DataSet"), sumdsi = sumResults.find("DataSet");
	      ! dsi->isNull(); dsi.nextSameName(), sumdsi.nextSameName() ) {
		mil_eval_renameAndReparentOptVars( *dsi, fname, *sumdsi );
		for ( segi = dsi->find("Segment"), sumsegi = sumdsi->find("Segment");
		      ! segi->isNull(); segi.nextSameName(), sumsegi.nextSameName() ) {
			mil_eval_renameAndReparentOptVars( *segi, fname, *sumsegi );
		}
	}
}


QUB_Tree run_mil_eval_tree( std::vector<QUB_Tree>& data, QUB_Tree cfg,
						    QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB )
{
	// find and remove all models
	std::vector<QUB_Tree> models;
	QUB_TreeIter modi = cfg.find("ModelFile");
	while ( ! modi->isNull() ) {
		models.push_back( *modi );
		modi.remove();
		if ( modi->name() != "ModelFile" )
			modi.next("ModelFile");
	}

	std::vector<QUB_Tree> results;
	for ( int mi = 0; mi < int(models.size()); ++mi ) {
		cfg.insertChild( models[mi] );
		mil_eval_tree eval( data, cfg, rptCB, rsltCB, pctCB );
		results.push_back( eval.execute() );
		cfg.removeChild( models[mi] );
	}

	QUB_Tree result = results[0];
	if ( results.size() > 1 ) {
		mil_eval_renameVars( result, models[0] );
		for ( int ri = 1; ri < int(results.size()); ++ri )
			mil_eval_renameAndReparentVars( results[ri], models[ri], result );
	}

	return result;
}

// *************************************************************************

extern "C" QUBOPT_API QTR_Impl *
miltreeiface( QTR_Impl *cfg, QTR_Impl **data,
			  QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB )
{
	QUB_Tree cfgNode( cfg );

	std::vector<QUB_Tree> dataNodes;
	QTR_Impl **di = data;
	while ( *di )
		dataNodes.push_back( QUB_Tree( *(di++) ) );

	QUB_Tree resultNode = run_mil_eval_tree( dataNodes, cfgNode, rptCB, rsltCB, pctCB );
	
	QTR_Impl *result = resultNode.getImpl();
	if ( result ) QTR_INCREF( result );
	return result;
}

extern "C" QUBOPT_API mil_export_model *
new_mil_model( QTR_Impl *model, int nchannel ){
	mil_export_model *mdl = new mil_export_model;
	mdl->model            = new mil_model( QUB_Tree(model), 1000.0 );
	mdl->metamodel        = new mil_metamodel( *(mdl->model), nchannel );

	mdl->nclass        = mdl->model->nclass;
	mdl->classOfState  = mdl->model->mi.clazz;

	mdl->nmetastate    = mdl->metamodel->nmetastate;
	mdl->nmetaclass    = mdl->metamodel->nmetaclass;
	mdl->nmetapath     = mdl->metamodel->nmetapath;
	mdl->maxnmetastate = mdl->metamodel->maxnmetastate;
	mdl->maxnmetapath  = mdl->metamodel->maxnmetapath;

	mdl->metaclass     = mdl->metamodel->metaclass;
	mdl->nmetagroup    = mdl->metamodel->nmetagroup;
	mdl->metapr        = mdl->metamodel->mprx;
	mdl->metastate     = mdl->metamodel->metastate;
	mdl->metapath      = mdl->metamodel->metapath;
	mdl->subpath       = mdl->metamodel->subpath;
	mdl->metaindex     = mdl->metamodel->metaindex;

	return mdl;
	}

/*
	int              nmetastate, nmetaclass, nmetapath;
	fq::vector<int>    metaclass, nmetagroup;
	fq::vector<double> mprx;
	matrix<int>      metastate, metapath, subpath, metaindex;
	int              maxnmetastate, maxnmetapath;
*/

extern "C" QUBOPT_API void
del_mil_model( mil_export_model *mdl ){
	delete mdl->metamodel;
	delete mdl->model;
	delete mdl;
	}


//-----
extern "C" QUBOPT_API
QTR_Impl * QUB_Idl_ReadDWT(const char *filename){
	QUB_Tree dwt = QUB_Tree::Create("DataSet");

	dwt["sampling"].setNumData(QTR_TYPE_DOUBLE, 1, 1, (void*)0);
	double& samp = * (double *) dwt["sampling"].data();
	samp = 1;
	bool noSamplingGiven = true;

	QUB_Tree seg;

	ifstream in;
	in.open(filename);
	if ( ! in )
		return NULL;

	int len;
	char * line;
	char * pchBuf, * pTok;

	while ( ! in.eof() ) {
		int i;
		int ndwell = 0;
		int nclass = 0;
		fq::vector<double> amp, sd;
		double start, at = 0.0;
		istream_linereader rdr;

		//----------------------------- Following section rewritten 8-2004 JTB
		//----- Find a 'Segment:'
		do { 
			line = rdr.readLine(in, len);

			if( memcmp("Sampling duration (ms):",line,23)== 0 	// ICE : sampling.
					&& atof(line+24) > 0.0 ) {
				samp=atof(line+24);
				noSamplingGiven = false;
				}

			} while( ! in.eof() && memcmp(line,"Segment:", 8)!=0  );

		//----- copy read line to buffer
		pchBuf = new char[strlen(line)+1];
		strcpy(pchBuf,line);

		//----- Parse 'tokens' from input line
		pTok=strtok(pchBuf," \t");
		nclass = 0;
		start = at;

		while( (pTok=strtok(NULL," \t")) ) { // Ok to lose 1st token - it is always Segment:
			if( strcmp(pTok, "Dwells:") == 0 )			// QUB & ICE : Dwells:
			  if( (pTok=strtok(NULL," ")) )
					ndwell=atoi(pTok);

			if ( strcmp(pTok, "Sampling(ms):") == 0 )	// QUB : Sampling(ms):
			  if( (pTok=strtok(NULL," ")) ) {
					samp=atof(pTok);
					noSamplingGiven = false;
					}

			if ( strcmp(pTok, "Start(ms):") == 0 )		// QUB : Start(ms):
			  if( (pTok=strtok(NULL," ")) )
					start=max(at,atof(pTok));

			if ( strcmp(pTok, "ClassCount:") == 0 )		// QUB : ClassCount: (list of pairs of amp,sd)
			  if( (pTok=strtok(NULL," ")) ) { 
					nclass=atoi(pTok);
					amp.resize( nclass );
					sd.resize( nclass );
					for ( i=0; i<nclass && pTok!=NULL ; ++i ) { 
					  if( (pTok=strtok(NULL," ")) )
							amp[i]=atof(pTok);
					  if( (pTok=strtok(NULL," ")) )
							sd[i]=atof(pTok);
						}
					}

			if( strcmp(pTok, "Starting") == 0			// ICE : Starting time (ms):
					&& (pTok=strtok(NULL," "))
					&& strcmp(pTok, "time") == 0
					&& (pTok=strtok(NULL," "))
					&& strcmp(pTok, "(ms):") == 0
					&& (pTok=strtok(NULL," ")) )
				start=max(at,atof(pTok));

			if( strcmp(pTok, "Amp:") == 0 )				// ICE : Amp: (list of amps)
				while ( (pTok=strtok(NULL," "))!=NULL && *pTok < 'A' ) { 
					++nclass;
					amp.resize( nclass );
					sd.resize( nclass );
					amp[nclass-1]=atof(pTok);
					sd[nclass-1]=0.0;
					}
			}
		delete[] pchBuf;

		//=========================================================================
		/*  Changed to code above 8-2004 JTB  - for ICE compatibility.
		int len;
		char *line = rdr.readLine(in, len);
		while ( (ndwell == 0) && ! in.eof() ) {
			istrstream lineIn(line);
			lineIn.get(s, 16, ' ') >> i >> ws;
			lineIn.get(s, 16, ' ');
			if ( (strcmp(s, "Dwells:") == 0) && ! lineIn.eof() )
				lineIn >> ndwell >> ws;

			nclass = 0;
			start = at;
			while ( (nclass == 0) && ! lineIn.eof() ) {
				lineIn.get(s, 16, ' ');
				if ( strcmp(s, "Sampling(ms):") == 0 ) {
					lineIn >> samp >> ws;
					noSamplingGiven = false;
				}
				else if ( strcmp(s, "Start(ms):") == 0 )
					lineIn >> start >> ws;
				else if ( strcmp(s, "ClassCount:") == 0 ) {
					lineIn >> nclass >> ws;

					amp.resize( nclass );
					sd.resize( nclass );
					for ( i=0; lineIn && i<nclass; ++i )
						lineIn >> amp[i] >> sd[i];
				}
				else
					lineIn >> d >> ws;
			}
			if ( start < at )
				start = at;

			if ( ndwell == 0 )
				line = rdr.readLine(in, len);
		}
		*/

		//-------------------------------------------------------------------
		if ( ndwell ) {
			seg = dwt.insertChild( seg, "Segment" );
			seg["start"].setData(QTR_TYPE_DOUBLE, start);
			seg["Classes"].setNumData(QTR_TYPE_INT, ndwell, 1, (void*)0);
			seg["Durations"].setNumData(QTR_TYPE_FLOAT, ndwell, 1, (void*)0);
			int *idwt = (int *) seg["Classes"].data();
			float *tdwt = (float *) seg["Durations"].data();

			set<int> classesPresent;

			for ( i=0; in && i<ndwell; ++i ) {
				in >> idwt[i] >> tdwt[i];
				classesPresent.insert(idwt[i]);
				at += tdwt[i];
				if ( noSamplingGiven && (tdwt[i] < samp) )
					samp = tdwt[i];
			}
			ndwell = i; // drop events that couldn't be read

			// tune up classes: add amp, sd for missing high-numbered; renumber negative as highest number, neg amp
			int minCls = 0, maxCls = nclass - 1;
			for ( set<int>::iterator ci = classesPresent.begin();
				  ci != classesPresent.end();
				  ++ci ) {
				if ( *ci < minCls )
					minCls = *ci;
				if ( *ci > maxCls )
					maxCls = *ci;
			}
			while ( nclass <= maxCls ) {
				amp.push_back( (double) nclass );
				sd.push_back( (double) 0.1 );
				++nclass;
			}
			if ( minCls < 0 ) {
				int offset = maxCls - minCls + 1;
				for ( int ij=0; ij<ndwell; ++ij ) {
					if ( idwt[ij] < 0 )
						idwt[ij] += offset;
				}
				while ( nclass < offset ) {
					amp.push_back( (double) (nclass - offset) );
					sd.push_back( (double) 0.1 );
					++nclass;
				}
			}

			seg["DwellCount"].setData(QTR_TYPE_INT, ndwell);
			seg["amp"].setNumData(QTR_TYPE_DOUBLE, nclass, 1, &(amp[0]));
			seg["sd"].setNumData(QTR_TYPE_DOUBLE, nclass, 1, &(sd[0]));
		}
	}
	
	// HEKA (and others?) give event durations in seconds not milliseconds, and
	// QUB doesn't deal well with gigahertz sampling
	if ( noSamplingGiven && (samp < 1e-3/*megahertz*/) ) {
		samp *= 1.0e3;

		for ( QUB_TreeMonoIter segi = dwt.find("Segment"); ! segi->isNull(); segi.nextSameName() ) {
			int ndwt = (*segi)["DwellCount"].dataAsInt();
			float *tdwt = (float *) (*segi)["Durations"].data();
			for ( int i=0; i<ndwt; ++i )
				tdwt[i] *= 1.0e3f;
		}
	}
	
	QTR_INCREF( dwt.getImpl() );
	return dwt.getImpl();
}


typedef struct
{
  // short epi; // clampex files only
  int start; // samples
  float amp; 
  int dur; // samples
  float std;
  short cls;
  short notes; // bitmap: {0: too short; 1: amp manually adjusted; 2: rejected}
  float baseline; // cls 0 only
} evl_event;

#define EVL_CLOSED_BYTES 24
#define EVL_OPEN_BYTES   20


QTR_Impl * QUB_Idl_ReadEVL(char *filename)
{
	FILE *f = fopen(filename, "rb");
	if ( ! f )
		return NULL;
  
	fseek(f, 0, SEEK_END);
	int evlen = (int) ftell(f);
  
	fseek(f, 158, SEEK_SET);
	float sampling;
	fread(&sampling, sizeof(float), 1, f);
	sampling /= 1000.0f;
  
	fseek(f, 186, SEEK_SET);
	short nclass;
	fread(&nclass, sizeof(short), 1, f);
	fq::vector<double> amp(nclass, 0.0), std(nclass, 0.0);

	evl_event evt;
	int ndwt = 0;
	fq::vector<int> idwt;
	fq::vector<float> tdwt;

	int at = 256;
	int cur_sample = 0;
	float start = 0.0;

	int err;
  
	while ( at < evlen ) {
		err = fseek(f, at, SEEK_SET);
		// cerr << "seeking " << at << " (" << err << ")" << "\r\n";
    
		if ( EVL_OPEN_BYTES == (err = (int) fread(&evt, 1, EVL_OPEN_BYTES, f)) ) {
			// cerr << "read " << err << " bytes" << "\r\n";
			// start amp dur std cls
			if ( ndwt == 0 )
			  start = ((float) evt.start) * sampling;
			  /*
			  cerr << "\tstart:\t" << evt.start << "\r\n";
			  cerr << "\tamp:\t" << evt.amp << "\r\n";
			  cerr << "\tdur:\t" << evt.dur << "\r\n";
			  cerr << "\tstd:\t" << evt.std << "\r\n";
			  cerr << "\tcls:\t" << evt.cls << "\r\n";
			  cerr << "\tnotes:\t" << evt.notes << "\r\n";
			  cerr << "\r\n";
			  */
			// QUB-style missed events: append "missing" time to prev.event:
			if ( ndwt && (evt.start > cur_sample) ) {
			  tdwt[ndwt-1] += sampling * (float) (evt.start - cur_sample);
				cur_sample = evt.start;
			}
			/*
			if ( ndwt && (evt.start > cur_sample) ) {
				cerr << "breaking segment (" << evt.start << " > " << cur_sample << ")" << "\r\n";
				printSeg(nclass, amp, std, ndwt, idwt, tdwt, sampling, start);
			}
			*/
			++ndwt;
			idwt.push_back( evt.cls );
			tdwt.push_back( sampling * (float) evt.dur );
			amp[ evt.cls ] += evt.amp;
			std[ evt.cls ] += evt.std;
			cur_sample = evt.start + evt.dur;

			at += evt.cls ? EVL_OPEN_BYTES : EVL_CLOSED_BYTES;
		}
		else {
			break;
		}
	}
  
	fclose(f);

	
	QUB_Tree dwt = QUB_Tree::Create("DataSet");

	if ( ndwt ) {
		for ( int i=0; i<nclass; ++i ) {
			amp[i] /= ndwt;
			std[i] /= ndwt;
		}
			
		dwt["sampling"].setData( QTR_TYPE_DOUBLE, (double) sampling );

		QUB_Tree seg = dwt["Segment"];
		seg["start"].setData(QTR_TYPE_DOUBLE, start);
		seg["amp"].setNumData(QTR_TYPE_DOUBLE, nclass, 1, &(amp[0]));
		seg["sd"].setNumData(QTR_TYPE_DOUBLE, nclass, 1, &(std[0]));
		seg["DwellCount"].setData(QTR_TYPE_INT, ndwt);
		seg["Classes"].setNumData(QTR_TYPE_INT, ndwt, 1, &(idwt[0]));
		seg["Durations"].setNumData(QTR_TYPE_FLOAT, ndwt, 1, &(tdwt[0]));
	}
  
	QTR_INCREF( dwt.getImpl() );
	return dwt.getImpl();
}


void ReadTAC_InsertSeg(QUB_Tree dwt, QUB_Tree& seg, double start_ms, std::vector<double>& amps, std::vector<int>& idwt, std::vector<float>& tdwt)
{
	seg = dwt.insertChild(seg, "Segment");
	seg["start"].setData(QTR_TYPE_DOUBLE, start_ms);
	seg["amp"].setNumData(QTR_TYPE_DOUBLE, (int) amps.size(), 1, &(amps[0]));
	seg["DwellCount"].setData(QTR_TYPE_INT, (int) idwt.size());
	seg["Classes"].setNumData(QTR_TYPE_INT, (int) idwt.size(), 1, &(idwt[0]));
	seg["Durations"].setNumData(QTR_TYPE_FLOAT, (int) tdwt.size(), 1, &(tdwt[0]));
}

QTR_Impl * QUB_Idl_ReadTAC(char *filename)
{
	ifstream in;
	in.open(filename);
	if ( ! in )
		return NULL;

	istream_linereader rdr;
	int len;
	char * line;

	std::vector<int> idwt;
	std::vector<float> tdwt;
	std::vector<double> amps;
	int seg_ix = -1;
	double last_secs = 0.0;
	int last_level = 0;
	double smallest_dt = 1.0e5;
	double start_ms = 0.0;
	double seg_ms = 0.0;

	QUB_Tree dwt = QUB_Tree::Create("DataSet");
	QUB_Tree seg;

	while ( ! in.eof() ) {
		line = rdr.readLine(in, len);
		if ( (! len) || (line[0] < '0') || (line[0] > '9') )
			continue;

		istringstream fromline(line);
		int sweep;
		int level = -1; // to ensure it's a full line, make sure this is overwritten with non-neg
		double secs, preamp, amp, tag;
		fromline >> sweep >> secs >> preamp >> amp >> level >> tag;
		if ( level < 0 )
			continue;

		if ( seg_ix != sweep ) {
			if ( idwt.size() ) {
				ReadTAC_InsertSeg(dwt, seg, start_ms, amps, idwt, tdwt);
				idwt.clear();
				tdwt.clear();
				start_ms += seg_ms;
				last_level = 0;
			}
			seg_ix = sweep;
		}

		if ( secs < last_secs )
			last_secs = 0.0;
		double dt = 1e3 * (secs - last_secs);
		last_secs = secs;
		seg_ms += dt;
		if ( dt < smallest_dt )
			smallest_dt = dt;

		idwt.push_back(last_level);
		tdwt.push_back((float) dt);
		if ( ((int) amps.size()) <= level )
			amps.resize(level+1);
		amps[level] = amp * 1e12;
		last_level = level;
	}
	if ( idwt.size() ) {
		idwt.push_back(last_level);
		tdwt.push_back((float) smallest_dt);
		ReadTAC_InsertSeg(dwt, seg, start_ms, amps, idwt, tdwt);
	}
	// first event has unknown duration -- set all segs' first duration to smallest encountered:
	for (QUB_TreeMonoIter segi = dwt.find("Segment"); ! segi->isNull(); segi.nextSameName())
		((float*)(*segi)["Durations"].data())[0] = (float) smallest_dt;
	smallest_dt /= 12;
	dwt["sampling"].setData(QTR_TYPE_DOUBLE, smallest_dt);

	QTR_INCREF( dwt.getImpl() );
	return dwt.getImpl();
}





extern "C" QUBOPT_API double MIL_Files(char *data_path, char *model_path, int max_iter, double tdead_ms, char *output_path )
{
	QUB_Tree data( QUB_Idl_ReadDWT(data_path) );
	if ( data.isNull() )
		return -11.0; // no data
	QTR_DECREF( data.getImpl() );
	data["ExpCond"]["DeadTime"].setData(QTR_TYPE_DOUBLE, tdead_ms);

	QTR_Impl* adata[2];
	adata[0] = data.getImpl();
	adata[1] = NULL;

	QUB_Tree config( QUB_Tree::Create("Properties") );
	QUB_Tree model( QUB_Tree::Open(model_path, true) );
	if ( model.isNull() )
		return -13.0; // no model
	config.children().insert( model.clone() );
	config["MaxIterations"].setData(QTR_TYPE_INT, max_iter);
	config["Mode"].setData((max_iter > 0) ? "optimize" : "evaluate");
	config["use segments"].setData("together");

	QUB_Tree result( miltreeiface(config.getImpl(), adata, NULL, NULL, NULL) );
	if ( result.isNull() )
		return -17.0; // no MIL result
	QTR_DECREF( result.getImpl() );

	char filename[2048];
	strncpy(filename, output_path, 2000);
	int baselen = (int) strlen(filename);

	QUB_Tree model_out = result["ModelFile"];
	strncpy(filename+baselen, "_model.qmf", 40);
	model_out.saveAs(filename);

	ofstream report;
	strncpy(filename+baselen, "_report.txt", 40);
	report.open(filename);
	if ( report ) {
		report << "MIL_Files( "
		       << data_path << "," << endl
			   << "\t\t" << model_path << "," << endl
			   << "\t\t" << max_iter << "," << endl
			   << "\t\t" << tdead_ms << "," << endl
			   << "\t\t" << output_path << " )" << endl << endl;
		report << "LL:\t" << result["LL"].dataAsDouble(-7.0) << endl;
		int n = result["Gradient"].dataCount();
		if ( n ) {
			double *gg = (double*) result["Gradient"].data();
			report << "Gradient:";
			for (int i=0; i<n; ++i)
				report << "\t" << gg[i];
			report << endl;
		}
		report << "Iterations:\t" << result["Iterations"].dataAsInt(-1);
		switch ( result["ErrorCode"].dataAsInt(-240) ) {
			case -666:
				report << "* invalid starting point" << endl;
				break;
			case -1:
				report << "* evaluation failed" << endl;
				break;
			case -2:
				report << "* check failed" << endl;
				break;
			case -3:
				report << "* exceeded max iterations" << endl;
				break;
			case -240:
				report << "* incomplete results" << endl;
				break;
		}
		QUB_TreeMonoIter tci;
		int tcx;
		for (tci=result["DataSet"]["HistogramGroup"].find("TimeConstants"), tcx=1;
			 ! tci->isNull(); ++tcx, tci.nextSameName()) {
			n = (*tci)["Tau"].dataCount();
			if ( n ) {
				report << endl << "Time constants for class " << tcx << ":" << endl;
				double *tau = (double*) (*tci)["Tau"].data();
				double *amp = (double*) (*tci)["Amp"].data();
				for (int i=0; i<n; ++i)
					report << "\tTau:\t" << tau[i] << "\tAmp:\t" << amp[i] << endl;
			}
		}
		report.flush();
	}

	return result["LL"].dataAsDouble(-7.0); // def "no LL"
}


