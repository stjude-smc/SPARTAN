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

#include <stdexcept>
#include <string.h>
#include "skm_eval_tree.h"
#include "stat_eval_tree.h"
using namespace fq;

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

class skm_workspace {	
	// can idealize segments for one file (dataset with same conditions)
	public:
		skm_workspace( QUB_Tree dataset, skm_model& mdl, double *aMatrixOverride,
						QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag );
		~skm_workspace(); 

		int runSKM( bool reestimate, int maxIter, double llConv, bool fixQ, int seg=-1 ); // default: all segs

		QUB_Tree         dataSet; // filled originals
		QUB_Tree         dataOut; // results and idealization
		std::vector<QUB_Tree> segments;
		std::vector<QUB_Tree> segOuts;

	protected:
		void setSegIdealization( int seg_ix, int ndata ); // changes stateSeq to classSeq
		void setSegStat( int seg_ix );
		//void setSegHist( int seg_ix, double *data, int ndata );
		void setErrIdlAndStat( int seg_ix );

		void setupData();
		void setupArrays();

		int              maxndata;
		double           dtMS;

		skm_metamodel    model;

		QTR_Callback     fillDataCB;
		std::ostream&    milerr;
		int*             Stopped;

		double *         initialAMatrix; // ignore if null, otherwise use it instead of model rates

		fq::vector<float>  mutia, mutipr, mutiamp, mutixms; // being reestimated
		fq::vector<int>    stateSeq;
	};

typedef CountedPtr<skm_workspace> skm_worksptr;

class skm_eval_tree{
	public:
		skm_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
						QTR_Callback fillDatCB, QTR_Callback rptCB, QTR_Callback percentCB );

		QUB_Tree execute();

	protected:
		void initResults();
		void finishResults();

		bool canRun();

		skm_model model;
		std::vector<skm_worksptr> workspaces;
		QUB_Tree cfg;
		QTR_Callback fillDataCB, reportCB, pctCB;

		mil_reportstream milerr;

		bool reestimate, fixQ;
		int maxIter;
		double llConv;

		int *Stopped;

		QUB_Tree result;
	};

void skm_mutimakv(int nchannel, int nstate, int npath, int *state, int *path, 
				  float *ratio, int& nmutistate, int& nmutipath, int& nmuticlass, 
				  int*& mutistate, int*& mutipath);
float argmax(int n, float* x, int& imaxpos);
void newa(int ndata, int* istate, int nstate, float* a);
void newpr(int nstate, int* istate, float* pr);
void newamp(int ndata, int* idata, int* istate, float* bin, int* state, int nclass, float* amp, int *fix);
void newxms(int ndata, int* idata, int* istate, float* bin, int* state, int nclass, float* amp, float* xms, int *fix);
void binning_dbl(int ndata, double* data, int nbin, float* bin, int* idata);
bool viterbi(int ndata, int* data, int nstate, float* pr, float* a, int nbin, float* b, int* istate, float& ll);
int skm_dbl(int ndata, double* data, int nstate, int nclass, int* state,
			float* pr, float* a, float* amp, float* xms,  
			int* istate, float& ll, int& it,
			int fFixA, int *fixAmp, int *fixSD, int MaxIter, double llConv, int Viterbi, int *Stopped);


// ---------------------------------------------------------------------------------
template< class T >

void reorderArray( T *arr, int n, int *order )	// same as in modelTree
{
	std::vector<T> cpy( n );
	memcpy( &(cpy[0]), arr, n * sizeof(T) );
	for ( int i=0; i<n; ++i )
		arr[i] = cpy[ order[i] ];
}

QUB_Tree skm_get_multiclasses(QUB_Tree model, int channelCount)
{
	skm_model skmodel(model, NULL, NULL);
	skm_metamodel mmodel(skmodel, channelCount, 0.1/*irrelevant dtMS for A matrix*/);
	QUB_Tree result = QUB_Tree::Create("");
	
	result["amp"].setNumData(QTR_TYPE_DOUBLE, mmodel.nmuticlass, 1, (void*)0);
	double *amp = (double *) result["amp"].data();
	result["sd"].setNumData(QTR_TYPE_DOUBLE, mmodel.nmuticlass, 1, (void*)0);
	double *sd = (double *) result["sd"].data();
	for ( int i=0; i<mmodel.nmuticlass; ++i ) {
		amp[i] = mmodel.mutiamp[i];
		sd[i] = mmodel.mutixms[i];
	}
	return result;
}

extern "C" QUBOPT_API QTR_Impl * skm_get_multiclasses_iface(QTR_Impl *model, int channelCount)
{
	QUB_Tree result = skm_get_multiclasses(QUB_Tree(model), channelCount);
	QTR_Impl *resultPtr = result.getImpl();
	QTR_INCREF(resultPtr);
	return resultPtr;
}

skm_model::skm_model( QUB_Tree modelNode, int *fFixAmp, int *fFixSD )
			: node( modelNode ), 
			tree( node ) {
	int uncondensedNClass = modelNode["Amps"].dataCount();

	tree.condenseClasses();
	tree.toMdlinf( mi );
	nclass = tree.classesPresent();
	nchannel = node["ChannelCount"].dataAsInt(1);

	int i;
	path.resize( mi.npath * 4 ); // a b 0 drug
	for( i=0; i<mi.npath; ++i ) {
		path[4*i  ] = mi.path[i][0];
		path[4*i+1] = mi.path[i][1];
		path[4*i+2] = 0;
		path[4*i+3] = mi.path[i][2];
		}

	x.resize( mi.x.size() );
	for( i=0; i<x.size(); ++i )
		x[i] = float(mi.x[i]);

	datinf di;
	tree.toDatinf(di);

	amp.resize( nclass );
	xms.resize( nclass );
	for ( i=0; i<nclass; ++i ) {
		amp[i] = float(di.i[i]);
		xms[i] = float(sqrt( di.r[i][0] ));
		}

	fixAmp.resize( uncondensedNClass );
	fixSD.resize( uncondensedNClass );
	for ( i=0; i<uncondensedNClass; ++i ) {
		fixAmp[i] = fFixAmp ? fFixAmp[i] : 0;
		fixSD[i]  = fFixSD  ? fFixSD[i]  : 0;
		}

	int * condensation = (int *) modelNode["Condensed"]["fromCondensed"].data();
	reorderArray( &(fixAmp[0]), uncondensedNClass, condensation );
	reorderArray( &(fixSD[0]), uncondensedNClass, condensation );

	q.resize( mi.nstate * mi.nstate );
	xtoq(mi.nstate, mi.npath, path, NULL, 0, x, q, 0, NULL, NULL);

	tree.restoreClasses();
	}

// ********************************************************************************************

skm_metamodel::skm_metamodel( skm_model& mdl, int nchan, double dtMS )
			: model( mdl ), 
			nchannel( nchan ) {
	mdlinf& mi = mdl.mi;
	int nstate = mi.nstate;
	int i;

	// ratio: same as amp but offset so ratio[0] = 0.0; used to classify mutistates
	fq::vector<float> ratio( model.nclass );
	// excessNoise: like ratio for std.dev; used to construct multichannel xms
	fq::vector<float> excessNoise( model.nclass );
	ratio[0] = 0.0;
	excessNoise[0] = 0.0;
	float var0 = model.xms[0];
	var0 *= var0;

	for ( i=1; i<model.nclass; ++i ) {
		ratio[i] = model.amp[i] - model.amp[0];
		excessNoise[i] = model.xms[i] * model.xms[i] - var0;
		}

	skm_mutimakv( nchannel, nstate, mi.npath, mi.clazz, model.path, ratio, 
					nmutistate, nmutipath, nmuticlass, mutistate, mutipath);

	cmutistate.resize( nmutistate );
	for ( i=0; i<nmutistate; ++i )
		cmutistate[i] = mutistate[i*(nstate+1)+0]; // mutistate[i][0] -- muticlass of mutistate i

	mutiq.resize( nmutistate * nmutistate );
	mutia.resize( nmutistate * nmutistate );

	try {
		qtomutiq(nstate, mi.npath, model.path, nmutistate, mutistate, nmutipath, 
				mutipath, model.q, mutiq, 0, NULL, NULL);
		mexp(nmutistate, float(dtMS*1.0e-3), mutiq, mutia, 0, NULL, NULL, 0); 
		}
	catch (...) {
		//unsigned int fstat = _clearfp(); // just in case
		//_fpreset();
		/*
		if ( fstat & _SW_OVERFLOW )
			milerr << "Numerical Overflow: generating A matrix" << endl;
		else if ( fstat & _SW_UNDERFLOW )
			milerr << "Numerical Underflow: generating A matrix" << (*segi+1) << endl;
		else if ( fstat & _SW_ZERODIVIDE )
			milerr << "Division by Zero: generating A matrix" << (*segi+1) << endl;
		else
			milerr << "Exception in SKM: generating A matrix" << endl;      */
		}

	mutipr.resize( nmutistate );
	for ( i=0; i<nmutistate; ++i )
		mutipr[i] = float(1.0/(float)nmutistate);

	mutiamp.resize( nmuticlass );
	mutixms.resize( nmuticlass );
	mutiFixAmp.resize( nmuticlass );
	mutiFixSD.resize( nmuticlass );

	for (i=0; i<nmuticlass; i++) {
		mutiamp[i] = model.amp[0]; // closed_amp + sum(open_amp === ratio)
		//mutixms[i] = 0.0;          // sqrt( sum( xms*xms ) )
		// now this should probably be // sqrt( closed_noise^2 + sum(excess_noise^2) ) ??
		mutixms[i] = var0; // will be sqrt( closed_noise^2 [var0] + sum(excessNoise^2) )
		//float sum=0.0; // # states contributing to this mutistate; should always be == nchannel
		bool allAmpFixed = true;
		bool allSDFixed = true;

                int k;
		for (k=0; mutistate[k*(nstate+1)+0]!=i; k++)    // mutistate[k][0]
			; // find the first mutistate in this muticlass.
		      // if (e.g. 2 sub-opens == 1 open 1 closed)
		      //   then the first example may not have representative s.d. and allXFixed
		
		for (int j=0; j<nstate; j++) {                      // mutistate[k][j+1]
			int countStateJ = mutistate[k*(nstate+1)+j+1];
			mutiamp[i] += ((float) countStateJ)*ratio[mi.clazz[j]];
			mutixms[i] += ((float) countStateJ)*excessNoise[mi.clazz[j]];
			//sum += countStateJ;
			                       // find if there is a contributing non-fixed class -- this one doesn't count unless
			if ( countStateJ &&    // it's a contributor (multiclass i involves one or more channel in class of j)
			     (mi.clazz[j]      // and either this one is not the baseline
				  || (i==0)  ) ) { // or we're in the baseline multiclass (since baseline contributes only to multi-baseline)
				if ( ! model.fixAmp[mi.clazz[j]] )  // all it takes is one unfixed contributor to make the multiclass unfixed.
					allAmpFixed = false;
				if ( ! model.fixSD[mi.clazz[j]] )
					allSDFixed = false;
			}
		}
		// mutixms[i] /= sum;
		mutixms[i] = float(sqrt( mutixms[i] ));
		mutiFixAmp[i] = allAmpFixed;
		mutiFixSD[i] = allSDFixed;
	}
}

skm_metamodel::~skm_metamodel(){
	delete [] mutistate;
	delete [] mutipath;
	}	

// ********************************************************************************************
skm_workspace::skm_workspace( QUB_Tree dataset, skm_model& mdl, double *aMatrixOverride,
							 QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag )
			: dataSet( dataset ),
			  dataOut( QUB_Tree::Create("DataSet") ), 
			  dtMS( dataset["sampling"].dataAsDouble( 1.0e-2 ) ),
			  model( mdl, (int) dataset["ExpCond"].find("ChannelCount")->dataAsDouble( mdl.nchannel ), dtMS ),
			  fillDataCB( fillDatCB ), 
			  milerr( errstream ), 
			  Stopped( StopFlag ),
			  initialAMatrix( aMatrixOverride ){
	setupData();
	setupArrays();
	}

skm_workspace::~skm_workspace(){
	}	


#pragma warning (disable: 4101)	// Disable unreferenced variables message 
int skm_workspace::runSKM( bool reestimate, int maxIter, double llConv, bool fixQ, int seg ) {
	int rtnVal = 0;
	
	std::vector<int> segs;
	if ( seg < 0 )
		for ( int i=0; i<int(segments.size()); ++i )
			segs.push_back( i );
		else if ( seg < int(segments.size() ) )
			segs.push_back( seg );
		
	for ( std::vector<int>::iterator segi = segs.begin(); segi != segs.end(); ++segi ) {
		try {
			QUB_Tree segIn = segments[ *segi ];
			if ( (*(segIn.find("Channel"))).isNull() ) {
				QTR_DoCallback(fillDataCB, segIn);
				segIn["SKMMustFreeData"];
				}
			QUB_Tree channel = segIn["Channel"];
			double *data = (double *) channel.data();
			int ndata = channel.dataCount();
			
			// get fresh model
			mutia = model.mutia;
			mutipr = model.mutipr;
			mutiamp = model.mutiamp;
			mutixms = model.mutixms;
			
			// if there's an a matrix override, overwrite mutia
			if ( initialAMatrix ) {
				for ( int i=model.nmutistate * model.nmutistate - 1; i>=0; --i )
					mutia[i] = float(initialAMatrix[i]);
				}
			
			QUB_Tree segOut = segOuts[ *segi ];
			
			float ll = 0.0;
			int it = 0;
			
			if ( ! ( Stopped && *Stopped ) ) {
				// run skm
				int skmErr = skm_dbl(ndata, data, model.nmutistate, model.nmuticlass, model.cmutistate, 
					mutipr, mutia, mutiamp, mutixms, 
					stateSeq, ll, it, fixQ, model.mutiFixAmp, model.mutiFixSD, maxIter, llConv, !reestimate, Stopped);
				segOuts[ *segi ]["ErrorCode"].setData( QTR_TYPE_INT, skmErr );
				
				if ( skmErr == 0 ) {
					setSegIdealization( *segi, ndata ); // from stateSeq
					setSegStat( *segi ); // from mutiamp, mutisd, mutia
					}
				else
					setErrIdlAndStat( *segi );
				
				segOut["LL"].setData( QTR_TYPE_DOUBLE, (double) ll );
				segOut["Iterations"].setData( QTR_TYPE_INT, it );
				}
			else {
				setErrIdlAndStat( *segi );
				segOuts[ *segi ]["ErrorCode"].setData( QTR_TYPE_INT, -2 );
				}
			
			QUB_TreeIter mustFree = segIn.find("SKMMustFreeData");
			if ( ! (*mustFree).isNull() ) {
				mustFree.remove();
				segIn.find("Channel").remove();
				}
			}
		catch (const std::overflow_error &oe) {
			milerr << "Numerical Overflow: segment " << (*segi+1) << endl;
			++rtnVal;
			}
		catch (const std::underflow_error &oe) {
			milerr << "Numerical Underflow: segment " << (*segi+1) << endl;
			++rtnVal;
			}
		catch (const std::exception &ee) {
			milerr << "Exception: " << ee.what() << ": segment " << (*segi+1) << endl;
			++rtnVal;
			}
		catch (...) {
			milerr << "Exception in SKM Core: segment " << (*segi+1) << endl;
			++rtnVal;
			}
		}
	return rtnVal;
	}
#pragma warning ( default : 4101)


QUB_Tree MakeIdealNode( int *byPoint, int ndata, float dtMS ){
	QUB_Tree result = QUB_Tree::Create("");
	
	int ii, di; // ideal, data indices
	int curIdl, points, ndwt;
	
	// count dwells
	curIdl = byPoint[0];
	ndwt = 1;
	
	for ( di=0; di<ndata; ++di ) {
		if ( byPoint[di] != curIdl ) {
			curIdl = byPoint[di];
			++ndwt;
			}
		}
	
	// setup storage
	result["DwellCount"].setData( QTR_TYPE_INT, ndwt );
	
	QUB_Tree classes = result["Classes"];
	classes.setNumData( QTR_TYPE_INT, ndwt, 1 );
	int *idwt = (int *) classes.data();
	
	QUB_Tree durs    = result["Durations"];
	durs.setNumData( QTR_TYPE_FLOAT, ndwt, 1 );
	float *tdwt = (float *) durs.data();
	
	// make dwells
	curIdl = byPoint[0];
	points = 1;
	ii = 0;
	
	for ( di=1; di<ndata; ++di ) {
		if ( byPoint[di] != curIdl ) {
			idwt[ii] = curIdl;
			tdwt[ii] = ((float) points) * dtMS;
			++ii;
			curIdl = byPoint[di];
			points = 1;
			}
		else {
			++points;
			}
		}
	idwt[ii] = curIdl;
	tdwt[ii] = ((float) points) * dtMS;
	
	return result;
	}

void skm_workspace::setSegIdealization( int seg_ix, int ndata ){
	QUB_Tree seg = segOuts[ seg_ix ];
	int i;
	
	QUB_Tree idlState = MakeIdealNode( stateSeq, ndata, float(dtMS) );
	QUB_Tree idwtState = idlState["Classes"];
	QUB_Tree tdwtState = idlState["Durations"];
	idwtState.setName( "States" );
	tdwtState.setName( "StateDurations" );
	// wait and append 'em at the end
	
	for ( i=0; i<ndata; ++i )
		stateSeq[i] = model.cmutistate[stateSeq[i]];
	
	QUB_Tree idl = MakeIdealNode( stateSeq, ndata, float(dtMS) );
	seg.appendChild( idl["DwellCount"] );
	seg.appendChild( idl["Classes"] );
	seg.appendChild( idl["Durations"] );
	
	seg.appendChild( idwtState );
	seg.appendChild( tdwtState );
	}

void skm_workspace::setSegStat( int seg_ix ){
	QUB_Tree seg = segOuts[ seg_ix ];
	QUB_Tree amp = seg["amp"];
	QUB_Tree sd  = seg["sd"];
	QUB_Tree A = seg["A"];
	
	amp.setNumData( QTR_TYPE_DOUBLE, model.nmuticlass, 1 );
	sd.setNumData( QTR_TYPE_DOUBLE, model.nmuticlass, 1 );
	A.setNumData( QTR_TYPE_DOUBLE, model.nmutistate, model.nmutistate );
	
	double *ampData = (double *) amp.data();
	double *sdData = (double *) sd.data();
	double *aData = (double *) A.data();
	
	for ( int i=0; i<model.nmuticlass; ++i ) {
		ampData[i] = mutiamp[i];
		sdData[i] = mutixms[i];
	}
	
	for ( int i=model.nmutistate * model.nmutistate - 1; i>=0; --i )
	  aData[i] = mutia[i];
}

void skm_workspace::setErrIdlAndStat( int seg_ix ){
	setSegStat( seg_ix );
	segOuts[ seg_ix ]["DwellCount"].setData( QTR_TYPE_INT, 0 );
	}

//-----
// replicate dataSet in dataOut, with the exception of DataSet["Segment"]["Channel"]
// get references to dataOut["Segment"]s in segOuts
// get size and pointers from dataSet["Segment"]["Channel"]s in ndata, datas, maxndata
void skm_workspace::setupData() {
	maxndata = 0;
	int ndat;
	
	QUB_TreeIter dsout = dataOut.end();
	
	for ( QUB_TreeMonoIter metai = dataSet.children(); ! (*metai).isNull(); ++metai ) {
		QUB_Tree meta = *metai;
		if ( meta.name() == "Segment" ) {
			segments.push_back( meta );
			
			ndat = meta.dataAs( 1, (int)0 ) - meta.dataAs( 0, (int)0 )+ 1;
			if ( ndat > maxndata )
				maxndata = ndat;
			
			QUB_Tree seg = meta.clone(false);
			segOuts.push_back( seg );
			dsout.insert( seg );
			++dsout;
			
			QUB_TreeIter segout = seg.end();
			for ( QUB_TreeMonoIter segmetai = meta.children(); ! (*segmetai).isNull(); ++segmetai ) {
				QUB_Tree segmeta = *segmetai;
				if ( segmeta.name() != "Channel" ) {
					segout.insert( segmeta.clone() );
					++segout;
					}
				}
			}
		else {
			dsout.insert( meta.clone() );
			++dsout;
			}
		}
	}

void skm_workspace::setupArrays(){
	stateSeq.resize( maxndata );
	}

//----------------------------------------------------------------------------
skm_eval_tree::skm_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
				  QTR_Callback fillDatCB, QTR_Callback rptCB, QTR_Callback percentCB  )
		  : model( config["ModelFile"].clone(), (int*) config["FixAmp"].data(), (int*) config["FixSD"].data() ),
			cfg( config ), 
			fillDataCB( fillDatCB ), 
			reportCB( rptCB ), 
			pctCB( percentCB ),
			milerr( rptCB ) {
	
	reestimate = (0!=config["Reestimate"].dataAsInt(0));
	maxIter    = config["MaxIterations"].dataAsInt(10);
	llConv     = config["ConvLL"].dataAsDouble(0.01);
	fixQ       = (0!=config["FixQ"].dataAsInt(0));

	double *amatOverride = 0;
	QUB_Tree amatOverNode = *( config.find("InitialAMatrix") );
	if ( ! amatOverNode.isNull() )
		amatOverride = (double *) amatOverNode.data();
	
	QUB_Tree stopFlag = * (config.find("StopFlag"));
	if ( stopFlag.isNull() )
		Stopped = 0;
	else
		Stopped = (int *) stopFlag.data();

	for ( std::vector<QUB_Tree>::iterator dsi = dataSets.begin(); dsi != dataSets.end(); ++dsi )
		workspaces.push_back( skm_worksptr( new skm_workspace( *dsi, model, amatOverride, fillDataCB, milerr, Stopped ) ) );
	}

void skm_eval_tree::initResults(){
	result = QUB_Tree::Create("SKM Result");
	result.appendClone( cfg );

	for ( int i=0; i<int(workspaces.size()); ++i ) {
		skm_workspace& wspace = * workspaces[i];
		result.appendChild( wspace.dataOut );
		}
	}

void skm_eval_tree::finishResults(){
	}

bool skm_eval_tree::canRun(){
	return true;
	}

QUB_Tree skm_eval_tree::execute(){
	initResults();
	
	QUB_Tree pctNode = QUB_Tree::Create("");
	pctNode.setData( QTR_TYPE_INT, 0 );
	int *pctData = (int *) pctNode.data();
	
	int nseg = 0;
	for ( int ws = 0; ws < int(workspaces.size()); ++ws )
	  nseg += (int) workspaces[ws]->segments.size();
	int nsegDone = 0;
	
	if ( canRun() ) {
		for ( int i=0; i<int(workspaces.size()); ++i ) {
			for ( int j=0; j<int(workspaces[i]->segments.size()); ++j, ++nsegDone ) {
				workspaces[i]->runSKM( reestimate, maxIter, llConv, fixQ, j );
				*pctData = (100 * nsegDone) / nseg;
				QTR_DoCallback(pctCB, pctNode);
				
				if ( Stopped && *Stopped )
					break;
				}
			
			if ( Stopped && *Stopped )
				break;
			}
		finishResults();
		}
	
	return result;
	}

//-----------------------------------------------------------------------------
#define SKM_MAX_POINTS 1000000

QUB_Tree skm_uldt_makeSubsegs( QUB_Tree dataSet, std::vector<QUB_Tree>& segs, int nstate, std::vector<QUB_Tree>& subsegs, std::vector<int>& subsegCount ){
	int maxpoint = SKM_MAX_POINTS / nstate;

	QUB_Tree subsegDS = QUB_Tree::Create("DataSet");
	subsegDS.appendClone( dataSet["FileName"] );
	subsegDS.appendClone( * dataSet.find("ChannelCount") );
	subsegDS.appendClone( dataSet["sampling"] );
	subsegDS.appendClone( dataSet["ADChannelCount"] );
	subsegDS.appendClone( dataSet["ActiveChannel"] );

	QUB_Tree procData = * ( dataSet.find("ProcessData") );
	if ( ! procData.isNull() )
		subsegDS.appendClone( procData );

	QUB_TreeIter subsegput = subsegDS.end();
	for ( QUB_TreeMonoIter segi = dataSet.find("Segment"); ! (*segi).isNull(); segi.nextSameName() ) {
		QUB_Tree seg_in = *segi;
		segs.push_back( seg_in );

		int *bounds = (int *) seg_in.data();
		int npoint = bounds[1] - bounds[0] + 1;
		int nsubs = npoint / maxpoint + ( (npoint % maxpoint) ? 1 : 0 );
		subsegCount.push_back( nsubs );

		int *subbounds = NULL; // bounds of most recent subseg
		for ( int isub = 0; isub < nsubs; ++isub ) {
			QUB_Tree subseg = seg_in.clone(); // yeah this copies wrong start time and length, but they will not be used or kept
			subsegDS.appendChild( subseg );
			subsegs.push_back( subseg );

			subbounds = (int *) subseg.data();
			subbounds[0] = bounds[0] + isub * maxpoint;
			subbounds[1] = subbounds[0] + maxpoint - 1;
			}
                if ( subbounds )
		  subbounds[1] = bounds[1]; // final subseg may be short
		}
	return subsegDS;
	}

void skm_ultd_avgAmpSD( QUB_Tree seg, std::vector<QUB_Tree>& subsegs ) {
	// and a-matrix weighted by # points
	QUB_Tree subamp = subsegs[0]["amp"];
	QUB_Tree subsd = subsegs[0]["sd"];
	QUB_Tree subA = subsegs[0]["A"];
	int nclass = subamp.dataCount();
	int nstate = subA.dataRows();
	std::vector<double> zeroPerClass( nclass, 0.0 );
	std::vector<double> zeroPerStates( nstate*nstate, 0.0 );
	double ll = 0.0;
	
	QUB_Tree amp = seg["amp"];
	amp.setNumData( QTR_TYPE_DOUBLE, nclass, 1, &(zeroPerClass[0]) );
	double *ampData = (double *) amp.data();
	
	QUB_Tree sd = seg["sd"];
	sd.setNumData( QTR_TYPE_DOUBLE, nclass, 1, &(zeroPerClass[0]) );
	double *sdData = (double *) sd.data();
	
	QUB_Tree A = seg["A"];
	A.setNumData( QTR_TYPE_DOUBLE, nstate, nstate, &(zeroPerStates[0]) );
	double *aData = (double *) A.data();
	
	std::vector<int> subsize( subsegs.size() );
	int sumsize = 0;
	
	for ( int i=0; i<int(subsegs.size()); ++i ) {
		QUB_Tree sub = subsegs[i];
		int *subbounds = (int *) sub.data();
		subsize[i] = subbounds[1] - subbounds[0] + 1;
		
		if ( ! (subsegs[i]["amp"].dataCount() && subsegs[i]["sd"].dataCount() && subsegs[i]["A"].dataCount()) )
		  continue;
		
		double *subampData = (double *) subsegs[i]["amp"].data();
		double *subsdData = (double *) subsegs[i]["sd"].data();
		double *subaData = (double *) subsegs[i]["A"].data();
		
		// careful of segs with errors and thus missing info?
		for ( int cls=0; cls<nclass; ++cls ) {
			ampData[cls] += subampData[cls] * subsize[i];
			sdData[cls]  += subsdData[cls] * subsdData[cls] * subsize[i];
			}
		for ( int j=nstate*nstate-1; j>=0; --j )
		  aData[j] += subaData[j] * subsize[i];
		
		ll += subsegs[i]["LL"].dataAs(0, 0, 0.0);
		sumsize += subsize[i];
		}
	
	if ( sumsize ) {
	  for ( int cls=0; cls<nclass; ++cls ) {
	    ampData[cls] = ampData[cls] / sumsize;
	    sdData[cls]  = sqrt( sdData[cls] / sumsize );
	  }
	  for ( int i=nstate*nstate-1; i>=0; --i )
	    aData[i] /= sumsize;
	}
	
	seg["LL"].setData(QTR_TYPE_DOUBLE, ll);
}

QUB_Tree skm_ultd_redoSubseg( QUB_Tree subseg, QUB_Tree fullseg, QUB_Tree config, 
							 QTR_Callback fillDataCB, QTR_Callback reportCB ) {
	QUB_Tree dataSet = fullseg.parent();

	QUB_Tree subsegDS = QUB_Tree::Create("DataSet");
	subsegDS.appendClone( dataSet["FileName"] );
	subsegDS.appendClone( * dataSet.find("ChannelCount") );
	subsegDS.appendClone( dataSet["sampling"] );
	subsegDS.appendClone( subseg );

	QUB_Tree procData = * ( dataSet.find("ProcessData") );
	if ( ! procData.isNull() )
		subsegDS.appendClone( procData );

	std::vector<QUB_Tree> datasets;
	datasets.push_back( subsegDS );

	QUB_Tree subsegConfig = config.clone();
	subsegConfig["Reestimate"].setData( QTR_TYPE_INT, 0 );

	QUB_Tree mdlAmp = subsegConfig["ModelFile"]["Amps"];
	QUB_Tree mdlSD  = subsegConfig["ModelFile"]["Stds"];
	double *mdlAmpData = (double *) mdlAmp.data();
	double *mdlSDData  = (double *) mdlSD.data();

	QUB_Tree avgAmp = fullseg["amp"];
	QUB_Tree avgSD  = fullseg["sd"];
	double *avgAmpData = (double *) avgAmp.data();
	double *avgSDData  = (double *) avgSD.data();

	int nclass = (mdlAmp.dataCount() < avgAmp.dataCount()) ? mdlAmp.dataCount() : avgAmp.dataCount();
	for ( int i=0; i<nclass; ++i ) {
		mdlAmpData[i] = avgAmpData[i];
		mdlSDData[i]  = avgSDData[i];
		}

	skm_eval_tree subeval( datasets, subsegConfig, fillDataCB, reportCB, QTR_NullCallback );
	return subeval.execute()["DataSet"]["Segment"];
	}

// returns # of unusual subsegs re-viterbi'd
// unusual: abs( seg.amp - avgamp ) > avgsd for some class
							 int skm_ultd_redoUnusual( QUB_Tree seg, std::vector<QUB_Tree>& subsegs, QUB_Tree config, 
						 QTR_Callback fillDataCB, QTR_Callback reportCB ) {
	
	if ( (int) seg.parent().find("ChannelCount")->dataAsDouble( config["ModelFile"].find("ChannelCount")->dataAsInt(1) ) > 1 )
		return 0; // until we can estimate single-channel amps from multi-channel unconstrained estimates
	
	int num_unusual = 0;
	
	double *avgamp = (double *) seg["amp"].data();
	double *avgsd  = (double *) seg["sd"].data();
	int nclass = seg["amp"].dataCount();
	
	for ( int i=0; i<int(subsegs.size()); ++i ) {
		double *segamp = (double *) subsegs[i]["amp"].data();
		bool isUnusual = false;
		for ( int cls=0; cls<nclass; ++cls ) {
			if ( fabs( segamp[cls] - avgamp[cls] ) > avgsd[cls] ) {
				isUnusual = true;
				break;
				}
			}
		
		if ( isUnusual ) {
			++num_unusual;
			subsegs[i] = skm_ultd_redoSubseg( subsegs[i], seg, config, fillDataCB, reportCB );
			}
		}
	
	return num_unusual;
	}

void skm_ultd_concatIdl( QUB_Tree seg, std::vector<QUB_Tree>& subsegs ){
	std::vector<int> subndwt;
	int sumndwt = 0;
	int i;

	for ( i=0; i<int(subsegs.size()); ++i ) {
		int ndwt = subsegs[i]["DwellCount"].dataAsInt(0);
		subndwt.push_back( ndwt );
		sumndwt += ndwt;
		}

	QUB_Tree clss = seg["Classes"];
	clss.setNumData( QTR_TYPE_INT, sumndwt, 1 );
	int *idwt = (int *) clss.data();

	QUB_Tree durs = seg["Durations"];
	durs.setNumData( QTR_TYPE_FLOAT, sumndwt, 1 );
	float *tdwt = (float *) durs.data();

	sumndwt = 0; // now == dwells copied so far
	int lastClass = -1; // class of last dwell copied, to merge adjacent dwells of same class
	for ( i=0; i<int(subsegs.size()); ++i ) {
		int *subidwt = (int *) subsegs[i]["Classes"].data();
		float *subtdwt = (float *) subsegs[i]["Durations"].data();
		if ( subidwt[0] == lastClass ) {
			subtdwt[0] += tdwt[sumndwt-1];
			--sumndwt;
			}
		memcpy( idwt + sumndwt, subidwt, subndwt[i] * sizeof(int) );
		memcpy( tdwt + sumndwt, subtdwt, subndwt[i] * sizeof(float) );
		sumndwt += subndwt[i];
		lastClass = idwt[sumndwt-1];
		}
	
	// careful of segs with errors and thus missing info?
	seg["DwellCount"].setData( QTR_TYPE_INT, sumndwt );
	clss.resizeData( sumndwt );
	durs.resizeData( sumndwt );
	}

void skm_ultd_concatStateIdl( QUB_Tree seg, std::vector<QUB_Tree>& subsegs ){
	std::vector<int> subndwt;
	int sumndwt = 0;
	int i;

	for ( i=0; i<int(subsegs.size()); ++i ) {
		int ndwt = subsegs[i]["StateDurations"].dataCount();
		subndwt.push_back( ndwt );
		sumndwt += ndwt;
	}

	QUB_Tree stts = seg["States"];
	stts.setNumData( QTR_TYPE_INT, sumndwt, 1 );
	int *idwt = (int *) stts.data();

	QUB_Tree durs = seg["StateDurations"];
	durs.setNumData( QTR_TYPE_FLOAT, sumndwt, 1 );
	float *tdwt = (float *) durs.data();

		// careful of segs with errors and thus missing info?

	sumndwt = 0; // now == dwells copied so far
	void * v1, * v2;
	for ( i=0; i<int(subsegs.size()); ++i ) {
		// Note : rolling up v1, v2 into memcpy causes internal compile error in vc7.
		// Try again after a service pack ... ?  {JB} 2/2005
		v1=subsegs[i]["States"].data();
		v2=subsegs[i]["StateDurations"].data();
		memcpy( idwt + sumndwt, v1, subndwt[i] * sizeof(int) );
		memcpy( tdwt + sumndwt, v2, subndwt[i] * sizeof(float) );
		sumndwt += subndwt[i];
	}
}

void skm_ultd_dropFirstLast( QUB_Tree seg, QUB_Tree config )
{
	if ( config["DropFirstLast"].dataAsInt(0) ) {
		int dropCls = -1;
		if ( config["DropIfClass"].dataAsInt(0) )
			dropCls = config["DropClass"].dataAsInt();

		double dt = seg.parent()["sampling"].dataAsDouble(0.01); // in millisec just like the dwells
		int ndwt = seg["DwellCount"].dataAsInt();
		QUB_Tree clss = seg["Classes"];
		int *idwt = (int *) clss.data();
		QUB_Tree durs = seg["Durations"];
		float *tdwt = (float *) durs.data();
		int *segBounds = (int *) seg.data();

		float dropFirstTm = 0.0;
		int dropFirstPts = 0;
		if ( (dropCls < 0) || (idwt[0] == dropCls) ) {
			dropFirstTm = tdwt[0];
			dropFirstPts = (int) floor(dropFirstTm / dt + 0.5);
		}

		float dropLastTm = 0.0;
		int dropLastPts = 0;
		if ( (dropCls < 0) || (idwt[ndwt-1] == dropCls) ) {
			dropLastTm = tdwt[ ndwt-1 ];
			dropLastPts = (int) floor(dropLastTm / dt + 0.5);
		}

		int eventsDropped = 0;
		if ( dropFirstPts ) {
			++ eventsDropped;

			memmove( idwt, idwt+1, (ndwt-1) * sizeof(int) );
			memmove( tdwt, tdwt+1, (ndwt-1) * sizeof(float) );
		}
		if ( dropLastPts ) {
			++eventsDropped;
		}
		if ( eventsDropped ) {
			ndwt -= eventsDropped;
			QUB_Tree stts = seg["States"];
			QUB_Tree stdurs = seg["StateDurations"];

			if ( ndwt <= 0 ) {
				seg["DwellCount"].setData( QTR_TYPE_INT, 0 );
				seg["length"].setData( QTR_TYPE_DOUBLE, 0.0 );
				clss.clearData();
				durs.clearData();
				stts.clearData();
				stdurs.clearData();
			}
			else {
				segBounds[0] += dropFirstPts;
				segBounds[1] -= dropLastPts;

				seg["DwellCount"].setData( QTR_TYPE_INT, ndwt );
				clss.resizeData( ndwt );
				durs.resizeData( ndwt );

				QUB_Tree startNode = seg["start"];
				startNode.setData( QTR_TYPE_DOUBLE, startNode.dataAsDouble() + dropFirstTm );

				QUB_Tree lengthNode = seg["length"];
				lengthNode.setData( QTR_TYPE_DOUBLE, lengthNode.dataAsDouble() - dropFirstTm - dropLastTm );

				int *sidwt = (int *) stts.data();
				int sndwt = stts.dataCount();
				float *stdwt = (float *) stdurs.data();

				int dropFirstStateDwells = 0;
				while ( dropFirstTm > 0.0 ) { // warning: destroying dropFirstTime from here on out
					dropFirstTm -= stdwt[ dropFirstStateDwells ];
					++dropFirstStateDwells;
				}
				int dropLastStateDwells = 0;
				while ( dropLastTm > 0.0 ) {
					dropLastTm -= stdwt[ sndwt - dropLastStateDwells - 1 ];
					++dropLastStateDwells;
				}
				
				if ( dropFirstStateDwells ) {
					memmove( sidwt, sidwt+dropFirstStateDwells, (sndwt - 1) * sizeof(int) );
					memmove( stdwt, stdwt+dropFirstStateDwells, (sndwt - 1) * sizeof(float) );
				}
				sndwt -= dropFirstStateDwells + dropLastStateDwells;
				stts.resizeData( sndwt );
				stdurs.resizeData( sndwt );
			}
		}
	}
}

QUB_Tree skm_ultd_rejoin( QUB_Tree resultDS, QUB_Tree origDS,
						 std::vector<QUB_Tree>& segs_in, QUBOPT_VAR_NOT_USED std::vector<QUB_Tree>& subsegs_in, std::vector<int>& subsegCount,
						  QUB_Tree config, QTR_Callback fillDataCB, QTR_Callback reportCB )
{
	QUB_Tree joinedDS = origDS.clone();

	QUB_TreeMonoIter segi = joinedDS.find("Segment");
	QUB_TreeMonoIter subsegi = resultDS.find("Segment");
	for ( int si = 0; si < int(segs_in.size()); ++si, segi.nextSameName() ) {
		QUB_Tree seg_out = *segi;

		std::vector<QUB_Tree> subsegs_out;
		for ( int ssi = 0; ssi < subsegCount[si]; ++ssi, subsegi.nextSameName() )
			subsegs_out.push_back( *subsegi );

		skm_ultd_avgAmpSD( seg_out, subsegs_out ); // and a-matrix

		if ( skm_ultd_redoUnusual( seg_out, subsegs_out, config, fillDataCB, reportCB ) )
			skm_ultd_avgAmpSD( seg_out, subsegs_out );

		skm_ultd_concatIdl( seg_out, subsegs_out );
		skm_ultd_concatStateIdl( seg_out, subsegs_out );

		skm_ultd_dropFirstLast( seg_out, config );
		}

	return joinedDS;
	}

void skm_ultd_mergeIdlAndStat( QUB_Tree idlDS, QUB_Tree statDS ){
	QUB_TreeMonoIter idli, stati;
	for ( idli = idlDS.find("Segment"), stati = statDS.find("Segment");
		  !( (*idli).isNull() || (*stati).isNull() );
		  idli.nextSameName(), stati.nextSameName() ) {
		QUB_Tree iseg = *idli;
		QUB_Tree statseg = *stati;
		QUB_TreeIter iseg_inserter = iseg.children(); // insert at beginning

		while ( ! statseg.child().isNull() ) {
			if ( iseg.find( statseg.child().name() )->isNull() )
				iseg_inserter.insert( statseg.child() );
			else
				statseg.removeChild( statseg.child() );
			}
		}
	}

QUB_Tree skm_ultd_oneDS( QUB_Tree dataSet, QUB_Tree config, QTR_Callback fillDataCB, 
						QTR_Callback reportCB, QTR_Callback pctCB, bool doStat ) {
	int nstate = 0;
	for ( QUB_TreeMonoIter statei = config["ModelFile"]["States"].find("State"); ! statei->isNull(); statei.nextSameName() )
		++nstate;

	std::vector<QUB_Tree> seg_in;
	std::vector<QUB_Tree> subseg_in;
	std::vector<int> subsegCount;
	QUB_Tree subsegDS = skm_uldt_makeSubsegs( dataSet, seg_in, nstate, subseg_in, subsegCount );

	std::vector<QUB_Tree> dataSets_in;
	dataSets_in.push_back( subsegDS );

	skm_eval_tree eval( dataSets_in, config, fillDataCB, reportCB, pctCB ); // maybe this stage should count for less than 100%?
	QUB_Tree subsegResult = eval.execute();

	QUB_Tree joinedDS;
	//QUB_Tree X = config["StopFlag"]; 
	if ( ! config["StopFlag"].dataAsInt() )
		joinedDS = skm_ultd_rejoin( subsegResult["DataSet"], dataSet, seg_in, subseg_in, subsegCount,
											 config, fillDataCB, reportCB );
	
	// let QUB handle stat unless...
	if ( doStat && (! config["StopFlag"].dataAsInt()) ) {
		QUB_Tree statDS = skm_ultd_runStat( joinedDS, config, fillDataCB, reportCB );
		statDS = statDS["DataSet"];
		skm_ultd_mergeIdlAndStat( joinedDS, statDS );
		}

	return joinedDS;
	}

QUB_Tree skm_unlimited( std::vector<QUB_Tree>& dataSets, QUB_Tree config, QTR_Callback fillDataCB, 
					   QTR_Callback reportCB, QTR_Callback pctCB, bool doStat ) {
	QUB_Tree result = QUB_Tree::Create("SKM Result");
	result.appendClone( config );

	for ( int i=0; i<int(dataSets.size()); ++i ) {
		QUB_Tree dsResult = skm_ultd_oneDS( dataSets[i], config, fillDataCB, reportCB, pctCB, doStat );
		result.appendChild( dsResult );
		}
	return result;
	}

//----------------------------------------------------------------------------
extern "C" QUBOPT_API QTR_Impl *
skmtreeiface( QTR_Impl *cfg, QTR_Impl **data, QTR_Callback fillDataCB, 
			  QTR_Callback rptCB, QTR_Callback pctCB ) {
	QUB_Tree cfgNode( cfg );

	std::vector<QUB_Tree> dataNodes;
	QTR_Impl **di = data;
	while ( *di )
		dataNodes.push_back( QUB_Tree( *(di++) ) );

	QUB_Tree resultNode = skm_unlimited( dataNodes, cfgNode, fillDataCB, rptCB, pctCB );
	
	QTR_Impl *result = resultNode.getImpl();
	QTR_INCREF( result );
	return result;
	}

//-----------------------------------------------------------------------------
// multi-channel Markov model creation : this function will make a multichannel 
// markov model out of a single channel markov model
// Input nchannel, nstate, npath, state, path, ratio
// Output nmutistate, nmutipath, nmuticlass, mutistate, mutipath
//-----------------------------------------------------------------------
void skm_mutimakv(int nchannel, int nstate, int npath, int *state, int *path, float *ratio, 
				  int& nmutistate, int& nmutipath, int& nmuticlass, int*& mutistate, int*& mutipath){
	int i,k,l,m,n,s,t, iPaths=0, iState;
	
    //-- Determine whether open is negative or positive and use to set sort order 
	// for muticlasses in eclass().
	// Note : The try catch is probably not effective here.  Accessing beyond the 
	// bounds of the array would likely run into more data and not throw an exception.
	// Better to ensure that the passed array ratio[] has an element [1] always.
	double openDirection;		
	try {			
		openDirection = ratio[1];
		}
	catch (...) {
		openDirection = 1.0;
		}

	//-- Precalc model size for memory allocation 
	// Calculates the number of Combinations of channels/states per row in mutistate.
	// For each state, add to temp log(state+nchann)-log(state) and then antilog on the sum
	float tmp=0.0;
	for( n=1; n<nstate; n++) 
		tmp += float(log(double(n+nchannel))-log(double(n))); 
	n = int(exp(tmp)+1);						//estimates nMutiState - may be 1 greater than necc.
	mutistate = new int[n * (nstate+1)];		//create empty mutistate array.
	float * aCond = new float[n];		// Conductance class sums
	
	// Create a 2d style array for accessing mutistate array 
	int ** aMutiState = new int * [n];
	for( i=0; i<n; ++i ) 
		aMutiState[i] = &(mutistate[i*(nstate+1)]);

	//----- Fill in mutistate , acond, iPaths.
	// mutistate is an integer array of (estimated nMutiState) * (nstate+1)
	// Each (nState+1) group is a 1-based array of number of channels in each state,
	// and [0] element is set to the eclass result class for the row
	// iterate through combinations of states / channels and set nMutiState
	// aCond is the effective conductance <aCond[]> = sum[all states] ( channels * ratio ) 
	// iPaths is the number of mutistate connections which have only one channel transition
	nmutistate=0;
	CChannelCombine oMeta(nchannel, nstate);
	do {
		aCond[nmutistate]=0.0;
		for ( i=0; i<nstate; ++i ) {
			iState=oMeta.State0(i);
			aMutiState[nmutistate][i+1]=iState;
			aCond[nmutistate] += ((float) iState) * ratio[state[i]];
			if( iState > 0 )								// pre calc # of paths 
				iPaths += (nstate-1);
			}
		++nmutistate;
		} while( oMeta.GetNextCombine() );

	//----- conductance classes 
	float * aClasses = new float[nmutistate];	// output buffer, discarded unused
	int *  ia = new int[nmutistate];			// mapping of data to classes
	
	eclass( aCond, nmutistate, nmuticlass, aClasses, ia, openDirection);		// Fill C,ia with classes
	
	for (n=0; n<nmutistate; n++) 
		aMutiState[n][0] = ia[n];
	
	delete [] aCond;
	delete [] aClasses;
	delete [] ia;

	//-----
	// precalc the number of transition paths for memory allocation
	// iterates through all combinations of mutistates and counts those 
	// which differ only by 1 transition.
	// New method : for each non zero state, add paths for the number of 
	//  places one of the channels could go.  This is precalced above.
	// This calculates a maximum possible, but the actual number used is less
	/*
	n=0;
	for(k=0; k<nmutistate; k++)  
		for (l=0; l<nmutistate; l++) { 
			s=0;
			for (i=0; i<nstate; i++)  
				s += int(fabs(aMutiState[k][i+1]-aMutiState[l][i+1])); 
			if (s==2) 
				n++; 
			} 
	mutipath = new int[n*4];
	*/
	mutipath = new int[iPaths*4];
	
	// transition paths - for each mutipath pair (m,n), if they differ 
	// by only 1 channel state change then find the original states (k,l)
	// and check if a path exists between k  & l -->> add a mutipath if so.
	//
	// Note : Consider a transition where  two channels change state in a 
	// single delta T - is this accounted for elsewhere ? {JB}
	// In a simple 2 state 3 channel model currently (C,C,C) we are only 
	// setting mutipaths for (C,C,O),(C,O,C),(O,C,C) but all states 
	// including (O,O,O) are reachable in a single delta T ?
	//
	nmutipath = 0; 
	for( m=0; m<nmutistate; m++ ) {
		for( n=0; n<nmutistate; n++ ) {
			s=0;
			for( i=0; i<nstate; i++) 
				s += abs(aMutiState[m][i+1]-aMutiState[n][i+1]);
			if (s==2) {	// if muti[m] can become muti[n] by only 1 channel changing state 
				// Find k,l the source and dest states 
				for (i=0; i<nstate; i++) {
					t = aMutiState[m][i+1]-aMutiState[n][i+1];
					if (t==1) 
						k = i;
					if (t==-1) 
						l = i;
					}
				// Check for a path from k to l in model 
				for (i=0; i<npath; i++) {
					if (path[i*4+0]==k && path[i*4+1]==l) {   
						// add a nutipath 
						mutipath[nmutipath*4+0] = m;
						mutipath[nmutipath*4+1] = n;
						mutipath[nmutipath*4+2] = k;
						mutipath[nmutipath*4+3] = l;
						nmutipath++;
						}
					}
				}
			}
		}

	delete [] aMutiState;
	}

//------------------------------------------------------------------------
// skm core routine
int skm_dbl(int ndata, double* data, int nstate, int nclass, int* state,
			float* pr, float* a, float* amp, float* xms,  
			int* istate, float& ll, int& it,
			int fFixA, int *fFixI, int *fFixSD, int MaxIter, double tol, 
			int Viterbi, int *Stopped){
    if ( Stopped && *Stopped )
		return -2;
	
	int i;
	
	// binning the data
	int nbin  = 100;
	// we compromise between std.dev resolution and number of bins
	double minxms = xms[0];
	for (i=1; i<nclass; ++i)
	  minxms = min(minxms, (double) xms[i]);
	if ( minxms == 0.0 ) minxms = 1e-3;
	nbin = max(10, min(10240, (int) (20 * fabs(amp[nclass-1] - amp[0]) / minxms)));
	minxms = 5 * fabs(amp[nclass-1] - amp[0]) / nbin;
	if ( minxms == 0.0 ) minxms = 1e-3;
	for (i=0; i<nclass; ++i)
	  xms[i] = (float) max(minxms, (double) xms[i]);

	fq::vector<float> bin( nbin );
	fq::vector<int> idata( ndata );
	binning_dbl(ndata,data,nbin,bin,idata);
	
	// skm loop for idealization and reestimation
	fq::vector<float> b( nstate * nbin );
	int fFixPr = 0;
	float  oldll = 0.0f, newll;
	
	for (it=0; it<MaxIter; it++, oldll=ll=newll) {
		for (i=0; i<nstate; i++) 
			for (int j=0; j<nbin; j++)    
				b[i*nbin+j] = gaussian(amp[state[i]],xms[state[i]],bin[j]);
		
		if ( ! viterbi(ndata,idata,nstate,pr,a,nbin,b,istate,newll) )
			return -1;
		
		if ( Stopped && *Stopped )
			return -2;
		
		if ( !Viterbi ) {
		  
			if ( !fFixPr )
				newpr(nstate,istate,pr);
			if ( !fFixA )
				newa(ndata,istate,nstate,a);
			newamp(ndata,idata,istate,bin,state,nclass,amp,fFixI);
			newxms(ndata,idata,istate,bin,state,nclass,amp,xms,fFixSD);
			}
		
		if (it>1 && fabs(newll-oldll)<tol){
			ll = newll;
			break;
			}
		
		if (Viterbi)  {
			ll = newll;
			break;
		}
	}
	
	// free memory (now done automatically by autodeleter)
	return 0;
}

//------------------------------------------------------------------------
// Viterbi detection
bool viterbi(int ndata, int* data, int nstate, float* pr, float* a,
			 int nbin, float* b, int* istate, float& ll) {
	// data[] : bin index for each data point
	// istate[] : state index for each data point
	// transform probabilities into log domain
	fq::vector<float> logBx( nstate * nbin );
	fq::vector<float> logAx( nstate * nstate );
	fq::vector<float> logPrx( nstate );
	float* logB  = logBx;
	float* logA  = logAx;
	float* logPr = logPrx;
	int i,j;

	for( i=0; i<nstate; i++) {
		for( j=0; j<nbin; j++)       // logB[i][j], b[i][j]
			logB[i*nbin+j] = float( (b[i*nbin+j]>0) ? log(b[i*nbin+j]) : -1.0e10 );
               
		for( j=0; j<nstate; j++)         // logA[i][j], a[i][j]
			logA[i*nstate+j] = float( (a[i*nstate+j]>0) ? log(a[i*nstate+j]) : -1.0e10 );
      
		logPr[i] = float( (pr[i]>0.0) ? log(pr[i]) : -1.0e10 );
		}
	
	// Viterbi recursion
	fq::vector<float> fx(nstate * 2), gx(nstate * ndata), ux(nstate);
	float *f = fx, *g = gx, *u = ux;
	int fdi; // f data index; alternating 0 or 1

	for( i=0; i<nstate; i++){	
		f[i*2] = logPr[i] + logB[i*nbin+data[0]];	// f[i][0], logB[i][data[0]]
					      //f[i*ndata+0] = logPr[i] + logB[i*nbin+data[0]];   ??
		g[i*ndata+0] = 0;       // g[i][0]
		}
	
	int t;
	for (t=1; t<ndata; t++) 
		for( j=0; j<nstate; j++) {
			fdi = (t - 1) % 2;
			for (int ij=0; ij<nstate; ij++)
				u[ij] = f[ij*2+fdi]+logA[ij*nstate+j]+logB[j*nbin+data[t]];
            //u[i] = f[i*ndata+t-1]+logA[i*nstate+j]+logB[j*nbin+data[t]];
			// f[i][t-1], logA[i][j], logB[j][data[t]]      
			int imax;
			f[j*2+(t%2)] = argmax(nstate,u,imax);  // f[j][t] // sets imax | u[imax] = argmax(...)
			//f[j*ndata+t] = argmax(nstate,u,i);  // f[j][t] // sets imax | u[imax] = argmax(...)
			g[j*ndata+t] = float(imax);                   // g[j][t]
			}
	
	// backtrack to recover the state sequence 
	for( i=0; i<nstate; i++)
		u[i] = f[i*2+((ndata-1)%2)];             // f[i][ndata-1]
	//u[i] = f[i*ndata+ndata-1];             // f[i][ndata-1]
	int imax;
	ll = argmax(nstate,u,imax);
	istate[ndata-1] = imax;
	
	for (t=ndata-2; t>=0; t--) {
		i = istate[t+1];
		istate[t] = int(g[i*ndata+t+1]);           // g[i][t+1]
		}
	
	return true;
	}

//-----------------------------------------------------------------------
// Returns max element in x[0..n-1] and sets imaxpos to position of its 
// first occurrence.
float argmax(int n, float* x, int& imaxpos){
	float xm = x[0];
	imaxpos = 0;
	for( int j=1; j<n; j++) {
		if (xm<x[j]) {
			xm = x[j];
			imaxpos = j;
			}
		}
	
	return xm;
	}	

//-----------------------------------------------------------------------
// k-means reestimation
void newa(int ndata, int* istate, int nstate, float* a){
	int i, j, s, t;

	// sizeof a is 4 since it is a ptr.  This will set the first member only to 0.
	// TODO - could this possibly be the intention ? 
	// If so, better choice would be sizeof(*a) or sizeof(float)
	// More likely it should be memset(a,0,nState*nState*sizeof(float))
	// Also check if the a array happens to always be 0's on entry to this function.
	// {JB}
	//memset(a, 0, sizeof(a));
	memset(a, 0, nstate*nstate*sizeof(float));
	
	for (t=0; t<ndata-1; t++) {
		i = istate[t];
		j = istate[t+1];
		a[i*nstate+j] += 1.0f;                  // a[i][j]
		}
	
	for (i=0; i<nstate; i++) {
		s=0;
		for(j=0; j<nstate; j++) 
			s += int(a[i*nstate+j]);                 // a[i][j]
		if (s > 0.0) {
			for (j=0; j<nstate; j++) 
			  a[i*nstate+j] /= (float) s;              // a[i][j]
			}
		}
	}

//-----------------------------------------------------------------------
void newpr(int nstate, int* istate, float* pr){
	/*
	for (int i=0; i<nstate; i++) 
		pr[i] = (istate[0]==i) ? float(1.0) : float(0.0);
	*/
	for( int i=0; i<nstate; i++) 
		pr[i] = float(0.0);
	if(istate[0]>=0 && istate[0]<nstate) //TODO - always true ? -- remove
		pr[istate[0]]=1.0;
	}

//-----------------------------------------------------------------------
void newamp(int ndata, int* idata, int* istate, float* bin, 
			int* state, int nclass, float* amp, int *fixAmp ) {
	float * newamp = new float[nclass];
	int * count = new int[nclass];
	int i, t;
	
	for( i=0; i<nclass; i++) {
		if ( ! fixAmp[i] ) {
			newamp[i] = 0.0;
			count[i] = 0;
			}
		}
  
	for( t=0; t<ndata; t++) {
		i = state[istate[t]];
		if ( ! fixAmp[i] ) {
			newamp[i] += bin[idata[t]];
			count[i] += 1;
			}
		}

	for (i=0; i<nclass; i++)
		if ( count[i] && ! fixAmp[i] )
		  amp[i] = newamp[i] / (float) count[i];
		
	delete [] newamp;
	delete [] count;
	}

//-----------------------------------------------------------------------
void newxms(int ndata, int* idata, int* istate, float* bin,
			int* state, int nclass, float* amp, float* xms, int *fix ) {
	float * newxms = new float[nclass];
	int * count = new int[nclass];
	int i;
	
	for( i=0; i<nclass; i++) {
		if ( ! fix[i] ) {
			newxms[i] = 0.0;
			count[i] = 0;
			}
		}
	
	for( int t=0; t<ndata; t++) {
		i = state[istate[t]];
		if ( ! fix[i] ) {
			newxms[i] += (bin[idata[t]]-amp[i])*(bin[idata[t]]-amp[i]);
			count[i] += 1;
			}
		}
	
	for( i=0; i<nclass; i++)
		if ( count[i] && ! fix[i] )
		  xms[i] = float(sqrt(newxms[i]/(float)count[i]));
		
	delete [] newxms;
	delete [] count;
	}

//-------------------------------------------------------------------------
// Given data[0..nData-1] and nBin number of bins, calculate bin[o..nBin-1]
// array of bins and idata[0..nData-1] array of bin #'s for data items. 
void binning_dbl(int ndata, double* data, int nbin, float* bin/*median*/, int* idata){
	int i;

	//-- Calculate min and max values 
	double nMax=data[0], nMin=data[0];
	for (i=1; i<ndata; i++){
		if(nMax<data[i])
			nMax = data[i];
		if(nMin>data[i])
			nMin = data[i];
		}

	//-- Calculate bins
	double binwidth = (nMax - nMin) / nbin;
	if ( binwidth == 0.0 ) 
		binwidth = 1.0; // avoid divide by 0
	double halfwidth = binwidth / 2.0;
	for ( i=0; i<nbin; i++ )
		bin[i] = float(nMin + halfwidth + i * binwidth);
	
	//-- Set bin # for each data 
	for ( i=0; i<ndata; i++ )
		idata[i] = max(0, min(nbin-1, int(floor((data[i]-nMin)/binwidth))  ));
	}
