/* Copyright 1998-2012 Research Foundation State University of New York   */

/* This file is part of QUB Express.                                      */

/* QUB Express is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by   */
/* the Free Software Foundation, either version 3 of the License, or      */
/* (at your option) any later version.                                    */

/* QUB Express is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          */
/* GNU General Public License for more details.                           */

/* You should have received a copy of the GNU General Public License,     */
/* named LICENSE.txt, in the QUB Express program directory.  If not, see  */
/* <http://www.gnu.org/licenses/>.                                        */

#include <stdexcept>
#include "amp_eval_tree.h"
#include "mil_eval_tree.h"
#include "mpl_eval_tree.h"
#include "qmatrixutil.h"
#include "qublib.h"
#include "skm_eval_tree.h"

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

class amp_workspace {
	// can evaluate segments for one file (dataset with same conditions)
	public:
		amp_workspace( QUB_Tree dataset, mil_model& mdl, QUB_Tree histoSettings,
						QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag );
		~amp_workspace(); // release metastate

		int segCount();
		int classCount();

		int runAmp( bool autoInit, bool autoInitUp, bool doIdealize, int seg=-1 ); // default: all segs

		fq::vector<double> z;      // constrained parameters
		bool             zValid; // false if the initial guess didn't meet constraints

		QUB_Tree         dataSet; // filled originals
		QUB_Tree         dataOut; // results and idealization
		std::vector<QUB_Tree> segments;
		std::vector<QUB_Tree> segOuts;

	protected:
		bool runViterbi( int seg_ix );
		//bool runStat( int seg_ix );

		bool resetZ();
		bool updateZ( fq::vector<double>& newz );

		void setupData();
		void setupArrays();

		void initStartingAmps( int seg_ix, bool autoInit, bool autoInitUp );
                void setAmpResult( QUB_Tree segOut, double ll, std::vector<int>& state_out );

		int              maxndata;

		mpl_metamodel    model;

		QUB_Tree         histSettings;
		QTR_Callback     fillDataCB;
		std::ostream&    milerr;
		int*             Stopped;

		fq::vector<double> C;
		double           V, dt;
		int              nfir;
		fq::vector<double> fir;

		// constraints: x = mtx * z + vct
		matrix<double>   mtx;
		fq::vector<double> vct;

		// from runAmp():
		int              nmetastate,**metastate;
		matrix<double>   q, qq, a, dx, dfa;
		tensor<double>   dq, dqq, da;
		matrix<double*>  dqAsPtr, dqqAsPtr, daAsPtr;
	};

typedef CountedPtr<amp_workspace> amp_worksptr;

class amp_eval_tree{
	public:
		amp_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
						QTR_Callback fillDatCB, QTR_Callback rptCB, QTR_Callback pctCB );

		QUB_Tree execute();

	protected:
		void initResults();
		void finishResults();

		bool canRun();

		mil_model model;
		std::vector<amp_worksptr> workspaces;
		QUB_Tree cfg;
		QTR_Callback fillDataCB, reportCB, percentCB;

		mil_reportstream milerr;

		bool autoInit, autoInitUp, doIdealize;

		int *Stopped;

		QUB_Tree result;
	};

//--------------------------------------------------------------------------------------

amp_workspace::amp_workspace( QUB_Tree dataset, mil_model& mdl, QUB_Tree histoSettings,
							 QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag )
			: dataSet( dataset ),
			  dataOut( QUB_Tree::Create("DataSet") ),
			  model( mdl, (int) dataset["ExpCond"].find("ChannelCount")->dataAsDouble( mdl.nchannel ) ),
			  histSettings( histoSettings ),
			  fillDataCB( fillDatCB ), 
			  milerr( errstream ), 
			  Stopped( StopFlag ),
			  mtx( 2 * mdl.mi.npath + 1, 2 * mdl.mi.npath + 1 ), 
			  vct( 2 * mdl.mi.npath + 1 ) {
	V    = dataset["ExpCond"]["Voltage"].dataAsDouble( 0.0 );
	dt   = 1.0e-3 * dataset["sampling"].dataAsDouble( 1.0e-2 );

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

	metastate = 0;

	QUB_Tree firNode = dataset["FIR"];
	nfir = firNode.dataCount();
	fir.resize( nfir + 1 );
	fir[0] = 1.0;
	for ( int i=0; i<nfir; ++i )
		fir[i] = firNode.dataAs( i, (double) 0.0 );

	setupData();
	setupArrays();
	}

amp_workspace::~amp_workspace(){
	if ( metastate )
		free2((char **) metastate);
	}

int amp_workspace::segCount(){
  return (int) segments.size();
	}

int amp_workspace::classCount(){
	return model.nmuticlass;
	}

void discardMBaumOutput(QUBOPT_VAR_NOT_USED char * msg) {
	}

#pragma warning (disable: 4101)	// Disable unreferenced variables message 
int amp_workspace::runAmp( bool autoInit, bool autoInitUp, bool doIdealize, int seg ){
	int rtnVal = 0;
	int processFlag = 2; // started
	
	std::vector<int> segs;
	if ( seg < 0 )
		for ( int i=0; i<int(segments.size()); ++i )
			segs.push_back( i );
		else if ( seg < int(segments.size()) )
			segs.push_back( seg );
		
	mdlinf& mi = model.model.mi;
	int i, j;
	
	  for ( std::vector<int>::iterator segi = segs.begin(); segi != segs.end(); ++segi ) {
		try {
			updateZ( z ); // reset a matrix and stuff
			
			QUB_Tree segIn = segments[ *segi ];
			if ( (*(segIn.find("Channel"))).isNull() ) {
				QTR_DoCallback(fillDataCB, segIn);
				segIn["AMPMustFreeData"];
				}
			QUB_Tree channel = segIn["Channel"];
			double *data = (double *) channel.data();
			int ndata = channel.dataCount();
			
			initStartingAmps( *segi, autoInit, autoInitUp );
			
			QUB_Tree segOut = segOuts[ *segi ];
			double *amp = (double *) segOut["amp"].data();
			double *xms = (double *) segOut["sd"].data();
			
			double ll = 0.0;
			
			for (i=0; i<model.nmutistate; i++) 
				for (j=0; j<model.nmutistate; j++) 
					if (a[i][j]>=1.0 || a[i][j]<0.0) {
						milerr << "Invalid A Matrix (" << a[i][j] << ")" << endl;
						ll += 1.0e10; // as in mil_eval...dunno why but it keeps it from getting stuck in lnsrch
						return ++rtnVal;
					}
					
			if ( Stopped && *Stopped )
			  return ++rtnVal;

			std::vector<int> state_out(ndata);
					
			mbaum(ndata, data, model.nmutistate, a, nmetastate, metastate, mi.pr,
			      model.nmuticlass, amp, xms, model.nar, model.ar, nfir, fir, &ll, & state_out[0],
			      discardMBaumOutput, &processFlag);
					
			setAmpResult( segOut, ll, state_out );
			
			if ( doIdealize ) {
			  runViterbi( *segi );
			  // runStat( *segi );
			}
			
			QUB_TreeIter mustFree = segIn.find("AMPMustFreeData");
			if ( ! (*mustFree).isNull() ) {
			  mustFree.remove();
			  segIn.find("Channel").remove();
			}
		}
		catch (const std::overflow_error &oe) {
			// reason number 0x80007 to hate microsoft:
#ifdef _WIN32
			unsigned int fstat = _clearfp(); // == 0x80007; should be 0
			// when you catch a floating-point overflow exception,
			// you have to clear the fp status register yourself,
			// or you will get an identical exception next time you use fp
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
				milerr << "Exception in AMP Core: segment " << (*segi+1) << endl;
			
			++rtnVal;
			}
		}

	return rtnVal;
	}
#pragma warning ( default : 4101)


template<class To, class From>
void copyvecmem( fq::vector<To>& to, From *from, int count ){
	to.resize( count );
	for ( int i=0; i<count; ++i )
		to[i] = (To) from[i];
	}

template<class To, class From>
void copyvec( fq::vector<To>& to, fq::vector<From>& from ){
	copyvecmem( to, &(from[0]), from.size() );
	}

bool amp_workspace::runViterbi( int seg_ix ){
	// temp. move segIn to its own DataSet for skm
	// hold pos in input DataSet using segpos
	// and replace it when done
	// (to pass in the already loaded data)

	QUB_Tree segIn = segments[ seg_ix ];

	QUB_TreeIter segpos = dataSet.children();
	while ( (! (*segpos).isNull()) && (! (*segpos).equals(segIn) ) )
	   ++segpos;
	segpos.remove();
	QUB_Tree ds = QUB_Tree::Create("DataSet");
	ds.appendClone( dataSet["FileName"] );

	ds.appendClone( dataSet["sampling"] );
	ds.appendClone( dataSet["ExpCond"] );
	ds.appendClone( dataSet["ADChannelCount"] );
	ds.appendClone( dataSet["ActiveChannel"] );
	ds.appendClone( dataSet["Channels"] );
	
	QUB_Tree procData = * ( dataSet.find("ProcessData") );
	if ( ! procData.isNull() )
		ds.appendClone( procData );

	ds.appendChild( segIn );
	std::vector<QUB_Tree> datasets;
	datasets.push_back( ds );

	QUB_Tree config = QUB_Tree::Create("Properties");
	config["Reestimate"].setData( QTR_TYPE_INT, 0 );
	config.appendClone( histSettings );
	
	// update model["Amps"], model["Stds"] with amp results -- careful of consolidation -- maybe need mdlTree.fromDatinf(di)?
	// this is most surely wrong with >1 channels, but that is not supported (yet) for AMP
	ModelTree skmodel( model.model.node.clone() );
	config.appendChild( skmodel.node );
	datinf di;
	skmodel.toDatinf( di );
	double *amp = (double *) segOuts[seg_ix]["amp"].data();
	double *xms = (double *) segOuts[seg_ix]["sd"].data();
	for ( int i=0; i<di.i.size(); ++i ) {
		di.i[i] = amp[i];
		di.r[i][0] = xms[i] * xms[i];
		}
	skmodel.fromDatinf( di );
	skmodel.restoreClasses();
	

	QUB_Tree result = skm_unlimited( datasets, config, fillDataCB, QTR_NullCallback, QTR_NullCallback, true );
	
	QUB_TreeIter segOut = segOuts[ seg_ix ].end();
	QUB_Tree skmseg = result["DataSet"]["Segment"];
	std::vector<QUB_Tree> setaside;
	while ( ! skmseg.child().isNull() ) {
		string name = skmseg.child().name();
		if ( name == "Channel" || name == "AMPMustFreeData" ) {
			setaside.push_back( skmseg.child() );
			skmseg.removeChild( skmseg.child() );
			}
		else if ( segOuts[ seg_ix ].find(name)->isNull() )
			segOut.insert( skmseg.child() );
		else
			skmseg.removeChild( skmseg.child() );
		}
	for ( std::vector<QUB_Tree>::iterator sai = setaside.begin(); sai != setaside.end(); ++sai )
		skmseg.appendChild( *sai );

	QUB_Tree histGroup = segOuts[ seg_ix ]["HistogramGroup"];
	ostringstream ost;
	ost << "Segment " << (seg_ix + 1);
	histGroup["Title"].setData( ost.str() );

	segpos.insert( segIn );
	
	return true;
	}

bool amp_workspace::resetZ(){
	if ( metastate )
		free2((char **) metastate);

	metamakv(model.nmutistate,model.nmuticlass,model.mutistate,model.nar+nfir,
			&nmetastate,&metastate);

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

bool amp_workspace::updateZ( fq::vector<double>& newz ){
	mdlinf& mi = model.model.mi;

	z = newz;

	ztox(z.size(), 2*mi.npath, mtx, vct, z, mi.x, dx);
	xtoq(mi.nstate,mi.npath,mi.path,C,V,mi.x,q,z.size(),dx,dqAsPtr); 
	qtomutiq(mi.nstate,mi.npath,mi.path,model.nmutistate,model.mutistate,model.nmutipath,
			model.mutipath,q,qq,z.size(),dqAsPtr,dqqAsPtr);
	mexp(model.nmutistate,dt,qq,a,z.size(),dqqAsPtr,daAsPtr,0); 

	return true;
	}

void amp_workspace::setupData() {
	// replicate dataSet in dataOut, with the exception of DataSet["Segment"]["Channel"]
	// get references to dataOut["Segment"]s in segOuts
	// get size and pointers from dataSet["Segment"]["Channel"]s in ndata, datas, maxndata
	maxndata = 0;
	int ndat;
	
	QUB_TreeIter dsout = dataOut.end();
	
	for ( QUB_TreeMonoIter metai = dataSet.children(); ! (*metai).isNull(); ++metai ) {
		QUB_Tree meta = *metai;
		if ( meta.name() == "Segment" ) {
			segments.push_back( meta );
			
			ndat = meta.dataAs(1, (int)0) - meta.dataAs(0, (int)0) + 1;
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

void amp_workspace::setupArrays(){
	mdlinf& mi = model.model.mi;
	zValid = resetZ();

	int nx = mi.x.size();
	int nz = z.size(); // why was this 1 before?
	int i;

	q.resize( mi.nstate, mi.nstate );
	qq.resize( nmetastate, nmetastate );
	a.resize( nmetastate+1, nmetastate+1 );
	dx.resize( nx, nz );
	dfa.resize( nmetastate, nmetastate );
	dq.resize( mi.nstate, mi.nstate, nz );
	dqq.resize( nmetastate, nmetastate, nz );
	da.resize( nmetastate, nmetastate, nz );

	dqAsPtr.resize( mi.nstate, mi.nstate );
	for ( i=0; i<mi.nstate; ++i )
		for ( int j=0; j<mi.nstate; ++j )
			dqAsPtr[i][j] = &(dq[i][j][0]);

	dqqAsPtr.resize( nmetastate, nmetastate );
	for ( i=0; i<nmetastate; ++i )
		for ( int j=0; j<nmetastate; ++j )
			dqqAsPtr[i][j] = &(dqq[i][j][0]);

	daAsPtr.resize( nmetastate, nmetastate );
	for ( i=0; i<nmetastate; ++i )
		for ( int j=0; j<nmetastate; ++j )
			daAsPtr[i][j] = &(da[i][j][0]);
	}

void amp_workspace::initStartingAmps( int seg_ix, bool autoInit, bool autoInitUp ) {
	// setup result amp, sd, ar, and maybe init
	// a seg's amps to evenly-spaced between lowest and highest points
	QUB_Tree segOut = segOuts[ seg_ix ];
	QUB_Tree ampNode = *(segOut.find("amp"));
	if ( ampNode.isNull() ) {
		ampNode = segOut.appendClone( model.model.tree.node["Amps"] );
		ampNode.setName( "amp" );
		ampNode.resizeData( model.nmuticlass );
		}

	QUB_Tree sdNode = *(segOut.find("sd"));
	if ( sdNode.isNull() ) {
		sdNode = segOut.appendClone( model.model.tree.node["Stds"] );
		sdNode.setName( "sd" );
		sdNode.resizeData( model.nmuticlass );
		}

	if ( autoInit ) {
		double *amp = (double *) ampNode.data();
		double *xms = (double *) sdNode.data();
		
		QUB_Tree channel = segments[seg_ix]["Channel"];
		double *data = (double *) channel.data();
		int ndata = channel.dataCount();
		
		double lo, hi, x, inc;
		int j, k;
		
		lo = hi = data[0];
		for ( j=1; j<ndata; ++j ) {
			x = data[j];
			if ( x < lo )
				lo = x;
			if ( x > hi )
				hi = x;
			}
		if ( model.nmuticlass == 1 ) {
			amp[0] = (hi + lo) / 2.0;
			inc = (hi - lo);
		}
		else if ( model.nmuticlass > 1 ) {
			inc = (hi - lo) / (model.nmuticlass - 1);
			x = lo;
			for ( k=0; k<model.nmuticlass; ++k, x += inc )
				amp[k] = x;
			}
		for ( k=0; k<model.nmuticlass; ++k )
			xms[k] = inc / 10.0;
		
		if ( ! autoInitUp ) {
			fq::vector<double> ampcopy( model.nmuticlass );
			fq::vector<double> xmscopy( model.nmuticlass );
			memcpy( ampcopy, amp, model.nmuticlass * sizeof(double) );
			memcpy( xmscopy, xms, model.nmuticlass * sizeof(double) );
			
			for ( int ki=0, kk=model.nmuticlass-1; ki<model.nmuticlass; ++ki, --kk ) {
				amp[ki] = ampcopy[kk];
				xms[ki] = xmscopy[kk];
			}
			}
		}
	}


void setAmpIdeal(QUB_Tree seg, std::vector<int>& state_out, int *clazz)
{
  std::vector<int> ff, ll, cc;
  int last_c = -1;
  int n = 0;
  int at = seg.dataAsInt(0);
  ff.push_back(at);
  for (std::vector<int>::iterator st=state_out.begin(); st!=state_out.end(); ++st, ++at) {
    int c = clazz[*st];
    if ( c != last_c ) {
      if ( n ) {
	cc.push_back(last_c);
	ll.push_back(at-1);
	ff.push_back(at);
      }
      last_c = c;
      n = 1;
    } else {
      ++n;
    }
  }
  ll.push_back(at-1);
  cc.push_back(last_c);

  seg["DwellCount"].setData(QTR_TYPE_INT, (int) ff.size());
  seg["Firsts"].setNumData(QTR_TYPE_INT, (int) ff.size(), 1, (void*) & ff[0]);
  seg["Lasts"].setNumData(QTR_TYPE_INT, (int) ll.size(), 1, (void*) & ll[0]);
  seg["Classes"].setNumData(QTR_TYPE_INT, (int) cc.size(), 1, (void*) & cc[0]);
}

void amp_workspace::setAmpResult( QUB_Tree segOut, double ll, std::vector<int>& state_out )
{
  int n = segOut["amp"].dataCount();
  double *amp = (double *) segOut["amp"].data();
  segOut["delta_amp"].setNumData(QTR_TYPE_DOUBLE, n, 1, (void*)0);
  double *delta_amp = (double *) segOut["delta_amp"].data();
  for (int i=0; i<n; ++i)
    delta_amp[i] = amp[i] - amp[0];
  
  segOut["LL"].setData( QTR_TYPE_DOUBLE, ll );

  setAmpIdeal(segOut, state_out, model.model.mi.clazz);
  // model.model.mi.clazz
  // clazz; segOut.dataAsInt(0); ff, ll, cc; DwellCount
}

// *************************************************************************
amp_eval_tree::amp_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
							  QTR_Callback fillDatCB, QTR_Callback rptCB, 
							  QTR_Callback pctCB  )
		: model( config["ModelFile"], 1000.0 ), // searchLimit (not used)
		  cfg( config ), 
		  fillDataCB( fillDatCB ), 
		  reportCB( rptCB ), 
		  percentCB( pctCB ),
		  milerr( rptCB ) {
	autoInit = (0!=config["AutoInit"].dataAsInt(0));
	autoInitUp = (0!=config["AutoInitUp"].dataAsInt(0));
	doIdealize = (0!=config["Idealize"].dataAsInt(0));

	QUB_Tree stopFlag = * (config.find("StopFlag"));
	if ( stopFlag.isNull() )
		Stopped = 0;
	else
		Stopped = (int *) stopFlag.data();

	QUB_Tree histSettings = config["Histograms"];

	for ( std::vector<QUB_Tree>::iterator dsi = dataSets.begin(); dsi != dataSets.end(); ++dsi )
		workspaces.push_back( amp_worksptr( new amp_workspace( *dsi, model, histSettings, fillDataCB, milerr, Stopped ) ) );
	}

void amp_eval_tree::initResults(){
	result = QUB_Tree::Create("AMP Result");
	result.appendClone( cfg );

	for ( int i=0; i<int(workspaces.size()); ++i ) {
		amp_workspace& wspace = * workspaces[i];
		result.appendChild( wspace.dataOut );
		}
	}

void amp_eval_tree::finishResults(){
	// calc global averages

	int nclass = workspaces[0]->classCount();

	QUB_Tree ampNode = result["amp"];
	ampNode.setNumData( QTR_TYPE_DOUBLE, nclass, 1 );
	double *amp = (double *) ampNode.data();

	QUB_Tree sdNode = result["sd"];
	sdNode.setNumData( QTR_TYPE_DOUBLE, nclass, 1 );
	double *sd = (double *) sdNode.data();

	int i, j, k;
	for ( k=0; k<nclass; ++k )
		amp[k] = sd[k] = 0.0;

	int nseg = 0;

	for ( i=0; i<int(workspaces.size()); ++i ) {
		for ( j=0; j<workspaces[i]->segCount(); ++j ) {
			double *segamp = (double *) workspaces[i]->segOuts[j]["amp"].data();
			double *segsd  = (double *) workspaces[i]->segOuts[j]["sd"].data();
			
			if ( segamp && segsd ) {
				++nseg;
				int nclassw = workspaces[i]->classCount();
				for ( k=0; k<nclassw; ++k ) {
					amp[k] += segamp[k];
					sd[k]  += segsd[k] * segsd[k];
					}
				}
			}
		}

	if ( nseg ) {
		for ( k=0; k<nclass; ++k ) {
			amp[k] /= nseg;
			sd[k]  = sqrt( sd[k] / nseg );
			}
		}
	
	// make final model with updated amp, sd
	QUB_Tree finalModel = result.appendClone( cfg["ModelFile"] );
	finalModel.setData("Final");

	ModelTree finalMTree( finalModel );
	finalMTree.condenseClasses();
	double *mdlamp = (double *) finalModel["Amps"].data();
	double *mdlsd  = (double *) finalModel["Stds"].data();
	for ( k=0; k<nclass; ++k ) {
		mdlamp[k] = amp[k];
		mdlsd[k]  = sd[k];
		}
	finalMTree.restoreClasses();
	}

bool amp_eval_tree::canRun(){
	bool can_run = true;
	
	// all of the following must be true:
	//   only one channel until we get meta-amp, -sd, -ar in mpl_metamodel
	for ( int i=0; i<int(workspaces.size()); ++i ) {
		int nchannel = (int) workspaces[i]->dataSet["ExpCond"].find("ChannelCount")->dataAsDouble( model.nchannel );
		if ( nchannel != 1 ) {
			milerr << "Can't AMP with multiple channels" << endl;
			can_run = false;
			}
		}
	
	// make sure loops in detailed balance are stimulus-friendly
	QUB_Tree loopMsg = DoCheckLoopStimuli(model.node);
	if ( loopMsg.dataCount() ) {
	  milerr << loopMsg.dataAsString();
	  can_run = false;
	}

	return can_run;
	}

QUB_Tree amp_eval_tree::execute(){
	initResults();
	
	if ( canRun() ) {
		int iSeg = 0, nSeg = 0;
		for ( int ii=0; ii<int(workspaces.size()); ++ii )
			nSeg += workspaces[ii]->segCount();
		
		QUB_Tree pctNode = QUB_Tree::Create("");
		pctNode.setData( QTR_TYPE_INT, 0 );
		int *pctData = (int *) pctNode.data();
		
		for ( int i=0; i<int(workspaces.size()); ++i ) {
			for ( int j=0; j<workspaces[i]->segCount(); ++j, ++iSeg ) {
				workspaces[i]->runAmp( autoInit, autoInitUp, doIdealize, j );
				
				*pctData = (100 * iSeg) / nSeg;
				QTR_DoCallback(percentCB, pctNode);
				
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

//----------------------------------------------------------------------------------

extern "C" QUBOPT_API QTR_Impl * amptreeiface( QTR_Impl *cfg, QTR_Impl **data,
			QTR_Callback fillDataCB, QTR_Callback rptCB, QTR_Callback pctCB ) {
	QUB_Tree cfgNode( cfg );

	std::vector<QUB_Tree> dataNodes;
	QTR_Impl **di = data;
	while ( *di )
		dataNodes.push_back( QUB_Tree( *(di++) ) );

	amp_eval_tree eval( dataNodes, cfgNode, fillDataCB, rptCB, pctCB );
	QUB_Tree resultNode = eval.execute();
	
	QTR_Impl *result = resultNode.getImpl();
	QTR_INCREF( result );
	return result;
	}

