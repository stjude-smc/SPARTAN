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
#include "mpl_eval_tree.h"
#include "qmatrixutil.h"
#include "qublib.h"

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

class mpl_workspace // can evaluate segments
                   // for one file (dataset with same conditions)
{
public:
	mpl_workspace( QUB_Tree dataset, mil_model& mdl,
					QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag );
	~mpl_workspace(); // release channel data

	void reset(); // back to original rates

	int segCount();

	// evaluate(...) is where the action is.
	// evaluate multi workspaces = sum( eachWorkspace.evaluate() )
	// call updateZ() first if z has changed
	// returns number of failed segments
	int evaluate( fq::vector<double>& gz, double& ll, int flagGrad, int seg=-1 ); // default: all segs

	bool resetZ( bool alsoResetModel=true ); // if false, you get the modified rates
	bool updateZ( fq::vector<double>& newz );
	void updateModel(); // update tree from most recent z
	void updateSD( matrix<double>& hessian ); // into the tree

	fq::vector<double> z;      // constrained parameters
	bool             zValid; // false if the initial guess didn't meet constraints

	QUB_Tree         dataSet; // filled originals
	QUB_Tree         dataOut; // results and idealization
	std::vector<QUB_Tree> segOuts;

protected:
	void setupData();
	void setupArrays();

	// data pointers
	int              maxndata;
	std::vector<int>      ndata;
	std::vector<double*>  datas;

	mpl_metamodel    model;

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

	// from evaluate:
	int              nmetastate,**metastate;
	fq::vector<double> scale;
	matrix<double>   q, qq, a, dx, dfa, alpha, beta, b;
	tensor<double>   dq, dqq, da;
	matrix<double*>  dqAsPtr, dqqAsPtr, daAsPtr;
};

typedef CountedPtr<mpl_workspace> mpl_worksptr;

class mpl_eval_tree : public dfp_optimizable
{
public:
	mpl_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
					QTR_Callback fillDatCB, QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB );

	QUB_Tree execute();

	virtual bool getStartingZ( fq::vector<double> &z );
	virtual bool evaluate( fq::vector<double> &z, fq::vector<double> &gz,
							double &ll, int &nf, int &ndf, bool flagGrad=true );
	virtual bool checkContinue( int iter, int nf, int ndf, double ll,
								fq::vector<double> &z, fq::vector<double> &gz );

protected:
	void initResults();
	void finishResults();
	void getMilHist();
	// QUB_Tree milHistToLogMS( QUB_Tree histGroup ); // returns arg

	void sendResult( QUB_Tree resNode, int iter, bool final=false, int errCode = 0 );

	bool canRun();
	bool evaluateAll( fq::vector<double> &z, fq::vector<double> &gz,
					  double &ll, int &nf, int &ndf, bool flagGrad=true );

	dfp_result runDFP();

	void blessMPLSegment( QUB_Tree seg );
	void blessMPLResult( QUB_Tree result );

	mil_model model;
	std::vector<mpl_worksptr> workspaces;
	QUB_Tree cfg;
	QTR_Callback fillDataCB, reportCB, resultCB, percentCB;

	mil_reportstream milerr;

	bool opt, together;
	int maxIter, maxRestarts;
	double maxStep, convLL, convGrad;
	
	int *Stopped;

	mpl_worksptr currSpace; // when together == false
	int currSeg, iSeg, nSeg;
	int iterations;
	
	QUB_Tree result;
	QUB_Tree pctNode;
	int *pctPtr;
};

//---------------------------------------------------------------------------------
mpl_metamodel::mpl_metamodel( mil_model& mdl, int nchan )
		: model( mdl ), 
		  nchannel( nchan ), 
		  idrug( mdl.mi.npath ) {
	model.tree.condenseClasses();
	
	mdlinf& mi = mdl.mi;
	int i, k, mc;
	
	for ( i=0; i<mi.npath; ++i )
		idrug[i] = mi.path[i][2];

	xlimit.resize( mi.x.size(), 2 );
	for ( i=0; i<mi.x.size(); ++i ) {
		xlimit[i][0] = mi.xlimit[0][i];
		xlimit[i][1] = mi.xlimit[1][i];
	}

	datinf di;
	model.tree.toDatinf(di);

	ratio.resize( model.nclass );
	for ( i=0; i<model.nclass; ++i )
		ratio[i] = di.i[i] - di.i[0];

	mutimakv( nchannel, mi.nstate, mi.npath, mi.clazz, mi.path, ratio,
		&nmutistate, &nmutipath, &nmuticlass, &mutistate, &mutipath);

	fq::vector<double> excessVar( model.nclass );
	double var0 = di.r[0][0];
	for ( i=0; i<model.nclass; ++i )
		excessVar[i] = di.r[i][0] - var0;

	if ( nchannel == 1 ) {
		amp = di.i;
		xms.resize( model.nclass );
		ar = di.r;
		nar = di.r.nc - 1;
		for ( i=0; i<model.nclass; ++i ) {
			xms[i] = sqrt( var0 + excessVar[i] );
			ar[i][0] = 1.0;
		}
		muticlazz = mi.clazz;
		mutipr = mi.pr;
	}
	else if ( mi.nstate ) {
		amp = fq::vector<double>(nmuticlass, di.i[0]);
		xms = fq::vector<double>(nmuticlass, var0);
		ar = matrix<double>(nmuticlass, 1, 1.0);
		nar = 0;

		muticlazz = fq::vector<int>(nmutistate);
		for ( k=0; k<nmutistate; ++k )
			muticlazz[k] = mutistate[k][0];

		for ( mc=0; mc<nmuticlass; ++mc ) {
			for ( k=0; mutistate[k][0] != mc; k++ )
				; // find the first mutistate in this muticlass

			for ( i=0; i<mi.nstate; ++i ) {
				amp[mc] += mutistate[k][i+1] * ratio[mi.clazz[i]];
				xms[mc] += mutistate[k][i+1] * excessVar[mi.clazz[i]];
			}
			xms[mc] = sqrt(xms[mc]);
		}
		mutipr.resize(nmutistate);
		CalcMultiPr_mpl( nchannel, mi.nstate, nmutistate, mutistate, mi.pr, mutipr );
	}

	model.tree.restoreClasses();
}

mpl_metamodel::~mpl_metamodel()
{
	free2( (char**) mutistate );
	free2( (char**) mutipath );
}

// ********************************************************************************************

mpl_workspace::mpl_workspace( QUB_Tree dataset, mil_model& mdl,
							 QTR_Callback fillDatCB, std::ostream& errstream, int *StopFlag )
: dataSet( dataset ), dataOut( QUB_Tree::Create("DataSet") ),
  model( mdl, (int) dataset["ExpCond"].find("ChannelCount")->dataAsDouble( mdl.nchannel ) ),
  fillDataCB( fillDatCB ), milerr( errstream ), Stopped( StopFlag ),
  mtx( 2 * mdl.mi.npath + 1, 2 * mdl.mi.npath + 1 ), vct( 2 * mdl.mi.npath + 1 )
{
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

	QUB_Tree firNode = * dataset.find("FIR");
	nfir = firNode.dataCount();
	fir.resize( nfir + 1 );
	fir[0] = 1.0;
	for ( int i=0; i<nfir; ++i )
		fir[i+1] = firNode.dataAs( i, (double) 0.0 );

	setupData();
	setupArrays();
}

mpl_workspace::~mpl_workspace()
{
	for ( QUB_TreeMonoIter segi = dataSet.find("Segment"); ! (*segi).isNull(); segi.nextSameName() ) {
		QUB_TreeIter mustFree = (*segi).find("MPLMustFreeData");
		if ( ! (*mustFree).isNull() ) {
			mustFree.remove();
			(*segi).find("Channel").remove();
		}
	}

	if ( metastate )
		free2((char **) metastate);
}

void mpl_workspace::reset()
{
	resetZ();
}

int mpl_workspace::segCount()
{
  return (int) ndata.size();
}

void GammaToUnblessedIdealNode( QUB_Tree seg, double **alpha, double **beta, int nstate, int ndat, QUBOPT_VAR_NOT_USED double dt, int *clazz )
{
	std::vector<int> idwt( ndat );
	std::vector<int> lasts( ndat );
	int i, j, maxj, count;
	double g, maxg;
	int at = seg.dataAs(0, 0, 0) - 1;

	// lazy two-pass: first idwt = max(alpha * beta)
	for ( i=0; i<ndat; ++i ) {
		maxj = 0;
		maxg = alpha[0][i] * beta[0][i];
		for ( j=1; j<nstate; ++j ) {
			g = alpha[j][i] * beta[j][i];
			if ( g > maxg ) {
				maxg = g;
				maxj = j;
			}
		}
		idwt[i] = clazz[maxj];
	}

	// now i is dwell index, j is data index, maxj is cur class
	i = 0;
	maxj = idwt[0];
	count = 1;
	for ( j=1; j<ndat; ++j ) {
		if ( idwt[j] != maxj ) {
			at += count;
			lasts[i] = at;

			idwt[i] = maxj;
			++i;
			maxj = idwt[j];
			count = 1;
		}
		else {
			++count;
		}
	}
	idwt[i] = maxj;
	lasts[i] = at + count;
	++i; // now == ndwt

	seg["TempDwellCount"].setData( QTR_TYPE_INT, i );
	seg["TempClasses"].setNumData( QTR_TYPE_INT, i, 1, &(idwt[0]) );
	seg["TempLasts"].setNumData( QTR_TYPE_INT, i, 1, &(lasts[0]) );
	
	QUB_Tree firstsNode = seg["TempFirsts"];
	firstsNode.setNumData( QTR_TYPE_INT, i, 1, (void*)0 );
	int *firsts = (int *) firstsNode.data();
	firsts[0] = seg.dataAs(0,0,0);
	for ( j=1; j<i; ++j )
		firsts[j] = lasts[j-1] + 1;
}

#pragma warning (disable: 4101)	// Disable unreferenced variables message 
int mpl_workspace::evaluate( fq::vector<double>& gz, double& ll, int flagGrad, int seg )
{
	int rtnVal = 0;
	
	std::vector<int> segs;
	if ( seg < 0 )
		for ( int i=0; i<int(datas.size()); ++i )
			segs.push_back( i );
	else if ( seg < int(datas.size()) )
		segs.push_back( seg );

	ll = 0.0;
	if( flagGrad )
		dzerov( gz.size(), &(gz[0]) );

	int i, j, k;
	mdlinf& mi = model.model.mi;
	int nz = gz.size();

	for ( std::vector<int>::iterator segi = segs.begin(); segi != segs.end(); ++segi ) {
		try {
			int ndat = ndata[ *segi ];
			double *data = datas[ *segi ];

			/* check the total starting pr of the states int the first dwell 
			for (sum=0, i=0; i<model.nmetagroup[idwell[0]]; i++)
				sum += model.mpr[model.metaindex[idwell[0]][i]];
			if (sum < 1.0e-5)   // discard this segment 
				continue;
			*/

			double segll = 0.0;
			fq::vector<double> seggz( gz.size(), 0.0 );

			for (i=0; i<model.nmutistate; i++) 
			 for (j=0; j<model.nmutistate; j++) 
			  if (a[i][j]>=1.0 || a[i][j]<0.0) {
				  milerr << "Invalid A Matrix (" << a[i][j] << ")" << endl;
				  ll += 1.0e10; // as in mil_eval...dunno why but it keeps it from getting stuck in lnsrch
				 return ++rtnVal;
			  }

			if ( Stopped && *Stopped )
				  return ++rtnVal;

			mvectb(ndat, data, model.amp, model.xms, model.nar, model.ar,
				   nfir, fir, nmetastate, metastate, b);

			mfward(model.nmutistate, model.nar, nfir, nmetastate, metastate,
				   ndat, model.mutipr, a, b, alpha, scale);

			if ( Stopped && *Stopped )
				  return ++rtnVal;

			mbward(model.nmutistate, model.nar, nfir, nmetastate, metastate,
			      ndat, a, b, scale, beta);

			if ( Stopped && *Stopped )
				  return ++rtnVal;

			for ( i=ndat-1; i>=0; --i )
			   if ( scale[i] <= 0.0 )
				   scale[i] = 0.001; // mlogl takes log of each...what's a replacement for negatives?

			segll = mlogl(model.nar, nfir, ndat, scale);
			ll -= segll;
			mdlogl(model.nmutistate, model.nar, nfir, nmetastate, metastate,
				   ndat, alpha, beta, b, dfa);

			for (k=0; k<nz; k++) {
				for (i=0; i<mi.nstate; i++) 
					for (j=0; j<mi.nstate; j++)
						seggz[k] -= dfa[i][j]*da[i][j][k];
			} 
			for (k=0; k<nz; ++k)
				gz[k] += seggz[k];

			UpdateResult( segOuts[*segi], segll, seggz );
			GammaToUnblessedIdealNode( segOuts[*segi], alpha, beta, model.nmutistate, ndat, dt, &(model.muticlazz[0]) );

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

	return rtnVal;
}
#pragma warning ( default : 4101)

bool mpl_workspace::resetZ( bool alsoResetModel )
{
	if ( alsoResetModel )
		model.model.reset();

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

bool mpl_workspace::updateZ( fq::vector<double>& newz )
{
	mdlinf& mi = model.model.mi;

	z = newz;

	ztox(z.size(), 2*mi.npath, mtx, vct, z, mi.x, dx);
	xtoq(mi.nstate,mi.npath,mi.path,C,V,mi.x,q,z.size(),dx,dqAsPtr); 
	qtomutiq(mi.nstate,mi.npath,mi.path,model.nmutistate,model.mutistate,model.nmutipath,
			model.mutipath,q,qq,z.size(),dqAsPtr,dqqAsPtr);
	mexp(model.nmutistate,dt,qq,a,z.size(),dqqAsPtr,daAsPtr,0); 

	return true;
}

void mpl_workspace::updateModel()
{
	model.model.updateRates();
}

void mpl_workspace::updateSD( matrix<double>& hessian )
{ // after updateZ(); writes to node
	mdlinf& mi = model.model.mi;

	int i, j;
	double da, db;
	
	// we already have x (from updateZ()), and mtx == dx_z

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

void mpl_workspace::setupData()
// replicate dataSet in dataOut, with the exception of DataSet["Segment"]["Channel"]
// get references to dataOut["Segment"]s in segOuts
// get size and pointers from dataSet["Segment"]["Channel"]s in ndata, datas, maxndata
{
	maxndata = 0;
	int segix = 0;
	int ndat;

	QUB_TreeIter dsout = dataOut.end();
	
	for ( QUB_TreeMonoIter metai = dataSet.children(); ! (*metai).isNull(); ++metai ) {
		QUB_Tree meta = *metai;
		if ( meta.name() == "Segment" ) {
			QUB_Tree seg = meta.clone(false);
			seg["Index"].setData( QTR_TYPE_INT, segix++ );
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

			if ( (*(meta.find("Channel"))).isNull() ) {
				QTR_DoCallback(fillDataCB, meta);
				meta["MPLMustFreeData"];
			}
			QUB_Tree channel = meta["Channel"];
			ndat = channel.dataCount();
			if ( ndat > maxndata )
				maxndata = ndat;
			ndata.push_back( ndat );
			datas.push_back( (double *) channel.data() );

			seg["amp"].setNumData(QTR_TYPE_DOUBLE, model.nmuticlass, 1, &(model.amp[0]));
			seg["sd"].setNumData(QTR_TYPE_DOUBLE, model.nmuticlass, 1, &(model.xms[0]));
		}
		else {
			dsout.insert( meta.clone() );
			++dsout;
		}
	}
}

void mpl_workspace::setupArrays()
{
	mdlinf& mi = model.model.mi;
	zValid = resetZ( false );

	int nx = mi.x.size();
	int nz = z.size();
	int i;

	scale.resize( maxndata );
	q.resize( mi.nstate, mi.nstate );
	qq.resize( nmetastate, nmetastate );
	a.resize( nmetastate, nmetastate );
	dx.resize( nx, nz );
	dfa.resize( nmetastate, nmetastate );
	alpha.resize( nmetastate, maxndata );
	beta.resize( nmetastate, maxndata );
	b.resize( nmetastate, maxndata );
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

// *************************************************************************

mpl_eval_tree::mpl_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
							  QTR_Callback fillDatCB, QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB )
  : model( config["ModelFile"], config["SearchLimit"].dataAsDouble(1000.0) ),
    cfg( config ), fillDataCB( fillDatCB ), reportCB( rptCB ), resultCB( rsltCB ), percentCB( pctCB ),
	milerr( rptCB )
{
	opt = config["Mode"].dataAsString() == "optimize";
	together = config["use segments"].dataAsString() == "together";
	
	QUB_Tree dfp = cfg["DFP"];
	maxIter = dfp["MaxIterations"].dataAsInt( 100 );
	maxRestarts = dfp["MaxRestarts"].dataAsInt( 0 );
	maxStep = dfp["MaxStep"].dataAsDouble( 1.0 );
	convLL = dfp["ConvLL"].dataAsDouble( 0.0001 );
	convGrad = dfp["ConvGrad"].dataAsDouble( 0.00005 );

	QUB_Tree stopFlag = * (config.find("StopFlag"));
	if ( stopFlag.isNull() )
		Stopped = 0;
	else
		Stopped = (int *) stopFlag.data();

	for ( std::vector<QUB_Tree>::iterator dsi = dataSets.begin(); dsi != dataSets.end(); ++dsi )
		workspaces.push_back( mpl_worksptr( new mpl_workspace( *dsi, model, fillDataCB, milerr, Stopped ) ) );

	pctNode = QUB_Tree::Create("");
	pctNode.setData( QTR_TYPE_INT, 0 );
	pctPtr = (int *) pctNode.data();

	nSeg = 0;
	for ( int i=0; i<int(workspaces.size()); ++i )
		nSeg += workspaces[i]->segCount();
}

void mpl_eval_tree::initResults()
{
	result = QUB_Tree::Create("MPL Result");
	result.appendClone( cfg );

	if ( together ) {
		result["LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
		result["Gradient"].setData( QTR_TYPE_DOUBLE, 0.0 );
		result["Iterations"].setData( QTR_TYPE_INT, 0 );
	}

	for ( int i=0; i<int(workspaces.size()); ++i ) {
		mpl_workspace& wspace = * workspaces[i];
		result.appendChild( wspace.dataOut );

		for ( int j=0; j<wspace.segCount(); ++j ) {
			QUB_Tree seg = wspace.segOuts[j];
			seg["LL"].setData( QTR_TYPE_DOUBLE, 0.0 );
			seg["Gradient"].setData( QTR_TYPE_DOUBLE, 0.0 );
			seg["Iterations"].setData( QTR_TYPE_INT, 0 );
		}
	}
}

void mpl_eval_tree::finishResults()
{
	for ( QUB_TreeMonoIter dsi = result.find("DataSet"); ! (*dsi).isNull(); dsi.nextSameName() ) {
		for ( QUB_TreeMonoIter segi = (*dsi).find("Segment"); ! (*segi).isNull(); segi.nextSameName() ) {
			QUB_Tree seg = *segi;
			seg.find("TempClasses").remove();
			seg.find("TempFirsts").remove();
			seg.find("TempLasts").remove();
			seg.find("TempDwellCount").remove();
			seg.find("TempLL").remove();
			seg.find("TempGradient").remove();
			seg.find("Index").remove();
		}
	}
	result.find("TempLL").remove();
	result.find("TempGradient").remove();

	getMilHist();
}

void mpl_eval_tree::getMilHist()
{
	QUB_TreeMonoIter dsi, mildsi, milsegi;
	QUB_Tree histGroup;

	std::vector<QUB_Tree> datasets;
	for ( dsi = result.find("DataSet"); ! (*dsi).isNull(); dsi.nextSameName() ) {
		QUB_Tree dataset = *dsi;
		dataset["DeadTime"].setData( QTR_TYPE_DOUBLE, 0.99 * dataset["sampling"].dataAsDouble( 1.0e-2 ) );
		datasets.push_back( dataset );

		double sampMS = (*dsi)["sampling"].dataAsDouble( 1.0e-2 );
		for ( QUB_TreeMonoIter segi = dsi->find("Segment"); ! segi->isNull(); segi.nextSameName() ) {
			int ndwt = (*segi)["DwellCount"].dataAsInt();
			int *firsts = (int *) (*segi)["Firsts"].data();
			int *lasts = (int *) (*segi)["Lasts"].data();

			QUB_Tree dursNode = (*segi)["Durations"];
			dursNode.setNumData(QTR_TYPE_FLOAT, ndwt, 1, (void *) 0 );
			float *tdwt = (float *) dursNode.data();
			for ( int i=0; i<ndwt; ++i )
				tdwt[i] = float(sampMS * (lasts[i] - firsts[i] + 1));
		}
	}

	QUB_Tree milConfig = QUB_Tree::Create("Properties");
	milConfig.appendClone( model.tree.node );
	milConfig["Mode"].setData("histograms");
	milConfig.appendClone( cfg["use segments"] );
	milConfig.appendClone( cfg["Histograms"] );
	milConfig.appendClone( cfg["DFP"] );

	try {
	  mil_eval_tree mil( datasets, milConfig, reportCB, QTR_NullCallback, QTR_NullCallback );
		QUB_Tree milResult = mil.execute();
		
		for ( dsi = result.find("DataSet"), mildsi = milResult.find("DataSet");
			  ! ( (*dsi).isNull() || (*mildsi).isNull() );
			  dsi.nextSameName(), mildsi.nextSameName() ) {
			
			histGroup = * ((*mildsi).find("HistogramGroup"));
			if ( ! histGroup.isNull() )
				(*dsi).appendChild( histGroup ); // milHistToLogMS( histGroup ) );

			QUB_TreeMonoIter segi;
			for ( segi = (*dsi).find("Segment"), milsegi = (*mildsi).find("Segment");
				  ! ( (*segi).isNull() || (*milsegi).isNull() );
				  segi.nextSameName(), milsegi.nextSameName() ) {
				
				histGroup = * ((*milsegi).find("HistogramGroup"));
				if ( ! histGroup.isNull() )
					(*segi).appendChild( histGroup ); // milHistToLogMS( histGroup ) );
			}
		}
	} catch (...) {}
}
/*
QUB_Tree mpl_eval_tree::milHistToLogMS( QUB_Tree histGroup )
{
	for ( QUB_TreeMonoIter histi = histGroup.find("Histogram"); ! (*histi).isNull(); histi.next("Histogram") ) {
		QUB_Tree hist = *histi;
		hist["XLabel"].setData( "duration [log10 ms]" );

		float *bins = (float *) hist["Bins"].data();
		for ( int i=hist["Bins"].dataCount()-1; i>=0; --i )
			bins[i] = log10( bins[i] * 1.0e3 );
	}
	return histGroup;
}
*/
void mpl_eval_tree::sendResult( QUB_Tree resNode, int iter, bool final, int errCode )
{
	if ( final ) {
		resNode["Final"];
		resNode["ErrorCode"].setData( QTR_TYPE_INT, errCode );
	}
	else {
		resNode["Iterations"].setData( QTR_TYPE_INT, iter );
	}
	QTR_DoCallback(resultCB, resNode);
}

bool mpl_eval_tree::canRun()
{
	bool can_run = true;

	std::vector<double> concs;
	
	// all of the following must be true:
	//   only one channel until we get meta-amp, -sd, -ar in mpl_metamodel
	//   (and collect the concentrations present for the next one)
	for ( int i=0; i<int(workspaces.size()); ++i ) {
		concs.push_back( workspaces[i]->dataSet["ExpCond"]["Ligand"].dataAsDouble( 1.0 ) );
		
		int nchannel = (int) workspaces[i]->dataSet["ExpCond"].find("ChannelCount")->dataAsDouble( model.nchannel );
		int nar = model.node["NAr"].dataAsInt(); // mpl requires all classes to have the same NAr anyway
		if ( (nchannel > 1) && (nar > 0) ) {
			milerr << "Can't MPL with multiple channels and correlated noise" << endl;
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

bool mpl_eval_tree::evaluateAll( fq::vector<double> &z, fq::vector<double> &gz,
								double &ll, int &nf, int &ndf, bool flagGrad )
{
	fq::vector<double> seggz( gz.size() );
	double segll;

	int i;

	for ( i=0; i<gz.size(); ++i )
		gz[i] = 0.0;
	ll = 0.0;

	for ( i=0; i<int(workspaces.size()); ++i ) {
		mpl_workspace& wspace = * workspaces[i];

		if ( ! wspace.updateZ( z ) ) {
			return false;
		}

		for ( int j=0; j<wspace.segCount(); ++j ) {
			if ( wspace.evaluate( seggz, segll, flagGrad, j ) != 0 )
				return false;

			UpdateResult( wspace.segOuts[j], segll, seggz );

			ll += segll;
			for ( int k=0; k<gz.size(); ++k )
				gz[k] += seggz[k];

			if ( Stopped && *Stopped )
				return false;
		}

		if ( Stopped && *Stopped )
			return false;
	}

	UpdateResult( result, ll, gz );

	++nf;
	if ( flagGrad )
		++ndf;

	return true;
}

dfp_result mpl_eval_tree::runDFP()
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
		case -1:
		case -2:
			break;
		case -3:
			milerr << "Exceeded Maximum Iterations" << endl;
			doRestart = true;
			break;
		default:
			milerr << "Unknown Error in MPL Core" << endl;
		}
	}

	return res;
}

QUB_Tree mpl_eval_tree::execute()
{
	initResults();
	int err = 123; // not even run

	if ( opt ) {
		if ( together ) {
			result.appendChild( model.node );
			if ( canRun() ) {
				dfp_result dfpres = runDFP();
				workspaces[0]->updateSD( dfpres.hessian );
				err = dfpres.err;
			}
			sendResult( result, -1, true, err );
		}
		else {
			bool can_run = canRun();
			iSeg = 0;
			for ( int i=0; i<int(workspaces.size()); ++i ) {
				currSpace = workspaces[i];
				for ( currSeg=0; currSeg<currSpace->segCount(); ++currSeg, ++iSeg ) {
					currSpace->reset();
					currSpace->segOuts[currSeg].appendChild( model.node );
					
					ostringstream mdlname;
					mdlname << "Final " << (currSpace->segOuts[currSeg]["Index"].dataAsInt() + 1);
					model.node.setData( mdlname.str() );

					if ( can_run ) {
						dfp_result dfpres = runDFP();
						currSpace->updateSD( dfpres.hessian );
						err = dfpres.err;
					}
					sendResult( currSpace->segOuts[currSeg], -1, true, err );

					if ( Stopped && *Stopped )
						break;
				}

				if ( Stopped && *Stopped )
					break;
			}
		}
	}
	else { // evaluate
		fq::vector<double>& z = workspaces[0]->z;
		fq::vector<double> gz( z.size() );
		double ll;
		int nf, ndf;

		if ( canRun() )
			evaluateAll( z, gz, ll, nf, ndf, true );
		blessMPLResult( result );
	}

	finishResults();
	return result;
}

bool mpl_eval_tree::getStartingZ( fq::vector<double> &z )
{
	z = workspaces[0]->z;
	return workspaces[0]->zValid;
}

bool mpl_eval_tree::evaluate( fq::vector<double> &z, fq::vector<double> &gz,
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

		UpdateResult( currSpace->segOuts[currSeg], ll, gz );
	}

	if ( rtnVal )
		rtnVal = ! (Stopped && *Stopped);

	return rtnVal;
}

bool mpl_eval_tree::checkContinue( int iter, QUBOPT_VAR_NOT_USED int nf, QUBOPT_VAR_NOT_USED int ndf, QUBOPT_VAR_NOT_USED double ll,
									QUBOPT_VAR_NOT_USED fq::vector<double> &z, QUBOPT_VAR_NOT_USED fq::vector<double> &gz )
{
	++iterations;

	workspaces[0]->updateModel();

	if ( together ) {
		blessMPLResult( result );
		sendResult( result, iterations );
	}
	else {
		blessMPLResult( currSpace->segOuts[currSeg] );
		blessMPLSegment( currSpace->segOuts[currSeg] );
		sendResult( currSpace->segOuts[currSeg], iterations );
	}

	*pctPtr = (100 * iter) / maxIter;
	if ( ! together )
		*pctPtr = (*pctPtr + 100 * iSeg) / nSeg;
	if ( *pctPtr > 100 )
		*pctPtr = 100;

	QTR_DoCallback(percentCB, pctNode);

	return ! (Stopped && *Stopped);
}

void BlessMPLChild( QUB_Tree parent, string from, string to )
{
	QUB_Tree fromNode = *(parent.find(from));
	if ( fromNode.isNull() )
		return;

	QUB_Tree toNode = parent[to];
	if ( fromNode.dataType() != toNode.dataType() || fromNode.dataCount() != toNode.dataCount() )
		toNode.setData( fromNode.dataType(), fromNode.dataSize(), fromNode.dataCount(), fromNode.data() );
	else
		memcpy( toNode.data(), fromNode.data(), fromNode.dataSize() * fromNode.dataRows() );
}

void mpl_eval_tree::blessMPLSegment( QUB_Tree seg )
{
	BlessMPLChild( seg, "TempClasses", "Classes" );
	BlessMPLChild( seg, "TempFirsts", "Firsts" );
	BlessMPLChild( seg, "TempLasts", "Lasts" );
	BlessMPLChild( seg, "TempDwellCount", "DwellCount" );
	
	if ( (*(seg.find("amp"))).isNull() ) {
		QUB_Tree amp = seg.appendClone( cfg["ModelFile"]["Amps"] );
		amp.setName("amp");
		amp.resizeData( model.nclass );

		QUB_Tree sd  = seg.appendClone( cfg["ModelFile"]["Stds"] );
		sd.setName("sd");
		sd.resizeData( model.nclass );
	}
}

void mpl_eval_tree::blessMPLResult( QUB_Tree res )
{
	BlessResult( result ); // Gradient, LL as in mil

	for (QUB_TreeMonoIter dsi = res.find("DataSet"); ! (*dsi).isNull(); dsi.nextSameName() )
		for ( QUB_TreeMonoIter segi = (*dsi).find("Segment"); ! (*segi).isNull(); segi.nextSameName() )
			blessMPLSegment( *segi );
}

// *************************************************************************

extern "C" QUBOPT_API QTR_Impl *
mpltreeiface( QTR_Impl *cfg, QTR_Impl **data,
			  QTR_Callback fillDataCB, QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB )
{
	QUB_Tree cfgNode( cfg );

	std::vector<QUB_Tree> dataNodes;
	QTR_Impl **di = data;
	while ( *di )
		dataNodes.push_back( QUB_Tree( *(di++) ) );

	mpl_eval_tree eval( dataNodes, cfgNode, fillDataCB, rptCB, rsltCB, pctCB );
	QUB_Tree resultNode = eval.execute();
	
	QTR_Impl *result = resultNode.getImpl();
	QTR_INCREF( result );
	return result;
}

