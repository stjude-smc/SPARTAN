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

#include <float.h>
#include <stdexcept>
#include <string.h>
#include "stat_eval_tree.h"
#include "CountedPtr.h"
#include "modelTree.h"
#include "qmatrixutil.h"

using namespace fq;

#include "LoadDataCB_Pager.h"

//#include <boost/math/special_functions/next.hpp> // compiler trouble; replacement in milutil
#define TDEAD_DELTA_FLOATS 5

#define STAT_MAX_POINTS 1000000

//----- Local Functions 
//class UniqueNamer;
//class stat_workspace;
//class stat_eval_tree;

void amppdf(int nclass, float* amp, float* xms, float* pe, int iBins, float* xh, float* yh);

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

//====================================================================================
class UniqueNamer {
	private:
		hash_map<string, int> basecount;
	public:
		string get(string base){
			int & count = basecount[base];
			if ( count++ ) {
				ostringstream ost;
				ost << base << ' ' << count;
				base = ost.str();
				}
			return base;
			}
	};

//===================================================================================
class stat_workspace { 
	// can idealize segments for one file (dataset with same conditions)
	public:
		stat_workspace( QUB_Tree dataset, int histoBinCount, QTR_Callback fillDatCB, 
						std::ostream& errstream, int *StopFlag );
		~stat_workspace(); 

		int segCount();
		int classCount(int segi=0); // == segment["amp"].dataCount()

		int runStat( int seg=-1 ); // default: all segs

		QUB_Tree         dataSet; // originals with idealization and filled data
		QUB_Tree         dataOut; // stats and histograms
		std::vector<QUB_Tree> segments;
		std::vector<QUB_Tree> segOuts;

		QUB_Tree MakeAmpHistNode( LoadDataCB_Pager& pager, int nclass, float *amp, float *xms, float *occProb );
			// int binCount, int nclass, 
			//, std::ostream& milerr );
		bool amphst_paged(LoadDataCB_Pager& pager, int nh, float* xh, float* yh );

	protected:
		void runIdlStat( int seg_ix );
		void runDataStat( int seg_ix );

		void setupData();
		void setupArrays();

		int              maxndata;
		double           dtMS;

		int              histBinCount;

		QTR_Callback     fillDataCB;
		std::ostream&    milerr;
		int *            Stopped;
	};

stat_workspace::stat_workspace( QUB_Tree dataset, int histoBinCount, QTR_Callback fillDatCB, 
								std::ostream& errstream, int *StopFlag )
		: dataSet( dataset ),
		  dataOut( QUB_Tree::Create("DataSet") ), 
		  dtMS( dataset["sampling"].dataAsDouble( 1.0e-2 ) ),
		  histBinCount( histoBinCount ),
		  fillDataCB( fillDatCB ), 
		  milerr( errstream ), 
		  Stopped( StopFlag ) {
	setupData();
	setupArrays();
	}

stat_workspace::~stat_workspace(){
	}

int stat_workspace::segCount(){
  return (int) segments.size();
	}

int stat_workspace::classCount(int segi){
	return segments.size() ? segments[segi]["amp"].dataCount() : 0 ;
	}

bool stat_workspace::amphst_paged( LoadDataCB_Pager& pager, int nh, float* xh, float* yh){
	int i, pi;
	double Imin, Imax, dataval;

	//----- Scan all data for min/max
	//todo - get this as an output from skm if available.
	pager.loadPageContaining(0);
	Imax = pager.data[0];
	Imin = pager.data[0];

	for (i=1, pi=1; i<pager.dataCount; ++i, ++pi){
		if ( pi == pager.pageSize ) {
			pager.loadPageContaining( i );
			pi = pager.indexWithinPage( i );
			}
		dataval = pager.data[pi];
		Imax = (Imax>dataval) ? Imax : dataval;
		Imin = (Imin<dataval) ? Imin : dataval;

		//-- Check for user cancelled 
		if( (i%1000 == 0) && (Stopped) && (*Stopped) ) { 
			return true;
			}
		}
	
	//----- Calculate bins midpoints xh[0..nh-1] and set count per bin to 0 
	int * aiCount = new int[nh];
	double binwidth = max((Imax - Imin) / nh, 0.01);

	for ( i=0; i<nh; i++ ) {
		xh[i] = float(Imin + binwidth/2.0 + i * binwidth);
		aiCount[i]=0;
		}

	//----- Rescan data and sum in bins
	int iBin;

	pager.loadPageContaining(0);
	for (i=0, pi=0; i<pager.dataCount; ++i, ++pi) {
		if ( pi == pager.pageSize ) {
			pager.loadPageContaining( i );
			pi = pager.indexWithinPage( i );
			}

		iBin = (int) floor( (pager.data[pi] - Imin) / binwidth );
		++aiCount[ min(nh-1,max(0,iBin)) ];

		//-- Check for user cancelled 
		if( (i%1000 == 0) && (Stopped) && (*Stopped) ) { 
			delete [] aiCount;
			return true;
			}
		}

	//----- Normalize 
	for (i=0; i<nh; i++) 
	  yh[i] = float(aiCount[i]) / float(pager.dataCount);

	delete[] aiCount;

	return true;
	}


//-----
#pragma warning (disable: 4101)
int stat_workspace::runStat( int seg ){

	int rtnVal = 0;

	std::vector<int> segs;
	if ( seg < 0 )
		for ( int i=0; i<int(segments.size()); ++i )
			segs.push_back( i );
	else if ( seg < int(segments.size()) )
		segs.push_back( seg );

	for ( std::vector<int>::iterator segi = segs.begin(); segi != segs.end(); ++segi ) {
		try {
			runIdlStat( *segi );
			runDataStat( *segi );
			}
		catch (const std::overflow_error &oe) {
#ifdef _WIN32
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
			//unsigned int fstat = _clearfp(); // just in case
			//_fpreset();
			
			milerr << "Exception: " << ee.what() << ": segment " << (*segi+1) << endl;
			++rtnVal;
			}
		catch (...) {
			milerr << "Exception in Stat Core: segment " << (*segi+1) << endl;
			++rtnVal;
			}
		}
	return rtnVal;
	}
#pragma warning ( default : 4101)

QUB_Tree stat_workspace :: MakeAmpHistNode( LoadDataCB_Pager& pager, int nclass, float *amp, float *xms, float *occProb ){
	QUB_Tree hist = QUB_Tree::Create("Histogram");
	hist["BinCount"].setData( QTR_TYPE_INT, histBinCount );
	hist["XLabel"].setData( "Amplitude" );
	hist["YLabel"].setData( "Count / Total" );

	QUB_Tree bins = hist["Bins"];
	bins.setNumData( QTR_TYPE_FLOAT, histBinCount, 1 );
	float *binData = (float *) bins.data();

	QUB_Tree bars = hist["Bars"];
	bars.setNumData( QTR_TYPE_FLOAT, histBinCount, 1 );

	fq::vector<float> pdf( (nclass + 1) * histBinCount );
	
	amphst_paged( pager, histBinCount, binData, (float *) bars.data());
	amppdf( nclass, amp, xms, occProb, histBinCount, binData, pdf );

	QUB_TreeIter histChildren = hist.end();

	// sum pdf
	QUB_Tree line = QUB_Tree::Create("Line");
	histChildren.insert( line );
	++histChildren;

	line.setNumData( QTR_TYPE_FLOAT, histBinCount, 1 );
	memcpy( line.data(), &(pdf[histBinCount * nclass]), histBinCount * sizeof(float) );
	line["Color"].setData( QTR_TYPE_INT, 0 );

	// each pdf
	for ( int cls=0; cls<nclass; ++cls ) {
		line = QUB_Tree::Create("Line");
		histChildren.insert( line );
		++histChildren;

		line.setNumData( QTR_TYPE_FLOAT, histBinCount, 1 );
		memcpy( line.data(), &(pdf[histBinCount * cls]), histBinCount * sizeof(float) );
		line["Color"].setData( QTR_TYPE_INT, cls );
		}

	return hist;
	}

void stat_workspace::runDataStat( int seg ){
	const double TM_EPS = 1.0e-3;
	int nclass = classCount(seg);
	int nchan = dataSet.find("ADChannelCount")->dataAsInt(1);
	
	QUB_Tree segIn = segments[ seg ];
	int ndwt = segIn["DwellCount"].dataAsInt();
	int *idwt = (int *) segIn["Classes"].data();
	float *tdwt = (float *) segIn["Durations"].data();
	
	QUB_Tree segout = segOuts[ seg ];
	double *ampIdl = (double *) segout["amp"].data();
	double *xmsIdl = (double *) segout["sd"].data();
	
	QUB_Tree histGroup = segout["HistogramGroup"];
	ostringstream ost;
	ost << "Segment " << (seg + 1);
	histGroup["Title"].setData( ost.str() );
	
	UniqueNamer namer;
	QUB_TreeMonoIter chanmeta = dataSet["Channels"].find("Channel");
	
	for ( int chan=0; chan<nchan; ++chan, chanmeta.nextSameName() ) {
		string channelPrefix = namer.get( chanmeta->find("Name")->dataAsString() );
		
		* (int *) dataSet["ActiveChannel"].data() = chan;
		LoadDataCB_Pager pager( segIn, fillDataCB, STAT_MAX_POINTS, milerr );
		
		std::vector<int> pointsOfClass( nclass, 0 );
		int ii, di, dpi, cls;
		double tm;
		
		QUB_Tree ampNode = segout.appendChild( channelPrefix + " Amp" );
		ampNode.setNumData(QTR_TYPE_DOUBLE, nclass, 1, (void *)0);
		double *amp = (double *) ampNode.data();
		QUB_Tree xmsNode = segout.appendChild( channelPrefix + " Std" );
		xmsNode.setNumData(QTR_TYPE_DOUBLE, nclass, 1, (void *)0);
		double *xms = (double *) xmsNode.data();
		
		for ( cls=0; cls<nclass; ++cls )
			amp[cls] = xms[cls] = 0.0;
		
		di = 0; // true index in data
		pager.loadPageContaining( di );
		dpi = pager.indexWithinPage( di );
		
		for ( ii=0; ii<ndwt; ++ii ) {
			// changed loop from for to while 7-28-2004 JB - Was crashing when dtMS = 0.2 , 
			// skipping over the +/- TM_EPS window.  see also 2nd occurrence below.
			// tm = 0.0;
			// while ( fabs(tm - tdwt[ii]) > TM_EPS ) 		// tm < tdwt[ii] && 
			for( tm=0.0; tm<tdwt[ii]-TM_EPS; tm+=dtMS ) { 
				amp[ idwt[ii] ] += pager.data[ dpi ];
				pointsOfClass[ idwt[ii] ] += 1;
				++di;
				++dpi;
				if ( dpi == pager.pageSize ) {
					pager.loadPageContaining( di );
					dpi = pager.indexWithinPage( di );
					}
				}
			}
		
		for ( cls=0; cls<nclass; ++cls )
			if ( pointsOfClass[cls] )
				amp[cls] /= pointsOfClass[cls];
			
		di = 0;
		pager.loadPageContaining( di );
		dpi = pager.indexWithinPage( di );
		for ( ii=0; ii<ndwt; ++ii ) {
			for( tm=0.0; tm<tdwt[ii]-TM_EPS; tm+=dtMS ) { 
				double diff = (pager.data[ dpi ] - amp[ idwt[ii] ]);
				xms[ idwt[ii] ] += diff * diff;
				++di;
				++dpi;
				if ( dpi == pager.pageSize ) {
					pager.loadPageContaining( di );
					dpi = pager.indexWithinPage( di );
					}
				}
			}
		
		for ( cls=0; cls<nclass; ++cls )
			if ( pointsOfClass[cls] )
				xms[cls] = sqrt( xms[cls] / pointsOfClass[cls] );
			
		double *occData = (double*) segout["occupancy"].data();
		fq::vector<float> occFloats( nclass ), ampFloats( nclass ), xmsFloats( nclass );
		for ( cls=0; cls<nclass; ++cls ) {
			occFloats[cls] = float(occData[cls]);
			ampFloats[cls] = float(ampIdl[cls]);
			xmsFloats[cls] = float(xmsIdl[cls]);
			}
		
		QUB_Tree hist = MakeAmpHistNode( pager, /*histBinCount,*/ 
										nclass,  ampFloats, xmsFloats, occFloats/*, milerr */ );
		hist["Title"].setData( channelPrefix );
		histGroup.appendChild( hist );
		}
	}

template<class T>
T median(std::vector<T>& lst) {
	sort(lst.begin(), lst.end());
	int h = ((int) lst.size()) / 2;
	if ( lst.size() == 0 )
		return 0.0;
	else if ( lst.size() % 2 )
		return lst[h];
	else
		return T((lst[h] + lst[h-1]) / 2.0);
}

void stat_workspace::runIdlStat( int seg ){
	int nclass = classCount(seg);

	QUB_Tree segIn = segments[ seg ];
	int ndwt = segIn["DwellCount"].dataAsInt();
	int *idwt = (int *) segIn["Classes"].data();
	float *tdwt = (float *) segIn["Durations"].data();

	QUB_Tree segOut = segOuts[ seg ];
	
	QUB_Tree lifNode = segOut["lifetime"];
	lifNode.setLineComment("ms");
	if ( int(lifNode.dataCount()) != nclass )
		lifNode.setNumData( QTR_TYPE_DOUBLE, nclass, 1 );
	double *lif = (double *) lifNode.data();

	QUB_Tree occNode = segOut["occupancy"];
	if ( int(occNode.dataCount()) != nclass )
		occNode.setNumData( QTR_TYPE_DOUBLE, nclass, 1 );
	double *occ = (double *) occNode.data();

	QUB_Tree neventNode = segOut["nevent"];
	if ( int(neventNode.dataCount()) != nclass )
		neventNode.setNumData( QTR_TYPE_INT, nclass, 1 );
	int *nevent = (int *) neventNode.data();
	for ( int lameMicrosoftForLoopScoping = 0; lameMicrosoftForLoopScoping<nclass; ++lameMicrosoftForLoopScoping )
		nevent[ lameMicrosoftForLoopScoping ] = 0;
	
	std::vector<double> sumDur( nclass, 0.0 );
	double sumAllDur = 0.0;

	QUB_Tree firstLatNode = segOut["first latency"];
	firstLatNode.setLineComment( "ms" );
	firstLatNode.setData( QTR_TYPE_DOUBLE, 0.0 );
	double& firstLat = * (double *) firstLatNode.data();
	bool firstLatNotFound = true;
	for ( int i=0; i<ndwt; ++i ) {
		nevent[ idwt[i] ] += 1;
		sumDur[ idwt[i] ] += tdwt[i];
		sumAllDur         += tdwt[i];

		if ( firstLatNotFound ) {
			if ( idwt[i] )
				firstLatNotFound = false;
			else
				firstLat += tdwt[i];
			}
		}

	for ( int cls=0; cls<nclass; ++cls ) {
		lif[cls] = nevent[cls] ? (sumDur[cls] / nevent[cls]) : 0.0;
		occ[cls] = (sumAllDur != 0.0) ? (sumDur[cls] / sumAllDur) : 0.0;
		}

	// compute each class's median lifetimes
	std::vector<double> medlif(nclass, 0.0);
	std::vector<double> lif1;
	for ( int cls=0; cls<nclass; ++cls ) {
		lif1.clear();
		for ( int i=0; i<ndwt; ++i )
			if ( idwt[i] == cls )
				lif1.push_back( tdwt[i] );
		medlif[cls] = median(lif1);
	}
	QUB_Tree medlifNode = segOut["median_lifetime"];
	medlifNode.setNumData(QTR_TYPE_DOUBLE, nclass, 1, (void*)&(medlif[0]));
}

void stat_workspace::setupData() {
	// replicate dataSet in dataOut, with the exception of DataSet["Segment"]["Channel"] and idealization
	// get references to dataOut["Segment"]s in segOuts
	// get size and pointers from dataSet["Segment"]["Channel"]s in ndata, datas, maxndata
	maxndata = 0;
	int ndat;

	QUB_TreeIter dsout = dataOut.end();
	
	for ( QUB_TreeMonoIter metai = dataSet.children(); ! (*metai).isNull(); ++metai ) {
		QUB_Tree meta = *metai;
		if ( meta.name() == "Segment" ) {
			segments.push_back( meta );

			ndat = meta.dataAs( 1, (int)0 ) - meta.dataAs( 0, (int)0 ) + 1;
			if ( ndat > maxndata )
				maxndata = ndat;

			QUB_Tree seg = meta.clone(false);
			segOuts.push_back( seg );
			dsout.insert( seg );
			++dsout;

			QUB_TreeIter segout = seg.end();
			for ( QUB_TreeMonoIter segmetai = meta.children(); ! (*segmetai).isNull(); ++segmetai ) {
				QUB_Tree segmeta = *segmetai;
				if ( (segmeta.name() != "Channel") &&
					 (segmeta.name() != "Classes") && (segmeta.name() != "Durations") &&
					 (segmeta.name() != "Firsts") &&  (segmeta.name() != "Lasts") &&
					 (segmeta.name() != "States") &&  (segmeta.name() != "StateDurations") ) {
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

void stat_workspace::setupArrays(){
	}

//===================================================================================
QUB_Tree ApplyDeadTime(QUB_Tree DataSet, float tdMS)
{
	double sampMS = DataSet["sampling"].dataAsDouble(0.1);
	QUB_Tree result = DataSet.clone();
	for (QUB_TreeIter segi = result.find("Segment"); ! segi->isNull(); segi.nextSameName()) {
		int ndwt = (*segi)["DwellCount"].dataAsInt();
		int *idwt = (int*) (*segi)["Classes"].data();
		float *tdwt = (float*) (*segi)["Durations"].data();

		// first, drop initial too-short events
		int i = 0;
		double t = 0.0;
		while ( i<ndwt && ((float_distance(tdwt[i], tdMS) >= TDEAD_DELTA_FLOATS) && (tdwt[i]<tdMS)) ) {
			++i;
			t += tdwt[i];
		}
		if (i == ndwt) {
			segi.remove();
			continue;
		} else if ( i > 0 ) {
			segi->dataAs(0,0,(int)0) += (int)((t+0.5*sampMS) / sampMS);
		}
		// then join subsequent too-shorts with the previous event
		int oi;
		for (oi=0; i<ndwt; ++i) {
			if ( oi>0 && idwt[i] == idwt[oi-1] ) // join same-class
				tdwt[oi-1] += tdwt[i];
			else if ( oi>0 && ((float_distance(tdwt[i], tdMS) >= TDEAD_DELTA_FLOATS) && (tdwt[i] <= tdMS)) ) // join too-short
				tdwt[oi-1] += tdwt[i];
			else {
				idwt[oi] = idwt[i];
				tdwt[oi] = tdwt[i];
				++oi;
			}
		}

		int *firsts = (int*) segi->find("Firsts")->data();
		int *lasts = (int*) segi->find("Lasts")->data();
		if ( firsts && lasts ) {
			firsts[0] = segi->dataAs(0,0,(int)0);
			for (i=0; i<oi; ++i) {
				lasts[i] = firsts[i] + (int)((tdwt[i]+0.5*sampMS) / sampMS) - 1;
				if ( (i+1) < oi )
					firsts[i+1] = lasts[i] + 1;
			}
			(*segi)["Firsts"].resizeData(oi);
			(*segi)["Lasts"].resizeData(oi);
		}

		(*segi)["DwellCount"].dataAs(0,0,(int)0) = oi;
		(*segi)["Classes"].resizeData(oi);
		(*segi)["Durations"].resizeData(oi);
	}
	return result;
}

//===================================================================================
typedef CountedPtr<stat_workspace> stat_worksptr;

class stat_eval_tree{
	public:
		stat_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
						QTR_Callback fillDatCB, QTR_Callback rptCB );

		QUB_Tree execute();

	protected:
		void initResults();
		void finishResults();

		bool canRun();

		std::vector<stat_worksptr> workspaces;
		QUB_Tree cfg;
		QTR_Callback fillDataCB, reportCB;

		mil_reportstream milerr;

		int *Stopped;

		QUB_Tree result;
	};

stat_eval_tree::stat_eval_tree( std::vector<QUB_Tree>& dataSets, QUB_Tree config,
								QTR_Callback fillDatCB, QTR_Callback rptCB  )
			  : cfg( config ), 
				fillDataCB( fillDatCB ), 
				reportCB( rptCB ),
				milerr( rptCB ) {
	QUB_Tree stopFlag = * (config.find("StopFlag"));
	if ( stopFlag.isNull() )
		Stopped = 0;
	else
		Stopped = (int *) stopFlag.data();

	int histBinCount = config["Histograms"]["BinCount"].dataAsInt(64);
	bool doDeadTime  = config["ApplyDeadTime"].dataAsInt();

	for ( std::vector<QUB_Tree>::iterator dsi = dataSets.begin(); dsi != dataSets.end(); ++dsi )
		if ( doDeadTime )
		  workspaces.push_back( stat_worksptr( new stat_workspace( ApplyDeadTime(*dsi, (float) (*dsi)["ExpCond"]["DeadTime"].dataAsDouble(1e-9)), histBinCount, fillDataCB, milerr, Stopped ) ) );
		else
			workspaces.push_back( stat_worksptr( new stat_workspace( *dsi, histBinCount, fillDataCB, milerr, Stopped ) ) );
	}

void stat_eval_tree::initResults(){
	result = QUB_Tree::Create("Stat Result");
	result.appendClone( cfg );

	for ( int i=0; i<int(workspaces.size()); ++i ) {
		stat_workspace& wspace = * workspaces[i];
		result.appendChild( wspace.dataOut );
		}
	}

void stat_eval_tree::finishResults(){
	// calc global averages

	int nclass = workspaces[0]->classCount(0);

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
				for ( k=0; k<nclass; ++k ) {
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
	}

bool stat_eval_tree::canRun(){
	bool can_run = true;
	return can_run;
	}

QUB_Tree stat_eval_tree::execute(){
	initResults();

	if ( canRun() ) {
		for ( int i=0; i<int(workspaces.size()); ++i ) {
			workspaces[i]->runStat();

			if ( Stopped && *Stopped )
				break;
			}
		finishResults();
		}

	return result;
	}

//====================================================================================

//--------------------------------------------------------------------------------
// nclass = number of classes
// amp = amplitude of each class
// xms = deviation of each class
// pe = occupancy of each class
// iBins = number of bins
// xh = median of each bin
// yh = an array of arrays, holding each class,
//      for each class, holds the expected number of points in that bin
//	   depending on the standard deviation and amplitude of that class
//	   the final array is a summation for all the other arrays
void amppdf(int nclass, float* amp, float* xms, float* pe, int iBins, float* xh, float* yh) {
	float sum=0.0, binsum;
	int iBin;

	// Sum per bin and overall
	for( iBin=0; iBin<iBins; iBin++) {
		binsum = 0.0;      
		for (int iClass=0; iClass<nclass; iClass++) {	// sum all classes for bin
			yh[iClass*iBins+iBin] = pe[iClass]*gaussian(amp[iClass],xms[iClass],xh[iBin]);   // yh[iClass][iBin]
			binsum += yh[iClass*iBins+iBin];     
			}
		sum += (yh[nclass*iBins+iBin]=binsum);
		}

	// Normalize 
	if( sum != 0.0 ) 
		for( iBin=0; iBin<iBins; iBin++)
			for( int iClass=0; iClass<=nclass; iClass++) 
				yh[iClass*iBins+iBin] /= sum; 
	}


//=================================================================================

// The main calling function From C++ for AMP
QUB_Tree skm_ultd_runStat( QUB_Tree dataSet, QUB_Tree config, QTR_Callback fillDataCB, QTR_Callback reportCB ){
	std::vector<QUB_Tree> datasets;
	datasets.push_back( dataSet );

	QUB_Tree statConfig = QUB_Tree::Create("Properties");
	statConfig.appendClone( config["Histograms"] );
	statConfig.appendClone( config["ApplyDeadTime"] );
	
	stat_eval_tree eval( datasets, statConfig, fillDataCB, reportCB );
	return eval.execute();
	}

// The main calling function From Delphi for SKM 
extern "C" QUBOPT_API QTR_Impl * stattreeiface( QTR_Impl *cfg, QTR_Impl **data,
			  QTR_Callback fillDataCB, QTR_Callback rptCB ) {
	QUB_Tree cfgNode( cfg );

	std::vector<QUB_Tree> dataNodes;
	QTR_Impl **di = data;
	while ( *di )
		dataNodes.push_back( QUB_Tree( *(di++) ) );

	stat_eval_tree eval( dataNodes, cfgNode, fillDataCB, rptCB );
	QUB_Tree resultNode = eval.execute();
	
	QTR_Impl *result = resultNode.getImpl();
	QTR_INCREF( result );
	return result;
	}

