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

#include <string.h>
#include "qubopt.h"
#include "LoadDataCB_Pager.h"
using namespace std;


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


QUB_Pager::QUB_Pager(int datacount, int pagesize, int pagestart, double *Data)
: dataCount(datacount),
  pageSize(pagesize),
  pageStart(pagestart),
  data(Data)
{}
QUB_Pager::~QUB_Pager()                        {}
void QUB_Pager::nextPage()                     { throw 1; }
void QUB_Pager::loadPageContaining( QUBOPT_VAR_NOT_USED int ix )   { throw 1; }


LoadDataCB_Pager::LoadDataCB_Pager( QUB_Tree segNode, QTR_Callback loadDataCB, int pagesize, std::ostream& milerr, int samplesOffset )
	: QUB_Pager(((int *) segNode.data())[1] - ((int *) segNode.data())[0] + 1,
				pagesize, - pagesize, NULL),
	  segFirst( ((int *) segNode.data())[0] ),
	  segLast(  ((int *) segNode.data())[1] ),
	  ds( QUB_Tree::Create("DataSet") ),
	  seg( ds["Segment"] ),
	  LoadData( loadDataCB ),
	  err( milerr ),
	  offset( samplesOffset )
{

	ds.appendClone( segNode.parent()["FileName"] );
	ds.appendClone( segNode.parent()["ADChannelCount"] );
	ds.appendClone( segNode.parent()["ActiveChannel"] );

	QUB_Tree procData = * ( segNode.parent().find("ProcessData") );
	if ( ! procData.isNull() )
		ds.appendClone( procData );

	seg.setNumData( QTR_TYPE_INT, 2, 1 );
	loadBounds = (int *) seg.data();
	if ( dataCount < pageSize )
		pageSize = dataCount;
	pageStart = - pageSize;

	QUB_TreeMonoIter ix = segNode.find("Index");
	if ( ! (*ix).isNull() )
		seg.appendClone( (*ix) );
}

LoadDataCB_Pager::~LoadDataCB_Pager()
{}

void LoadDataCB_Pager::nextPage() {
	loadPageContaining(pageStart + pageSize);
}

void LoadDataCB_Pager::loadPageContaining( int ix ) { // 0 <= ix < dataCount
	//err << "loading " << ix << endl;
	if ( ix >= dataCount )
		return;
	//err << ix << " < " << dataCount << endl;

	if ( data && (ix >= loadBounds[0] - segFirst) && (ix <= loadBounds[1] - segFirst) )
		return;
	//err << ix << " < (" << loadBounds[0] << " - " << segFirst << ") or " << ix << " > " << loadBounds[1] << " - " << segFirst << ")" << endl;

	seg.find("Channel").remove();
	int rqPage = ix / pageSize;
	pageStart = pageSize * rqPage;
	
	// don't go outside seg bounds: fill before with first, after with last
	int actualLoad[2];
	actualLoad[0] = segFirst + offset + pageSize * rqPage;
	actualLoad[1] = segFirst + offset + pageSize * (rqPage + 1) - 1;
	loadBounds[0] = max(segFirst, actualLoad[0]);
	loadBounds[1] = min(segLast, actualLoad[1]);
	int slops[2];
	slops[0] = min(loadBounds[0] - actualLoad[0], pageSize);
	slops[1] = min(actualLoad[1] - loadBounds[1], pageSize);

	if ( loadBounds[1] < loadBounds[0] ) {
		QUB_Tree chan = seg["Channel"];
		chan.setNumData(QTR_TYPE_DOUBLE, pageSize, 1, (void*)0);
		data = (double *) seg["Channel"].data();
		for ( int i=0; i<pageSize; ++i )
			data[i] = 0.0;
	}
	else {
		QTR_DoCallback(LoadData, seg);
		QUB_Tree chan = seg["Channel"];
		chan.resizeData( pageSize );
		data = (double *) seg["Channel"].data();
		if ( slops[0] )
			memmove(data+slops[0], data, (pageSize-slops[0]-slops[1])*sizeof(double));
		double firstVal = (slops[0] == pageSize) ? 0.0 : data[slops[0]];
		double lastVal = (slops[1] == pageSize) ? 0.0 : data[pageSize-slops[1]-1];
		for ( int i=0; i<slops[0]; ++i )
			data[i] = firstVal;
		for ( int i=0; i<slops[1]; ++i )
			data[pageSize-i-1] = lastVal;
	}

	loadBounds[0] = actualLoad[0];
	loadBounds[1] = actualLoad[1];
}
