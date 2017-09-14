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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include <string.h>
#include <float.h>

#include "modelTree.h"
using namespace fq;
/*
#include <fstream>
ofstream logstream;

static std::ofstream& logg()
{
	if ( ! logstream.is_open() )
		logstream.open("C:\\qublog.txt");
	return logstream;
}
*/

#ifndef _WIN32
  #define _isnan isnan
  #define _finite finite
#endif

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

//-------------------------------------------------------
ModelTree::ModelTree() : node( QUB_Tree::Create( "ModelFile" ) ){
	// setup basic structures
	node["States"];
	node["Rates"];
	node["Constraints"];
	node["ChannelCount"].setData( QTR_TYPE_INT, 1 );
	
	std::vector<double> amps( 10 ), stds( 10 ), ar( 10 );
	for( int i=0; i<10; ++i ) {
		amps[i] = i;
		stds[i] = 0.5 + 0.5 * i;
		ar[i] = 0.0;
		}
	std::vector<int> nar( 10, 0 );

	node["Amps"].setNumData( QTR_TYPE_DOUBLE, 10, 1, &(amps[0]) );
	node["Stds"].setNumData( QTR_TYPE_DOUBLE, 10, 1, &(stds[0]) );
	node["NAr"].setNumData( QTR_TYPE_INT, 10, 1, &(nar[0]) );

	QUB_TreeIter arIter = node["Ars"].children();
	for ( int j=0; j<10; ++j ) {
		arIter.insert("Ar");
		arIter->setNumData( QTR_TYPE_DOUBLE, 10, 1, &(ar[0]) );
		++arIter;
		/*
		QUB_Tree arNode = QUB_Tree::Create("Ar");
		arNode.setNumData( QTR_TYPE_DOUBLE, 10, 1, &(ar[0]) );
		arIter.insert( arNode );
		++arIter;
		*/
		}
	}

ModelTree::ModelTree( QUB_Tree mdlNode ) : node( mdlNode ), nparam( 0 ) {
	updateParams();
	verifyArs();
	}

ModelTree::ModelTree( const ModelTree& copy )
			: node( copy.node ), indexOfParam( copy.indexOfParam ), nparam( copy.nparam ){
	}

ModelTree::~ModelTree(){
	}

void ModelTree::updateParams(){
	for ( QUB_TreeMonoIter rate = node["Rates"].find("Rate"); ! rate->isNull(); rate.nextSameName() ) {
		int *P = (int *) (*rate)["P"].data();
		QUB_TreeMonoIter pname = rate->find("PNames")->find("PName");
		for ( int i=0; i<2; ++i, pname.nextSameName() ) {
			if ( P[i] ) {
				int& ix = indexOfParam[ pname->dataAsString() ];
				if ( ix == 0 )
					ix = 2 + nparam++;
				}
			}

		int *Q = (int *) (*rate)["Q"].data();
		QUB_TreeMonoIter qname = rate->find("QNames")->find("QName");
		for ( int ii=0; ii<2; ++ii, qname.nextSameName() ) {
			if ( Q[ii] ) {
				int& ix = indexOfParam[ qname->dataAsString() ];
				if ( ix == 0 )
					ix = 2 + nparam++;
				}
			}
		}
	}

void ModelTree::verifyArs(){
	QUB_Tree ars = node["Ars"];
	QUB_Tree firstAr = ars["Ar"];

	if ( firstAr.dataCount() < 10 ) {
		std::vector<double> zeros(10, 0.0);
		firstAr.setNumData( QTR_TYPE_DOUBLE, 10, 1, &(zeros[0]) );
		}
	
	QUB_TreeIter ari = ars.find("Ar");
	ari.next("Ar");
	for ( int i=1; i<10; ++i, ari.next("Ar") ) {
		if ( ari->isNull() )
			ari.insert( firstAr.clone() );
		}
	}

void StatesToMdlinf( QUB_Tree states, mdlinf& mi ){
	mi.clazz.clear();
	mi.pr.clear();

	for ( QUB_TreeMonoIter sti = states.find("State"); ! (*sti).isNull(); sti.nextSameName() ) {
		QUB_Tree state = *sti;
		mi.clazz.push_back( state["Class"].dataAsInt() );
		mi.pr.push_back( state["Pr"].dataAsDouble() );
		}
	mi.nstate = mi.clazz.size();
	}

// note apply search limit to mdlinf not this thing
void RatesToMdlinf( QUB_Tree rates, hash_map<string, int>& indexOfParam, mdlinf& mi ){
	mi.path.clear();
	mi.x.clear();

	for ( QUB_TreeMonoIter ri = rates.find("Rate"); ! (*ri).isNull(); ri.nextSameName() ) {
		QUB_Tree rate = *ri;

		fq::vector<int> pathRow(4);
		pathRow[0] = ( (int *) rate["States"].data() )[0];
		pathRow[1] = ( (int *) rate["States"].data() )[1];
		//pathRow[2] = ( (int *) rate["P"].data() )[0];
		if ( ( (int *) rate["P"].data() )[0] )
			pathRow[2] = indexOfParam[ rate["PNames"].child().dataAsString() ];
		else
			pathRow[2] = 0;
		if ( ( (int *) rate["Q"].data() )[0] )
			pathRow[3] = indexOfParam[ rate["QNames"].child().dataAsString() ];
		else
			pathRow[3] = 1;
		mi.path.push_back( pathRow );

		pathRow[0] = ( (int *) rate["States"].data() )[1];
		pathRow[1] = ( (int *) rate["States"].data() )[0];
		//pathRow[2] = ( (int *) rate["P"].data() )[1];
		if ( ( (int *) rate["P"].data() )[1] )
			pathRow[2] = indexOfParam[ rate["PNames"].child().sibling().dataAsString() ];
		else
			pathRow[2] = 0;
		if ( ( (int *) rate["Q"].data() )[1] )
			pathRow[3] = indexOfParam[ rate["QNames"].child().sibling().dataAsString() ];
		else
			pathRow[3] = 1;
		mi.path.push_back( pathRow );

		if ( rate["k0"].dataAsDouble() <= 0.0 )
		  mi.x.push_back( log( 1e-10 ) );
		else
		  mi.x.push_back( log( ( (double *) rate["k0"].data() )[0] ) );
		mi.x.push_back( ( (double *) rate["k1"].data() )[0] );
		if ( rate["k0"].dataAs(1, 0, (double)0.0) <= 0.0 )
		  mi.x.push_back( log( 1e-10 ) );
		else
		  mi.x.push_back( log( ( (double *) rate["k0"].data() )[1] ) );
		mi.x.push_back( ( (double *) rate["k1"].data() )[1] );

		for ( int i=0; i<4; ++i ) {
			mi.xlimit[0].push_back( 0.0 );
			mi.xlimit[1].push_back( 1.0e15 );
			}
		}
	mi.npath = mi.path.nr;
	}

int findPathIndex( matrix<int> &path, int start, int end ){
	for ( int i=0; i<path.nr; ++i )
		if ( path[i][0] == start && path[i][1] == end )
			return i;
	return -1;
	}

//-----
class rate_ref {
	public:
		double *valPtr;
		int start, end;
		bool protect;

		rate_ref( fq::vector<double> &x, matrix<int> &path, int st, int en )
				: valPtr( &(x[ 2*findPathIndex(path,st,en) ]) ),
				  start( st ), end( en ), protect( false ){
			}
		rate_ref()
				: valPtr( 0 ), start( 0 ), end( 0 ), protect( false ){
			}
		rate_ref( const rate_ref &r )
			: valPtr( r.valPtr ), start( r.start ), end( r.end ), protect( r.protect ){
			}
		rate_ref& operator= (const rate_ref &r) {
			valPtr = r.valPtr;
			start = r.start;
			end = r.end;
			protect = r.protect;
			return *this;
			}

		bool operator== (const rate_ref &r) {
			return (start == r.start) && (end == r.end);
			}

		double &val() { 
			return *valPtr; 
			}
	};


void markRateProtected( std::vector<rate_ref> &refs, int start, int end ){
	for ( std::vector<rate_ref>::iterator ref = refs.begin(); ref != refs.end(); ++ref )
		if ( ref->start == start && ref->end == end )
			ref->protect = true;
	}

void BalanceLoopsFromConstraintTree( QUB_Tree constraints, mdlinf& mi ){
	std::vector<int> loopStates, otherStates;
	int loopSize;
	std::vector<rate_ref> fwd, bwd;
	std::vector<rate_ref> alreadyBalanced; // do not modify to balance another loop
	int i;
	double fprod, bprod;
	std::vector<rate_ref>::iterator ref;

	for ( QUB_TreeMonoIter loop = constraints.find("LoopBal"); ! (*loop).isNull(); loop.nextSameName() ) {
		loopStates.resize( loop->dataCount() );
		loop->readDataInto( &(loopStates[0]), 0, (int) loopStates.size() - 1 );
		loopSize = (int) loopStates.size();
		if ( loopStates[0] == loopStates[ loopSize - 1 ] )
			--loopSize;

		for ( i=1; i<loopSize; ++i ) {
			fwd.push_back( rate_ref( mi.x, mi.path, loopStates[i-1], loopStates[i] ) );
			bwd.push_back( rate_ref( mi.x, mi.path, loopStates[loopSize-i], loopStates[loopSize-(i+1)] ) );
		}
		fwd.push_back( rate_ref( mi.x, mi.path, loopStates[loopSize-1], loopStates[0] ) );
		bwd.push_back( rate_ref( mi.x, mi.path, loopStates[0], loopStates[loopSize-1] ) );

		for ( ref=alreadyBalanced.begin(); ref != alreadyBalanced.end(); ++ref ) {
			markRateProtected( fwd, ref->start, ref->end );
			markRateProtected( bwd, ref->start, ref->end );
		}

		for ( QUB_TreeMonoIter fixs = constraints.find("FixRate"); ! (*fixs).isNull(); fixs.nextSameName() ) {
			otherStates.resize( fixs->dataCount() );
			fixs->readDataInto( &(otherStates[0]), 0, ((int) otherStates.size())-1 );
			markRateProtected( fwd, otherStates[0], otherStates[1] );
			markRateProtected( bwd, otherStates[0], otherStates[1] );
		}
		for ( QUB_TreeMonoIter scales = constraints.find("ScaleRate"); ! (*scales).isNull(); scales.nextSameName() ) {
			otherStates.resize( scales->dataCount() );
			scales->readDataInto( &(otherStates[0]), 0, ((int) otherStates.size())-1 );
			markRateProtected( fwd, otherStates[0], otherStates[1] );
			markRateProtected( fwd, otherStates[2], otherStates[3] );
			markRateProtected( bwd, otherStates[0], otherStates[1] );
			markRateProtected( bwd, otherStates[2], otherStates[3] );
		}
		//for ( QUB_TreeMonoIter imbals = constraints.find("LoopImbal"); ! (*scales).isNull(); scales.nextSameName() ) {
		// Changed 2/2005 {JB}
		for ( QUB_TreeMonoIter imbals = constraints.find("LoopImbal"); ! (*imbals).isNull(); imbals.nextSameName() ) {
			otherStates.resize( imbals->dataCount() );
			imbals->readDataInto( &(otherStates[0]), 0, ((int) otherStates.size())-1 );
			if ( otherStates[ otherStates.size()-1 ] != otherStates[0] )
				otherStates.push_back( otherStates[0] );
			for ( i=1; i<int(otherStates.size()); ++i ) {
				markRateProtected( fwd, otherStates[i-1], otherStates[i] );
				markRateProtected( fwd, otherStates[i], otherStates[i-1] );
				markRateProtected( bwd, otherStates[i-1], otherStates[i] );
				markRateProtected( bwd, otherStates[i], otherStates[i-1] );
			}
		}
		
		fprod = bprod = 0.0;
		for ( ref = fwd.begin(); ref != fwd.end(); ++ref )
			fprod += ref->val();
		for ( ref = bwd.begin(); ref != bwd.end(); ++ref )
			bprod += ref->val();

		if ( fprod != bprod ) {
			for ( ref = fwd.begin(); ref != fwd.end(); ++ref )
				if ( ! ref->protect ) {
					fprod -= ref->val();
					ref->val() = bprod - fprod;
					fprod = bprod;
					break;
				}
		}
		if ( fprod != bprod ) { // still, if all the last ones were protected
			for ( ref = bwd.begin(); ref != bwd.end(); ++ref )
				if ( ! ref->protect ) {
					bprod -= ref->val();
					ref->val() = fprod - bprod;
					bprod = fprod;
					break;
				}
		}
    
		for ( ref = fwd.begin(); ref != fwd.end(); ++ref )
			alreadyBalanced.push_back( *ref );
		for ( ref = bwd.begin(); ref != bwd.end(); ++ref )
			alreadyBalanced.push_back( *ref );

		fwd.clear();
		bwd.clear();
	}
}

// returns number of constraints (including voltage-rates of 0.0 that are auto-fixed)
int ConstraintsToMdlinf( QUB_Tree constraints, mdlinf& mi, std::vector<int>& cns_src_ix ){
	BalanceLoopsFromConstraintTree( constraints, mi );

	int ncns = 0;
	int npath = mi.path.m();
	int nx = 2 * npath;

	fq::vector<double>& x = mi.x;
	matrix<int>& path = mi.path;
	matrix<double>& mtx = mi.mtx;
	fq::vector<double>& vct = mi.vct;

	mtx.clear();
	vct.clear();
	cns_src_ix.clear();

	fq::vector<double> u( nx ), v( nx );
	int k, k1, k2, cns_ix=0;

	for ( QUB_TreeMonoIter cnsi = constraints.children(); ! (*cnsi).isNull(); ++cnsi, ++cns_ix ) {
		string cname = (*cnsi).name();
		std::vector<int> cnstates( cnsi->dataCount() );
		cnsi->readDataInto( &(cnstates[0]), 0, ((int) cnstates.size())-1 );

		if ( cname == "FixRate" ) {
			k = findPathIndex( path, cnstates[0], cnstates[1] );
			v = fq::vector<double>( nx, 0.0 );
			v[2*k] = 1.0;
			mtx.push_back( v );
			vct.push_back( x[2*k] );
			cns_src_ix.push_back(cns_ix);
			++ncns;
		}
		else if ( cname == "FixExp" ) {
			k = findPathIndex( path, cnstates[0], cnstates[1] );
			v = fq::vector<double>( nx, 0.0 );
			v[2*k+1] = 1.0;
			mtx.push_back( v );
			vct.push_back( x[2*k+1] );
			cns_src_ix.push_back(cns_ix);
			++ncns;
		}
		else if ( cname == "ScaleRate" ) {
			k1 = findPathIndex( path, cnstates[0], cnstates[1] );
			k2 = findPathIndex( path, cnstates[2], cnstates[3] );
			v = fq::vector<double>( nx, 0.0 );
			v[2*k1] = -1.0;
			v[2*k2] = 1.0;
			mtx.push_back( v );
			vct.push_back( x[2*k2] - x[2*k1] );
			cns_src_ix.push_back(cns_ix);
			++ncns;
		}
		else if ( cname == "ScaleExp" ) {
			k1 = findPathIndex( path, cnstates[0], cnstates[1] );
			k2 = findPathIndex( path, cnstates[2], cnstates[3] );
			if ( fabs(x[2*k1+1]) > 1.0e-6 && fabs(x[2*k2+1]) > 1.0e-6 ) {
				v = fq::vector<double>( nx, 0.0 );
				v[2*k1+1] = 1.0;
				v[2*k2+1] = - x[2*k1+1] / x[2*k2+1];
				mtx.push_back( v );
				vct.push_back( 0.0 );
				cns_src_ix.push_back(cns_ix);
				++ncns;
			}
		}
		else if ( cname == "Generalized" ) {
			v = fq::vector<double>( nx, 0.0 );
			QUB_Tree coeffs = (*cnsi)["Coeffs"];
			memcpy(v, coeffs.data(), nx*sizeof(double));
			mtx.push_back( v );
			vct.push_back( (*cnsi)["Value"].dataAsDouble() );
			cns_src_ix.push_back(cns_ix);
			++ncns;
		}
		else if ( cname == "LoopBal" ) {
			bool closed = cnstates[ cnstates.size() - 1 ] == cnstates[0];
			int nloop = ((int) cnstates.size()) + (closed ? 0 : 1); // # states in loop + 1 overlap
			std::vector<int> loop( nloop );
			int i;
			for ( i=0; i<nloop-1; ++i )
				loop[i] = cnstates[i];
			loop[nloop-1] = loop[0];

			bool flag = false;
			u = v = fq::vector<double>( nx, 0.0 );
			for ( int j=0; j<nloop-1; ++j ) {
				k1 = findPathIndex( path, loop[j], loop[j+1] );
				if ( fabs( x[2*k1+1]) > 1.0e-6 ) {
					flag = true;
					u[2*k1+1] = 1.0;
				}
				v[2*k1] = 1.0;
				
				k2 = findPathIndex( path, loop[j+1], loop[j] );
				if ( fabs( x[2*k2+1]) > 1.0e-6 ) {
					flag = true;
					u[2*k2+1] = -1.0;
				}
				v[2*k2] = -1.0;
				//u[2*k1+1] = 1.0;   u[2*k2+1] = -1.0;
				//v[2*k1]   = 1.0;   v[2*k2]   = -1.0;
			}
			mtx.push_back( v );
			vct.push_back( 0.0 );
			cns_src_ix.push_back(cns_ix);
			++ncns;

			if ( flag ) {
				mtx.push_back( u );
				vct.push_back( 0.0 );
				cns_src_ix.push_back(-1); // not removable by pareconstraints
				++ncns;
			}
		}
		else if ( cname == "LoopImbal" ) {
			bool closed = cnstates[ cnstates.size() - 1 ] == cnstates[0];
			int nloop = ((int) cnstates.size()) + (closed ? 0 : 1); // # states in loop + 1 overlap
			std::vector<int> loop( nloop );
			int i;
			for ( i=0; i<nloop-1; ++i )
				loop[i] = cnstates[i];
			loop[nloop-1] = loop[0];

			bool flag = false;
			u = v = fq::vector<double>( nx, 0.0 );
			double uSum = 0.0;
			double vSum = 0.0;
			for ( int j=0; j<nloop-1; ++j ) {
				k1 = findPathIndex( path, loop[j], loop[j+1] );
				if ( fabs( x[2*k1+1]) > 1.0e-6 ) {
					flag = true;
					u[2*k1+1] = 1.0;
				}
				v[2*k1] = 1.0;
				
				k2 = findPathIndex( path, loop[j+1], loop[j] );
				if ( fabs( x[2*k2+1]) > 1.0e-6 ) {
					flag = true;
					u[2*k2+1] = -1.0;
				}
				v[2*k2] = -1.0;
				//u[2*k1+1] = 1.0;   u[2*k2+1] = -1.0;
				//v[2*k1]   = 1.0;   v[2*k2]   = -1.0;
				vSum += x[2*k1] - x[2*k2];
				uSum += x[2*k1+1] - x[2*k2+1];
			}
			mtx.push_back( v );
			vct.push_back( vSum );
			cns_src_ix.push_back(cns_ix);
			++ncns;

			if ( flag ) {
				mtx.push_back( u );
				vct.push_back( 0.0 ); //  usum?
				cns_src_ix.push_back(-1); // not removable by pareconstraints
				++ncns;
			}
		}

	}

	// fix exp when voltage constant is 0
	for ( k=0; k<npath; ++k ) {
		if ( fabs( x[2*k+1] ) < 1.0e-6 ) {
			v = fq::vector<double>( nx, 0.0 );
			v[2*k+1] = 1.0;
			mtx.push_back( v );
			vct.push_back( 0.0 );
			cns_src_ix.push_back(-1); // not removable by pareconstraints
			++ncns;
		}
	}

	return ncns;
}

// returns # constraints
int ModelTree::toMdlinf( mdlinf& mi ){
	if ( (* node.find("Condensed")).isNull() )
		condenseClasses();

	StatesToMdlinf( node["States"], mi );
	RatesToMdlinf( node["Rates"], indexOfParam, mi );
	return ConstraintsToMdlinf( node["Constraints"], mi, cns_src_ix );
	}

template< class T >
void reorderArray( T *arr, int n, int *order ){	// same as in SKM_EVAL_TREE
  if ( ! arr ) return;
	std::vector<T> cpy( n );
	memcpy( &(cpy[0]), arr, n * sizeof(T) );
	for ( int i=0; i<n; ++i )
		arr[i] = cpy[ order[i] ];
	}

void reorderNodes( QUB_Tree parent, int n, int *order ){
	std::vector<QUB_Tree> nodes( n );
	QUB_TreeIter ii = parent.children();
	int i;
	for ( i=0; i<n; ++i )
		nodes[i] = ii.remove();

	for( i=0, ii = parent.children(); i<n; ++i, ++ii )
		ii.insert( nodes[ order[i] ] );
	}

void ModelTree::condenseClasses(){
	int nclass = node["Amps"].dataCount();
	double *amp = (double *) node["Amps"].data();
	double *std = (double *) node["Stds"].data();
	int    *nar = (int    *) node["NAr"].data();

	QUB_Tree ars = node["Ars"];

	int newix = 0, cls;
	std::vector<int> occupied(nclass, 0);
	QUB_TreeMonoIter state;
	int noccupied = 0;
	for ( state = node["States"].children(); !(*state).isNull(); ++state )
		occupied[ (*state)["Class"].dataAsInt() ] = 1;
	for ( cls = 0; cls < nclass; ++cls )
		noccupied += occupied[ cls ];

	QUB_Tree condensed = node["Condensed"];
	QUB_Tree toCondNode = condensed["toCondensed"];
	QUB_Tree fromCondNode = condensed["fromCondensed"];
	condensed["ClassCount"].setData( QTR_TYPE_INT, noccupied );
	
	toCondNode.setNumData( QTR_TYPE_INT, nclass, 1 );
	int *toCond = (int *) toCondNode.data();
	fromCondNode.setNumData( QTR_TYPE_INT, nclass, 1 );
	int *fromCond = (int *) fromCondNode.data();

	for ( cls=0; cls<nclass; ++cls ) {
		if ( occupied[cls] ) {
			toCond[cls] = newix;
			fromCond[newix] = cls;
			++newix;
			}
		}
	for ( cls=0; cls<nclass; ++cls ) {
		if ( ! occupied[cls] ) {
			toCond[cls] = newix;
			fromCond[newix] = cls;
			++newix;
			}
		}

	reorderArray( amp, nclass, fromCond );
	reorderArray( std, nclass, fromCond );
	reorderArray( nar, nclass, fromCond );
	reorderNodes( ars, nclass, fromCond );

	for ( state = node["States"].children(); ! (*state).isNull(); ++state ) {
		int& clz = ( (int *) (*state)["Class"].data() )[0];
		clz = toCond[ clz ];
		}
	}

void ModelTree::restoreClasses(){
	QUB_Tree condensed = node["Condensed"];
	int *toCond = (int *) condensed["toCondensed"].data();
	int *fromCond = (int *) condensed["fromCondensed"].data();

	if ( ! (toCond && fromCond) )
		return;

	int nclass = node["Amps"].dataCount();
	reorderArray( (double *) node["Amps"].data(), nclass, toCond );
	reorderArray( (double *) node["Stds"].data(), nclass, toCond );
	reorderArray( (int    *) node["NAr"].data(),  nclass, toCond );
	reorderNodes( node["Ars"], nclass, toCond );

	for ( QUB_TreeMonoIter state = node["States"].children(); ! (*state).isNull(); ++state ) {
		int& clz = ( (int *) (*state)["Class"].data() )[0];
		clz = fromCond[ clz ];
		}

	node.find("Condensed").remove();
	}

int ModelTree::classesPresent(){
	return node["Condensed"]["ClassCount"].dataAsInt();
	}

// just nchannel, amp, sd, ar
void ModelTree::toDatinf( datinf& di ){
	if ( (* node.find("Condensed")).isNull() )
		condenseClasses();

	di.nchannel = node["ChannelCount"].dataAsInt();

	di.i.clear();
	di.r.clear();
	int nclass = node["Condensed"]["ClassCount"].dataAsInt();
	double *amp = (double *) node["Amps"].data();
	double *std = (double *) node["Stds"].data();
	int    *nar = (int    *) node["NAr"].data();
	int    maxnar = 0;
	int    i;
	for ( i=0; i<nclass; ++i )
		if ( nar[i] > maxnar )
			maxnar = nar[i];

	QUB_TreeMonoIter ari = node["Ars"].children();
	for ( i=0; i<nclass; ++i, ++ari ) {
		di.i.push_back( amp[i] );
		fq::vector<double> noiseRow( maxnar+1, 0.0 );
		noiseRow[0] = std[i] * std[i];

		double *ar = (double *) (*ari).data();
		for ( int j=0; j<nar[i]; ++j )
			noiseRow[j+1] = ar[j];
		di.r.push_back( noiseRow );
		}
	}

void ModelTree::fromDatinf( datinf& di ){
	if ( (* node.find("Condensed")).isNull() )
		condenseClasses();

	* ( (int *) node["ChannelCount"].data() ) = di.nchannel;

	int nclass = node["Condensed"]["ClassCount"].dataAsInt();
	double *amp = (double *) node["Amps"].data();
	double *std = (double *) node["Stds"].data();
	int    *nar = (int    *) node["NAr"].data();

	QUB_TreeMonoIter ari = node["Ars"].children();
	for ( int i=0; i<nclass; ++i, ++ari ) {
		amp[i] = di.i[i];
		std[i] = sqrt( di.r[i][0] );

		double *ar = (double *) (*ari).data();
		for ( int j=0; j<nar[i]; ++j )
			ar[j] = di.r[i][j+1];
		}
	}

void ModelTree::updateRates( mdlinf& mi ){ 
	// assume that the rates in question came from ModelTree::toMdlinf( mi ), so they're in order
	int i = 0;
	double *k0, *k1;

	for ( QUB_TreeMonoIter rate_i = node["Rates"].children(); ! (*rate_i).isNull(); ++rate_i, i += 4 ) {
		k0 = (double *) (*rate_i)["k0"].data();
		k1 = (double *) (*rate_i)["k1"].data();

		if ( _isnan(mi.x[i+0]) || mi.x[i+0] == DBL_MAX || mi.x[i+0] > 10e20)
			k0[0] = 0.0;
		else
			k0[0] = exp( mi.x[i+0] );
		if ( _isnan(mi.x[i+1]) || mi.x[i+1] == DBL_MAX || mi.x[i+1] > 10e20)
			k1[0] = 0.0;
		else
			k1[0] = mi.x[i+1];
		if ( _isnan(mi.x[i+2]) || mi.x[i+2] == DBL_MAX || mi.x[i+2] > 10e20)
			k0[1] = 0.0;
		else
			k0[1] = exp( mi.x[i+2] );
		if ( _isnan(mi.x[i+3]) || mi.x[i+3] == DBL_MAX || mi.x[i+3] > 10e20)
			k1[1] = 0.0;
		else
			k1[1] = mi.x[i+3];
		}
	}

void ModelTree::updateRates( mdlinf& mi, fq::vector<double>& xsd ){
	int i = 0;
	double *k0, *k1;

	for ( QUB_TreeMonoIter rate_i = node["Rates"].children(); ! (*rate_i).isNull(); ++rate_i, i += 4 ) {
		k0 = (double *) (*rate_i)["k0"].data();
		k1 = (double *) (*rate_i)["k1"].data();

		if ( _isnan(mi.x[i+0]) || mi.x[i+0] == DBL_MAX || mi.x[i+0] > 10e20)
			k0[0] = 0.0;
		else
			k0[0] = exp( mi.x[i+0] );
		if ( _isnan(mi.x[i+1]) || mi.x[i+1] == DBL_MAX || mi.x[i+1] > 10e20)
			k1[0] = 0.0;
		else
			k1[0] = mi.x[i+1];
		if ( _isnan(mi.x[i+2]) || mi.x[i+2] == DBL_MAX || mi.x[i+2] > 10e20)
			k0[1] = 0.0;
		else
			k0[1] = exp( mi.x[i+2] );
		if ( _isnan(mi.x[i+3]) || mi.x[i+3] == DBL_MAX || mi.x[i+3] > 10e20)
			k1[1] = 0.0;
		else
			k1[1] = mi.x[i+3];
		
		if ( ! (*rate_i)["dk0"].data() ) {
		  (*rate_i)["dk0"].setNumData(QTR_TYPE_DOUBLE, 2, 1, (void*) 0);
		  (*rate_i)["dk1"].setNumData(QTR_TYPE_DOUBLE, 2, 1, (void*) 0);
		}
		k0 = (double *) (*rate_i)["dk0"].data();
		k1 = (double *) (*rate_i)["dk1"].data();

		if ( _isnan(xsd[i+0]) || xsd[i+0] == DBL_MAX || xsd[i+0] > 10e20)
			k0[0] = 0.0;
		else
			k0[0] = xsd[i+0];
		if ( _isnan(xsd[i+1]) || xsd[i+1] == DBL_MAX || xsd[i+1] > 10e20)
			k1[0] = 0.0;
		else
			k1[0] = xsd[i+1];
		if ( _isnan(xsd[i+2]) || xsd[i+2] == DBL_MAX || xsd[i+2] > 10e20)
			k0[1] = 0.0;
		else
			k0[1] = xsd[i+2];
		if ( _isnan(xsd[i+3]) || xsd[i+3] == DBL_MAX || xsd[i+3] > 10e20)
			k1[1] = 0.0;
		else
			k1[1] = xsd[i+3];
		}
	}


bool FindPathInTree(map<int, list<int> >& edges, std::vector<int>& path, int from, int to, set<int>& visited){
	visited.insert(from);
	for (list<int>::iterator ei = edges[from].begin(); ei != edges[from].end(); ++ei) {
		path.push_back( *ei );
		if ( *ei == to )
			return true;
		if ( (visited.find(*ei) == visited.end()) && FindPathInTree(edges, path, *ei, to, visited) )
			return true;
		path.pop_back();
		}
	visited.erase(from);
	return false;
	}

bool FindPathInTree(map<int, list<int> >& edges, std::vector<int>& path, int from, int to){
	path.push_back(from);
	set<int> visited;
	return FindPathInTree(edges, path, from, to, visited);
	}

extern "C" QUBOPT_API QTR_Impl* BalanceAllLoops( QTR_Impl *mdl ){
	using namespace boost;

	QUB_Tree outModel = QUB_Tree(mdl).clone();

	// remove existing loop constraints
	QUB_Tree constraints = outModel["Constraints"];
	QUB_TreeIter ci;
	while ( ! (ci = constraints.find("LoopBal"))->isNull() )
		ci.remove();

	// build boost graph
	typedef adjacency_list < vecS, vecS, undirectedS,
		property<vertex_distance_t, int>, property < edge_weight_t, int > > Graph;
	typedef std::pair < int, int >E;

	int num_nodes = 0;
	for ( QUB_TreeMonoIter si = outModel["States"].find("State"); ! si->isNull(); si.nextSameName() )
		++num_nodes;

	std::vector<E> edges;
	std::vector<int> weights;
	for ( QUB_TreeMonoIter ei = outModel["Rates"].find("Rate"); ! ei->isNull(); ei.nextSameName() ) {
		QUB_Tree states = * (ei->find("States"));
		edges.push_back( E(states.dataAs(0,0,0), states.dataAs(1,0,0)) );
		weights.push_back(1);
	}

	Graph g(num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g); 
	for (std::size_t j = 0; j < edges.size(); ++j) {
		graph_traits<Graph>::edge_descriptor e; bool inserted;
		tie(e, inserted) = add_edge(edges[j].first, edges[j].second, g);
		weightmap[e] = weights[j];
	}
	std::vector < graph_traits < Graph >::vertex_descriptor > p(num_vertices(g));

	// find minimum spanning tree
/*
	property_map<Graph, vertex_distance_t>::type distance = get(vertex_distance, g);
	property_map<Graph, vertex_index_t>::type indexmap = get(vertex_index, g);
	prim_minimum_spanning_tree
	(g, *vertices(g).first, &p[0], distance, weightmap, indexmap, 
		default_dijkstra_visitor());
*/
	prim_minimum_spanning_tree(g, &p[0]);

	// the loops are found by observing which edges were removed,
	// and for each of those (from a to b), listing the states along the remaining path from a to b
	map<int, list<int> > keptEdges;
	for (std::size_t i = 0; i != p.size(); ++i) {
		if ( i == p[i] )
			continue;

		keptEdges[(int)i].push_back( (int) p[(int)i] );
		keptEdges[(int) p[i]].push_back( (int)i );

		for (int j = ((int) edges.size())-1; j >= 0; --j) {
		  if ( (edges[j].first == (int)p[(int)i] && edges[j].second == (int)i)
		       || (edges[j].first == ((int)i) && edges[j].second == (int)p[(int)i]) ) {
				edges.erase( edges.begin() + j );
				break;
			}
		}
	}
	for (std::vector<E>::iterator ei=edges.begin(); ei!=edges.end(); ++ei) {
		// find path in tree from ei->first to ei->second
		std::vector<int> path;
		if ( FindPathInTree(keptEdges, path, ei->first, ei->second) ) {
			if ( path.size() > 2 ) {
				QUB_Tree constraint = constraints.appendChild("LoopBal");
				constraint.setNumData(QTR_TYPE_INT, (int) path.size(), 1, &path[0]);
			}
		}
		else {
			// complain
		}
	}

	outModel = DoPareLoopConstraints( outModel );
	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
	}


QUB_Tree DoBalanceAllLoops(QUB_Tree model){
	QTR_Impl *outModelPtr = BalanceAllLoops(model.getImpl());
	QUB_Tree outModel(outModelPtr);
	QTR_DECREF(outModelPtr);
	return outModel;
	}

template<class T>
T svdprecis() { return T(1.0e-6); }

template<>
double svdprecis<double>() { return 1.0e-12; }

template<class T>
int rank(matrix<T>& mtx)
{
	matrix<T> U(mtx.nr+1, mtx.nc+1);
	for ( int i=0; i<mtx.nr; ++i )
		for ( int j=0; j<mtx.nc; ++j )
			U[i+1][j+1] = mtx[i][j];
	fq::vector<T> W(U.nc);
	matrix<T> Vt(U.nc, U.nc);
	svdcmp(U, mtx.nr, mtx.nc, W, Vt, cerr);
	T precis = svdprecis<T>();
	int result = 0;
	for ( int i=0; i<W.n; ++i )
		if ( W[i] > precis )
			++result;
	return result;
}
/*
extern "C" QUBOPT_API QTR_Impl* PareConstraints_rank( QTR_Impl *mdl )
{
	QUB_Tree outModel = QUB_Tree(mdl).clone();

	ModelTree mtree(outModel);
	mdlinf mi;
	mtree.toMdlinf(mi);
	int ncon = mi.mtx.nr;
	int nourcon = 0;
	for (QUB_TreeMonoIter cnsmi = outModel.find("Constraints")->children();
		! cnsmi->isNull(); ++cnsmi)
		++nourcon;

	matrix<double> mtxk0only(nourcon, mi.mtx.nc/2+1); // should this be transposed? is vct appropriate as last column?
	for (int i=0; i<nourcon; ++i) {
		for (int j=0; j<mtxk0only.nc-1; j++)
			mtxk0only[i][j] = mi.mtx[i][2*j];
		mtxk0only[i][j] = mi.vct[i];
	}
	matrix<double> mtx;
	//for (int i=nourcon; i<ncon; ++i) skip the auto-generated fix-k1 constraints
	//	mtx.push_back( mi.mtx[i] );
	std::vector<int> keep;
	int lastrank = rank(mtx);
	for (int i=0; i<nourcon; ++i) {
		matrix<double> trymtx = mtx;
		trymtx.push_back( mtxk0only[i] );
		int tryrank = rank(trymtx);
		if ( tryrank > lastrank ) {
			mtx = trymtx;
			lastrank = tryrank;
			keep.push_back(i);
		}
	}


	QUB_TreeIter cnsi = outModel.find("Constraints")->children();
	int at = 0;
	for (int i=0; i<keep.size(); ++i, ++cnsi) {
		while (at++ < keep[i])
			cnsi.remove();
	}
	while ( ! cnsi->isNull() )
		cnsi.remove();

	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
	}
*/
template<class T>
std::ostream& operator<< (std::ostream& out, fq::vector<T>& v)
{
	out << '[';
	for (int i=0; i<v.n; ++i)
		out << (i?", ":"") << v[i];
	out << ']';
	return out;
}

template<class T>
std::ostream& operator<< (std::ostream& out, matrix<T>& m)
{
	out << '[';
	for (int i=0; i<m.nr; ++i)
		out << (i?"\n ":"") << m[i];
	out << ']';
	return out;
}

int svdsetupsolve(matrix<double>& A, int nr, int nc, matrix<double>& VWiUt, std::ostream& err)
// returns rank(A)
{
	const double TOL = 1.0e-10;
	matrix<double> U = A;
	fq::vector<double> W(nc+1);
	matrix<double> VWi(nc+1, nc+1);
	int rank = 0;
	svdcmp(U, nr, nc, W, VWi, err); // does VWi contain V or Vt?
	//if ( nc == 49 ) {
	//	err << "A:" << endl << A << endl;
	//	err << "U:" << endl << U << endl;
	//	err << "W:" << endl << W << endl;
	//	err << "V:" << endl << VWi << endl;
	//}
	for (int c=1; c<=nc; ++c) { // Vt * Wi
		double w = W[c];
		if ( abs(w) < TOL )
			w = 0.0;
		else {
			++rank;
			if ( w < 0.0 ) {
				for (int r=1; r<nc; ++r)
					VWi[c][r] *= -1; // fix negative sv by also negating its col in V (row in Vt)
				w = - 1.0 / w;
			}
			else
				w = 1.0 / w;
		}
		for (int r=1; r<=nc; ++r)
			VWi[r][c] *= w;
	}
	//if ( nc == 49 )
	//	err << "VWi:" << endl << VWi << endl;
	for (int m=1; m<=nc; ++m) { // VWi * Ut
		for (int n=1; n<=nr; ++n) {
			double sum = 0.0;
			for (int k=1; k<=nc; ++k)
				sum += VWi[m][k] * U[n][k];
			VWiUt[m][n] = sum;
		}
	}
	//if ( nc == 49 )
	//	err << "VWiUt:" << endl << VWiUt << endl;
	return rank;
}

bool svdsolve(matrix<double>& A, int nr, int nc, matrix<double>& VWiUt, fq::vector<double>& B, fq::vector<double>& X, QUBOPT_VAR_NOT_USED std::ostream& err)
{
	const double TOL = 1.0e-10;

	//X = VWiUt * B;
	for (int r=1; r<=nc; ++r) {
		double sum = 0.0;
		for (int k=1; k<=nr; ++k)
			sum += VWiUt[r][k] * B[k];
		X[r] = sum;
	}
	//if ( nc == 49 ) {
	//	err << "B:" << endl << B << endl;
	//	err << "X:" << endl << X << endl;
	//}
	//fq::vector<double> BB = A * X;
	//for (int i=1; i<B.n; ++i)
	//	if ( abs(B[i] - BB[i]) > TOL )
	//		return false;
	for (int r=1; r<=nr; ++r) {
		double sum = 0.0;
		for (int k=1; k<=nc; ++k)
			sum += A[r][k] * X[k];
		if ( abs(B[r] - sum) > TOL )
			return false;
	}
	return true;
}

extern "C" QUBOPT_API QTR_Impl* PareConstraints( QTR_Impl *mdl, QTR_Impl *removed )
{
	QUB_Tree outModel = QUB_Tree(mdl).clone();
	QUB_Tree remCns = QUB_Tree(removed);
	ModelTree mtree(outModel);
	mdlinf mi;
	int nboth = mtree.toMdlinf(mi);

	if ( nboth > 0 ) {
		int nk = mi.mtx.nc;
		int npad = max(nboth, nk+1);
		std::vector<int> keep(nboth, 1);
		matrix<double> mtxvTpad(npad+1, npad+1), VWiUt(npad+1, npad+1);
		fq::vector<double> nextRow(npad+1), X(npad+1);
		int nkept = 0;
		int lastrank = 0;
		for (int ic=0; ic<nboth; ++ic) {
			if ( keep[ic] ) {
				int i;
				for (i=0; i<nk; ++i)
					mtxvTpad[i+1][nkept+1] = mi.mtx[ic][i];
				mtxvTpad[i+1][nkept+1] = mi.vct[ic];
				++nkept;
				int thisrank = svdsetupsolve(mtxvTpad, nk+1, nkept, VWiUt, cerr);
				if ( thisrank == lastrank ) {
					// err << "dropped sneaky " << ic << endl;
					keep[ic] = 0;
					--nkept;
					continue;
				}
				else {
					lastrank = thisrank;
				}
				for (int jc=ic+1; jc<nboth; ++jc) {
					if ( keep[jc] ) {
						int ij;
						for (ij=0; ij<nk; ++ij)
							nextRow[ij+1] = mi.mtx[jc][ij];
						nextRow[ij+1] = mi.vct[jc];
						if ( svdsolve(mtxvTpad, nk+1, nkept, VWiUt, nextRow, X, cerr) ) {
							keep[jc] = 0;
							// err << "drop " << jc << " at " << ic << ": " << nextRow << endl;
						}
					}
				}
			}
		}

		QUB_TreeIter cnsi = outModel.find("Constraints")->children();
		for (size_t i=0; i<keep.size(); ++i) {
			if ( keep[i] )
				++cnsi;
			else {
			  QUB_Tree rem = *cnsi;
			  cnsi.remove();
			  remCns.appendChild(rem);
			}
		}
	}
	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
}

extern "C" QUBOPT_API QTR_Impl* PareConstraints_bak( QTR_Impl *mdl )
{
	QUB_Tree outModel = QUB_Tree(mdl).clone();
	ModelTree mtree(outModel);
	mdlinf mi;

	QUB_TreeMonoIter cnsmi = outModel.find("Constraints")->children();
	int nboth = 0;
	for ( ; ! cnsmi->isNull(); ++cnsmi )
		++nboth;

	if ( nboth > 0 ) {
		int nk0 = mi.mtx.nc / 2;
		int npad = max(nboth, nk0+1);
		std::vector<int> keep(nboth, 1);
		std::vector<int> remove;
		matrix<double> mtxvTpad(npad+1, npad+1), VWiUt(npad+1, npad+1);
		fq::vector<double> nextRow(npad+1), X(npad+1);
		int nkept = 0;
		int lastrank = 0;
		for (int ic=0; ic<nboth; ++ic) {
			if ( keep[ic] ) {
				int i;
				for (i=0; i<nk0; ++i)
					mtxvTpad[i+1][nkept+1] = mi.mtx[ic][2*i];
				mtxvTpad[i+1][nkept+1] = mi.vct[ic];
				++nkept;
				int thisrank = svdsetupsolve(mtxvTpad, nk0+1, nkept, VWiUt, cerr);
				if ( thisrank == lastrank ) {
					// err << "dropped sneaky " << ic << endl;
					keep[ic] = 0;
					--nkept;
					if ( mtree.cns_src_ix[ic] >= 0 )
					  remove.push_back( mtree.cns_src_ix[ic] );
					continue;
				}
				else {
					lastrank = thisrank;
				}
				for (int jc=ic+1; jc<nboth; ++jc) {
					if ( keep[jc] ) {
						for (i=0; i<nk0; ++i)
							nextRow[i+1] = mi.mtx[jc][2*i];
						nextRow[i+1] = mi.vct[jc];
						if ( svdsolve(mtxvTpad, nk0+1, nkept, VWiUt, nextRow, X, cerr) ) {
							keep[jc] = 0;
							if ( mtree.cns_src_ix[jc] >= 0 )
							  remove.push_back( mtree.cns_src_ix[jc] );
							// err << "drop " << jc << " at " << ic << ": " << nextRow << endl;
						}
					}
				}
			}
		}

		QUB_TreeIter cnsi = outModel.find("Constraints")->children();
		int icns=0;
		size_t irem=0;
		sort(remove.begin(), remove.end());
		for ( ; irem<remove.size(); ++irem ) {
		  while ( icns < remove[irem] ) {
		    ++icns;
		    ++cnsi;
		  }
		  cnsi.remove();
		  ++icns;
		}
	}
	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
}

extern "C" QUBOPT_API QTR_Impl* PareLoopConstraints( QTR_Impl *mdl )
{
  if ( true )
    return PareConstraints( mdl, NULL );
  // issues with the following:
  //   doesn't count constraints correctly (loop constraints can be 1 or 2 rows)
  //   doesn't deal with loop imbalance

	QUB_Tree outModel = QUB_Tree(mdl).clone();
	ModelTree mtree(outModel);
	mdlinf mi;
	int nboth = mtree.toMdlinf(mi);

	// assuming it's non-loops followed by loops, balanced only, no exp anything
	QUB_TreeMonoIter cnsmi = outModel.find("Constraints")->children();
	int nother = 0;
	for ( ; (! cnsmi->isNull()) && (cnsmi->name() != "LoopBal"); ++cnsmi )
		++nother;
	int nloop = 0;
	for ( ; ! cnsmi->isNull(); ++cnsmi )
		++nloop;
	// int nboth = nother + nloop;

	if ( nboth > 0 ) {
		int nk0 = mi.mtx.nc / 2;
		int npad = max(nboth, nk0+1);
		std::vector<int> keep(nboth, 1);
		matrix<double> mtxvTpad(npad+1, npad+1), VWiUt(npad+1, npad+1);
		fq::vector<double> nextRow(npad+1), X(npad+1);
		int nkept = nother; // keep all the non-loop constraints
		for (int ic=0; ic<nother-1; ++ic) {
			int i;
			for (i=0; i<nk0; ++i)
				mtxvTpad[i+1][ic+1] = mi.mtx[ic][2*i];
			mtxvTpad[i+1][ic+1] = mi.vct[ic];
		}
		//int nkept = 0;
		int lastrank = 0;
		for (int ic=/*0*/max(nother-1,0); ic<nboth; ++ic) {
			if ( keep[ic] ) {
				int i;
				for (i=0; i<nk0; ++i)
					mtxvTpad[i+1][nkept+1] = mi.mtx[ic][2*i];
				mtxvTpad[i+1][nkept+1] = mi.vct[ic];
				++nkept;
				int thisrank = svdsetupsolve(mtxvTpad, nk0+1, nkept, VWiUt, cerr);
				if ( thisrank == lastrank ) {
					// err << "dropped sneaky " << ic << endl;
					keep[ic] = 0;
					--nkept;
					continue;
				}
				else {
					lastrank = thisrank;
				}
				for (int jc=ic+1; jc<nboth; ++jc) {
					if ( keep[jc] ) {
						int ij;
						for (ij=0; ij<nk0; ++ij)
							nextRow[ij+1] = mi.mtx[jc][2*ij];
						nextRow[ij+1] = mi.vct[jc];
						if ( svdsolve(mtxvTpad, nk0+1, nkept, VWiUt, nextRow, X, cerr) ) {
							keep[jc] = 0;
							// err << "drop " << jc << " at " << ic << ": " << nextRow << endl;
						}
					}
				}
			}
		}

		QUB_TreeIter cnsi = outModel.find("Constraints")->children();
		for (int i=0; i<(int) keep.size(); ++i) {
			if ( keep[i] )
				++cnsi;
			else
				cnsi.remove();
		}
	}
	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
}

/*
extern "C" QUBOPT_API QTR_Impl* PareConstraints_sparse( QTR_Impl *mdl )
{
	const int NMAX = 8192;
	const double THRESH = 1.0e-3;
	const int NREP = 24;
	const int ITOL = 1;
	const double TOL = 1.0e-3;

	QUB_Tree outModel = QUB_Tree(mdl).clone();
	ModelTree mtree(outModel);
	mdlinf mi;
	mtree.toMdlinf(mi);

	// assuming it's non-loops followed by loops, balanced only, no exp anything
	// maybe this should reorder constraints just in case
	QUB_TreeMonoIter cnsmi = outModel.find("Constraints")->children();
	int nother = 0;
	for ( ; (! cnsmi->isNull()) && (cnsmi->name() != "LoopBal"); ++cnsmi )
		++nother;
	int nloop = 0;
	for ( ; ! cnsmi->isNull(); ++cnsmi )
		++nloop;
	int nboth = nother + nloop;

	if ( nboth > 0 ) {
		ofstream err;
		err.open("pareconstraints_log.txt");
		int nk0 = mi.mtx.nc / 2;
		int npad = max(nboth, nk0+1);
		vector<int> keep(nboth, 1);
		CSparseSolver solver(err, NMAX);
		matrix<double> mtxvTpad(npad+1, npad+1);
		fq::vector<double> nextRow(npad+1);
		int nkept = nother; // keep all the non-loop constraints
		for (int ic=0; ic<nother-1; ++ic) {
			for (int i=0; i<nk0; ++i)
				mtxvTpad[i+1][ic+1] = mi.mtx[ic][2*i];
			mtxvTpad[i+1][ic+1] = mi.vct[ic];
		}
		for (int ic=nother-1; ic<nboth; ++ic) {
			if ( keep[ic] ) {
				for (int i=0; i<nk0; ++i)
					mtxvTpad[i+1][nkept+1] = mi.mtx[ic][2*i];
				mtxvTpad[i+1][nkept+1] = mi.vct[ic];
				++nkept;
				solver.SetA(mtxvTpad, npad, THRESH);
				for (int jc=ic+1; jc<nboth; ++jc) {
					if ( keep[jc] ) {
						for (int i=0; i<nk0; ++i)
							nextRow[i+1] = mi.mtx[jc][2*i];
						nextRow[i+1] = mi.vct[jc];
						int iter;
						double tol;
						fq::vector<double> X(npad+1);
						for (int rep=0; rep<NREP; ++rep) {
							for (int i=1; i<=npad; ++i)
								X[i] = (i>nkept) ? 0.0 : (rand() * 0.1 / INT_MAX);
							iter = 0;
							solver.linbcg(nextRow, X, ITOL, TOL, npad, iter, tol);
							if ( iter < npad+1 )
								break;
						}
						if ( tol<TOL ) {
							keep[jc] = 0;
							err << "drop " << jc << " at " << ic << ": " << nextRow << endl;
						}
					}
				}
			}
		}

		QUB_TreeIter cnsi = outModel.find("Constraints")->children();
		for (int i=0; i<keep.size(); ++i) {
			if ( keep[i] )
				++cnsi;
			else
				cnsi.remove();
		}
	}
	QTR_Impl *outModelPtr = outModel.getImpl();
	QTR_INCREF(outModelPtr);
	return outModelPtr;
}
*/

QUB_Tree DoPareConstraints(QUB_Tree model, QUB_Tree removedContainer){
  QTR_Impl *outModelPtr = PareConstraints(model.getImpl(), removedContainer.getImpl());
	QUB_Tree outModel(outModelPtr);
	QTR_DECREF(outModelPtr);
	return outModel;
	}

QUB_Tree DoPareLoopConstraints(QUB_Tree model){
	QTR_Impl *outModelPtr = PareLoopConstraints(model.getImpl());
	QUB_Tree outModel(outModelPtr);
	QTR_DECREF(outModelPtr);
	return outModel;
	}


typedef pair<int, int> CheckLoopStates;

typedef struct {
  string pname, qname;
  double k1;
} CheckLoopRate;

typedef map<CheckLoopStates, CheckLoopRate> CheckLoopRates;

extern "C" QUBOPT_API QTR_Impl* CheckLoopStimuli( QTR_Impl *mdl )
{
  QUB_Tree outModel = QUB_Tree(mdl).clone();
  ModelTree mtree(outModel);
  mdlinf mi;
  mtree.toMdlinf(mi);
  
  ostringstream msg;
  
  // detailed balance requires:
  // * same # of (each) ligand in each direction
  // * more than one (of each) voltage, with either opposite signs or opposite directions

  // get the pname and qname for each rate
  CheckLoopRates rates;
  QUB_TreeMonoIter ri = outModel.find("Rates")->find("Rate");
  for ( ; ! ri->isNull(); ri.nextSameName() ) {
    int *states = (int*) (*ri)["States"].data();
    int *p = (int*) (*ri)["P"].data();
    int *q = (int*) (*ri)["Q"].data();
	double *k1 = (double*) (*ri)["k1"].data();
    QUB_Tree pname = (*ri)["PNames"].child();
    QUB_Tree qname = (*ri)["QNames"].child();
    CheckLoopRate& rate = rates[CheckLoopStates(states[0], states[1])];
    if ( p[0] )
      rate.pname = pname.dataAsString();
    if ( q[0] )
      rate.qname = qname.dataAsString();
	rate.k1 = k1[0];
    pname = pname.sibling();
    qname = qname.sibling();
    CheckLoopRate& rate2 = rates[CheckLoopStates(states[1], states[0])];
    if ( p[1] )
      rate2.pname = pname.dataAsString();
    if ( q[1] )
      rate2.qname = qname.dataAsString();
	rate2.k1 = k1[1];
  }

  QUB_TreeMonoIter cnsmi = outModel.find("Constraints")->children();
  for ( ; ! cnsmi->isNull(); ++cnsmi ) {
    if ( cnsmi->name() == "LoopBal" || cnsmi->name() == "LoopImbal" ) {
      std::vector<int> cnstates( cnsmi->dataCount() );
      cnsmi->readDataInto( &(cnstates[0]), 0, ((int) cnstates.size())-1 );
      bool closed = cnstates[ cnstates.size() - 1 ] == cnstates[0];
      int nloop = ((int) cnstates.size()) + (closed ? 0 : 1); // # states in loop + 1 overlap
      std::vector<int> loop( nloop );
	  ostringstream loopTxt;
	  for (int i=0; i<nloop-1; ++i ) {
          loop[i] = cnstates[i];
		  loopTxt << (i ? " " : "") << (loop[i]+1);
	  }
      loop[nloop-1] = loop[0];

	  map<string, int> lig, ligFwd, ligBwd, volt, voltPlus, voltMinus; // name -> count
      for (int i=0; i<nloop-1; ++i) {
	CheckLoopRate& rateF = rates[CheckLoopStates(loop[i], loop[i+1])];
	if ( rateF.pname.size() ) {
	  lig[rateF.pname] = 1;
	  ++ligFwd[rateF.pname];
	}
	if ( rateF.qname.size() ) {
		++volt[rateF.qname];
	  if ( rateF.k1 > 0 )
		  ++voltPlus[rateF.qname];
	  else if ( rateF.k1 < 0 )
		  ++voltMinus[rateF.qname];
	}
	CheckLoopRate& rateB = rates[CheckLoopStates(loop[i+1], loop[i])];
	if ( rateB.pname.size() ) {
	  lig[rateB.pname] = 1;
	  ++ligBwd[rateB.pname];
	}
	if ( rateB.qname.size() ) {
		++volt[rateB.qname];
	  if ( rateB.k1 > 0 )
		  ++voltMinus[rateB.qname];
	  else if ( rateB.k1 < 0 )
		  ++voltPlus[rateB.qname];
	}
      }
      for (map<string,int>::iterator ligi=lig.begin(); ligi!=lig.end(); ++ligi) {
	if ( ligFwd[ligi->first] != ligBwd[ligi->first] )
		msg << "In loop: " << loopTxt.str() << ": " << ligi->first << " must influence the same number of rates forward and backward." << endl;
      }
      for (map<string,int>::iterator voli=volt.begin(); voli!=volt.end(); ++voli) {
	if ( voli->second == 1 )
	  msg << "In loop: " << loopTxt.str() << ": " << voli->first << " must influence more than one rate." << endl;
	else if ( (! voltPlus[voli->first]) || (! voltMinus[voli->first]) )
		msg << "In loop: " << loopTxt.str() << ": at least one k1 must have opposite sign or direction." << endl;
      }
    }
  }
  QUB_Tree msgTree = QUB_Tree::Create("");
  msgTree.setData(msg.str());
  QTR_Impl *outPtr = msgTree.getImpl();
  QTR_INCREF(outPtr);
  return outPtr;
}


QUB_Tree DoCheckLoopStimuli(QUB_Tree model){
	QTR_Impl *outPtr = CheckLoopStimuli(model.getImpl());
	QUB_Tree out(outPtr);
	QTR_DECREF(outPtr);
	return out;
	}


typedef pair<int, int> Rate;
typedef set<Rate> RateSet;
typedef map<Rate, RateSet> RateMMap;

void GetScaledRates_AddRelations(RateMMap& direct, RateSet& related, const Rate& start){
	related.insert(start);
	for (RateSet::iterator di=direct[start].begin(); di!=direct[start].end(); ++di)
		if (related.find(*di) == related.end())
			GetScaledRates_AddRelations(direct, related, *di);
	}

QUB_Tree GetScaledRates(QUB_Tree constraints, int from, int to, bool exp){
	string cname = exp ? "ScaleExp" : "ScaleRate";

	RateMMap directScaled;
	for (QUB_TreeMonoIter ci = constraints.find(cname); ! ci->isNull(); ci.nextSameName()) {
		Rate a(ci->dataAs(0,0,0), ci->dataAs(1,0,0));
		Rate b(ci->dataAs(2,0,0), ci->dataAs(3,0,0));
		directScaled[a].insert(b);
		directScaled[b].insert(a);
		}
	RateSet related;
	Rate start(from, to);
	GetScaledRates_AddRelations(directScaled, related, start);
	related.erase(start);

	QUB_Tree result = QUB_Tree::Create("ScaledRates");
	QUB_TreeIter result_ins = result.end();
	for (RateSet::iterator ri = related.begin(); ri != related.end(); ++ri) {
		QUB_Tree rateNode = QUB_Tree::Create("Rate");
		rateNode.setNumData(QTR_TYPE_INT, 2, 1, (void*)0);
		rateNode.dataAs(0,0,0) = ri->first;
		rateNode.dataAs(1,0,0) = ri->second;
		result_ins.insert(rateNode);
		++result_ins;
		}
	return result;
	}

QUB_Tree GetScaledRates(QUB_Tree constraints, int from, int to){
	return GetScaledRates(constraints, from, to, false);
	}

QUB_Tree GetScaledExpRates(QUB_Tree constraints, int from, int to){
	return GetScaledRates(constraints, from, to, true);
	}

extern "C" QUBOPT_API QTR_Impl* _GetScaledRates(QTR_Impl* constraints, int from, int to){
	QUB_Tree result = GetScaledRates(QUB_Tree(constraints), from, to);
	QTR_Impl *resultPtr = result.getImpl();
	QTR_INCREF(resultPtr);
	return resultPtr;
	}

extern "C" QUBOPT_API QTR_Impl* _GetScaledExpRates(QTR_Impl* constraints, int from, int to){
	QUB_Tree result = GetScaledExpRates(QUB_Tree(constraints), from, to);
	QTR_Impl *resultPtr = result.getImpl();
	QTR_INCREF(resultPtr);
	return resultPtr;
	}

//=============================================================================================
// Object array - support class for Model Merge 
class CItem {
	public :
		virtual ~CItem() { 
			}
	};

class CItemArray { 
	public :
		CItem ** m_aItems;
		int m_iItems;
		CItemArray();						// Constructor 
		~CItemArray();						// Destructor
		void Allocate( int iItems ) ;		// Set size of array - does not preserve contents
		void DeleteItems( ) ;				// delete all items in array, utilizing virtual destructor
		int GetCount(){ return m_iItems; }
	};

CItemArray::CItemArray() {						// Constructor 
	m_iItems=0;
	m_aItems=NULL;
	} 

CItemArray::~CItemArray() {						// Destructor
	Allocate(0);
	}

typedef CItem * PItem ;
void CItemArray::Allocate( int iItems ) {		// Set size of array - does not preserve contents
	CItem ** pNew=NULL;							// allocate new array 
	if( iItems>0 ) {
		pNew = new PItem[iItems];
		if( pNew==NULL ) 
			iItems=0;
		}

	for( int i=0; i<iItems || i<m_iItems; ++i ) // copy or delete items
		if( i<iItems ){
			if( i<m_iItems ) 
				pNew[i]=m_aItems[i];
			else 
				pNew[i]=NULL;
			}
		else
			delete m_aItems[i];

	m_iItems=iItems;							// Clean up 
	if( m_aItems!=NULL ) 
		delete[]m_aItems;
	m_aItems=pNew;
	}

void CItemArray::DeleteItems( ) {				// delete all items in array, utilizing virtual destructor
	for( int iItem=0; iItem<m_iItems; ++iItem )
		if( m_aItems[iItem]!=NULL ) 
			delete m_aItems[iItem];
	}

//=============================================================================================
// Model Merge, unmerge, and support methods.
// Main model class for merge, CMrgMdlFmt contains arrays of CMrgMdlPath, CMrgMdlState,
//	CMrgMdlConstraint, CMrgMdlClass.
const double AMPWIDTH = 0.0001; // minimal overlap allowed
const int CONSTR_FIXA = 0;
const int CONSTR_SCALEA = 1;
const int CONSTR_LOOP = 2;
const int CONSTR_FIXB = 3;
const int CONSTR_SCALEB = 4;

class CMrgMdlConstraint : public CItem { 
	public :
		int m_iType;
		int m_iStates;
		int * m_aiStates;
		CMrgMdlConstraint(){
			m_aiStates=NULL;
			m_iStates=0;
			}
		~CMrgMdlConstraint(){
			if( m_aiStates != NULL ) 
				delete [] m_aiStates;
			}
		void Set( int iType, int iStates, int i0=0, int i1=0, int i2=0, int i3=0 ) { 
			m_iType=iType;
			m_iStates=iStates;
			m_aiStates=new int[iStates];
			if( iStates>0 ) m_aiStates[0]=i0;
			if( iStates>1 ) m_aiStates[1]=i1;
			if( iStates>2 ) m_aiStates[2]=i2;
			if( iStates>3 ) m_aiStates[3]=i3;
			}
	};

class CMrgMdlPath : public CItem {
	public : 
		int m_aiState[2];
		int m_aiPIndex[2];
		int m_aiQIndex[2];
		double m_anK0[2];
		double m_anK1[2];
		string m_acPName[2];
		string m_acQName[2];
		void Set( int iFrom, int iTo, int iFwdDrug, int iBwdDrug, int iFwdVolt, int iBwdVolt,
					double nFwdK0, double nFwdK1, double nBwdK0, double nBwdK1, 
					string cPName0, string cPName1, string cQName0, string cQName1) {
			m_aiState[0]=iFrom;
			m_aiState[1]=iTo;
			m_aiPIndex[0]=iFwdDrug;
			m_aiPIndex[1]=iBwdDrug;
			m_aiQIndex[0]=iFwdVolt;
			m_aiQIndex[1]=iBwdVolt;
			m_anK0[0]=nFwdK0;
			m_anK0[1]=nBwdK0;
			m_anK1[0]=nFwdK1;
			m_anK1[1]=nBwdK1;
			m_acPName[0]=cPName0;
			m_acPName[1]=cPName1;
			m_acQName[0]=cQName0;
			m_acQName[1]=cQName1;
			}
	};

class CMrgMdlClass : public CItem { 
	public : 
		CMrgMdlClass() { 
			m_iNAR=0;
			}
		double m_nAmp;
		double m_nDev;
		int m_iNAR;
		int m_iMerged;		// track number of states merged with this one 
		void Merge( CMrgMdlClass * pCl ) { 
            m_nAmp = (m_iMerged * m_nAmp + pCl->m_nAmp)/(m_iMerged+1);
			m_nDev = (m_iMerged * m_nDev + pCl->m_nDev)/(m_iMerged+1);
			m_iMerged++;
			}

	};

class CMrgMdlState : public CItem { 
	public : 
		int m_iClass;
		double m_nProb;
		double m_nXPos;
		double m_nYPos;
		void Set( int iClass, double nProb, double nXPos, double nYPos ) { 
			m_iClass=iClass;
			m_nProb=nProb;
			m_nXPos=nXPos;
			m_nYPos=nYPos;
			}
	};

class CMrgMdlFmt { 
	public :
		//----- Constructor / Destructor
		CMrgMdlFmt() { 
			nchannel=1;
			}

		~CMrgMdlFmt() { 
			}

		//----- Paths 
		CItemArray m_aPaths;
		CMrgMdlPath * Path( int iPath ) { 
			return (CMrgMdlPath*) m_aPaths.m_aItems[iPath];
			}
		int Paths() { 
			return m_aPaths.GetCount();
			}
		void AllocatePaths( int iPaths ) { 
			int iPrev=m_aPaths.GetCount();
			m_aPaths.Allocate(iPaths);
			for( int i=iPrev; i<iPaths; ++i ) 
				m_aPaths.m_aItems[i] = new CMrgMdlPath;
			}

		//----- Classes
		CItemArray m_aClasses;
		CMrgMdlClass * Class( int iClass) {
			return (CMrgMdlClass *) m_aClasses.m_aItems[iClass];
			}
		int Classes() { 
			return m_aClasses.GetCount();
			}
		bool Overlap( int i, int j ) {
            return fabs(Class(i)->m_nAmp-Class(j)->m_nAmp) 
				< (Class(i)->m_nDev+Class(j)->m_nDev)*AMPWIDTH;
			}
		void AllocateClasses( int iClasses ) { 
			int iPrev=m_aClasses.GetCount();
			m_aClasses.Allocate(iClasses);
			for( int i=iPrev; i<iClasses; ++i ) 
				m_aClasses.m_aItems[i] = new CMrgMdlClass;
			}

		//----- States 
		CItemArray m_aStates;
		CMrgMdlState * State( int iState) {
			return (CMrgMdlState *) m_aStates.m_aItems[iState];
			}
		int States() { 
			return m_aStates.GetCount();
			}
		void AllocateStates( int iStates ) { 
			int iPrev=m_aStates.GetCount();
			m_aStates.Allocate(iStates);
			for( int i=iPrev; i<iStates; ++i ) 
				m_aStates.m_aItems[i] = new CMrgMdlState;
			}

		//----- Constraints
		CItemArray m_aConstraints;
		CMrgMdlConstraint * Constraint( int iConstraint ) { 
			return (CMrgMdlConstraint *) m_aConstraints.m_aItems[iConstraint];
			}
		int Constraints() { 
			return m_aConstraints.GetCount();
			}
		void AllocateConstraints( int iConstraints ) { 
			int iPrev=m_aConstraints.GetCount();
			m_aConstraints.Allocate(iConstraints);
			for( int i=iPrev; i<iConstraints; ++i ) 
				m_aConstraints.m_aItems[i] = new CMrgMdlConstraint;
			}
		
		//----- Misc 
		int nchannel;		

	};

//--------------------------------------------------------------
// Combine 2 Models pMdl1 and pMdl2 into a merged model pMrgMdl.
// Mrg->States = M1->States * M2->States
// Mrg->Paths = M2->Paths * M1->States + M1->Paths * M2->States
// Constraints are original constraints plus a consraint for each copy of 
//    a path to scale it to the original path.
// Mrg->Classes <= M1->Classes - M2->Classes - 1.  New classes which overlap are 'reduced' to a single class.
void combine( CMrgMdlFmt * pMdl1, CMrgMdlFmt * pMdl2, CMrgMdlFmt * pMrgMdl, QUBOPT_VAR_NOT_USED std::ostream& quberr ){
    int j, k, m1, m2, iConstraint, iClass;

    pMrgMdl->nchannel = 1;

	//----- Build merge model class list 
	// Note : Mrg base amp is same as model 0;  other models contribute amplitude and noise 
	// relative to their base classes so their base classes should have the minimum noise.
	pMrgMdl->AllocateClasses( pMdl1->Classes() * pMdl2->Classes() );
	int * aiMap = new int[pMrgMdl->Classes()];		// map from merged full class list to reduced class list
	for( m2=0; m2<pMdl2->Classes(); m2++ )  
		for( m1=0; m1<pMdl1->Classes(); m1++ ) {
			iClass=m2*pMdl1->Classes()+m1;
			//-- Mrg(m1,m2).Amp = m1.amp + m2.amp - m2.amp(0)
			//-- Mrg(m1,m2).Dev = m1.Dev^2 + m2.Dev^2 - m2.Dev(0)^2 
			//-- Todo - consider m1.Dev^2 + (m2.Dev- m2.Dev(0))^2 
			//  For combining noise, assume all amplitudes are normally distributed.
			pMrgMdl->Class(iClass)->m_nAmp = pMdl1->Class(m1)->m_nAmp + pMdl2->Class(m2)->m_nAmp -pMdl2->Class(0)->m_nAmp;
			pMrgMdl->Class(iClass)->m_nDev = sqrt( pow(pMdl1->Class(m1)->m_nDev,2) + pow(pMdl2->Class(m2)->m_nDev,2) - pow(pMdl2->Class(0)->m_nDev,2) );
			aiMap[iClass] = iClass;

			//  Use aiMap to collect overlapping classes.  Maximize distinct classes... 
			for( j=0; j<iClass; j++ ) {
				if( pMrgMdl->Overlap(iClass,j) ) {
					aiMap[iClass] = j;							
					// if iClass overlaps j, and j overlaps some k=aiMap[j],  
					// then set iClass to k if it overlaps k also, or revert to iClass 
					if( aiMap[j] < j )
						aiMap[iClass] = ( pMrgMdl->Overlap(iClass,aiMap[j]) ? aiMap[j] : iClass );
					//  if iClass overlaps j, and some class k < iClass overlaps j 
					//  then revert i to itself if it does not overlap k 
					for( k=j+1; k<iClass; k++ ) {
						if( aiMap[k]==j && ! pMrgMdl->Overlap(iClass,k) ) {
							aiMap[iClass] = iClass;
							break;
							}
						}
					}
				}
			}

	//-- Reduce full class list and update aiMap;  
	int iDelta=0;							// accumulated removed class count
	for( iClass=1; iClass<pMrgMdl->Classes(); iClass++ ){
		if( aiMap[iClass]!=iClass ){
			aiMap[iClass]=aiMap[aiMap[iClass]];	// if this was mapped to another class, find where that class has gone.
			++iDelta;
			}
		else {
			aiMap[iClass]-=iDelta;
			pMrgMdl->Class(aiMap[iClass])->m_iMerged=0;
			}
		pMrgMdl->Class(aiMap[iClass])->Merge( pMrgMdl->Class(iClass) );
		}
	if( iDelta > 0 ) 
		pMrgMdl->AllocateClasses( pMrgMdl->m_aClasses.GetCount()-iDelta );

	//----- Build states list : prob=prob1*prob2; xpos=evenly spaced m1 states; ypos=evenly spaced m2 states
	int iState=0;
	pMrgMdl->AllocateStates( pMdl1->States()*pMdl2->States() );
	for( m2=0; m2<pMdl2->States(); ++m2 ) 
		for( m1=0; m1<pMdl1->States(); ++m1 ) 
			pMrgMdl->State(iState++)->Set(
						aiMap[pMdl2->State(m2)->m_iClass*pMdl1->Classes()+pMdl1->State(m1)->m_iClass], 
						pMdl1->State(m1)->m_nProb * pMdl2->State(m2)->m_nProb, 
						100*(m1+1)/(pMdl1->States()+1), 
						100*(m2+1)/(pMdl2->States()+1) );
	delete [] aiMap;

	//----- Build paths list 
	int iPath=0;
	pMrgMdl->AllocatePaths( pMdl1->Paths()*pMdl2->States() + pMdl2->Paths()*pMdl1->States() );
	for( m2=0; m2<pMdl2->States(); m2++ )  
		for( m1=0; m1<pMdl1->Paths(); m1++ )
			pMrgMdl->Path(iPath++)->Set(
						m2*pMdl1->States()+pMdl1->Path(m1)->m_aiState[0], 
						m2*pMdl1->States()+pMdl1->Path(m1)->m_aiState[1],
						pMdl1->Path(m1)->m_aiPIndex[0], 
						pMdl1->Path(m1)->m_aiPIndex[1], 
						pMdl1->Path(m1)->m_aiQIndex[0], 
						pMdl1->Path(m1)->m_aiQIndex[1], 
						pMdl1->Path(m1)->m_anK0[0], 
						pMdl1->Path(m1)->m_anK1[0], 
						pMdl1->Path(m1)->m_anK0[1], 
						pMdl1->Path(m1)->m_anK1[1], 
						pMdl1->Path(m1)->m_acPName[0], 
						pMdl1->Path(m1)->m_acPName[1],
						pMdl1->Path(m1)->m_acQName[0],
						pMdl1->Path(m1)->m_acQName[1]);

	for( m1=0; m1<pMdl1->States(); m1++ )  
		for( m2=0; m2<pMdl2->Paths(); m2++ )
			pMrgMdl->Path(iPath++)->Set(
						pMdl2->Path(m2)->m_aiState[0]*pMdl1->States()+m1,
						pMdl2->Path(m2)->m_aiState[1]*pMdl1->States()+m1,
						pMdl2->Path(m2)->m_aiPIndex[0],
						pMdl2->Path(m2)->m_aiPIndex[1],
						pMdl2->Path(m2)->m_aiQIndex[0], 
						pMdl2->Path(m2)->m_aiQIndex[1], 
						pMdl2->Path(m2)->m_anK0[0],
						pMdl2->Path(m2)->m_anK1[0],
						pMdl2->Path(m2)->m_anK0[1],
						pMdl2->Path(m2)->m_anK1[1],
						pMdl2->Path(m2)->m_acPName[0], 
						pMdl2->Path(m2)->m_acPName[1],
						pMdl2->Path(m2)->m_acQName[0],
						pMdl2->Path(m2)->m_acQName[1]);
    
	//----- Build constraints list 
    //  For the direct product constraints on the paths, use the
    //  same indexing scheme as for the direct product paths.
    //  Put the second factor internal constraints just before
    //  the corresponding direct product constraints.
	int iFrom, iTo, iOffset;
	CMrgMdlConstraint * pMrgCon, * pSrcCon;
	pMrgMdl->AllocateConstraints( pMdl1->Constraints() + pMdl2->Constraints() + 2*(pMrgMdl->Paths()-pMdl1->Paths()-pMdl2->Paths()) );
	iConstraint=0;

	//-- copy base constraints from Mdl1 to MrgMdl
	for( int iCon=0; iCon<pMdl1->Constraints(); iCon++ ){
		pSrcCon = pMdl1->Constraint(iCon);
		pMrgCon = pMrgMdl->Constraint(iConstraint++);
		pMrgCon->Set( pSrcCon->m_iType, pSrcCon->m_iStates );
		for( j=0; j<pMdl1->Constraint(iCon)->m_iStates; ++j ) 
			pMrgCon->m_aiStates[j]=pSrcCon->m_aiStates[j];
       }
	iConstraint=pMdl1->Constraints(); // redundant, but safe.

	//-- add constraints to scale paths for each m2.[1..n] copy of M1 with m2.[0] path
	for( m2=1; m2<pMdl2->States(); m2++)  
		for( m1=0; m1<pMdl1->Paths(); m1++ ) {
			iFrom = pMdl1->Path(m1)->m_aiState[0];
			iTo = pMdl1->Path(m1)->m_aiState[1];
			iOffset = m2*pMdl1->States();
			pMrgMdl->Constraint(iConstraint++)->Set( CONSTR_SCALEA, 4, iFrom, iTo, iOffset+iFrom, iOffset+iTo );
			pMrgMdl->Constraint(iConstraint++)->Set( CONSTR_SCALEA, 4, iTo, iFrom, iOffset+iTo, iOffset+iFrom );
			}

	//-- copy base constraints from mdl2 to MrgMdl
	for( int iCon=0; iCon<pMdl2->Constraints(); iCon++ ){
		pSrcCon = pMdl2->Constraint(iCon);
		pMrgCon = pMrgMdl->Constraint(iConstraint++);
		pMrgCon->Set( pSrcCon->m_iType, pSrcCon->m_iStates );
		for( j=0; j<pMdl2->Constraint(iCon)->m_iStates; ++j ) 
			pMrgCon->m_aiStates[j]=pSrcCon->m_aiStates[j];
       }

	//-- add constraints to scale paths for each m1.[1..n] copy of M2 with m1.[0] path 
	for( m1=1; m1<pMdl1->States(); m1++)  
		for( m2=0; m2<pMdl2->Paths(); m2++ ) {  
			iFrom = pMdl2->Path(m2)->m_aiState[0];
			iTo = pMdl2->Path(m2)->m_aiState[1];
			pMrgMdl->Constraint(iConstraint++)->Set( CONSTR_SCALEA, 4, iFrom*pMdl1->States(),iTo*pMdl1->States(),iFrom*pMdl1->States()+m1,iTo*pMdl1->States()+m1);
			pMrgMdl->Constraint(iConstraint++)->Set( CONSTR_SCALEA, 4, iTo*pMdl1->States(),iFrom*pMdl1->States(),iTo*pMdl1->States()+m1,iFrom*pMdl1->States()+m1);
			}

	}

//----- Convert QUB_Tree format to CMrgMdlFmt
bool ModelTreeToMdlfinf( ModelTree mt, CMrgMdlFmt * pmif, QUBOPT_VAR_NOT_USED std::ostream& quberr  ){
	if ( (* mt.node.find("Condensed")).isNull() )
		mt.condenseClasses();

	//-- Count States
	QUB_Tree qtStates = mt.node["States"];
	QUB_TreeMonoIter sti;
	int iStates=0; 
	for( sti=qtStates.find("State"); ! (*sti).isNull(); sti.nextSameName() )
		++iStates;
	//-- Copy states 
	pmif->AllocateStates( iStates );
	iStates=0;
	for( sti = qtStates.find("State"); ! (*sti).isNull(); sti.nextSameName() ) {
		QUB_Tree state = *sti;
		int ii = state["Class"].dataAsInt();
		CMrgMdlState * ps = pmif->State(iStates);
		ps->m_iClass = ii;
		pmif->State(iStates)->m_nProb = state["Pr"].dataAsDouble();
		++iStates;
		}

	//-- Count Rates
	QUB_Tree qtRates = mt.node["Rates"];
	QUB_TreeMonoIter ri;
	int iPaths=0;
	for( ri=qtRates.find("Rate"); ! (*ri).isNull(); ri.nextSameName() )
		++iPaths;
	//-- Copy Rates
	pmif->AllocatePaths(iPaths);
	iPaths=0;
	for( ri=qtRates.find("Rate"); ! (*ri).isNull(); ri.nextSameName() ) {
		QUB_Tree rate = *ri;
		pmif->Path(iPaths++)->Set(
				((int *) rate["States"].data())[0], 
				((int *) rate["States"].data())[1], 
				((int *) rate["P"].data())[0] ? mt.indexOfParam[ rate["PNames"].child().dataAsString() ] : 0, 
				((int *) rate["P"].data())[1] ? mt.indexOfParam[ rate["PNames"].child().sibling().dataAsString() ] : 0, 
				((int *) rate["Q"].data())[0] ? mt.indexOfParam[ rate["QNames"].child().dataAsString() ] : 0, 
				((int *) rate["Q"].data())[1] ? mt.indexOfParam[ rate["QNames"].child().sibling().dataAsString() ] : 0, 
				((double *) rate["k0"].data())[0],
				((double *) rate["k1"].data())[0],
				((double *) rate["k0"].data())[1],
				((double *) rate["k1"].data())[1], 
				rate["PNames"].child().dataAsString(), 
				rate["PNames"].child().sibling().dataAsString(),
				rate["QNames"].child().dataAsString(),
				rate["QNames"].child().sibling().dataAsString()
				);
		}

	//----- Count Constraints
	QUB_Tree qtConstraints = mt.node["Constraints"];
	QUB_TreeMonoIter cnsi;
	int iConstraints=0;				// count constraints
	for( cnsi = qtConstraints.children(); ! (*cnsi).isNull(); ++cnsi ) 
		++iConstraints;
	//----- Copy Constraints
	pmif->AllocateConstraints(iConstraints);		// allocate constraints
	iConstraints=0;
	for( cnsi = qtConstraints.children(); ! (*cnsi).isNull(); ++cnsi ) {
		CMrgMdlConstraint * pConst = pmif->Constraint(iConstraints++);	// ptr to next constraint 
		string cname = (*cnsi).name();									// string rep of constraint type in node
		int iStatesC = cnsi->dataCount();								// # of states in list for constraint

		//----- setup constraint object 
		if ( cname.find("FixRate") != string::npos && iStatesC==2 )
			pConst->Set(CONSTR_FIXA, 2);
		else if ( cname.find("FixExp") != string::npos && iStatesC==2 )
			pConst->Set(CONSTR_FIXB, 2);
		else if ( cname.find("ScaleRate") != string::npos && iStatesC==4 )
			pConst->Set(CONSTR_SCALEA, 4);
		else if ( cname.find("ScaleExp") != string::npos && iStatesC==4 )
			pConst->Set(CONSTR_SCALEB, 4);
		else if ( cname.find("LoopBal") != string::npos && iStatesC>2 )
			pConst->Set(CONSTR_LOOP, iStatesC );
		else
			--iConstraints;	// unsupported constraint; backtrack.
		
		//----- Copy list of states 
		if( pConst->m_iStates > 0 ) { // i.e. supported constraint type
			cnsi->readDataInto( pConst->m_aiStates, 0, pConst->m_iStates-1 );
			if( pConst->m_iType==CONSTR_LOOP && pConst->m_aiStates[pConst->m_iStates-1]==pConst->m_aiStates[0] ) 
				pConst->m_iStates--;			// remove closure on closed loops.  i.e.  1,2,3,1 => 1,2,3
			}
		}
	
	//-- Channel count 
	pmif->nchannel = mt.node["ChannelCount"].dataAsInt();

	//-- condensed list of classes : amp,dev
	pmif->AllocateClasses( mt.node["Condensed"]["ClassCount"].dataAsInt() );
	double * anAmp = (double *) mt.node["Amps"].data();
	double * anDev = (double *) mt.node["Stds"].data();
	for( int i=0; i<pmif->Classes(); ++i ) {
		pmif->Class(i)->m_nAmp=anAmp[i];
		pmif->Class(i)->m_nDev=anDev[i];
		}

	return true;
	}

//----- Convert CMrgMdlFmt format to QUB_Tree 
QUB_Tree MdlFinfToTree( CMrgMdlFmt * pmif ){
	QUB_Tree product = QUB_Tree::Create("ModelFile");
	int i;

	//--States 
	QUB_Tree states = product["States"];
	for( i=0; i<pmif->States(); ++i ) {
		QUB_Tree state = states.appendChild("State");
		state["x"].setData( QTR_TYPE_DOUBLE, pmif->State(i)->m_nXPos );
		state["y"].setData( QTR_TYPE_DOUBLE, pmif->State(i)->m_nYPos );
		state["Class"].setData( QTR_TYPE_INT, pmif->State(i)->m_iClass );
		state["Pr"].setData( QTR_TYPE_DOUBLE, pmif->State(i)->m_nProb );
		}

	//--Rates 
	QUB_Tree rates = product["Rates"];
	string cx;
	for ( i=0; i<pmif->Paths(); ++i ) {
		QUB_Tree rate = rates.appendChild("Rate");
		rate["States"].setNumData( QTR_TYPE_INT, 2, 1, pmif->Path(i)->m_aiState );
		rate["k0"].setNumData( QTR_TYPE_DOUBLE, 2, 1, pmif->Path(i)->m_anK0 );
		rate["k1"].setNumData( QTR_TYPE_DOUBLE, 2, 1, pmif->Path(i)->m_anK1 );

		rate["PNames"].appendChild("PName").setData( pmif->Path(i)->m_acPName[0] );
		rate["PNames"].appendChild("PName").setData( pmif->Path(i)->m_acPName[1] );
		
		rate["QNames"].appendChild("QName").setData( pmif->Path(i)->m_acQName[0] );
		rate["QNames"].appendChild("QName").setData( pmif->Path(i)->m_acQName[1] );

		rate["P"].setNumData( QTR_TYPE_INT, 2, 1, pmif->Path(i)->m_aiPIndex );
		rate["Q"].setNumData( QTR_TYPE_INT, 2, 1, pmif->Path(i)->m_aiQIndex );
		}

	//-- Constraints
	std::vector<string> constraintNames;
	constraintNames.push_back("FixRate");
	constraintNames.push_back("ScaleRate");
	constraintNames.push_back("LoopBal");
	constraintNames.push_back("FixExp");
	constraintNames.push_back("ScaleExp");
	// will merge ever return a LoopImbal constraint?

	QUB_Tree constraints = product["Constraints"];
	for ( i=0; i<pmif->Constraints(); ++i ) {
		QUB_Tree constraint = constraints.appendChild( constraintNames[ pmif->Constraint(i)->m_iType ] );
		constraint.setNumData( QTR_TYPE_INT, pmif->Constraint(i)->m_iStates, 1, &(pmif->Constraint(i)->m_aiStates[0]) );
		}

	//-- Misc 
	product["ChannelCount"].setData( QTR_TYPE_INT, 1 );

	//-- create Amp, Std, NAr; Ars
	int tenclass = max (10, pmif->Classes());
	std::vector<double> tenamps( tenclass, 0.0 ), tenstds( tenclass, 0.0 ), tenar( tenclass, 0.0 );
	for( i=0; i<pmif->Classes(); ++i ) {
		tenamps[i] = pmif->Class(i)->m_nAmp;
		tenstds[i] = pmif->Class(i)->m_nDev;
		}
	std::vector<int> tennar( tenclass, 0 );

	QUB_Tree amps = product.appendChild("Amps");
	QUB_Tree stds = product.appendChild("Stds");
	QUB_Tree nar = product.appendChild("NAr");

	amps.setNumData(QTR_TYPE_DOUBLE, tenclass, 1, &(tenamps[0]));
	stds.setNumData(QTR_TYPE_DOUBLE, tenclass, 1, &(tenstds[0]));
	nar.setNumData(QTR_TYPE_INT, tenclass, 1, &(tennar[0]));

	QUB_TreeIter arIter = product["Ars"].children();
	for ( i=0; i<tenclass; ++i ) {
		arIter.insert("Ar");
		arIter->setNumData( QTR_TYPE_DOUBLE, 10, 1, &(tenar[0]) );
		++arIter;	
		}

	return product;
	}

//----  Convert from Model trees to older models, merge, convert back to result model tree.
QUB_Tree MergeModels( QUB_Tree factors, std::ostream& quberr ) {
	CMrgMdlFmt * pM0=NULL, * pM1=NULL, * pM2=NULL;
	QUB_TreeMonoIter modi=factors.children(); 

	if( ! (*modi).isNull() && (pM0=new CMrgMdlFmt) != NULL  ) {
		ModelTreeToMdlfinf( ModelTree(*modi), pM0, quberr );
		++modi;
		}

	while( ! (*modi).isNull() 
			&& (pM1=new CMrgMdlFmt)!=NULL
			&& (pM2=new CMrgMdlFmt)!=NULL ) {
		ModelTreeToMdlfinf( ModelTree(*modi), pM1, quberr );
		combine( pM0, pM1, pM2, quberr );
		delete pM0;
		delete pM1;
		pM0=pM2;
		++modi;
		}
	
	QUB_Tree product = DoPareConstraints(MdlFinfToTree( pM0 ), QUB_Tree()/*discard removed constraints*/);
	delete pM0;
	
	return product;
	}

int UnmergeModels( QUB_Tree product, QUB_Tree factors, std::ostream& quberr ){
	//-- create Mdlfinf M for QUB_Tree product 
	CMrgMdlFmt M;
	ModelTreeToMdlfinf( ModelTree(product), &M, quberr );

	//-- create model inf for first factor
	QUB_TreeMonoIter modi = factors.children();
	QUB_TreeMonoIter path_i, state_i;
	double * pK0, * pK1;
	int  iPaths, iStates, iMergePaths=0, iMergeStates=1;
	while( ! (*modi).isNull() ) { 
		ModelTree mtree( *modi );

		//-- Count States 
		iStates=0;
		for( state_i = mtree.node["States"].children(); ! (*state_i).isNull(); ++state_i ) 
			++iStates;

		iPaths=0;
		for( path_i = mtree.node["Rates"].children(); ! (*path_i).isNull(); ++path_i ) {
			pK0 = (double *) (*path_i)["k0"].data();
			pK1 = (double *) (*path_i)["k1"].data();

			pK0[0] = M.Path( iMergePaths*iStates + iPaths )->m_anK0[0];
			pK1[0] = M.Path( iMergePaths*iStates + iPaths )->m_anK1[0];
			pK0[1] = M.Path( iMergePaths*iStates + iPaths )->m_anK0[1];
			pK1[1] = M.Path( iMergePaths*iStates + iPaths )->m_anK1[1];
			++iPaths;
			}

		mtree.restoreClasses();

		iMergePaths = iMergePaths * iStates  + iMergeStates * iPaths ; 
		iMergeStates = iMergeStates * iStates;

		++modi;
		}

	return 0; // negative for error; 0 for ok, positive for warning
	}

extern "C" QUBOPT_API QTR_Impl* DoMergeModels( QTR_Impl *factors, QTR_Callback reportCB ){
	QUB_Tree factorsTree( factors );
	mil_reportstream quberr( reportCB );

	QUB_Tree product = MergeModels( factorsTree, quberr );
	QTR_Impl *prodImpl = product.getImpl();
	if ( prodImpl )
		QTR_INCREF( prodImpl );

	return prodImpl;
	}

extern "C" QUBOPT_API int DoUnmergeModels( QTR_Impl *product, QTR_Impl *factors, QTR_Callback reportCB ){
	QUB_Tree productTree( product );
	QUB_Tree factorsTree( factors );

	mil_reportstream quberr( reportCB );

	return UnmergeModels( productTree, factorsTree, quberr );
	}

