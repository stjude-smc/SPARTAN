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

#ifndef SKM_EVAL_TREE_H
#define SKM_EVAL_TREE_H

#include "qubopt.h"
#include "CountedPtr.h"
#include "modelTree.h"
#include "qmatrixutil.h"

using namespace std;

class skm_model{
	public:
		skm_model( QUB_Tree modelNode, int *fFixAmp=0, int *fFixSD=0 );
		
		QUB_Tree        node;
		ModelTree       tree;
		mdlinf          mi;

		int             nclass, nchannel;
		fq::vector<int>   path, fixAmp, fixSD;
		fq::vector<float> x, q, amp, xms;
	};

class skm_metamodel{
	public:
		skm_metamodel( skm_model& mdl, int nchan, double dtMS );
		~skm_metamodel();

		skm_model& model;
		int nchannel;

		int             nmutistate, nmuticlass, nmutipath;
		int             *mutistate, *mutipath; // allocated in skm_mutimakv() [skm-utils.cpp]
		fq::vector<int>   cmutistate, mutiFixAmp, mutiFixSD;
		fq::vector<float> mutiq, mutia, mutipr, mutiamp, mutixms;
	};

QUB_Tree MakeIdealNode( int *byPoint, int ndata, float dtMS );

QUB_Tree skm_unlimited( vector<QUB_Tree>& dataSets, QUB_Tree config,
					    QTR_Callback fillDataCB, QTR_Callback reportCB, QTR_Callback pctCB, bool doStat = false );

#ifdef __cplusplus
extern "C" {
#endif


QUBOPT_API QTR_Impl *
skmtreeiface( QTR_Impl *cfg, QTR_Impl **data,
			  QTR_Callback fillDataCB, QTR_Callback rptCB, QTR_Callback pctCB );

QUBOPT_API QTR_Impl * skm_get_multiclasses_iface(QTR_Impl *model, int channelCount);


#ifdef __cplusplus
}
#endif


#endif
