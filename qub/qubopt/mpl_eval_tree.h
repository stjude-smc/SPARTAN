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

#ifndef MPL_EVAL_TREE_H
#define MPL_EVAL_TREE_H

#include "qubopt.h"
#include "mil_eval_tree.h"

	class mpl_metamodel{
		public:
			mpl_metamodel( mil_model& mdl, int nchan );
			~mpl_metamodel();

			mil_model&       model;
			int              nchannel;

			// generated on construction:
			fq::vector<int>    idrug; // idrug[i] = mi.path[i][2]
			fq::matrix<double>   xlimit; // mi.xlimit transposed
			fq::vector<double> ratio; // ratio[i] == i

			fq::vector<double> amp, xms;
			fq::matrix<double>   ar;
			int              nar;

			int              nmutistate, nmuticlass, nmutipath;
			int              **mutistate, **mutipath;
			fq::vector<int>    muticlazz;
			fq::vector<double> mutipr;
		};

#ifdef __cplusplus
extern "C" {
#endif

	QUBOPT_API QTR_Impl *
	mpltreeiface( QTR_Impl *cfg, QTR_Impl **data, QTR_Callback fillDataCB,
				  QTR_Callback rptCB, QTR_Callback rsltCB, QTR_Callback pctCB );

#ifdef __cplusplus
}
#endif


#endif
