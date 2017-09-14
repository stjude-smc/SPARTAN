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

#ifndef STAT_EVAL_TREE_H
#define STAT_EVAL_TREE_H

#include "qubopt.h"

	QUB_Tree skm_ultd_runStat( QUB_Tree dataSet, QUB_Tree config, QTR_Callback fillDataCB, QTR_Callback reportCB );
	QUB_Tree ApplyDeadTime(QUB_Tree DataSet, float tdMS);

#ifdef __cplusplus
extern "C" {
#endif


//----- From Stat_eval_tree
	QUBOPT_API QTR_Impl * stattreeiface( QTR_Impl *cfg, QTR_Impl **data, QTR_Callback fillDataCB, QTR_Callback rptCB );


#ifdef __cplusplus
}
#endif


#endif
