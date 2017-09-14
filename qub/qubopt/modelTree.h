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

#ifndef MODELTREE_H
#define MODELTREE_H

#include "qubopt.h"

#ifdef _WIN32
#include <hash_map>
using stdext::hash_map;
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;

namespace __gnu_cxx {
/**
        Explicit template specialization of hash of a string class,
        which just uses the internal char* representation as a wrapper.
 */
template <>
struct hash<std::string> {
        size_t operator() (const std::string& x) const {
                return hash<const char*>()(x.c_str());
	// hash<const char*> already exists
        }
};
}
#endif

	class ModelTree	{
		public:
			ModelTree();
			ModelTree( QUB_Tree mdlNode );
			ModelTree( const ModelTree& copy );
			virtual ~ModelTree();

			QUB_Tree node;
			hash_map<string, int> indexOfParam;
			int nparam;
			std::vector<int> cns_src_ix;

			void updateParams(); // if a rate's ligand or voltage could have changed since construction
			void verifyArs(); // something in qub likes to remove unused classes but a lot of old code here depends on 'em

			int toMdlinf( mdlinf& mi );
			// just nchannel, amp, sd, ar
			void toDatinf( datinf& di );
			void fromDatinf( datinf& di ); // same nclass, nar (meant for amp updates)

			// old algos don't work if there are no states of class 0 (or any class < max used class)
			// do this before .toMdlinf(), restore after committing changes
			void condenseClasses();
			void restoreClasses();
			int  classesPresent();

			void updateRates( mdlinf& mi );
			void updateRates( mdlinf& mi, fq::vector<double>& xsd );

			// update amp, sd, ar

			// generate metamodel?
		};
	QUB_Tree DoCheckLoopStimuli(QUB_Tree model);

	QUB_Tree GetScaledRates(QUB_Tree constraints, int from, int to);
	QUB_Tree GetScaledExpRates(QUB_Tree constraints, int from, int to);

	QUB_Tree DoBalanceAllLoops(QUB_Tree model);
	QUB_Tree DoPareConstraints(QUB_Tree model, QUB_Tree removed);
	QUB_Tree DoPareLoopConstraints(QUB_Tree model);



#ifdef __cplusplus
extern "C" {
#endif

	QUBOPT_API QTR_Impl* DoMergeModels( QTR_Impl *factors, QTR_Callback reportCB );
	QUBOPT_API int DoUnmergeModels( QTR_Impl *product, QTR_Impl *factors, QTR_Callback reportCB );

	QUBOPT_API QTR_Impl* BalanceAllLoops( QTR_Impl *model );

	QUBOPT_API QTR_Impl* PareConstraints( QTR_Impl *model, QTR_Impl *removed );
	QUBOPT_API QTR_Impl* PareLoopConstraints( QTR_Impl *model );

	QUBOPT_API QTR_Impl* CheckLoopStimuli( QTR_Impl *model );

	QUBOPT_API QTR_Impl* _GetScaledRates(QTR_Impl* constraints, int from, int to);
	QUBOPT_API QTR_Impl* _GetScaledExpRates(QTR_Impl* constraints, int from, int to);

#ifdef __cplusplus
}
#endif


#endif
