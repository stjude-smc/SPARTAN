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

#ifndef DFP_OPTIMIZE_H
#define DFP_OPTIMIZE_H

#include "qubopt.h"

	class dfp_optimizable{
		public:
			virtual ~dfp_optimizable() {}
			
			virtual bool getStartingZ( fq::vector<double> &z ) = 0;
			virtual bool evaluate( fq::vector<double> &z, fq::vector<double> &gz,
							double &ll, int &nf, int &ndf, bool flagGrad=true ) = 0;
			virtual bool checkContinue( int iter, int nf, int ndf, double ll,
							fq::vector<double> &z, fq::vector<double> &gz )=0;
		};

	class dfp_result{
		public:
			fq::vector<double> z;
			double ll;
			int iter;
			int nfcn, ndfcn;
			fq::matrix<double> hessian;
			int err; // ? 
		};

	dfp_result dfp_optimize( dfp_optimizable *optable, int maxIter, double convLL, 
							double convGrad, double stepMax = 1.0 );

	// dll exported for delphi.	Both return 0 for success
	typedef int (*dfp_dll_func)(void*, double*, double*, double*, int*, int*, int);
	//                          caller, params,  grads,   ll,      nf,   ndf,  (bool) computeGrads
	typedef int (*dfp_dll_check)(void*, int, int, int, double, double*, double*);
	//                          caller, iters, nf, ndf,   ll     params,  grads



#ifdef __cplusplus
extern "C" {
#endif

	QUBOPT_API
	double dfp_opt_dll( void *call_obj, int npar, double *initPar, double **hessianOut,
						dfp_dll_func func, dfp_dll_check check,
						int maxIter, double convLL, double convGrad, double stepMax,
						int *errOut );
		// errOut:  0   success
		//          -1  func returned error
		//          -2  check returned error
		//          -3  exceeded max iter

#ifdef __cplusplus
}
#endif

#endif
