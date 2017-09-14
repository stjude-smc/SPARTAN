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

#ifndef QUBOPT_H
#define QUBOPT_H

#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1

#ifdef _USRDLL
#define _LIBMIL_LIB_
#endif 

//#if defined(_WIN32)
//	#define QUBOPT_DLLEXPORT __declspec(dllexport)
//  #define QUBOPT_DLLIMPORT __declspec(dllimport)
//#else
	#define QUBOPT_DLLEXPORT
    #define QUBOPT_DLLIMPORT
//#endif

#if !defined(_LIBMIL_LIB_)
	#define QUBOPT_API QUBOPT_DLLIMPORT
#else
	#define QUBOPT_API QUBOPT_DLLEXPORT
#endif

#ifdef __GNUC__
#define QUBOPT_VAR_NOT_USED __attribute__ ((unused))
#else
#define QUBOPT_VAR_NOT_USED
#endif


// disable error related to dll export of STL method 
#pragma warning (disable: 4251)

//#include <FLOAT.h>
#include <map>
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iostream>
#include <list>
#include <fstream>
#include <set>
#include <iomanip>

#include <QUB_Tree.h>

#include "qub_reportstream.h"

#include "milutil.h"
#include "matrixutil.h"

//----- Model Structures   (Formerly ch_x.h)
	typedef struct tag_mdlinf{
		fq::vector<int> clazz; // class of state i		
		fq::matrix<int> path;    // [from, to, drug] of path i		
		fq::vector<double> x ; // [k0, k1] for each path		
		fq::vector<double> pr;  // start prob of state i		
		
		int nstate;  // derived		
		int npath;   // derived // one-directional		
		
		fq::vector<double> xlimit[2];		
		fq::matrix<double> mtx;		
		fq::vector<double> vct;		
		} mdlinf;

	class datinf {
		public:
			int nchannel;		
			fq::vector<double> i;       // single amplitudes		
			fq::matrix<double> r;       // single correlations		
		};

#endif
