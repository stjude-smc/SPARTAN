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

#ifndef MATRIXUTIL_H
#define MATRIXUTIL_H

#include <iostream>
#include "qubopt.h"
#include "milutil.h"   
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

QUBOPT_API int svdcmp_dll( double **a, int m, int n, double *w, double **v );
int svdcmp(double **a, int m, int n, double *w, double **v, std::ostream &msgOut);
int gaussj0(double** a, int n, double **b, int m, std::ostream &msgOut);
int gaussj_invert(fq::matrix<double>& mat, int n, std::ostream& msgOut);

#ifdef __cplusplus
}
#endif

#endif
