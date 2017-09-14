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

#ifndef MDLREP_H
#define MDLREP_H

#include <iostream>
#include "matrix.h"

	void mil_expQ_full(int n, double* w, double** p, double** p1, double t, double** a);



bool mdlrep_AMatrix( double** q, int n, double sampling, double** a, std::ostream& milerr );
// requires a[n+1][n+1]

class mdlrep_Eigen_entry
{
public:
	double val;
	fq::vector<double> vec;

	mdlrep_Eigen_entry(double* w, double** evect, int n, int i)
		: val(w[i]), vec(n)
	{
		for ( int j=0; j<n; ++j )
			vec[j] = evect[j][i];
	}
};

bool mdlrep_Eigen_entry_cmp(const mdlrep_Eigen_entry& a, const mdlrep_Eigen_entry& b);

bool mdlrep_Eigen( double** q, int n, double* w, double** evect );
// returns sorted eigenvalues, column eigenvectors

bool mdlrep_Spectrum( double *w, double **evect, int n, fq::tensor<double>& Ak, std::ostream& milerr );

bool mdlrep_Current( double* w, fq::tensor<double>& Ak, fq::vector<double>& p0, double* stateAmp, int n, double* B, double* T );

bool mdlrep_Peq( double** q, int n, fq::vector<double>& Peq, std::ostream& milerr );

double mdlrep_Ieq(double* Peq, double* stateAmp, int n);

bool mdlrep_FirstLatency( double** q, double* p0, int* clazz, int n, double tMax, int tBins,
						  double* flBins, double* flVals, double *Amp, double *Tau, std::ostream& milerr );


#endif
