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

#include "mdlrep.h"
#include "qmatrixutil.h"
#include "qublib.h"

using namespace std;
using namespace fq;

/* Lorin's methods are in:
     QUB_QFS_QMF_QMatrix.pas:
	   MakeAMatrixPEq
	 
	 QUB_QFS_QMF:
	   CalculateMath
	 
	 QUB_QVO_QVM:
	   prWriteModelStats
*/


bool mdlrep_AMatrix( double** q, int n, double sampling, double** a, std::ostream& milerr )
// requires a[n+1][n+1]
{
	// using the spectral method of mil

	fq::vector<double> w(n+1);
	matrix<double> p(n+1, n+1), p1(n+1, n+1);
	int posieigencount = 0;

	if ( qspctrm(n, q, w, p, p1, posieigencount, milerr) == 1 ) {
		mil_expQ_full(n, w, p, p1, sampling, a);
		return true;
	}
	return false;
}

bool mdlrep_Eigen_entry_cmp(const mdlrep_Eigen_entry& a, const mdlrep_Eigen_entry& b)
{
	return a.val < b.val;
}

bool mdlrep_Eigen( double** q, int n, double* w, double** evect )
// returns sorted eigenvalues, column eigenvectors
{
	// using the qublib-based method of Lorin

	matrix<double> negQ(n, n);
	fq::vector<double> wi(n);
	int i, j;
	double tmp;

	for ( i=0; i<n; ++i )
		for ( j=0; j<n; ++j )
			negQ[i][j] = - q[i][j];
	eigen(negQ, n, w, wi, evect);

	std::vector<mdlrep_Eigen_entry> eigens;
	for ( i=0; i<n; ++i )
		eigens.push_back( mdlrep_Eigen_entry(w, evect, n, i) );
	sort(eigens.begin(), eigens.end(), mdlrep_Eigen_entry_cmp);

	for ( i=0; i<n; ++i ) {
		tmp = eigens[i].val;
		if ( tmp < 1e-10 )
			w[i] = 0.0;
		else
			w[i] = - tmp;

		for ( j=0; j<n; ++j )
			evect[j][i] = eigens[i].vec[j];
	}

	return true;
}

bool mdlrep_Spectrum( QUBOPT_VAR_NOT_USED double * w, double **evect, int n, tensor<double>& Ak, std::ostream& milerr )
{
	// starting from mdlrep_Eigen's sorted outputs
	// following Lorin's method

	int i, j, k;
	matrix<double> evectI(n, n);

	for ( i=0; i<n; ++i )
		for ( j=0; j<n; ++j )
			evectI[i][j] = evect[i][j];
	gaussj_invert(evectI, n, milerr);

	for ( k=0; k<n; ++k )
		for ( i=0; i<n; ++i )
			for ( j=0; j<n; ++j )
				Ak[k][i][j] = evect[i][k] * evectI[k][j];

	return true;
}

bool mdlrep_Current( double* w, tensor<double>& Ak, fq::vector<double>& p0, double* stateAmp, int n, double* B, double* T )
{
	// using outputs from mdlrep_Eigen and mdlrep_Spectrum
	// following Lorin's method

	int i, j;
	double tmp;

	for ( i=1; i<n; ++i ) {
		fq::vector<double> v = p0 * Ak[i];
		for ( j=0, tmp=0.; j<n; ++j )
			tmp += v[j] * stateAmp[j];
		B[i-1] = tmp;
		T[i-1] = - 1.0e3 / w[i];
	}

	return true;
}

bool mdlrep_Peq( double** q, int n, fq::vector<double>& Peq, std::ostream& milerr )
{
	// following method from Yu's Maple file:
	// S = [Q | 1]; U = unit row; Peq = U * inv(S * ST)

	matrix<double> S(n, n+1, 1.0), ST(n+1, n, 1.0);
	fq::vector<double> U(n, 1.0);
	int i, j;
	for ( i=0; i<n; ++i ) {
		for ( j=0; j<n; ++j ) {
			S[i][j] = ST[j][i] = q[i][j];
		}
	}
	matrix<double> SST = S * ST;
	gaussj_invert(SST, n, milerr); // now SST is (S * ST)-1
	Peq = U * SST;

	return true;
}

double mdlrep_Ieq(double* Peq, double* stateAmp, int n)
{
	// Ieq = sum( Peq[i] * amp[clazz[i]] )
	// Lorin's method

	double Ieq = 0.0;
	int i;

	for ( i=0; i<n; ++i )
		Ieq += Peq[i] * stateAmp[i];
	return Ieq;
}

bool mdlrep_FirstLatency( double** q, double* p0, int* clazz, int n, double tMax, int tBins,
						  double* flBins, double* flVals, double *Amp, double *Tau, std::ostream& milerr )
{
	// as per Yu:
	// sum( P0(closed) * exp(Qcc*t) * (-Qcc) )
	int i, j, k, l, posieigencount = 0;

	fq::vector<int> closed;
	for ( i=0; i<n; ++i )
		if ( 0 == clazz[i] )
			closed.push_back(i);
	int nc = closed.size();

	fq::vector<double> P0c( nc );
	matrix<double> Qcc(nc, nc), negQcc(nc, nc);
	for ( i=0; i<nc; ++i ) {
		P0c[i] = p0[ closed[i] ];
		for ( j=0; j<nc; ++j ) {
			Qcc[i][j] = q[closed[i]][closed[j]];
			negQcc[i][j] = - Qcc[i][j];
		}
	}

	fq::vector<double> wcc(nc+1);
	matrix<double> pcc(nc+1, nc+1), p1cc(nc+1, nc+1);
	matrix<double> acc(nc, nc);
	qspctrm(nc, Qcc, wcc, pcc, p1cc, posieigencount, milerr);

	double deltaT = tMax / (tBins - 1);
	double t;
	
	for ( i=0, t=0.0; i<tBins; ++i, t += deltaT ) {
		mil_expQ_full(nc, wcc, pcc, p1cc, t, acc);
		fq::vector<double> V = P0c * acc * negQcc;
		flVals[i] = dsumv(V, nc);
		flBins[i] = t;
	}


	// firstlat(t) = sum( B[l] * exp(-t / T[l]) )
	for ( l=0; l<nc; ++l ) {
		t = 0.0;
		for ( i=0; i<nc; ++i )
			for ( j=0; j<nc; ++j )
				for ( k=0; k<nc; ++k )
					t += Qcc[i][j] * P0c[k] * pcc[k][l] * p1cc[l][i];
		Amp[l] = - t;
		Tau[l] = - 1.0 / wcc[l];
	}


	return true;
}


