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

#include "qublib.h"

#define HAVE_INLINE
#include <gsl/gsl_eigen.h>

/* 
 * Find the eigenvalues and eigenvectors of an n-D matrix A.
 * The eigenvalues are stored in the array wr[] (real par) 
 * and wi[] (imaginary part), and the corresponding vectors
 * are stored in the columns of matrix V. So a, w and v are
 * related by AV=VW.
 */
 
#include <vector>
using std::vector;

bypass_eigen_f bypass_eigen = 0;

extern "C" QUBOPT_API void set_bypass_eigen(bypass_eigen_f eigen_f)
{
	bypass_eigen = eigen_f;
}


//-------------------------------------------------------------------------
extern "C" void QUBOPT_API eigen(double **a, int n, double *wr, double *wi, double **v) {

	if ( bypass_eigen ) {
		double *by_a, *by_v, *by_v_inv;
		vector<double> stor_a, stor_v, stor_v_inv(n*n);
		by_v_inv = & stor_v_inv[0];
		if ( (n == 1) || (a[1] - a[0]) == n ) {
			by_a = a[0];
		} else {
			stor_a.resize(n*n);
			by_a = & stor_a[0];
			for (int i=0; i<n; ++i)
				for (int j=0; j<n; ++j)
					by_a[i*n+j] = a[i][j];
		}
		if ( (n == 1) || (v[1] - v[0]) == n ) {
			by_v = v[0];
		} else {
			stor_v.resize(n*n+1);
			by_v = & stor_v[0];
		}

		bypass_eigen(n, by_a, wr, wi, by_v, by_v_inv);

		if ( by_v != v[0] ) {
			for (int i=0; i<n; ++i)
				for (int j=0; j<n; ++j)
					v[i][j] = *(by_v++);
		}
		return;
	}

  gsl_matrix *Q = gsl_matrix_alloc(n,n);
  for (int i=0; i<n; ++i)
	  for (int j=0; j<n; ++j)
		  gsl_matrix_set(Q,i,j, a[i][j]);
  gsl_matrix_complex *V = gsl_matrix_complex_alloc(n,n);
  gsl_vector_complex *L = gsl_vector_complex_alloc(n);
  gsl_eigen_nonsymmv_workspace *ws = gsl_eigen_nonsymmv_alloc(n);
  gsl_eigen_nonsymmv(Q, L, V, ws);
  for (int i=0; i<n; ++i) {
    wr[i] = GSL_REAL(gsl_vector_complex_get(L, i));
    wi[i] = 0.0; // GSL_IMAG(gsl_vector_complex_get(L, i));
    for (int j=0; j<n; ++j)
      v[i][j] = GSL_REAL(gsl_matrix_complex_get(V,i,j));
  }
  gsl_matrix_free(Q);
  gsl_matrix_complex_free(V);
  gsl_vector_complex_free(L);
  gsl_eigen_nonsymmv_free(ws);
}
