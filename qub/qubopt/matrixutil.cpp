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

#include "matrixutil.h"
#include "qublib.h"

#ifdef _WIN32
  #include <windows.h>
#else
  #include <stdlib.h>
  #define BOOL int
  #define TRUE 1
  #define FALSE 0
#endif

#include <math.h>
#include <iostream>

using namespace std;
using namespace fq;

//----- Local functions ( optional declaration ) 
void balanc(double** a, int n);
void elmhes0(double** a, int n);		// TODO - not yet used 
void elmhes(double** a, int n);
void hmges(int n, double* w, double** v, double* x);
BOOL hqr(double** a, int n, double* wr, double* wi);
void swap(double *a, double *b, int iItems);


// ----------------------------------------------------------
// Swap two doubles - swap(&n,&n) function or SWAP(n,n) macro
double g_nTemp;
#define SWAP(g,h) {g_nTemp=(g);(g)=(h);(h)=g_nTemp;}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void swap(double *a, double *b, int iItems){
	double tmp;
	for( int i=0; i<iItems; ++i ) {
		tmp=a[i];
		a[i]=b[i];
		b[i]=tmp;
		}
	}

/*	---------------------------------------------------------------
	Gauss-Jordan linear equation solver
	a[0..n-1][0..n-1] is the input matrix
	b[0..n-1][0..m-1] is input containing the m right-hand side vectors
	on output, a is replaced by its matrix inverse,
	and b is replaced by the corresponding set of solution vectors.
*/

//--  Uses 0 based arrays.
//--  Added comments.  More meaningful variable names.  Better Swap usage.
BOOL gaussj0(double** a, int n, double **b, int m, ostream &msgOut){
	int iIter, iRow, iCol, iMaxCol, iMaxRow;	
	double 	big, dum, pivinv;
	int * indxc = int_alloc1D(n);	// [i] is the column of the i'th pivot element
	int * indxr = int_alloc1D(n);	// [i] is the original row for the i'th pivot 
	int * ipiv = int_alloc1D(n);
	
	//-- Set all pivots to 0 ( false ) 
	for( iRow=0; iRow<n; iRow++) 
		ipiv[iRow]=0;

	iMaxRow = iMaxCol = -1;

	//-- For each iteration, find max and {divide, swap}.
	for (iIter=0; iIter<n; iIter++) {
		//-- For i,j not yet pivoted, find where max abs(a[i][j]) occurs ==> Pivot 
		big=0.0;
		for (iRow=0; iRow<n; iRow++) {
			if (ipiv[iRow] != 1) {
				for (iCol=0; iCol<n; iCol++) {
					if (ipiv[iCol] == 0) { 
						if (fabs(a[iRow][iCol]) > big)
							big=fabs(a[iMaxRow=iRow][iMaxCol=iCol]);
						}
					else if (ipiv[iCol]>1) {
						msgOut << "Fatal err: Singular matrix in gauss-j" << endl;
						free(ipiv);
						free(indxr);
						free(indxc);
						return FALSE;
						}
					}
				}
			}
		if ( iMaxRow == -1 || iMaxCol == -1 ) {
			msgOut << "gauss-j found no pivot" << endl;
			free(ipiv);
			free(indxr);
			free(indxc);
			return FALSE;
			}

		//-- Got the pivot element, interchange rows to force the pivot element to 
		//-- be on the diagonal.  
		++(ipiv[iMaxCol]);
		if (iMaxRow != iMaxCol) {
			swap( a[iMaxRow], a[iMaxCol], n );
			if( b != NULL ) 
				swap( b[iMaxRow], b[iMaxCol], m );
			}

		indxr[iIter]=iMaxRow;
		indxc[iIter]=iMaxCol;
		pivinv = a[iMaxCol][iMaxCol];
		a[iMaxCol][iMaxCol]=1.0;		//-- this is unusual, but is fixed by the next if stmt.
		//-- divide row iCol by big in a and b
		for (iCol=0; iCol<n; iCol++) 
			a[iMaxCol][iCol] /= pivinv;
		for (iCol=0; iCol<m; iCol++) 
			b[iMaxCol][iCol] /= pivinv;

		//-- For each row, subtract a multiple of the maxrow which makes maxcol element 0
		for (iRow=0; iRow<n; iRow++) {
			if (iRow != iMaxCol) {
				dum=a[iRow][iMaxCol];
				a[iRow][iMaxCol]=0.0;
				for (iCol=0; iCol<n; iCol++) 
					a[iRow][iCol] -= a[iMaxCol][iCol]*dum;
				for (iCol=0; iCol<m; iCol++) 
					b[iRow][iCol] -= b[iMaxCol][iCol]*dum;
				}
			}
		}
	
	// unscramble the solution w.r.t. the column interchanges
	// by interchanging columns in the reverse order 
	for( iIter=n-1; iIter>=0; iIter--)
		if (indxr[iIter] != indxc[iIter])
			for (iRow=0; iRow<n; iRow++)
				SWAP( (a[iRow][indxr[iIter]]), (a[iRow][indxc[iIter]]) );

	free(ipiv);
	free(indxr);
	free(indxc);
	return TRUE;
	}

BOOL gaussj_invert(matrix<double>& mat, int n, std::ostream& msgOut){
	return gaussj0(mat, n, NULL, 0, msgOut );
	}



extern "C" QUBOPT_API int svdcmp_dll( double **a, int m, int n, double *w, double **v )
{
	return svdcmp( a, m, n, w, v, cerr );
}

int svdcmp(double **a, int m, int n, double *w, double **v, ostream &msgOut){
	int rtn = 0;

	fq::vector<double*> zrows_a, zrows_v;
	for (int i=1; i<=m; ++i)
		zrows_a.push_back(a[i]+1);
	for (int i=1; i<=n; ++i)
		zrows_v.push_back(v[i]+1);
	
	svdecomp0(&zrows_a[0], m, n, w+1, &zrows_v[0]);
	return rtn;
	}


