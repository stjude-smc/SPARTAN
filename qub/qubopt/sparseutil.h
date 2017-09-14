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

#ifndef SPARSEUTIL_H
#define SPARSEUTIL_H

#include <fstream>
#include "matrix.h"

// numerical recipes 2-7 except iindexx (nr 4-8)
// refashioned into CSparseSolver

void dsprsin(double **a, int n, double thresh, unsigned long nmax, double sa[],unsigned long ija[], std::ostream& errout);
// Converts a square matrix a[1..n][1..n] into row-indexed sparse storage mode. Only elements
// of a with magnitude =thresh are retained. Output is in two linear arrays with dimension
// nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
// number of elements filled of sa and ija on output are both ija[ija[1]-1]-1 (see text).

void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n, std::ostream& errout);
// Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x[1..n], giving
// a vector b[1..n].

void dsprstx(double sa[], unsigned long ija[], double x[], double b[],unsigned long n, std::ostream& errout);
// Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a vector
// x[1..n], giving a vector b[1..n].

void iindexx(unsigned long n, long arr[], unsigned long indx[], std::ostream& errout);
// Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is
// in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.

void dsprstp(double sa[], unsigned long ija[], double sb[], unsigned long ijb[], std::ostream& errout);
// Construct the transpose of a sparse square matrix, from row-index sparse storage arrays sa and
// ija into arrays sb and ijb.

void sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
			float sc[], unsigned long ijc[], std::ostream& errout);
// Matrix multiply A · BT where A and B are two sparse matrices in row-index storage mode, and
// BT is the transpose of B. Here, sa and ija store the matrix A; sb and ijb store the matrix B.
// This routine computes only those components of the matrix product that are pre-specified by the
// input index array ijc, which is not modified. On output, the arrays sc and ijc give the product
// matrix in row-index storage mode. For sparse matrix multiplication, this routine will often be
// preceded by a call to sprstp, so as to construct the transpose of a known matrix into sb, ijb.

void sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[], float thresh,
			unsigned long nmax, float sc[], unsigned long ijc[], std::ostream& errout);
// Matrix multiply A · BT where A and B are two sparse matrices in row-index storage mode, and
// BT is the transpose of B. Here, sa and ija store the matrix A; sb and ijb store the matrix
// B. This routine computes all components of the matrix product (which may be non-sparse!),
// but stores only those whose magnitude exceeds thresh. On output, the arrays sc and ijc
// (whose maximum size is input as nmax) give the product matrix in row-index storage mode.
// For sparse matrix multiplication, this routine will often be preceded by a call to sprstp, so as
// to construct the transpose of a known matrix into sb, ijb.

double snrm(unsigned long n, double sx[], int itol);
// Compute one of two norms for a vector sx[1..n], as signaled by itol. Used by linbcg.

class CSparseSolver
// wrapper on linbcg (to better manage storage)
{
public:
	CSparseSolver(std::ostream& err, int nmax);

	void SetA(double **a, int n, double thresh); // a 1-based square matrix

	void linbcg(double *b, double *x, int itol, double tol, int itmax, int& iter, double& err);
// Solves A · x = b for x[1..n], given b[1..n], by the iterative biconjugate gradient method.
// On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
// or 4, specifying which convergence test is applied (see text); itmax is the maximum number
// of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
// reset to the improved solution, iter is the number of iterations actually taken, and err is the
// estimated error. The matrix A is referenced only through the following routines atimes,
// which computes the product of either A or its transpose on a vector; and asolve, which solves
// A·x = b or AT
// ·x = b for some preconditioner matrix A (possibly the trivial diagonal part of A).

	void atimes(unsigned long n, double x[], double r[], int itrnsp);
	void asolve(unsigned long n, double b[], double x[], int itrnsp);

private:
	int nmax, n;
	fq::vector<unsigned long> ija;
	fq::vector<double> sa;
	std::ostream& errout;
};

#endif
