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

#include "sparseutil.h"
#include "qubopt.h"
#include <math.h>

using namespace std;

void dsprsin(double **a, int n, double thresh, unsigned long nmax, double sa[],unsigned long ija[], std::ostream& errout)
// Converts a square matrix a[1..n][1..n] into row-indexed sparse storage mode. Only elements
// of a with magnitude =thresh are retained. Output is in two linear arrays with dimension
// nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
// number of elements filled of sa and ija on output are both ija[ija[1]-1]-1 (see text).
{
	int i,j;
	unsigned long k;
	for (j=1;j<=n;j++) sa[j]=a[j][j];  // Store diagonal elements.
	ija[1]=n+2;                        // Index to 1st rowoff- diagonal element, if any.
	k=n+1;
	for (i=1;i<=n;i++) {               // Loop over rows.
		for (j=1;j<=n;j++) {           // Loop over columns.
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) errout << "sprsin: nmax too small" << endl;
				sa[k]=a[i][j];         // Store off-diagonal elements and their columns.
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;                  // As each rowi s completed, store index to next.
	}
}

void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n, std::ostream& errout)
// Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x[1..n], giving
// a vector b[1..n].
{
	unsigned long i,k;
	if (ija[1] != n+2) errout << "sprsax: mismatched vector and matrix" << endl;
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];                    // Start with diagonal term.
		for (k=ija[i];k<=ija[i+1]-1;k++)    // Loop over off-diagonal terms.
			b[i] += sa[k]*x[ija[k]];
	}
}

void dsprstx(double sa[], unsigned long ija[], double x[], double b[],unsigned long n, std::ostream& errout)
// Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a vector
// x[1..n], giving a vector b[1..n].
{
	unsigned long i,j,k;
	if (ija[1] != n+2) errout << "mismatched vector and matrix in sprstx" << endl;
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];     // Start with diagonal terms.
	for (i=1;i<=n;i++) {                    // Loop over off-diagonal terms.
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}

#define M 7
#define NSTACK 50

template<class T>
void SWAP(T& x, T& y)
{
	T tmp = x;
	x = y;
	y = tmp;
}

void iindexx(unsigned long n, long arr[], unsigned long indx[], std::ostream& errout)
// Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is
// in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.
{
	unsigned long i,indxt,ir=n,j,k,l=1;
	int jstack=0;
	long *istack;
	long a;
	fq::vector<long> _istack(1+NSTACK);  istack = _istack;
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir]);
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1]);
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j]);
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) errout << "NSTACK too small in indexx.";
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}

#undef M
#undef NSTACK

void dsprstp(double sa[], unsigned long ija[], double sb[], unsigned long ijb[], std::ostream& errout)
// Construct the transpose of a sparse square matrix, from row-index sparse storage arrays sa and
// ija into arrays sb and ijb.
{
	unsigned long j,jl(0),jm,jp,ju,k,m,n2,noff,inc,iv;
	double v;
	n2=ija[1];                               // Linear size of !matrix plus 2.
	for (j=1;j<=n2-2;j++) sb[j]=sa[j];       // Diagonal elements.
	iindexx(ija[n2-1]-ija[1],(long *)&ija[n2-1],&ijb[n2-1], errout);
	                          // Index all off-diagonal elements by their columns.
	jp=0;
	for (k=ija[1];k<=ija[n2-1]-1;k++) {      // Loop over output off-diagonal elements.
		m=ijb[k]+n2-1;                       // Use index table to store by (former) columns.
		sb[k]=sa[m];
		for (j=jp+1;j<=ija[m];j++) ijb[j]=k; // Fill in the index to any omitted rows.
		jp=ija[m];                           // Use bisection to find which row element
		                                     // m is in and put that into jl=1; ijb[k].
		ju=n2-1;
		while (ju-jl > 1) {
			jm=(ju+jl)/2;
			if (ija[jm] > m) ju=jm; else jl=jm;
		}
		ijb[k]=jl;
	}
	for (j=jp+1;j<n2;j++) ijb[j]=ija[n2-1];
	for (j=1;j<=n2-2;j++) {                  // Make a final pass to sort each rowb y
		jl=ijb[j+1]-ijb[j];                  // Shell sort algorithm.
		noff=ijb[j]-1;
		inc=1;
		do {
			inc *= 3;
			inc++;
		} while (inc <= jl);
		do {
			inc /= 3;
			for (k=noff+inc+1;k<=noff+jl;k++) {
				iv=ijb[k];
				v=sb[k];
				m=k;
				while (ijb[m-inc] > iv) {
					ijb[m]=ijb[m-inc];
					sb[m]=sb[m-inc];
					m -= inc;
					if (m-noff <= inc) break;
				}
				ijb[m]=iv;
				sb[m]=v;
			}
		} while (inc > 1);
	}
}

void sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
			float sc[], unsigned long ijc[], std::ostream& errout)
// Matrix multiply A · BT where A and B are two sparse matrices in row-index storage mode, and
// BT is the transpose of B. Here, sa and ija store the matrix A; sb and ijb store the matrix B.
// This routine computes only those components of the matrix product that are pre-specified by the
// input index array ijc, which is not modified. On output, the arrays sc and ijc give the product
// matrix in row-index storage mode. For sparse matrix multiplication, this routine will often be
// preceded by a call to sprstp, so as to construct the transpose of a known matrix into sb, ijb.
{
	unsigned long i,ijma,ijmb,j,m,ma,mb,mbb,mn(0);
	float sum;
	if (ija[1] != ijb[1] || ija[1] != ijc[1])
	errout << "sprspm: sizes do not match";
	for (i=1;i<=ijc[1]-2;i++) {          // Loop over rows.
		j=m=i;                           // Set up so that first pass through loop does the
		                                 // diagonal mn=ijc[i]; component.
		sum=sa[i]*sb[i];
		for (;;) {                       // Main loop over each component to be output.
			mb=ijb[j];
			for (ma=ija[i];ma<=ija[i+1]-1;ma++) {
				// Loop through elements in A’s row. Convoluted logic, following, accounts for the
				// various combinations of diagonal and off-diagonal elements.
				ijma=ija[ma];
				if (ijma == j) sum += sa[ma]*sb[j];
				else {
					while (mb < ijb[j+1]) {
						ijmb=ijb[mb];
						if (ijmb == i) {
							sum += sa[i]*sb[mb++];
							continue;
						} else if (ijmb < ijma) {
							mb++;
							continue;
						} else if (ijmb == ijma) {
							sum += sa[ma]*sb[mb++];
							continue;
						}
						break;
					}
				}
			}
			for (mbb=mb;mbb<=ijb[j+1]-1;mbb++) { // Exhaust the remainder of B’s row.
				if (ijb[mbb] == i) sum += sa[i]*sb[mbb];
			}
			sc[m]=sum;
			sum=0.0;                             // Reset indices for next pass through loop.
			if (mn >= ijc[i+1]) break;
			j=ijc[m=mn++];
		}
	}
}

void sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[], float thresh,
			unsigned long nmax, float sc[], unsigned long ijc[], std::ostream& errout)
// Matrix multiply A · BT where A and B are two sparse matrices in row-index storage mode, and
// BT is the transpose of B. Here, sa and ija store the matrix A; sb and ijb store the matrix
// B. This routine computes all components of the matrix product (which may be non-sparse!),
// but stores only those whose magnitude exceeds thresh. On output, the arrays sc and ijc
// (whose maximum size is input as nmax) give the product matrix in row-index storage mode.
// For sparse matrix multiplication, this routine will often be preceded by a call to sprstp, so as
// to construct the transpose of a known matrix into sb, ijb.
{
	unsigned long i,ijma,ijmb,j,k,ma,mb,mbb;
	float sum;
	if (ija[1] != ijb[1]) errout << "sprstm: sizes do not match";
	ijc[1]=k=ija[1];
	for (i=1;i<=ija[1]-2;i++) {                  // Loop over rows of A,
		for (j=1;j<=ijb[1]-2;j++) {              // and rows of B.
			if (i == j) sum=sa[i]*sb[j]; else sum=0.0e0;
			mb=ijb[j];
			for (ma=ija[i];ma<=ija[i+1]-1;ma++) {
				// Loop through elements in A’s row. Convoluted logic, following, accounts for the
				// various combinations of diagonal and off-diagonal elements.
				ijma=ija[ma];
				if (ijma == j) sum += sa[ma]*sb[j];
				else {
					while (mb < ijb[j+1]) {
						ijmb=ijb[mb];
						if (ijmb == i) {
							sum += sa[i]*sb[mb++];
							continue;
						} else if (ijmb < ijma) {
							mb++;
							continue;
						} else if (ijmb == ijma) {
							sum += sa[ma]*sb[mb++];
							continue;
						}
						break;
					}
				}
			}
			for (mbb=mb;mbb<=ijb[j+1]-1;mbb++) { // Exhaust the remainder of B’s row.
				if (ijb[mbb] == i) sum += sa[i]*sb[mbb];
			}
			if (i == j) sc[i]=sum;               // Where to put the answer...
			else if (fabs(sum) > thresh) {
				if (k > nmax) errout << "sprstm: nmax too small";
				sc[k]=sum;
				ijc[k++]=j;
			}
		}
		ijc[i+1]=k;
	}
}

double snrm(unsigned long n, double sx[], int itol)
// Compute one of two norms for a vector sx[1..n], as signaled by itol. Used by linbcg.
{
	unsigned long i,isamax;
	double ans;
	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i]; // Vector magnitude norm.
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {                   // Largest component norm.
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}

CSparseSolver::CSparseSolver(std::ostream& err, int maxn)
: nmax(maxn), n(0), ija(maxn+1), sa(maxn+1), errout(err)
{}

void CSparseSolver::SetA(double **a, int n_, double thresh)
{
	this->n = n_;
	dsprsin(a, n_, thresh, nmax, sa, ija, errout);
}

void CSparseSolver::atimes(unsigned long n_, double x[], double r[], int itrnsp)
{
	if (itrnsp) dsprstx(sa,ija,x,r,n_,errout);
	else dsprsax(sa,ija,x,r,n_,errout);
}

void CSparseSolver::asolve(unsigned long n_, double b[], double x[], QUBOPT_VAR_NOT_USED int itrnsp)
{
	unsigned long i;
	for(i=1;i<=n_;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
	// The matrix A is the diagonal part of A, stored in the first n elements of sa. Since the
	// transpose matrix has the same diagonal, the flag itrnsp is not used.
}

void CSparseSolver::linbcg(double *b, double *x, int itol, double tol, int itmax, int& iter, double& err)
// Solves A · x = b for x[1..n], given b[1..n], by the iterative biconjugate gradient method.
// On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
// or 4, specifying which convergence test is applied (see text); itmax is the maximum number
// of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
// reset to the improved solution, iter is the number of iterations actually taken, and err is the
// estimated error. The matrix A is referenced only through the user-supplied routines atimes,
// which computes the product of either A or its transpose on a vector; and asolve, which solves
// A·x = b or AT
// ·x = b for some preconditioner matrix A (possibly the trivial diagonal part of A).
{
	const double EPS = 1.0e-14;

	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;   // Double precision is a good idea in this routine.
	fq::vector<double> _p (1+n); p  = _p;
	fq::vector<double> _pp(1+n); pp = _pp;
	fq::vector<double> _r (1+n); r  = _r;
	fq::vector<double> _rr(1+n); rr = _rr;
	fq::vector<double> _z (1+n); z  = _z;
	fq::vector<double> _zz(1+n); zz = _zz;
	                               // Calculate initial residual.
	iter=0;
	atimes(n,x,r,0);               // Input to atimes is x[1..n], output is r[1..n];
	                               // the final 0 indicates that the matrix (not its
	                               // transpose) is to be used.
	for (j=1;j<=(unsigned long)n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	/* atimes(n,r,rr,0); */        // Uncomment this line to get the “minimum residual”
	                               // variant of the algorithm.
	if (itol == 1) {
		bnrm=snrm(n,b,itol);
		asolve(n,r,z,0);           // Input to asolve is r[1..n], output is z[1..n];
		                           // the final 0 indicates that the matrix A (not
		                           // its transpose) is to be used.
	}
	else if (itol == 2) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0);
		znrm=snrm(n,z,itol);
	} else errout << "illegal itol in linbcg";
	while (iter <= itmax) {       // Main loop.
		++(iter);
		asolve(n,rr,zz,1);         // Final 1 indicates use of transpose matrix AT.
		for (bknum=0.0,j=1;j<=(unsigned long)n;j++) bknum += z[j]*rr[j];
		                        // Calculate coefficient bk and direction vectors p and pp.
		if (iter == 1) {
			for (j=1;j<=(unsigned long)n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=(unsigned long)n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;            // Calculate coefficient ak, newi terate x, and new
		atimes(n,p,z,0);        // residuals r and rr.
		for (akden=0.0,j=1;j<=(unsigned long)n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(n,pp,zz,1);
		for (j=1;j<=(unsigned long)n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(n,r,z,0);        // Solve A · z = r and check stopping criterion.
		if (itol == 1)
			err=snrm(n,r,itol)/bnrm;
		else if (itol == 2)
			err=snrm(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				err=znrm/bnrm;              // Error may not be accurate, so loop again.
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (err <= 0.5*xnrm) err /= xnrm;
			else {
				err=znrm/bnrm;              // Error may not be accurate, so loop again.
				continue;
			}
		}
		errout << "iter=" << iter << " err=" << err << endl;
		if (err <= tol) break;
	}
}

