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


#define  PROCESS_STOPPED    0
#define  PROCESS_STOPPING	1
#define  PROCESS_STARTED	2
#define  PROCESS_STARTING	3

//----- Local functions 
void oupbaum(int iter, double ll, int nstate, double **a, int namp, double *amp, double *xms, 
				 int nar, double **ar, MSGOUTPUT msgoutput);

void mnewpr(int nar, int nfir, int nmetastate, double **gama, double *pr);
void mnewa(int nstate, int nar, int nfir, int nmetastate,
		   int **metastate, int ndata, double **alpha,
		   double **beta, double **b, double **a);

void mnewar(int nar, int nfir, int nmetastate, int **metastate,
			int ndata, double *data, double **gama, int namp,
			double *amp, double *xms, double **ar, double *fir);

void mnewamp1(int nar, int nfir, int nmetastate, int **metastate, 
			  int ndata, double *data, double **gama, double *xms, 
			  double **ar, double *fir, int namp, double *amp);

void gaussj(double **a, int n);
              
void mexasol(int nar, int nfir, int nmetastate, int **metastate, int ndata, double *data, double **gama, double *xms, double **ar, double *fir, int namp, double *amp);


//--------------------------------------------------------------------------
void oupbaum(int iter, double ll, QUBOPT_VAR_NOT_USED int nstate, QUBOPT_VAR_NOT_USED double ** a, int namp, double *amp, double *xms, QUBOPT_VAR_NOT_USED int nar, QUBOPT_VAR_NOT_USED double ** ar, MSGOUTPUT msgoutput){
	int i;
	char buffer[256];

	sprintf(buffer, "\titer: %d  logl: %.2f\r\n",iter,ll);
	msgoutput(buffer);

	for( i=0; i<namp; i++) {
		sprintf(buffer, "\t\tmulticlass %d: %.2f\t%.2f\r\n", i, amp[i], xms[i]);
		msgoutput(buffer);
		}
	}

//--------------------------------------------------------------------------
void mgamma(int nar, int nfir, int nmetastate, int ndata, double **alpha, double **beta, double *scale){
	int t0=nar+nfir;
	for( int i=0; i<nmetastate; i++)
		for( int t=t0; t<ndata; t++)
			beta[i][t] = alpha[i][t]*beta[i][t]/scale[t];
	}

//--------------------------------------------------------------------------
extern "C" void mfward(int nstate, int nar, int nfir, int nmetastate, 
			int **metastate, int ndata, double *pr, double **a,
			double **b, double **alpha, double *scale){
    int		I, J, i, j, k, n, t, t0 = nar + nfir;
    double  sum;

    for (scale[t0] = 0.0, J = 0; J < nmetastate; J++) {
		alpha[J][t0] = pr[J]*b[J][t0];
		scale[t0] += alpha[J][t0];
		}

    scale[t0] = 1.0/scale[t0];
	for (J = 0; J < nmetastate; J++) {
		alpha[J][t0] *= scale[t0];
		}
        
    for (t = t0 + 1; t < ndata; t++) {
		for (scale[t] = 0.0, J = 0; J < nmetastate; J++) {
			n = metastate[J][nar + nfir + 3 + nstate];
			j = metastate[J][0];
			for (sum = 0.0, k = 0; k < n; k++) {
				I = metastate[J][nar + nfir + 4 + nstate + k];
				i = metastate[I][0];
				sum += a[i][j]*alpha[I][t - 1];
				}
			alpha[J][t] = sum*b[J][t];
			scale[t] += alpha[J][t];
			}
		scale[t] = 1.0/scale[t];
		for (J = 0; J < nmetastate; J++) {
			alpha[J][t] *= scale[t];
			}
		}

	}

//--------------------------------------------------------------------------
extern "C" void mbward(QUBOPT_VAR_NOT_USED int nstate, int nar, int nfir, int nmetastate,
			int **metastate, int ndata, double **a, double **b,
			double *scale, double **beta){
    int		I, J, i, j, k, n, t, t0 = nar + nfir;
    double  sum;         

    for( I=0; I<nmetastate; I++ )
		beta[I][ndata-1] = scale[ndata-1];

    for( t=ndata-2; t>=t0; t-- ) {
		for( I=0; I<nmetastate; I++ ) {
			n = metastate[I][nar + nfir + 2];
			i = metastate[I][0];
			sum=0.0;
			for( k=0; k<n; k++ ) {
				J = metastate[I][nar + nfir + 3 + k];
				j = metastate[J][0];
				sum += a[i][j]*beta[J][t+1]*b[J][t+1];
				}
			beta[I][t] = sum*scale[t];
			}
		}
	}


//--------------------------------------------------------------------------
void 		mstateout(int nstate, int ndata, double **gama, int *state_out)
{
  if ( ! state_out )
    return;

  for (int d=0; d<ndata; ++d) {
    int iMax = 0;
    double gMax = gama[0][d];
    for (int i=1; i<nstate; ++i) {
      if ( gama[i][d] > gMax ) {
	gMax = gama[i][d];
	iMax = i;
      }
    }
    state_out[d] = iMax;
  }
}

//--------------------------------------------------------------------------
extern "C" void mbaum(int ndata, double *data, int nstate, double **a, int nmetastate, int **metastate, 
				 double *pr, int namp, double *amp, double *xms, int nar, double **ar, 
				 int nfir, double *fir, double *fret, int *state_out, MSGOUTPUT msgoutput, int *processflag) {
    int		it;     
	const int ITMAX = 100;
    double test, ftmp = 0;
	const double TOL = 1.0e-4;
     
    double ** alpha = dalloc2(nmetastate, ndata);
    double ** beta = dalloc2(nmetastate, ndata);
    double * scale = dalloc1(ndata); 
    double ** b	= dalloc2(nmetastate, ndata);   
    
    for( it=1; it<=ITMAX; it++, ftmp = *fret) { 
		if (*processflag != PROCESS_STARTED)
			break;

		mvectb(ndata, data, amp, xms, nar, ar, nfir, fir, nmetastate, metastate, b);
		mfward(nstate, nar, nfir, nmetastate, metastate, ndata, pr, a, b, alpha, scale);
		mbward(nstate, nar, nfir, nmetastate, metastate, ndata, a, b, scale, beta);
		*fret = mlogl(nar, nfir, ndata, scale);
		oupbaum(it, *fret, nstate, a, namp, amp, xms, nar, ar, msgoutput); 
        test = fabs(*fret - ftmp)/fabs(*fret);

        if (test < TOL) //fabs(*fret - ftmp) < 0.01
			break;

		mnewa(nstate, nar, nfir, nmetastate, metastate, ndata, alpha, beta, b, a);
		mgamma(nar, nfir, nmetastate, ndata, alpha, beta, scale);
		mnewpr(nar, nfir, nmetastate, beta, pr);
		mstateout(nstate, ndata, beta, state_out);

        //mnewamp(nar, nfir, nmetastate, metastate, ndata, data, beta, namp, amp);
        //mnewar(nar, nfir, nmetastate, metastate, ndata, data, beta, namp, amp, xms, ar, fir);
		mexasol(nar, nfir, nmetastate, metastate, ndata, data, beta, xms, ar, fir, namp, amp);
		}

    if( it>ITMAX ) 
		fprintf(stderr, "\tToo many iterations in baum\n");

	free2((char**)alpha);
	free2((char**)beta);
	free(scale);
	free2((char**)b);
	}

void mexasol(int nar, int nfir, int nmetastate, int **metastate, 
			int ndata, double *data, double **gama, double *xms, 
			double **ar, double *fir, int namp, double *amp){
    int		i, it, itmax = 10;
    double  sum1, sum2, test, tol = 1.0e-4;
    
    double * amp1 = dalloc1(namp);
    
    for( it=1; it<=itmax; it++) {
		mnewamp1(nar, nfir, nmetastate, metastate, ndata, data, gama, xms, ar, fir, namp, amp1);
		for (sum1 = sum2 = 0.0, i = 0; i < namp; i++) {
			sum1 += fabs(amp1[i] - amp[i]);
			sum2 += fabs(amp1[i] + amp[i]);
			amp[i] = amp1[i];
			}
		mnewar(nar, nfir, nmetastate, metastate, ndata, data, gama, namp, amp, xms, ar, fir);
		test = sum1 / sum2;
		if (test < tol) {
			free(amp1);
			return;
			}
		}
	free(amp1);
	}         


void mnewar(int nar, int nfir, int nmetastate, int **metastate,
			int ndata, double *data, double **gama, int namp,
			double *amp, double *xms, double **ar, double *fir) 
{       
    int		I, i, j, k, l, m, t, t0 = nar + nfir;
    double  sum, tmp1, tmp2, **r;
    
    r = dalloc2(nar + 2, nar + 2);

    for (m = 0; m < namp; m++) {
		for (k = 0; k <= nar; k++) {
			for (l = 0; l <= nar; l++) {
				r[k][l] = 0.0;
				}
			}

		for (sum = 0.0, I = 0; I < nmetastate; I++) {
			if (metastate[I][1] == m) {
				for (t = t0; t < ndata; t++) {
					for (k = 0; k <= nar; k++) {
						tmp1 = data[t - k];
						for (i = 0; i <= nfir; i++) {
							tmp1 -= fir[i]*amp[metastate[I][k + i + 1]];
							}
						for (l = 0; l <= nar; l++) {
							tmp2 = data[t - l];
							for (j = 0; j <= nfir; j++) {
								tmp2 -= fir[j]*amp[metastate[I][l + j + 1]];
								}
							r[k][l] += gama[I][t]*tmp1*tmp2;
							}
						}
					sum += gama[I][t];
					}
				}
			}

		for (k = 0; k <= nar; k++) {
			for (l = 0; l <= nar; l++) {
				r[k][l] /= sum;
				}
			}
		for (k = 0; k <= nar; k++) {
			r[k][nar + 1] =- r[k][0];
			if (k == 0) {
				r[k][0] =- 1;
			} else {
				r[k][0] = 0;
				}
			}
		gaussj(r, nar + 1);
		for (k = 1; k <= nar; k++) {
			ar[m][k] = r[k][nar + 1];
			}
		ar[m][0] = 1.0;
		xms[m] = sqrt(r[0][nar + 1]);
		}

	free2((char**)r);

	}


void mnewamp1(int nar, int nfir, int nmetastate, int **metastate, 
			  int ndata, double *data, double **gama, double *xms, 
			  double **ar, double *fir, int namp, double *amp)
{
    int		I, i, j, k, l, m, i1, j1, t, t0 = nar + nfir;
    double  tmp, sum, *r, *y, **a;                     
    
    r = dalloc1(nmetastate);
    y = dalloc1(nmetastate); 
    a = dalloc2(namp,namp + 1);
    
    for (I = 0; I < nmetastate; I++) {
		r[I] = y[I] = 0.0;
		m = metastate[I][1];
		for (t = t0; t < ndata; t++) {
			r[I] += gama[I][t];
			for (sum = 0.0, i = 0; i <= nar; i++) {
				sum += ar[m][i] * data[t - i];
				}
			y[I] += gama[I][t] * sum;
			}
		y[I] /= r[I];
		}
    
    for (m = 0; m < namp; m++) {
		for (l = 0; l <= namp; l++) {
			a[m][l] = 0.0;
			}
		for (i = 0; i <= nar; i++) {
			for (j = 0; j <= nfir; j++) {
				for (I = 0; I < nmetastate; I++) {
					if (metastate[I][i + j + 1] == m) {
						k = metastate[I][1];
						tmp = ar[k][i] * fir[j] / sqr(xms[k]) * r[I];
						a[m][namp] += tmp * y[I];
						for (i1 = 0; i1 <= nar; i1++) {
							for (j1 = 0; j1 <= nfir; j1++) {
								l = metastate[I][i1 + j1 + 1];
								a[m][l] += tmp * ar[k][i1] * fir[j1];
								}
							}
						}
					}
				}
			}
		}
    
    
	gaussj(a, namp);

    for (m = 0; m < namp; m++) {
		amp[m] = a[m][namp];
		}
    
    free(r);
    free(y);
    free2((char**)a);

	}
    

void mnewa(int nstate, int nar, int nfir, int nmetastate,
		   int **metastate, int ndata, double **alpha,
		   double **beta, double **b, double **a)
{       
    int		I, J, i, j, k, n, t, t0 = nar + nfir;
    double  chi, **a1, *a2;
    
    a1 = dalloc2(nstate, nstate);
    a2 = dalloc1(nstate);
    
    for (i = 0; i < nstate; i++) {
		for (j = 0; j < nstate; j++) {
			a1[i][j]=0.0;
			}
		a2[i] = 0.0;
		}             
    
    for (I = 0; I < nmetastate; I++) {

		i = metastate[I][0];
		n = metastate[I][nar + nfir + 2];
		for (k = 0; k < n; k++) {
			J = metastate[I][nar + nfir + 3 + k];
			j = metastate[J][0];
			for (t = t0; t < ndata - 1; t++) {
				chi = alpha[I][t]*a[i][j]*beta[J][t + 1]*b[J][t + 1];
				a1[i][j] += chi;
				a2[i] += chi;
				}
			}
		}
	
    for (i = 0; i < nstate; i++) {
		for (j = 0; j < nstate; j++) {
			a[i][j] = a1[i][j] / a2[i];
			}
		a[i][nstate] = a2[i] / (ndata - t0);
		}
    
    free2((char**)a1);
    free(a2);
	}

void mnewpr(int nar, int nfir, int nmetastate, double **gama, double *pr) {
    int	I, t0 = nar + nfir;
    for( I=0; I<nmetastate; I++ )
		pr[I] = gama[I][t0];
	}
    

/*
 * Solve linear equations "Ax=b" using the Gauss-Jordan elimination method. 
 * n is the dim of matrix a, and b is stored in the (n+1)-th column of a, 
 * i.e., b=a[][n]. On output, b is replaced by the solution.
 */

void gaussj(double **a, int n) {         
    int		*indxc, *indxr, *ipiv;
    int		icol, irow, i, j, k, l, ll;
    double	big, dum, pivinv, tmpa;

    indxc = ialloc1(n);
    indxr = ialloc1(n);
    ipiv  = ialloc1(n);
    
    for (j = 1; j <= n; j++) ipiv[j - 1] = 0;
    for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++) {
			if (ipiv[j - 1] != 1) {
				for (k = 1; k <= n; k++) {
					if (ipiv[k - 1] == 0) {
						if (fabs(a[j - 1][k - 1]) >= big) {
							big = fabs(a[j - 1][k - 1]);
							irow = j;
							icol = k;
							}
						}
					/*
					if (ipiv[k - 1] > 1) {
						//fprintf(stderr,"gauss-j: Singular Matrix-1\n");
						}
					*/
					}
				}
			}
		++(ipiv[icol - 1]);

		//-- Swap a[iRow-1][0..n] and a[iCol-1][0..n]
		if (irow != icol) {
			for( l=0; l<=n; l++) {
				tmpa = a[irow-1][l];
				a[irow-1][l] = a[icol-1][l];
				a[icol-1][l] = tmpa;
			}
		}

		indxr[i - 1] = irow;
		indxc[i - 1] = icol;
		/*
		if (a[icol - 1][icol - 1] == 0.0) {
			//fprintf(stderr,"gauss-j: Singular Matrix-2\n");
			}
		*/
		pivinv = 1.0 / a[icol - 1][icol - 1];
		a[icol - 1][icol - 1] = 1.0;
		for (l = 1; l <= n; l++) a[icol - 1][l - 1] *= pivinv;
		a[icol - 1][n] *= pivinv;
		for (ll = 1; ll <= n; ll++) {
			if (ll != icol) {
				dum = a[ll - 1][icol - 1];
				a[ll - 1][icol - 1] = 0.0;
				for (l = 1; l <= n; l++) 
					a[ll - 1][l - 1] -= a[icol - 1][l - 1]*dum;
				a[ll - 1][n] -= a[icol - 1][n]*dum;
				}
			}
		}
    
    for (l = n; l >= 1; l--) {
		if (indxr[l - 1] != indxc[l - 1]) {
			for (k = 1; k <= n; k++) {
				tmpa = a[k - 1][indxr[l - 1] - 1];
				a[k - 1][indxr[l - 1] - 1] = a[k - 1][indxc[l - 1] - 1];
				a[k - 1][indxc[l - 1] - 1] = tmpa;
				}
			}
		}
        
    free(ipiv);
    free(indxr);
    free(indxc);
	}

