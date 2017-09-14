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

/*
 * Calculate A=I*t^i0+Qt^(i0+1)/(i0+1)+Q^2*t^(i0+2)/(i0+2)+... 
 * and its derivatives.
 */

  
int nterm0(double **q, int n, double t){      
    int  i,j,m;
    double s=0.0,eps=1.0e-50;
           
    for (i=0; i<n; i++)   
	    for (j=0; j<n; j++) 
			s += t*fabs(q[i][j]);
	
	double r=1.0;
    for (m=0; r/(1.0-r/(m+2))>=eps; m++) 
		r*=s/(m+1);
	// Note : Qubopt version returns max(m,50) - only difference from this 
    return max(m, 30);
	}

extern "C" QUBOPT_API void mexp(int n, double t, double **q, double **a, int nz, 
					double ***dq, double ***da, int i0) {
    int  i,j,k,l,it,itmax;
    double  s,**u,**u1,***v,***v1;
    
    u  = dalloc2(n,n);
    u1 = dalloc2(n,n);
    v  = dalloc3(n,n,nz);
    v1 = dalloc3(n,n,nz);
    
	for (i=0; i<n; i++)    
		for (j=0; j<n; j++) {
			a[i][j]=0.0;
			u[i][j]=dirac(i,j)*pow(t,i0);
			for (k=0; k<nz; k++) 
				da[i][j][k]=v[i][j][k]=0.0;
			}
    
    itmax = nterm0(q,n,t);   
        
	for (it=0; it<=itmax; it++) {
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
				a[i][j] += u[i][j];
				for (k=0; k<nz; k++) 
					da[i][j][k]+=v[i][j][k];
				}
		for (i=0; i<n; i++)
			for (j=0; j<n; j++) { 
				for (s=0.,l=0; l<n; l++) 
					s+=u[i][l]*q[l][j]; 
				u1[i][j]=s; 
				} 
		for (k=0; k<nz; k++)  
			for (i=0; i<n; i++)         
				for (j=0; j<n; j++) {
					for (s=0.,l=0; l<n; l++) 
						s+=q[i][l]*v[l][j][k]+dq[i][l][k]*u[l][j];
					v1[i][j][k]=s;
					}
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
				u[i][j]=u1[i][j]*t/(it+1+i0);
				for (k=0; k<nz; k++) 
					v[i][j][k]=v1[i][j][k]*t/(it+1+i0);
				}
		}
        
	free2((char**)u);
	free2((char**)u1);
	free3((char***)v);
	free3((char***)v1);   
	}
