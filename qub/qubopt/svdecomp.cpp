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


void svdecomp0(double **a, int m, int n, double *w, double **v) {
    int  flag,i,its,j,jj,k,l,nm;
    double  anorm,c,f,g,h,s,scale,x,y,z,*rv1,ww,*aa,*vv;
    const int MAXITER = 30;
    
    rv1 = dalloc1(n);
    aa  = dalloc1(m);
    vv  = dalloc1(n);
    
    g=scale=anorm=0.0;
    for (i=1; i<=n; i++) {
		l=i+1;
		rv1[i-1]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i; k<=m; k++) scale += fabs(a[k-1][i-1]);
			if (scale) {
				for (k=i; k<=m; k++) {
					a[k-1][i-1] /= scale;
					s += a[k-1][i-1]*a[k-1][i-1];
					}
				f=a[i-1][i-1];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i-1][i-1]=f-g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=i ;k<=m; k++) 
						s += a[k-1][i-1]*a[k-1][j-1];
					f=s/h;
					for (k=i; k<=m; k++) 
						a[k-1][j-1] += f*a[k-1][i-1];
					}
				for (k=i; k<=m; k++) a[k-1][i-1] *= scale;
				}
			}
		w[i-1]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l; k<=n; k++) 
				scale += fabs(a[i-1][k-1]);
			if (scale) {
				for (k=l; k<=n; k++) {
					a[i-1][k-1] /= scale;
					s += a[i-1][k-1]*a[i-1][k-1];
					}
				f=a[i-1][l-1];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i-1][l-1]=f-g;
				for (k=l; k<=n; k++) 
					rv1[k-1]=a[i-1][k-1]/h;
				for (j=l; j<=m; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[j-1][k-1]*a[i-1][k-1];
					for (k=l; k<=n; k++) 
						a[j-1][k-1] += s*rv1[k-1];
					}
				for (k=l; k<=n; k++) 
					a[i-1][k-1] *= scale;
				}
			}
		anorm=max(anorm,(fabs(w[i-1])+fabs(rv1[i-1])));
		}
    for (i=n; i>=1; i--) {
		if (i < n) {
			if (g) {
				for (j=l; j<=n; j++) 
					v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[i-1][k-1]*v[k-1][j-1];
					for (k=l; k<=n; k++) 
						v[k-1][j-1] += s*v[k-1][i-1];
					}
				}
			for (j=l; j<=n; j++) 
				v[i-1][j-1]=v[j-1][i-1]=0.0;
			}
		v[i-1][i-1]=1.0;
		g=rv1[i-1];
		l=i;
		}
    for (i=min(m,n); i>=1; i--) {
		l=i+1;
		g=w[i-1];
		for (j=l; j<=n; j++) 
			a[i-1][j-1]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l; j<=n; j++) {
				for (s=0.0,k=l; k<=m; k++) 
					s += a[k-1][i-1]*a[k-1][j-1];
				f=(s/a[i-1][i-1])*g;
				for (k=i; k<=m; k++) 
					a[k-1][j-1] += f*a[k-1][i-1];
				}
			for (j=i; j<=m; j++) 
				a[j-1][i-1] *= g;
			} 
		else {
			for (j=i; j<=m; j++) a[j-1][i-1]=0.0;
			}
		++a[i-1][i-1];
		}
    for (k=n; k>=1; k--) {
		for (its=1; its<=MAXITER; its++) {
			flag=1;
			for (l=k; l>=1; l--) {
				nm=l-1;
				if (fabs(rv1[l-1])+anorm == anorm) {
					flag=0;
					break;
					}
				if (fabs(w[nm-1])+anorm == anorm) 
					break;
				}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l; i<=k; i++) {
					f=s*rv1[i-1];
					rv1[i-1]=c*rv1[i-1];
					if (fabs(f)+anorm == anorm) 
						break;
					g=w[i-1];
					h=pythag(f,g);
					w[i-1]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1; j<=m; j++) {
						y=a[j-1][nm-1];
						z=a[j-1][i-1];
						a[j-1][nm-1]=y*c+z*s;
						a[j-1][i-1]=z*c-y*s;
						}
					}
				}
			z=w[k-1];
			if (l == k) {
				if (z < 0.0) {
					w[k-1] = -z;
					for (j=1; j<=n; j++) 
						v[j-1][k-1] = -v[j-1][k-1];
					}
				break;
				}
			if (its==MAXITER)
				fprintf(stderr,"No convergence in svdecomp0\n");
			x=w[l-1];
			nm=k-1;
			y=w[nm-1];
			g=rv1[nm-1];
			h=rv1[k-1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c=s=1.0;
			for (j=l; j<=nm; j++) {
				i=j+1;
				g=rv1[i-1];
				y=w[i-1];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j-1]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1; jj<=n; jj++) {
					x=v[jj-1][j-1];
					z=v[jj-1][i-1];
					v[jj-1][j-1]=x*c+z*s;
					v[jj-1][i-1]=z*c-x*s;
					}
				z=pythag(f,h);
				w[j-1]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
					}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1; jj<=m; jj++) {
					y=a[jj-1][j-1];
					z=a[jj-1][i-1];
					a[jj-1][j-1]=y*c+z*s;
					a[jj-1][i-1]=z*c-y*s;
					}
				}
			rv1[l-1]=0.0;
			rv1[k-1]=f;
			w[k-1]=x;
			}
		}
	
 	// Reorders columns of a,v matrices and w vector ordered by abs(w) descending.
	// This appears to be an inefficient swap sort but check more thoroughly.
	// More efficient : Quick sort the W array and a parallel vector of column 
	// destinations then move the a and v columns efficiently.
   for (j=2; j<=n; j++) {
		ww=w[j-1];
		for (i=1; i<=m; i++) 
			aa[i-1]=a[i-1][j-1];
		for (i=1; i<=n; i++) 
			vv[i-1]=v[i-1][j-1];
		k=j-1;
		while ( k > 0 && -fabs(w[k-1])>-fabs(ww) ) {
			w[k+1-1]=w[k-1];
			for (i=1; i<=m; i++) 
				a[i-1][k+1-1]=a[i-1][k-1];
			for (i=1; i<=n; i++) 
				v[i-1][k+1-1]=v[i-1][k-1];
			k--;
			}
		w[k+1-1]=ww;
		for (i=1; i<=m; i++) 
			a[i-1][k+1-1]=aa[i-1];
		for (i=1; i<=n; i++) 
			v[i-1][k+1-1]=vv[i-1];
		}
    
    free(rv1);
    free(aa);
    free(vv);

	}             

void svdecomp0_mem(double **a, int m, int n, double *w, double **v, double *rv1/*n*/, double *aa/*m*/, double *vv/*n*/) {
    int  flag,i,its,j,jj,k,l,nm;
    double  anorm,c,f,g,h,s,scale,x,y,z,ww;

    g=scale=anorm=0.0;
    for (i=1; i<=n; i++) {
		l=i+1;
		rv1[i-1]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i; k<=m; k++) scale += fabs(a[k-1][i-1]);
			if (scale) {
				for (k=i; k<=m; k++) {
					a[k-1][i-1] /= scale;
					s += a[k-1][i-1]*a[k-1][i-1];
					}
				f=a[i-1][i-1];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i-1][i-1]=f-g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=i ;k<=m; k++) 
						s += a[k-1][i-1]*a[k-1][j-1];
					f=s/h;
					for (k=i; k<=m; k++) 
						a[k-1][j-1] += f*a[k-1][i-1];
					}
				for (k=i; k<=m; k++) a[k-1][i-1] *= scale;
				}
			}
		w[i-1]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l; k<=n; k++) 
				scale += fabs(a[i-1][k-1]);
			if (scale) {
				for (k=l; k<=n; k++) {
					a[i-1][k-1] /= scale;
					s += a[i-1][k-1]*a[i-1][k-1];
					}
				f=a[i-1][l-1];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i-1][l-1]=f-g;
				for (k=l; k<=n; k++) 
					rv1[k-1]=a[i-1][k-1]/h;
				for (j=l; j<=m; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[j-1][k-1]*a[i-1][k-1];
					for (k=l; k<=n; k++)
						a[j-1][k-1] += s*rv1[k-1];
					}
				for (k=l; k<=n; k++) 
					a[i-1][k-1] *= scale;
				}
			}
		anorm=max(anorm,(fabs(w[i-1])+fabs(rv1[i-1])));
		}
    for (i=n; i>=1; i--) {
		if (i < n) {
			if (g) {
				for (j=l; j<=n; j++) 
					v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[i-1][k-1]*v[k-1][j-1];
					for (k=l; k<=n; k++) 
						v[k-1][j-1] += s*v[k-1][i-1];
					}
				}
			for (j=l; j<=n; j++) 
				v[i-1][j-1]=v[j-1][i-1]=0.0;
			}
		v[i-1][i-1]=1.0;
		g=rv1[i-1];
		l=i;
		}
    for (i=min(m,n); i>=1; i--) {
		l=i+1;
		g=w[i-1];
		for (j=l; j<=n; j++) 
			a[i-1][j-1]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l; j<=n; j++) {
				for (s=0.0,k=l; k<=m; k++) 
					s += a[k-1][i-1]*a[k-1][j-1];
				f=(s/a[i-1][i-1])*g;
				for (k=i; k<=m; k++) 
					a[k-1][j-1] += f*a[k-1][i-1];
				}
			for (j=i; j<=m; j++)
				a[j-1][i-1] *= g;
			} 
		else {
			for (j=i; j<=m; j++) a[j-1][i-1]=0.0;
			}
		++a[i-1][i-1];
		}
    for (k=n; k>=1; k--) {
		for (its=1; its<=30; its++) {
			flag=1;
			for (l=k; l>=1; l--) {
				nm=l-1;
				if (fabs(rv1[l-1])+anorm == anorm) {
					flag=0;
					break;
					}
				if (fabs(w[nm-1])+anorm == anorm) 
					break;
				}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l; i<=k; i++) {
					f=s*rv1[i-1];
					rv1[i-1]=c*rv1[i-1];
					if (fabs(f)+anorm == anorm) 
						break;
					g=w[i-1];
					h=pythag(f,g);
					w[i-1]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1; j<=m; j++) {
						y=a[j-1][nm-1];
						z=a[j-1][i-1];
						a[j-1][nm-1]=y*c+z*s;
						a[j-1][i-1]=z*c-y*s;
						}
					}
				}
			z=w[k-1];
			if (l == k) {
				if (z < 0.0) {
					w[k-1] = -z;
					for (j=1; j<=n; j++) 
						v[j-1][k-1] = -v[j-1][k-1];
					}
				break;
				}
			if (its==30) 
				fprintf(stderr,"No convergence in svdecomp0_mem\n");
			x=w[l-1];
			nm=k-1;
			y=w[nm-1];
			g=rv1[nm-1];
			h=rv1[k-1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c=s=1.0;
			for (j=l; j<=nm; j++) {
				i=j+1;
				g=rv1[i-1];
				y=w[i-1];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j-1]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1; jj<=n; jj++) {
					x=v[jj-1][j-1];
					z=v[jj-1][i-1];
					v[jj-1][j-1]=x*c+z*s;
					v[jj-1][i-1]=z*c-x*s;
					}
				z=pythag(f,h);
				w[j-1]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
					}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1; jj<=m; jj++) {
					y=a[jj-1][j-1];
					z=a[jj-1][i-1];
					a[jj-1][j-1]=y*c+z*s;
					a[jj-1][i-1]=z*c-y*s;
					}
				}
			rv1[l-1]=0.0;
			rv1[k-1]=f;
			w[k-1]=x;
			}
		}
	
    for (j=2; j<=n; j++) {
		ww=w[j-1];
		for (i=1; i<=m; i++) 
			aa[i-1]=a[i-1][j-1];
		for (i=1; i<=n; i++) 
			vv[i-1]=v[i-1][j-1];
		k=j-1;
		while ( k > 0 && -fabs(w[k-1])>-fabs(ww) ) {
			w[k+1-1]=w[k-1];
			for (i=1; i<=m; i++)
				a[i-1][k+1-1]=a[i-1][k-1];
			for (i=1; i<=n; i++)
				v[i-1][k+1-1]=v[i-1][k-1];
			k--;
			}
		w[k+1-1]=ww;
		for (i=1; i<=m; i++)
			a[i-1][k+1-1]=aa[i-1];
		for (i=1; i<=n; i++)
			v[i-1][k+1-1]=vv[i-1];
		}
    
}

void svdecomp0_x_mem(long double **a, int m, int n, long double *w, long double **v, long double *rv1/*n*/, long double *aa/*m*/, long double *vv/*n*/) {
    int  flag,i,its,j,jj,k,l,nm;
    long double  anorm,c,f,g,h,s,scale,x,y,z,ww;

    g=scale=anorm=0.0;
    for (i=1; i<=n; i++) {
		l=i+1;
		rv1[i-1]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i; k<=m; k++) scale += fabsl(a[k-1][i-1]);
			if (scale) {
				for (k=i; k<=m; k++) {
					a[k-1][i-1] /= scale;
					s += a[k-1][i-1]*a[k-1][i-1];
					}
				f=a[i-1][i-1];
				g = -signl(sqrtl(s),f);
				h=f*g-s;
				a[i-1][i-1]=f-g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=i ;k<=m; k++) 
						s += a[k-1][i-1]*a[k-1][j-1];
					f=s/h;
					for (k=i; k<=m; k++)
						a[k-1][j-1] += f*a[k-1][i-1];
					}
				for (k=i; k<=m; k++) a[k-1][i-1] *= scale;
				}
			}
		w[i-1]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l; k<=n; k++) 
				scale += fabsl(a[i-1][k-1]);
			if (scale) {
				for (k=l; k<=n; k++) {
					a[i-1][k-1] /= scale;
					s += a[i-1][k-1]*a[i-1][k-1];
					}
				f=a[i-1][l-1];
				g = -signl(sqrtl(s),f);
				h=f*g-s;
				a[i-1][l-1]=f-g;
				for (k=l; k<=n; k++) 
					rv1[k-1]=a[i-1][k-1]/h;
				for (j=l; j<=m; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[j-1][k-1]*a[i-1][k-1];
					for (k=l; k<=n; k++)
						a[j-1][k-1] += s*rv1[k-1];
					}
				for (k=l; k<=n; k++) 
					a[i-1][k-1] *= scale;
				}
			}
		anorm=max(anorm,(fabsl(w[i-1])+fabsl(rv1[i-1])));
		}
    for (i=n; i>=1; i--) {
		if (i < n) {
			if (g) {
				for (j=l; j<=n; j++) 
					v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=l; k<=n; k++) 
						s += a[i-1][k-1]*v[k-1][j-1];
					for (k=l; k<=n; k++)
						v[k-1][j-1] += s*v[k-1][i-1];
					}
				}
			for (j=l; j<=n; j++) 
				v[i-1][j-1]=v[j-1][i-1]=0.0;
			}
		v[i-1][i-1]=1.0;
		g=rv1[i-1];
		l=i;
		}
    for (i=min(m,n); i>=1; i--) {
		l=i+1;
		g=w[i-1];
		for (j=l; j<=n; j++) 
			a[i-1][j-1]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l; j<=n; j++) {
				for (s=0.0,k=l; k<=m; k++) 
					s += a[k-1][i-1]*a[k-1][j-1];
				f=(s/a[i-1][i-1])*g;
				for (k=i; k<=m; k++) 
					a[k-1][j-1] += f*a[k-1][i-1];
				}
			for (j=i; j<=m; j++)
				a[j-1][i-1] *= g;
			} 
		else {
			for (j=i; j<=m; j++) a[j-1][i-1]=0.0;
			}
		++a[i-1][i-1];
		}
    for (k=n; k>=1; k--) {
		for (its=1; its<=30; its++) {
			flag=1;
			for (l=k; l>=1; l--) {
				nm=l-1;
				if (fabsl(rv1[l-1])+anorm == anorm) {
					flag=0;
					break;
					}
				if (fabsl(w[nm-1])+anorm == anorm)
					break;
				}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l; i<=k; i++) {
					f=s*rv1[i-1];
					rv1[i-1]=c*rv1[i-1];
					if (fabsl(f)+anorm == anorm) 
						break;
					g=w[i-1];
					h=pythagl(f,g);
					w[i-1]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1; j<=m; j++) {
						y=a[j-1][nm-1];
						z=a[j-1][i-1];
						a[j-1][nm-1]=y*c+z*s;
						a[j-1][i-1]=z*c-y*s;
						}
					}
				}
			z=w[k-1];
			if (l == k) {
				if (z < 0.0) {
					w[k-1] = -z;
					for (j=1; j<=n; j++) 
						v[j-1][k-1] = -v[j-1][k-1];
					}
				break;
				}
			if (its==30) 
				fprintf(stderr,"No convergence in svdecomp0_xmem\n");
			x=w[l-1];
			nm=k-1;
			y=w[nm-1];
			g=rv1[nm-1];
			h=rv1[k-1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythagl(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+signl(g,f)))-h))/x;
			c=s=1.0;
			for (j=l; j<=nm; j++) {
				i=j+1;
				g=rv1[i-1];
				y=w[i-1];
				h=s*g;
				g=c*g;
				z=pythagl(f,h);
				rv1[j-1]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1; jj<=n; jj++) {
					x=v[jj-1][j-1];
					z=v[jj-1][i-1];
					v[jj-1][j-1]=x*c+z*s;
					v[jj-1][i-1]=z*c-x*s;
					}
				z=pythagl(f,h);
				w[j-1]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
					}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1; jj<=m; jj++) {
					y=a[jj-1][j-1];
					z=a[jj-1][i-1];
					a[jj-1][j-1]=y*c+z*s;
					a[jj-1][i-1]=z*c-y*s;
					}
				}
			rv1[l-1]=0.0;
			rv1[k-1]=f;
			w[k-1]=x;
			}
		}
	
    for (j=2; j<=n; j++) {
		ww=w[j-1];
		for (i=1; i<=m; i++) 
			aa[i-1]=a[i-1][j-1];
		for (i=1; i<=n; i++) 
			vv[i-1]=v[i-1][j-1];
		k=j-1;
		while ( k > 0 && -fabsl(w[k-1])>-fabsl(ww) ) {
			w[k+1-1]=w[k-1];
			for (i=1; i<=m; i++) 
				a[i-1][k+1-1]=a[i-1][k-1];
			for (i=1; i<=n; i++) 
				v[i-1][k+1-1]=v[i-1][k-1];
			k--;
			}
		w[k+1-1]=ww;
		for (i=1; i<=m; i++) 
			a[i-1][k+1-1]=aa[i-1];
		for (i=1; i<=n; i++) 
			v[i-1][k+1-1]=vv[i-1];
		}
    
}




