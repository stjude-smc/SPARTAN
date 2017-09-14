/* Copyright 1998-2014 Research Foundation State University of New York */

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

#ifndef MILUTIL_H
#define MILUTIL_H

inline double DSQR( double a ) {	return (a==0.0?0.0:a*a);	}
inline double DMAX( double a, double b ) {	return (a>b?a:b);	}
inline double DMIN( double a, double b ) {	return (a<b?a:b);	}	// UNUSED
inline int IMIN( int a, int b ) {	return (a<b?a:b);	}
inline int IMAX( int a, int b ) {	return (a>b?a:b);	}

int *int_alloc1D(int size);
float *float_alloc1D(int size);
double *double_alloc1D(int size);
double **double_alloc2D(int row, int col);
double ***double_alloc3D(int row, int col, int depth);
void free_2D(char **ptr);
void free_3D(char ***ptr);
void mxm(int m, int n, int p, double **a, double **b, double **c);
void mxv(int m, int n, double **a, double *beta, double *alpha);
void vxm(int m, int n, double **a, double *beta, double *alpha);
int imaxv(int *d, int n);
double dsumv(double *d, int n);
void dzerom(int m, int n, double **a);
void dzerov(int n, double *v);

int float_distance(float A, float B);

#endif
