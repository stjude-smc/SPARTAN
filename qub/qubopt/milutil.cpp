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

#include "milutil.h"

#ifdef _WIN32
  #include <windows.h>
#else
  #include <stdlib.h>
  #define BOOL int
  #define TRUE 1
  #define FALSE 0
#endif

#include "milutil.h"

/*	--------------------------
	memory allocation routines

	One dimensional arrays are normal.  Two dimensional arrays are block allocated, 
then a 1D array of ptrs to the rows of the array is constructed and returned.  Three
dimensional arrays are similarly block allocated then returned as an array of ptrs
to an array of ptrs to rows.

	int *int_alloc1D(int size)								Allocate 1-D array of ints
	float *float_alloc1D(int size)							Allocate 1-D array of floats
	double *double_alloc1D(int size)						Allocate 1-D array of doubles
	double **double_alloc2D(int row, int col)				Allocate 2-D array of doubles
	double ***double_alloc3D(int row, int col, int depth)	Allocate 3-D array of doubles
	void free_2D(char **ptr)								Free 2-d array 
	void free_3D(char ***ptr)  								Free 3-d array

	void mxv(int m, int n, double **a, double *beta, double *alpha)
	void vxm(int m, int n, double **a, double *beta, double *alpha)
	int imaxv(int *d, int n)
	void mxm(int m, int n, int p, double **a, double **b, double **c)

	--------------------------
  
10/2004 - Unused functions commented out so we don't waste time trying to optimize them.
		- changed constructs like : 
				double **ptr = (double **)malloc(row*sizeof(double));	
		  to :	double **ptr = (double **)malloc(row*sizeof(double *));	
		  Need to check if the 2x allocation size was being used anywhere for padding.

void Message( char *szFormat, ... );
void Message( char *szFormat, ... ){
	char * pArgs;
	pArgs = (char *) &szFormat + sizeof( szFormat );
	
	char msg[ 2048 ];
	vsprintf( msg, szFormat, pArgs );
	}
		  
*/

#ifdef __GNUC__
#define MILUTIL_VAR_NOT_USED __attribute__ ((unused))
#else
#define MILUTIL_VAR_NOT_USED
#endif

void Message( MILUTIL_VAR_NOT_USED char * pchMessage ){
	// 9/05 - function was not doing anything...
	// left this shell as a place to put a breakpoint if desired.
	}

int *int_alloc1D(int size){
	int *ptr = (int *)malloc(size*sizeof(int));
	if( ptr == NULL ) 
		Message("Memory allocation error\n");
	return ptr;
	}

float *float_alloc1D(int size){
	float *ptr = (float *)malloc(size*sizeof(float));
	if( ptr == NULL ) 
		Message("Memory allocation error\n");
	return ptr;
	}

double *double_alloc1D(int size){
	double *ptr = (double *)malloc(size*sizeof(double));
	if( ptr == NULL ) 
		Message("Memory allocation error\n");
	return ptr;
	}

double **double_alloc2D(int row, int col){
	int i;
	double **ptr = (double **)malloc(row*sizeof(double *));
	double *ptr0 = (double *)malloc(row*col*sizeof(double));
	if( ptr==NULL || ptr0==NULL ) {
		free(ptr);
		Message("Memory allocation error\n");   
		return NULL;
		}  
	for(i=0;i<row;i++) 
		ptr[i]=&(ptr0[i*col]);
	return ptr;
	}

double ***double_alloc3D(int row, int col, int depth){
	int i,j;

	double ***ptr=(double ***)malloc(row*sizeof(double**));		// Should this be sizeof(double **) ?
	double **ptr0=(double **)malloc(row*col*sizeof(double*));	// Should this be sizeof(double *) ?
	double *ptr00=(double *)malloc(row*col*depth*sizeof(double));

	if( ptr==NULL || ptr0==NULL || ptr00==NULL ) { 
		free(ptr0);
		free(ptr);
		Message("Memory allocation error\n");  
		return NULL;
		}   

	for(i=0;i<row;i++) 
		ptr[i]=&(ptr0[i*col]);

	for(i=0;i<row;i++)
		for(j=0;j<col;j++) 
			ptr[i][j]=&(ptr00[(i*col+j)*depth]);

	return ptr;
	}

void free_2D(char **ptr){
	free(ptr[0]);
	free(ptr);
	}

void free_3D(char ***ptr){
	free(ptr[0][0]);
	free(ptr[0]);
	free(ptr);
	}

//----- Matrix * Matrix 
void mxm(int m, int n, int p, double **a, double **b, double **c){
	int	j,k,l;
	double s;
	
	for (k=0; k<m; k++) {
		for (l=0; l<p; l++) {
			s=0.0;
			for (j=0; j<n; j++) 
				s=s+a[k][j]*b[j][l];
			c[k][l]=s;
			}
		}
	}

/* Faster ?? Test - may provide hint to optimizer.
void mxm(int m, int n, int p, double **a, double **b, double **c){
	int	j,k,l;
	double s;
	double * a_k, * c_k;
	
	for( k=0; k<m; k++ ) {
		a_k_j=&(a[k][0]);
		c_k=c[k];
		for( l=0; l<p; l++ ) {
			s=0.0;
			for( j=0; j<n; j++ ) {
				s=s+(*a_k_j)*b[j][l];
				++a_k_j;
				}
			c_k[l]=s;
			}
		}
	}
*/

//----- Matrix * vector
void mxv(int m, int n, double **a, double *beta, double *alpha){
	int	j,k;
	double s;
	
	for (k=0; k<m; k++) {
		s=0.0;
		for (j=0; j<n; j++) 
			s=s+a[k][j]*beta[j];
		alpha[k]=s;
		}
	}

//----- Vector * Matrix 
void vxm(int m, int n, double **a, double *beta, double *alpha){
	int	j,k;
	double s;
	
	for (k=0; k<n; k++) {
		s=0.0;
		for (j=0; j<m; j++) 
			s=s+a[j][k]*beta[j];
		alpha[k]=s;
		}
	}

//----- max( element ) in integer vector 
int imaxv(int *d, int n){
	int	i,dm;
	dm=d[0];
	for (i=0; i<n; i++) 
		if (d[i]>dm) 
			dm=d[i];
		return dm;
	}

//----- Sum() double vector 
double dsumv(double *d, int n){
	int	i;
	double	s;
	for (s=0.0,i=0; i<n; i++) 
		s+=d[i];
	return s;
	}

//----- Initialize matrix to 0's 
void dzerom(int m, int n, double **a){
	int	i,j;
	for (i=0; i<m; i++) 
		for (j=0; j<n; j++) 
			a[i][j]=0.0;
	}

//----- Initialize vector to 0's 
void dzerov(int n, double *v){
	int	j;
	for (j=0; j<n; j++) 
		v[j]=0.0;
	//memset( v, 0, n*sizeof(double)
	}

int float_distance(float A, float B)
{
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    int aInt = *(int*)&A;
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;
    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;
    return (int) abs(aInt - bInt);
}

