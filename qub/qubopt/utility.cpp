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


#define EXIT_ERROR /*fprintf(stderr,"Memory allocation error\n"); exit(1);*/return NULL;

void swapn(double *a, double *b, int iItems){
	double tmp;
	for( int i=0; i<iItems; ++i ) {
		tmp=a[i];
		a[i]=b[i];
		b[i]=tmp;
		}
	}

void swapnl(long double *a, long double *b, int iItems){
	long double tmp;
	for( int i=0; i<iItems; ++i ) {
		tmp=a[i];
		a[i]=b[i];
		b[i]=tmp;
		}
	}

void swap(double * p1, double *p2 ){
	double n;
	n=*p1;
	*p1=*p2;
	*p2=n;
	}

void swapl(long double * p1, long double *p2 ){
	long double n;
	n=*p1;
	*p1=*p2;
	*p2=n;
	}

/* standard error handler */
void errexit(char *msg) {
    fprintf(stderr,"\n%s\n",msg);
    exit(1);
	}

/* memory allocation */
int *ialloc1(int size) {
    int * ptr = (int *)malloc(size*sizeof(int));
    if (ptr==NULL)
		EXIT_ERROR
    return ptr;
	}

double* dalloc1(int size) {
    double * ptr = (double *)malloc(size*sizeof(double));
    if (ptr==NULL)
       EXIT_ERROR
    return ptr;
	}

int **ialloc2(int row, int col) {
    int i;
    int **ptr = (int **)malloc(row*sizeof(int *));
    if (ptr==NULL)
       EXIT_ERROR
      
    ptr[0] = (int *)malloc(row*col*sizeof(int));
    if (ptr[0]==NULL)
       EXIT_ERROR
     
    for( i=1; i<row; i++) 
		ptr[i]=&(ptr[0][i*col]);
    return ptr;
	}

long double ** xalloc2(int row, int col) {
	int i;
	long double **ptr = (long double **)malloc(row*sizeof(long double *));
	if (ptr == NULL)
		EXIT_ERROR
	ptr[0] = (long double *)malloc(row*col*sizeof(long double));
	if (ptr[0] == NULL)
		EXIT_ERROR
	for (i=1; i<row; i++)
		ptr[i] = &(ptr[0][i*col]);
	return ptr;
}

double  ** dalloc2(int row, int col) {
    int i;
    double **ptr = (double **)malloc(row*sizeof(double *));
    if (ptr==NULL)
       EXIT_ERROR
    
    ptr[0] = (double *)malloc(row*col*sizeof(double));
    if (ptr[0]==NULL)
       EXIT_ERROR
    
    for (i=1; i<row; i++) 
		ptr[i]=&(ptr[0][i*col]);
    return ptr;
	}
        
void free2(char **ptr) {
    free(ptr[0]);
    free(ptr);
	}
       
double  *** dalloc3(int row, int col, int depth) {
    int i,j;
    double *** ptr = (double ***)malloc(row*sizeof(double **));
    if (ptr==NULL)
       EXIT_ERROR
    
    ptr[0] = (double **)malloc(row*col*sizeof(double *));
    if (ptr[0]==NULL)
       EXIT_ERROR
    
    for (i=1; i<row; i++) 
		ptr[i]=&(ptr[0][i*col]);
    
    ptr[0][0] = (double *)malloc(row*col*depth*sizeof(double));
    if (ptr[0][0]==NULL)
       EXIT_ERROR
    
    for (i=0; i<row; i++)
		for (j=0; j<col; j++) 
			ptr[i][j]=&(ptr[0][0][(i*col+j)*depth]);
    return ptr;
	}

void free3(char ***ptr) {
    free(ptr[0][0]);   
    free(ptr[0]);
    free(ptr);
	}

#define DSQR(a) ((a)==0.0 ? 0.0 : (a)*(a))

double pythag(double a, double b){
	a=fabs(a);
	b=fabs(b);
	if( a>b )
		return a*sqrt(1.0+DSQR(b/a));
	if ( b==0.0 )		
		return 0.0;
	return b*sqrt(1.0+DSQR(a/b));
	}

long double pythagl(long double a, long double b){
	a=fabsl(a);
	b=fabsl(b);
	if( a>b )
		return a*sqrtl(1.0+DSQR(b/a));
	if ( b==0.0 )
		return 0.0;
	return b*sqrtl(1.0+DSQR(a/b));
	}

#undef EXIT_ERROR

