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

void eclass(double *a, int n, int *nc, double *c, int *ia);

void nexcom0(int n, int k, int *r, int *mtc, int *t, int *h);


extern "C" QUBOPT_API void mutimakv(int nchannel, int nstate, int npath, 
			  int *state, int **path, double *amp, 
			  int *nmutistate, int *nmutipath, int *nmuticlass, 
			  int ***mutistate, int ***mutipath){
	int	i, k, l, m, n, s, t, iContinue, iStateIndex, *z, *ia, nn;
    double tmp, *a, *c;

	//model size for memory allocation
	tmp = 0.0;
    for( n=1; n<nstate; n++)
		tmp += log(double(n+nchannel)) - log(double(n));

	nn = (int)(exp(tmp) + 1);
    *mutistate = ialloc2(nn, nstate + 1);

	//muti-channel states
	int temp;
    z = ialloc1(nstate + 1);
    *nmutistate = iContinue = 0; 
    do {
		nexcom0(nchannel, nstate, z, &iContinue, &temp, &iStateIndex);
		for (i = 0; i < nstate; i++)
			(*mutistate)[*nmutistate][i + 1] = z[i + 1];
		(*nmutistate)++;
		} while (iContinue == 1);
    free(z);
    
	//conductance classes
    a  = dalloc1(*nmutistate);
    ia = ialloc1(*nmutistate);
    c  = dalloc1(*nmutistate);
        
    for (n = 0; n < *nmutistate; n++) {
		for (a[n] = 0., i = 0; i < nstate; i++) {
			k = state[i];
			a[n] += (*mutistate)[n][i + 1]*amp[k];
			}
		}

    eclass(a, *nmutistate, nmuticlass, c, ia);
    for( n=0; n<*nmutistate; n++ )
		(*mutistate)[n][0] = ia[n];
   
	//number of paths for memory allocation
	n=0;
	for( k=0; k<*nmutistate; k++ ) {
		for( l=0; l<*nmutistate; l++ ) {
			s=0;
			for( i=0; i<nstate; i++)
				s += abs((*mutistate)[k][i + 1] - (*mutistate)[l][i + 1]);
			if (s == 2)
				n++;
			}
		}
    
	*mutipath = ialloc2(n, 4); 
    
	//transition paths        
    for (*nmutipath = m = 0; m < *nmutistate; m++) {
		for( n=0; n<*nmutistate; n++) {
			s=0;
			for( i=0; i<nstate; i++)
				s += abs((*mutistate)[m][i + 1] - (*mutistate)[n][i + 1]);
			
			if (s == 2) {
				for( i=0; i<nstate; i++) {
					t = (*mutistate)[m][i + 1] - (*mutistate)[n][i + 1];
					if (t == 1) 
						k = i;
					if (t == -1) 
						l = i;
					}	
				
				for( i=0; i<npath; i++) {
					if (path[i][0] == k && path[i][1] == l) {
						(*mutipath)[*nmutipath][0] = m;
						(*mutipath)[*nmutipath][1] = n;
						(*mutipath)[*nmutipath][2] = k;
						(*mutipath)[*nmutipath][3] = l;
						(*nmutipath)++;
						}
					}
				}
			}
		}
    
    free(a);
    free(ia);  
    free(c);
	}

// See also CChannelCombine in QubOpt for description 
void nexcom0(int iChannels, int iStates, int *paStates, int *piContinue, int *pTemp, int *piStateIndex) {
	if( *piContinue==0 ) {
		paStates[1] = iChannels;
        *pTemp = iChannels;
        *piStateIndex = 0;
        for( int iState=2; iState<=iStates; iState++) 
			paStates[iState]=0;
		}
	else {
		if(*pTemp > 1) 
			*piStateIndex=0;
		(*piStateIndex)++;
		*pTemp = paStates[*piStateIndex];
		paStates[*piStateIndex] = 0;
		paStates[1] = *pTemp - 1;
		paStates[*piStateIndex + 1] = paStates[*piStateIndex + 1] + 1;
		} 
    
	*piContinue=(paStates[iStates] != iChannels) ? 1 : 0 ; 
	}

/* DZERO is the largest number in double precision which when added  to one gives one. */  
/* FZERO is the largest number in single precision which when added  to one gives one. */      
// #define FZERO   1.0e-8
// For notes - see eclass description in QubOpt 
#define DZERO   1.0e-15 

void eclass(double *a, int n, int *nc, double *c, int *ia) {
    int  i,k,m,nl,*il,*ii;
    double  amin;
	double openDirection = (n > 1) ? (a[1] - a[0]) : 1.0;
	
    il = ialloc1(n);
    ii = ialloc1(n);
	
    for (i=0; i<n; i++) il[i]=i;
    nl=n;
    *nc=0;
    while (nl>0) {
		k=il[0];
		amin=a[k];
		for (i=1; i<nl; i++) {
			k=il[i];
			if ( ((openDirection > 0.0) && (a[k]-amin<DZERO)) || ((openDirection < 0.0) && (amin-a[k]<DZERO)) )
				amin=a[k];
			}
		c[*nc]=amin;
		for (m=i=0; i<nl; i++) {
			k=il[i]; 
			if (fabs(a[k]-amin)<DZERO) ia[k]=*nc;
			else {
				ii[m]=k;
				m++;
				}
			}
		for (i=0; i<m; i++) il[i]=ii[i];
		nl=m;
		(*nc)++;
		}   
	
    free(il);
    free(ii);
	}


