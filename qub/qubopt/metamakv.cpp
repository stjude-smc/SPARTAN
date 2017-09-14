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

extern "C" QUBOPT_API void metamakv(int nstate, int namp, int **state, int nar, 
			  int *nmetastate, int ***metastate){
    int  I, J, i, j, k, l, n, s;
    
	//model size for memory allocation 
    k = nstate*(int)pow((double)namp, nar);   
    l = (nar + 3 + nstate*(namp + 1)) + 1; // +1 because it wasn't enough!
    *metastate = ialloc2(k, l);

	//meta-states      
	*nmetastate = 0;
    for( j=0; j<nstate; j++) {
		(*metastate)[*nmetastate][0] = j;
		for( k=1; k<=nar; k++)
			(*metastate)[*nmetastate][k] = j;
		(*nmetastate)++;
		}   

    n = *nmetastate;
    for( l=1; l<=nar; l++ ) {
		for( i=0; i<*nmetastate; i++ ) {
			for( j=1; j<namp; j++ ) {
				for( k=0; k<=nar; k++ ) {
					if (k == l)
						(*metastate)[n][k] = j;
					else
						(*metastate)[n][k]=(*metastate)[i][k];
					}
				n++;
				}
			}
		*nmetastate = n;
    	}

	//set the 2nd column to be conductance index
	for( n=0; n<*nmetastate; n++ ) {
		for( k=nar+1; k>=1; k-- )
			(*metastate)[n][k] = (*metastate)[n][k-1];
		i = (*metastate)[n][1];
		(*metastate)[n][1] = state[i][0];
		}

	//allowable transitions leaving meta-state I 
    for( I=0; I<*nmetastate; I++ ) {
		n=0;
		for( J=0; J<*nmetastate; J++) {
			s=0;
			for( k=2; k<=nar+1; k++)	// TODO - instead test expr. and break since s is not needed.
				s += abs((*metastate)[J][k] - (*metastate)[I][k - 1]);
			if( s==0 )
				(*metastate)[I][nar+3+(n++)] = J;
			}
		(*metastate)[I][nar + 2] = n;
    	}

	//allowable transitions entering meta-state J
    for( J=0; J<*nmetastate; J++ ) {
		n=0;
		for( I=0; I<*nmetastate; I++ ) {
			s=0;
			for( k=2; k<=nar+1; k++)
				s += abs((*metastate)[J][k] - (*metastate)[I][k - 1]);
			if (s == 0)
				(*metastate)[J][nar+4+nstate+(n++)] = I;
			}
		(*metastate)[J][nar+3+nstate] = n;
		}
	}
