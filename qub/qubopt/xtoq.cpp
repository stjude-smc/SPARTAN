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

// rate constants Q[0..nstate][0..nstate] and dq[nstate][nstate][nz]
void xtoq(int nstate, int npath, int **path, double *C, 
					double V, double *x, double **q, int nz, double **dx, double ***dq) {
    int i, j, k, m;
	double nQ;
     
	//First zero the Q and dq matrices
    for( i=0; i<nstate; i++) {
		for( j=0; j<nstate; j++) {
			q[i][j] = 0.0;
		    for( k=0; k<nz; k++) 
				dq[i][j][k] = 0.0;
			}
		}

	//Calculate the rates for the connected states and the diagonal
	for( m=0; m<npath; m++ ) {
		i = path[m][0];	//from state 
		j = path[m][1];	//to state
		k = path[m][2];	//drug index 
		nQ = ((k==-1) ? 1 : C[k]) * exp(x[2*m] + V*x[2*m + 1]);

		//the q matrix - the rate expression.  diagonal = -sum(row)
		q[i][i] -=  ( q[i][j] = nQ );	
		// the dq matrix 
		for( k=0; k<nz; k++)
			dq[i][i][k] -=  ( dq[i][j][k] = nQ *(dx[2*m][k] + V*dx[2*m + 1][k]) );
		}
	}
