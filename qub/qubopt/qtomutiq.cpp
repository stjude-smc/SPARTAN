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

extern "C" QUBOPT_API void qtomutiq(QUBOPT_VAR_NOT_USED int nstate, QUBOPT_VAR_NOT_USED int npath, QUBOPT_VAR_NOT_USED int ** path, int nmutistate,
			  int **mutistate, int nmutipath, int **mutipath, double **q,
			  double **qq, int nz, double ***dq, double ***dqq){   
    int	I, J, i, j, k, m;
    double s;                
    
 	//Rate constants
    for( I=0; I<nmutistate; I++)
		for( J=0; J<nmutistate; J++)
			qq[I][J]=0.0;

    for( m=0; m<nmutipath; m++) {
		I = mutipath[m][0];
		J = mutipath[m][1];
		i = mutipath[m][2];
		j = mutipath[m][3];
		qq[I][J] = mutistate[I][i + 1]*q[i][j];
		}

    for( I=0; I<nmutistate; I++) {
		s = 0.0;
		for( J=0; J<nmutistate; J++)
			if (J != I)
				s -= qq[I][J];
		qq[I][I] = s;
		}
    

	//derivatives to z-variables
    for( I=0; I<nmutistate; I++ )
		for( J=0; J<nmutistate; J++ )
			for( k=0; k<nz; k++ )
				dqq[I][J][k] = 0.0;
    
    for( m=0; m<nmutipath; m++ ) {
		I = mutipath[m][0];
		J = mutipath[m][1];
		i = mutipath[m][2];
		j = mutipath[m][3];
		for( k=0; k<nz; k++ )
			dqq[I][J][k] = mutistate[I][i + 1]*dq[i][j][k];
		} 
    
    for( k=0; k<nz; k++ ) {
		for( I=0; I<nmutistate; I++) {
			s = 0.0;
			for( J=0; J<nmutistate; J++) 
				if (J != I)
					s -= dqq[I][J][k];
			dqq[I][I][k] = s;
			}
		}
	}
