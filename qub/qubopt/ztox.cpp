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

// dx = a
// x = b + az	(matrix multiply)
void ztox(int nz, int nx, double **a, double *b, double *z, double *x, double **dx){
    int  i, j;
	double xsum;
                  
    for( i=0; i<nx; i++ ) {
		xsum = b[i];
		for( j=0; j<nz; j++ ) 
			xsum += (dx[i][j]=a[i][j]) * z[j];
		x[i]=xsum;
		}
	}
