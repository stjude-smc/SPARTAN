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

void pdftbl(int npd, double *pd){   
    int	k;
    double tmp;

    pd[npd] = 6.0/npd;
    for (tmp = 0.0, k = 0; k < npd; k++) {
		pd[k] = exp(-0.5*tmp*tmp)/sqrt(2.0*PI);
		tmp += pd[npd];
		}
	} 

extern "C" QUBOPT_API void mvectb(int ndata, double *data, double *amp, double *xms,
			int nar, double **ar, int nfir, double *fir,
			int nmetastate, int **metastate, double **b) {       
    int		I, i, k, m, n, t, t0 = nar + nfir, npd = 1000;
    double  tmp, sum, *pd;
	
    pd = dalloc1(npd + 1);
    pdftbl(npd, pd);
    for( I=0; I<nmetastate; I++) {
		m = metastate[I][1];
		for (t = t0; t < ndata; t++) {
			sum = 0.0;
			for (k = 0; k <= nar; k++) {
				tmp = 0.0;
				for (i = 0; i <= nfir; i++) {
					n = metastate[I][k + i + 1];
					tmp += fir[i]*amp[n];
					}
				sum += ar[m][k]*(data[t - k] - tmp);
				}
			sum /= xms[m];
			sum /= pd[npd];
			
			k = (int)fabs(sum);
			if (k >= npd)
				k = npd - 1;
			
			b[I][t] = pd[k] / xms[m];
			}
		}
	
	free(pd);
	}     


