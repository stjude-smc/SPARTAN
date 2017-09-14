/* Copyright 1998-2011 Research Foundation State University of New York */

/* This file is part of QuB.                                            */

/* QuB is free software; you can redistribute it and/or modify          */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* QuB is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License,   */
/* named LICENSE.txt, in the QuB program directory.  If not, see        */
/* <http://www.gnu.org/licenses/>.                                      */

#include "qublib.h"

int rdHeader(FILE *fp, int *nlead, double *dt, int *count) {
	unsigned short  sampling, scaling;
    int	leader;
	
	if (fread(&leader, sizeof(leader), 1, fp) != 1
		|| fread(&sampling, sizeof(sampling), 1, fp) != 1
		|| fread(&scaling, sizeof(scaling), 1, fp) != 1 )  
		return FAILURE;

	*count = scaling;
    *nlead = leader/sampling;                 
    *dt = sampling*1.0e-6; 
    return SUCCESS;
	}

int rdSegHeader(FILE *fp, int *tstart, int *ndata) {
	int start, length;

    if( fread(&start, sizeof(start), 1, fp) != 1
		|| fread(&length, sizeof(length), 1, fp) != 1) 
		return FAILURE;

	*tstart = start;
    *ndata = length;
    return SUCCESS;
	}

int rdSegData(FILE *fp, int count, int ndata, double *data) {
    int	i;
    short	* buffer = new short[ndata];
    
	if (buffer==NULL 
		|| fread(buffer, ndata, sizeof(short), fp) != sizeof(short)) {
		delete [] buffer;
        return FAILURE;
		}

    try {
		for( i=0; i<ndata; i++)
			data[i] = (double)buffer[i] / (double)count*100.0;
	    } 
    catch (...) {}
	
	delete [] buffer;
    return SUCCESS;
	}
    

