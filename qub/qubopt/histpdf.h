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

#ifndef HISTPDF_H
#define HISTPDF_H

#include "qubopt.h"

	void binning(float ts, int iSegs, int* ndwell, float** tdwell, int* nbin, float* bin, float tdead, double& rOut);
	void calhist(float ts, int ic, int iSegs, int* ndwell, int **idwell, float** tdwell, int iBins, float* bin, float* hst, float tdead=0.0f);
	void calhist_smooth(int ic, int iSegs, int* ndwell, int **idwell, float** tdwell, int iBins, float* bin, float* hst, float tdead, float sampling, float r);

#endif
