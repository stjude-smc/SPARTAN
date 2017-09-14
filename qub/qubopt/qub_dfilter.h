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

#ifndef DFILTER_H
#define DFILTER_H

#include "qubopt.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void * Filter32_p;
typedef void * Filter64_p;

QUBOPT_API Filter32_p DFII32_New();
QUBOPT_API void DFII32_Free(Filter32_p F);
QUBOPT_API int DFII32_Init(const char * FileName, Filter32_p F, float IC); //  IC = 0.0
QUBOPT_API float DFII32(float X, Filter32_p F);
QUBOPT_API float DFIIT32(float X, Filter32_p F);

QUBOPT_API Filter32_p DFII64_New();
QUBOPT_API void DFII64_Free(Filter32_p F);
QUBOPT_API int DFII64_Init(const char * FileName, Filter64_p F, double IC); // IC = 0.0
QUBOPT_API double DFII64(double X, Filter64_p F);
QUBOPT_API double DFIIT64(double X, Filter64_p F);

#ifdef __cplusplus
}
#endif

#endif
