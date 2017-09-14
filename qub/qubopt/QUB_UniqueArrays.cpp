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

#include "QUB_UniqueArrays.h"
#include <algorithm>
using namespace std;



template<class T>
bool operator< (const QUB_UniqueArray<T>& arr1, const QUB_UniqueArray<T>& arr2)
{
	int N = min(arr1.len, arr2.len);
	T *a1 = arr1.arr;
	T *a2 = arr2.arr;
	for ( int i=0; i<N; ++i, ++a1, ++a2 )
		if ( *a1 < *a2 )
			return true;
		else if ( *a1 > *a2 )
			return false;
	return false;
}


extern "C" QUBOPT_API  QUB_UniqueArrays_Integer_iface *   QUB_UniqueArrays_Integer_Create( int len )
{
	return new QUB_UniqueArrays_Integer(len);
}


