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

#ifndef QUB_UNIQUEARRAYS_H
#define QUB_UNIQUEARRAYS_H

#include <map>
#include <vector>
#include "qubopt.h"

using std::map;
using std::vector;

template<class T>
class QUB_UniqueArray
{
public:
	T*  arr;
	int len;

	QUB_UniqueArray(T *zarr, int zlen)
		: arr(zarr), len(zlen)
	{}
	QUB_UniqueArray(const QUB_UniqueArray<T>& copy)
		: arr(copy.arr), len(copy.len)
	{}
};

template<class T>
bool operator< (const QUB_UniqueArray<T>& arr1, const QUB_UniqueArray<T>& arr2);

template<class T>
class QUB_UniqueArrays
{
private:
	typedef map< QUB_UniqueArray<T>, int > ArrayMap;
	typedef vector< vector<T>* > ArrayList;
	// Had to store vector pointers so the underlying arrays would not be surprise-reallocated.
	// That way, the Map keys can point to the vector*'s underlying arrays.

	ArrayMap           arrIndex;
	ArrayList          arrs;
	int                len;

public:
	QUB_UniqueArrays(int arrlen)
		: len(arrlen)
	{}
	virtual ~QUB_UniqueArrays()
	{
		for ( int i=0; i<int(arrs.size()); ++i )
			delete arrs[i];
	}
	
	int insert(T* arr)
	{
		int ix;

		typename ArrayMap::iterator ai = arrIndex.find( QUB_UniqueArray<T>(arr, len) );
		if ( ai == arrIndex.end() ) {
			ix = (int) arrs.size();
			arrs.push_back( new vector<int> );
			arrs[ix]->insert(arrs[ix]->begin(), arr, arr+len);
			arrIndex[ QUB_UniqueArray<T>(&((*arrs[ix])[0]), len) ] = ix;
		}
		else {
			ix = ai->second;
		}
		return ix;
	}

	int count()    { return (int) arrs.size(); }
	T*  get(int i) { return &((*arrs[i])[0]); }
};
			



class QUB_UniqueArrays_Integer_iface
{
public:
  virtual ~QUB_UniqueArrays_Integer_iface() {}
	virtual void   Free() = 0;

	virtual int    insert(int *arr) = 0;

	virtual int    count() = 0;
	virtual int*   get(int i) = 0;
};

class QUB_UniqueArrays_Integer : public QUB_UniqueArrays_Integer_iface
{
private:
	QUB_UniqueArrays<int>  arrs;

public:
	QUB_UniqueArrays_Integer(int len)
		: arrs(len)
	{}
	virtual ~QUB_UniqueArrays_Integer() {}

	void   Free()            { delete this; }
	int    insert(int *arr)  { return arrs.insert(arr); }
	int    count()           { return arrs.count(); }
	int*   get(int i)        { return arrs.get(i); }
};

#ifdef __cplusplus
extern "C" {
#endif

QUBOPT_API  QUB_UniqueArrays_Integer_iface *   QUB_UniqueArrays_Integer_Create( int len );

#ifdef __cplusplus
}
#endif




#endif
