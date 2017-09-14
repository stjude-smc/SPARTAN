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

#ifndef QUB_HEAP_H
#define QUB_HEAP_H

#include <vector>
#include <algorithm>
using namespace std;

// first() returns biggest [using operator<(T&, T&)]

template<class T>
class QUB_Heap{
private:
	std::vector<T> elts;

public:
	typedef typename std::vector<T>::iterator iterator;

	void push(T& x) {
		elts.push_back(x);
		push_heap(elts.begin(), elts.end());
	}

	void resort() {
		make_heap(elts.begin(), elts.end());
	}

	T& first() {
		return elts[0];
	}

	void pop() {
		pop_heap(elts.begin(), elts.end());
		elts.pop_back();
	}

	int size() {
		return elts.size();
	}

	void erase(T& x) {
		iterator ii = find( elts.begin(), elts.end(), x );
		erase(ii);
	}

	void erase(iterator& ii) {
		int i = &(*ii) - &(elts[0]);
		if ( i >= int(elts.size()) )
			return;
		while ( i ) {
			int p = (i+1) / 2 - 1;
			elts[i] = elts[p];
			i = p;
		}
		// elts[0] = x;
		pop();
	}
	/*
		erase( find( elts.begin(), elts.end(), x ) );
		for ( vector<T>::iterator i = elts.begin(); i != elts.end(); ++i )
			if ( *i == x ) {
				elts.erase(i);
				break;
			}
		make_heap(elts.begin(), elts.end());
	*/

	iterator begin() { return elts.begin(); }
	iterator end()   { return elts.end(); }
};


#endif
