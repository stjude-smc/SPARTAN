#ifndef QUB_HEAP_H
#define QUB_HEAP_H

#include <vector>
#include <algorithm>
using namespace std;

// first() returns biggest [using operator<(T&, T&)]

template<class T>
class QUB_Heap{
private:
	vector<T> elts;

public:
	typedef typename vector<T>::iterator iterator;

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
