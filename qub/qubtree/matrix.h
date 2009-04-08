
/**
 * int size() const: if no const, then when calling b.size() from assgnmt
 * operator, because the argument b is const ref, the compiler will 
 * complain.
 *
 * One senario: CString::SetCString(const CString&) - the arg must be 
 * either const ref or CString, and cannot be CString&. Consider
 *   LPCSTR lp;
 *   obj.SetCString(lp);
 * Suppose there is a type conversion function. To compile, the compiler 
 * needs to create a temp variable and the code becomes:
 *   LPCSTR;
 *   CString temp(lp); 
 *   obj.SetCString(temp);
 * With CString&, the compiler must allow the function to change the 
 * contents of lp. But such change cannot propagate back through temp.
 * See Periodic journal of MSDN under key word: const reference.
 *
 * operator []: the return type should be just a ref, not const ref; 
 * otherwise, if a is a matrix, then a[i]=v would not be allowed -
 * in fact compiler generates a warning message:
 * 'passing `const fqvector<double>' as `this' argument of `const class 
 * fqvector<double> & fqvector<double>::operator =<double>(const class 
 * fqvector<double> &)' discards const.
 *
 * subscript[i]: sometimes, i may be some type other than int, Microsoft
 * compiler will complain. So we need to add overloaded subscript function
 * for each type.
 *
 * 4/18/00: xmatrix.cpp doesn't behave properly after adding the type
 * conversion function T* in fqvector. The bugs disappers after 
 * rearrange the const-ness, but the reason is unknown.
 */

#ifndef _MATRIX_
#define _MATRIX_
#include <math.h>
#include <assert.h>
#include <ostream>

#define true 1
#define false 0

// -------------------------------------------------------------------------
// fqvector class

template<class T> class fqvector {
	public:
		int n;             
		T* v; 
		int nbuf;      
		
		// empty Constructor
		fqvector() : n(0), v(0), nbuf(0) { 
			}            
		
		// copy Constructor
		fqvector(const fqvector<T>& b){
			n=nbuf=0, v=0; // v=0 is new 9/20/01 Chris 
			*this = b;                     // call assgn operator
			}
		
		fqvector(int m, T c) {
			n=nbuf=m;
			v=0;
			if (m > 0)
				v = new T[m];
			*this = c;                     // call assgn operator
			}

		fqvector(int m) { 
			n=nbuf=m;
			v=0; // v=0 is new 9/20/01 Chris
			if( m>0 ) 
				v = new T[m];
			for( int i=0; i<m; i++) 
				v[i]=(T)0;
			}
		
		// destructor
		~fqvector() {
			if (v != 0) 
				delete[] v;
			n=nbuf=0;
			}
		
		T& operator [] (int i) { 
			return v[i]; 
			} 
		T& operator [] (long i) { 
			return v[i]; 
			} 
		T& operator [] (unsigned long i) { 
			return v[i]; 
			} 
		
		// assgmt operators
		const fqvector<T>& operator =(const fqvector<T>& b){  // default assgnmt
			if (this == &b)
				return *this;
			
			if (n == b.size()) {
				for (int i=0; i<n; i++) 
					v[i] = b[i];
				return *this;
				}
			
			if (nbuf > 0) 
				delete[] v;
			n = nbuf = b.size();
			v = 0;
			if (n > 0) 
				v = new T[n];
			for (int i=0; i<n; i++)
				v[i] = b[i];          // recall that b[i] == *(b + i) == *(vector + scalar [add scalar to all elements]) == *vector = value of first element of vector after adding i.  Unfortunately the vector+T call itself calls vector=vector and the stack overflows.  Pretty sure this happens only with vector<int> (and maybe vector<long> if we had any) 
			return *this;
			}
		
		// set all members to a specified value 
		const fqvector<T>& operator =(const T c) {
			for (int i=0; i<n; i++)
				v[i] = c;
			return *this;
			}
		
		// cast (type conversion)
		operator T* () const { 
			return (T*)v; 
			}
		
		// resize while preserving
		void resize(int m){ // const fqvector<T>& resize(int m)                
			if ( m == n )
				return;
			
			fqvector<T> b(m);
			for (int i=0; i< (m<n ? m:n); i++)
				b[i] = v[i];
			*this = b;
			}
		
		// clear    
		void clear() {
			if (v!=0)
				delete[] v;   
			n=nbuf=0;
			}

		// append an element
		void push_back (const T c){ // const fqvector<T>& push_back (const T c)             
			if (n == nbuf){
				//----- If no room in preallocated space, grow buffer.
				T* b = new T[nbuf+nbuf+100];		//T* b = new T[nbuf+10000];		// Changed 11/2004 {JB}  generally these are growing from 0 to 1
				for ( int i=0; i<n; i++)
					b[i] = v[i];
				
				if (nbuf>0) 
					delete[] v;
				v = b;
				nbuf=nbuf+nbuf+100;		//nbuf += 10000;
				}
			v[n] = c;
			n++;
			}
		
		// empty, size
		int size() const { 
			return n; 
			} 
	};

template<class T> fqvector<T> operator -(const fqvector<T>& u, const fqvector<T>& v){
   return u+(-v);
	}

// ------------------------------------------------------------------
// matrix class
//
// Chris 9/20/01 fix:
// if you repeatedly 'T **m = (T**)theMatrix;'
// then operator T** would repeatedly new[] without delete[]
// * added delete[] where appropriate.
template<class T> class matrix {
	public:
		int nr;             // number of rows
		int nc;             // number of columns
		fqvector<T>* a;       // array of fqvectors
    
		T** ptr;            
		// used for type conversion from matrix to T**.
		// ptr[i] refers to &a[i][0], i.e. the elements 
		// of the fqvector a[i] instead of the fqvector itself.
		// this can be eliminated by changing a to T**.
		
		// empty Constructor 
		matrix() : nr(0), nc(0), a(0), ptr(0) { 
			}     
		
		matrix(const matrix<T>& b) {          // copy
			nr=nc=0; 
			ptr=0;
			a = 0; // new 9/20/01 Chris
			*this = b;   // call assgn operator
			}
		
		matrix(int m, int n, T c){
			nr=m, nc=n, a=0, ptr=0;
			if (m > 0){
				a = new fqvector<T>[m];
				for (int i=0; i<m; i++) 
					a[i].resize(n);
				}
			*this = c;   // call assgn operator
			}
		
		matrix(int m, int n) { 
			nr=m, nc=n, a=0, ptr=0;
			if (m > 0){
				a = new fqvector<T>[m];
				for (int i=0; i<m; i++) {
					a[i].resize(n);
					for (int j=0; j<n; j++)
						a[i][j] = (T)0;
					}
				}
			}
		
		// destructor
		~matrix() {
			if (a != 0) 
				delete[] a;
			if (ptr!=0) 
				delete[] ptr;
			}
		
		// assgmt operators
		const matrix<T>& operator = (const matrix<T>& b) { // default assgnmt
			if (this == &b)
				return *this;
			
			if (nr == b.m() && nc == b.n()){
				for (int i=0; i<nr; i++)
					for (int j=0; j<nc; j++)
						a[i][j] = b[i][j];
				return *this;
				}
			
			if (nr > 0) 
				delete[] a;   // invoke fqvector destructor???
			nr = b.m();
			nc = b.n();
			a = 0;
			if (nr > 0) 
				a = new fqvector<T>[nr];
			for (int i=0; i<nr; i++){
				a[i].resize(nc);
				for (int j=0; j<nc; j++)
					a[i][j] = b[i][j];
				}
			return *this;
			}
		
		const matrix<T>& operator = (const T c) {
			for (int i=0; i<nr; i++)
				for (int j=0; j<nc; j++)
					a[i][j] = c;
				return *this;
			}
		
		// other operators
		fqvector<T>& operator [] (int i) const {	// subscripting
			return a[i]; 
			} 
		
		// cast (type conversion)
		operator T** () { 
			if ( ptr ) 
				delete [] ptr; // new Chris 9/20/01
			ptr = new T*[nr];
			for (int i=0; i<nr; i++) 
				ptr[i] = a[i];
			return ptr; 
			}
		
		void operator *= (const matrix<T>& b){ // matrix<T>& operator *= (const matrix<T>& b)
			assert(nc == b.nr);
			int mm = nr, pp = nc, nn = b.nc;
			int i, j, k;
			matrix<T> c(mm, nn);
			for (i=0; i<mm; ++i)
				for (j=0; j<nn; ++j)
					for (k=0; k<pp; ++k)
						c[i][j] += a[i][k]*b[k][j];
					*this = c;
			}
		
		// functions
		void resize(int m, int n){ // matrix<T>& resize(int m, int n)                // resize, no copy ???
			if ( (m == nr) && (n == nc) )
				return;
			
			matrix<T> b(m,n);
			*this = b;
			}
		
		void clear(){                                  // clear
			if (a!=0)
				delete[] a;   // invoke fqvector destructor ???
			nr=nc=0;
			}
		
		void push_back(const fqvector<T>& v){ // fqvector<T>& push_back(fqvector<T>& v)             // append a new row
			assert (nc==0 || nc==v.size());
			nc = v.size();          
			matrix<T> b(nr+1,nc);
			int i,j;
			for (i=0; i<nr; i++)
				for (j=0; j<nc; j++)
					b[i][j] = a[i][j];
				for (j=0; j<nc; j++)
					b[nr][j] = v[j];
				*this = b;   // invoke fqvector destructor ???
			}
		
		int m() const { 
			return nr; 
			}                   // retrieve dimensions
		int n() const { 
			return nc; 
			}
	};

template<class T> matrix<T> operator *(const matrix<T>& a, const T b) {
	matrix<T> c = a;
	c *= b;
	return c;
	}

template<class T> matrix<T> operator *(const T a, const matrix<T>& b) {
	return b*a;
	}

template<class T> matrix<T> operator *(const matrix<T>& a, const matrix<T>& b) {
	matrix<T> c = a;
	c *= b;
	return c;
	}


// fqvector * matrix
template<class T> fqvector<T> operator *(const fqvector<T>& v, const matrix<T>& a) {
	int m = a.m();
	int n = a.n();
	assert(v.size()==m);
	fqvector<T> u(n);
	for (int j=0; j<n; j++){
		u[j] = 0.0;
		for (int i=0; i<m; i++)
			u[j] += v[i]*a[i][j];
		}
	return u;
	}

template<class T> fqvector<T> operator *(const matrix<T>& a, const fqvector<T>& v){
	int m = a.m();
	int n = a.n();
	assert(v.size()==n);
	fqvector<T> u(m);
	for (int i=0; i<m; i++) {
		u[i] = 0.0;
		for (int j=0; j<n; j++)
			u[i] += a[i][j]*v[j];
		}
	return u;
	}

// -----------------------------------------------------
// Tensor class
//
template<class T> class tensor {
	public:
		int nr;             // number of rows
		int nc;             // number of columns
		int nd;             // number of depth
		matrix<T>* a;       
		
		// constructors
		tensor() : nr(0), nc(0), nd(0), a(0) { 
			}          // default empty
		
		tensor(const tensor<T>& b){  // default copy
			nr=nc=nd=0; 
			*this = b;   // call assgn operator
			}
		
		tensor(int d, int m, int n, const T& c=(T)0){
			nd = d; 
			nr = m;
			nc = n;
			a = 0;
			if (nd>0) {
				a = new matrix<T>[nd];
				for (int i=0; i<nd; i++) 
					a[i].resize(nr,nc);
				}
			*this = c;   // call assgn operator
			}
		
		~tensor()  { 
			if (a != 0) 
				delete[] a;      // invoke matrix destructor???
			}
		
		// assgmt operators
		const tensor<T>& operator = (const tensor<T>& b){     // default assgnmt
			if (this == &b)
				return *this;
			
			if (nd > 0) 
				delete[] a;   // invoke matrix destructor???
			nr = b.m();
			nc = b.n();
			nd = b.d();
			if (nd > 0) 
				a = new matrix<T>[nd];
			for (int i=0; i<nd; i++)
				a[i] = b[i];
			return *this;
			}
		
		const tensor<T>& operator =(const T c) {
			for (int d=0; d<nd; d++)
				for (int i=0; i<nr; i++)
					for (int j=0; j<nc; j++)
						a[d][i][j] = c;
					return *this;
			}
		
		// other operators
		matrix<T>& operator [] (int i) const {       // subscripting
			return a[i]; 
			}
		
		// template<class T> inline T*** tensor<T>::operator T*** () const { return a; }
		
		// functions
		void resize(int d, int m, int n){ // tensor<T>& resize(int d, int m, int n)      // resize
			tensor<T> b(d,m,n);
			*this = b;
			//return *this; 
			}
		
		void clear(){                                // clear
			if (a != 0)
				delete[] a;   // invoke matrix destructor
			nd=nr=nc=0;
			a = 0;
			}
		
		int m() const { 
			return nr; 
			}                // retrieve dimensions
		int n() const { 
			return nc; 
			}
		int d() const { 
			return nd; 
			}
	};

#endif // _MATRIX_


/* todo - unused ? 
// unary plus    
   
   fqvector<T> operator +() const { return *this; }      
    
// unary minus    
   
   fqvector<T> operator - () const               
   {
      fqvector<T> b(n);
      for (int i=0; i<n; i++)
         b[i] = -v[i];
      return b;
   }
   
// equality

   bool operator == (const fqvector<T>& b)
   {
      if (n != b.size())
         return false;

      for (int i=0; i<n; i++) 
         if (! (v[i] == b[i]))
            return false;

      return true;
   }


// add and assign    
   
   void operator += (const T c)           
   {
      for (int i=0; i<n; i++)
         v[i] += c;
      //return *this;
   }

   void operator += (const fqvector<T>& b) //const fqvector<T>& operator += (const fqvector<T>& b)
   {
      for (int i=0; i<n; i++)
         v[i] += b[i];
      //return *this;
   }
    
// subtract and assign    
   
   void operator -= (const T c) // const fqvector<T>& operator -= (const T c)           
   {
      for (int i=0; i<n; i++)
         v[i] -= c;
      //return *this;
   }

   void operator -= (const fqvector<T>& b) // const fqvector<T>& operator -= (const fqvector<T>& b)
   {
      for (int i=0; i<n; i++)
         v[i] -= b[i];
      //return *this;
   }
    
// multiply and assign    
   
   void operator *= (const T c) // const fqvector<T>& operator *= (const T c)             
   {
      for (int i=0; i<n; i++)
         v[i] *= c;
      // return *this;
   }
*/
//   bool empty() const { return n==0; } 
// multiplication
/*
template<class T> fqvector<T> operator *(const fqvector<T>& v, const T c) 
{
   fqvector<T> u = v;
   u *= c;
   return u;
}

template<class T> fqvector<T> operator *(const T c, const fqvector<T>& v) 
{
   return v*c;
}
*/
// dot product
/* todo - unused ? 
template<class T> T operator *(const fqvector<T>& u, const fqvector<T>& v) 
{
   assert(u.size() == v.size());       
   T sum = 0.0;
   for (int i=0; i<v.size(); i++)
      sum += u[i]*v[i];
   return sum;
}

// division

template<class T> fqvector<T> operator / (const fqvector<T>& v, const T c) 
{
   fqvector<T> u(v);
   for (int i=0; i<u.size(); i++)
      u[i] /= c;
   return u;
}

// absolute values

template<class T> fqvector<T> abs(const fqvector<T>& v)   // complex ???
{
   fqvector<T> u(v);
   for (int i=0; i<u.size(); i++)
      u[i] = fabs(u[i]);
   return u;
}

// maximal element

template<class T> T fqv_max(const fqvector<T>& v) 
{
   T vm = v[0];
   for (int i=0; i<v.size(); i++)
      vm = (vm>v[i]) ? vm : v[i];
   return vm;
}

// summation

template<class T> T sum(const fqvector<T>& v)
{
   T sm = 0.0;
   for (int i=0; i<v.size(); i++)
      sm += v[i];
   return sm;
}

template<class T>
std::ostream& operator<< (std::ostream& out, const fqvector<T>& vec)
{
	out << '[';
	for ( int i=0; i<vec.size(); ++i )
		out << (i ? ", " : "") << vec[i];
	out << ']';
	return out;
}
*/
// addition
/* wtf does anyone use this?
template<class T> fqvector<T> operator +(const fqvector<T>& v, const T c)
{
   fqvector<T> u = v;
   u += c;
   return u;
}

template<class T> fqvector<T> operator +(const T c, const fqvector<T>& v)  
{
   return v+c;
}
*/
/* todo - unused ? 
template<class T> fqvector<T> operator +(const fqvector<T>& u, const fqvector<T>& v)  
{
   assert(u.size() == v.size());
   fqvector<T> z = u;
   z += v;
   return z;
}

// subtraction
/*
template<class T> fqvector<T> operator -(const fqvector<T>& v, const T c)  
{
   return v+(-c);
}

template<class T> fqvector<T> operator - (const T c, const fqvector<T>& v)  
{
   return -(v-c);
}
*/

/* TODO - unused ??? 
// append a fqvector
   void push_back (const fqvector<T>& u) // const fqvector<T>& push_back (const fqvector<T>& u)             
   {
      fqvector<T> b(n+u.size());
      int i;

      for (i=0; i<n; i++)
         b[i] = v[i];
      for (i=0; i<u.size(); i++)
         b[n+i] = u[i];
      *this = b;
      //return *this;
   }

// extract a sub-fqvector
   fqvector<T> slice (int k, int m)             
   {
      assert (k>=0 && k+m<=n);
      fqvector<T> b(m);
      for (int i=0; i<m; i++)
         b[i] = v[k+i];
      return b;
   }
*/

/* todo - unused 
    matrix<T> operator + () const { return *this; }      // unary plus
    
    matrix<T> operator - () const               // unary minus
    {
       matrix<T> b(nr,nc);
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             b[i][j] = -a[i][j];
       return b;
    }
    
    void operator += (const T c) // matrix<T>& operator += (const T c) // add and assign
    {
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             a[i][j] += c;
       //return *this;
    }

    void operator += (const matrix<T>& b) // matrix<T>& operator += (const matrix<T>& b)
    {
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             a[i][j] += b[i][j];
       // return *this;
    }
    
    void operator -= (const T c) // matrix<T>& operator -= (const T c)           // subtract and assign
    {
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             a[i][j] -= c;
       // return *this;
    }

    void operator -= (const matrix<T>& b) // matrix<T>& operator -= (const matrix<T>& b)
    {
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             a[i][j] -= b[i][j];
       // return *this;
    }
    
    void operator *= (const T c) // matrix<T>& operator *= (const T c)              // multiply and assign
    {
       for (int i=0; i<nr; i++)
          for (int j=0; j<nc; j++)
             a[i][j] *= c;
       // return *this;
    }
*/
/* todo - unused ? 
    matrix<T> t() const                            // transpose
    {
       matrix<T> b(nc,nr);
       for (int i=0; i<nc; i++)
          for (int j=0; j<nr; j++)
             b[i][j] = a[j][i];
       return b;
    }

    matrix<T> sub(int k, int m, int l, int n)      // submatrix
    {
       assert(k+m<=nr && l+n<=nc);
       matrix<T> b(m,n);
       for (int i=0; i<m; i++)
          for (int j=0; j<n; j++)
             b[i][j] = a[k+i][l+j];
       return b;
    }

    fqvector<T> diag()                        // diagonal fqvector
    {
       assert (nr==nc);
       fqvector<T> d(nr);
       for (int i=0; i<nr; i++)
          d[i] = a[i][i];
       return d;
    }
*/
/* todo - unused ?
    int size() const { assert(nr==nc); return nr; }
*/
/* todo - unused ? 
// operators as global functions
template<class T> matrix<T> operator +(const matrix<T>& a, const T b) 
{
   matrix<T> c = a;
   for (int i=0; i<c.m(); i++)
      for (int j=0; j<c.n(); j++)
         c[i][j] += b;
   return c;
}

template<class T> matrix<T> operator +(const T a, const matrix<T>& b) 
{
   return b+a;
}

template<class T> matrix<T> operator +(const matrix<T>& a, const matrix<T>& b) 
{
   assert(a.m()==b.m() && a.n()==b.n());
    matrix<T> c = a;
    for (int i=0; i<c.m(); i++)
      for (int j=0; j<c.n(); j++)
         c[i][j] += b[i][j];
   return c;
}

template<class T> matrix<T> operator -(const matrix<T>& a, const T b) 
{
   return a+(-b);
}

template<class T> matrix<T> operator -(const T a, const matrix<T>& b) {
   return -(b-a);
}

template<class T> matrix<T> operator -(const matrix<T>& a, const matrix<T>& b) {
   return a+(-b);
}
*/

/* todo - unused ?
// col fqvector * row fqvector
template<class T> matrix<T> mul(const fqvector<T>& u, const fqvector<T>& v) 
{
   int m = u.size();
   int n = v.size();
   matrix<T> a(m,n);
   for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
         a[i][j] = u[i]*v[j];
   return a;
}

template<class T> fqvector<T> diag(const matrix<T>& a)
{
   assert(a.m()==a.n());
   fqvector<T> v(a.n());
   for (int i=0; i<a.n(); i++)
      v[i] = a[i][i];
   return v;
}

template<class T> matrix<T> diag(const fqvector<T>& v)
{
   int n = v.size();
   matrix<T> a(n,n);
   for (int i=0; i<n; i++)
      a[i][i] = v[i];
   return a;
}

template<class T> matrix<T> diag(int n, const T c)
{
    matrix<T> a(n,n);
    for (int i=0; i<n; i++)
       a[i][i] = c;
    return a;
}

template<class T>
std::ostream& operator<< (std::ostream& out, const matrix<T>& mat)
{
	out << '[';
	for ( int i=0; i<mat.nr; ++i )
		out << (i ? "\n " : "") << mat[i];
	out << ']';
	return out;
}
*/

