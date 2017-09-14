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

#ifndef COUNTED_PTR_H
#define COUNTED_PTR_H

/* class for counted reference semantics
 * - deletes the object to which it refers when the last CountedPtr
 * that refers to it is destroyed.
 *
 * The C++ Standard Library, Nicolai Josuttis, 1999, p222
 * "Many thanks to Greg Colvin and Beman Dawes for feedback"
 *
 * < and == added by cnicolai 3/01 for use in sets and maps
 */

template <class T>
class CountedPtr {
 private:
  T* pptr;        // pointer to the value
  long *count;   // shared number of owners
  
 public:
  // initialize pointer with existing pointer
  // - requires that the pointer p is a return value of new
  explicit CountedPtr( T* p = 0 )
    : pptr(p), count(new long(1)) {
  }
  
  // copy pointer (one more owner)
  CountedPtr( const CountedPtr<T> &p ) throw()
    : pptr(p.ptr()), count(p.count) {
    ++*count;
  }
  
  // destructor (delete value if this was the last owner)
  ~CountedPtr() throw() {
    dispose();
  }
  
  // assignment (unshare old and share new value)
  CountedPtr<T>& operator= (const CountedPtr<T>& p) throw() {
    if ( this != &p ) {
      dispose();
      pptr = p.ptr();
      count = p.count;
      ++*count;
    }
    return *this;
  }
  
  // access the value to which the pointer refers
  T& operator*() const throw() {
    return *pptr;
  }
  T* operator->() const throw() {
    return pptr;
  }
  T* getPtr() const {
    return pptr;
  }
  T* ptr() const {
    return pptr;
  }
  
  
  bool operator< (const CountedPtr<T>& p) const {
    return ( (long)ptr() < (long)p.ptr() );
  }
  
  bool operator== (const CountedPtr<T>& p) const {
    return ( (long)ptr() == (long)p.ptr() );
  }
  
 private:
  void dispose() {
    if ( --*count == 0 ) {
      delete count;
      if ( pptr ) // 4/19/01 chris--because 0 is def.constr.arg.
	delete pptr;
    }
  }
  
};

// {
//   CountedPtr< MdlClass > mdlPtr( new MdlClass );
//   storeModel( mdlPtr );
// }
// CountedPtr< MdlClass > storage;
// void storeModel( CountedPtr< MdlClass > mdlPtr ) {
//     storage = mdlPtr;
// }
// void unstoreModel() {
//    storage = CountedPtr< MdlClass >( NULL );
// }

template <class T>
class CountedArrayPtr {
 private:
  T* pptr;          // pointer to the value
  long *count;     //  shared number of owners
  int len;        //   size of array
  bool weDelete; //    so stuff can return a pointer that is counted or not
  
 public:
  // initialize pointer with existing pointer
  // - requires that the pointer p is a return value of new []
  explicit CountedArrayPtr( T* p = 0, int length = 0, bool shouldDelete = true )
    : pptr(p), count(new long(1)), len( length ), weDelete( shouldDelete )
  {
  }
  
  // copy pointer (one more owner)
  CountedArrayPtr( const CountedArrayPtr<T> &p ) throw()
    : pptr(p.ptr()), count(p.count), len( p.len ), weDelete( p.weDelete )
  {
    ++*count;
  }
  
  // destructor (delete value if this was the last owner)
  ~CountedArrayPtr() throw() {
    dispose();
  }
  
  // assignment (unshare old and share new value)
  CountedArrayPtr<T>& operator= (const CountedArrayPtr<T>& p) throw() {
    if ( this != &p ) {
      dispose();
      pptr = p.ptr();
      count = p.count;
      len = p.len;
      weDelete = p.weDelete;
      ++*count;
    }
    return *this;
  }
  
  // access the value to which the pointer refers
  T& operator*() const throw() {
    return *pptr;
  }
  T* operator->() const throw() {
    return pptr;
  }
  T* getPtr() const {
    return pptr;
  }
  T* ptr() const {
    return pptr;
  }
  T& operator[] (int offset) const {
    return pptr[offset];
  }
  
  int size() const {
    return len;
  }
  
  bool operator< (const CountedArrayPtr<T>& p) const {
    return ( (long)ptr() < (long)p.ptr() );
  }
  
  bool operator== (const CountedArrayPtr<T>& p) const {
    return ( (long)ptr() == (long)p.ptr() );
  }
  
 private:
  void dispose() {
    if ( --*count == 0 ) {
      delete count;
      if ( pptr && weDelete ) // 4/19/01 chris--because 0 is def.constr.arg.
	delete [] pptr;
    }
  }
  
};

#endif /* COUNTED_PTR_H */

/*
Copyright 1999 by Addison Wesley Longman, Inc. and Nicolai
M. Josuttis. All rights reserved. Permission to use, copy, modify and
distribute this software for personal and educational use is hereby
granted without fee, provided that the above copyright notice appears
in all copies and that both that copyright notice and this permission
notice appear in supporting documentation, and that the names of
Addison Wesley Longman or the author are not used in advertising or
publicity pertaining to distribution of the software without specific,
written prior permission. Addison Wesley Longman and the author make
no representations about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.

ADDISON WESLEY LONGMAN AND THE AUTHOR DISCLAIM ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL ADDISON WESLEY LONGMAN
OR THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.  
*/
