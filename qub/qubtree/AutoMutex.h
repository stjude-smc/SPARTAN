#ifndef AUTO_MUTEX_H
#define AUTO_MUTEX_H

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

#include <vector>

#if defined(WIN32)

#include <windows.h>
#include <winbase.h>

// this would be useful to implement timeouts with critical sections, but TryEnterCriticalSection is only available in win2kpro or newer
//#include <Stopwatch.h>

// try critical section instead of mutex for speed boost
/*
#define AUTOMUTEX_WAIT_GRAIN 100

class Lock{
	protected:
		CRITICAL_SECTION crit;
	public:
		Lock() {
			InitializeCriticalSection(&crit);
		}
		~Lock() {
			DeleteCriticalSection(&crit);
		}
		bool lock( int timeoutMS = INFINITE ) {
			//if ( timeoutMS < 0 )
			//	timeoutMS = INFINITE;

			//if ( timeoutMS == INFINITE ) {
				EnterCriticalSection(&crit);
				return true;
			//}
			//else {
			//	Stopwatch sw;
			//	if ( 0 != TryEnterCriticalSection(&crit) ) {
			//		return true;
			//	}
			//	else {
			//		do {
			//			sw.waitMS( AUTOMUTEX_WAIT_GRAIN );
			//			if ( 0 != TryEnterCriticalSection(&crit) ) {
			//				return true;
			//			}
			//		} while ( sw.ms() < timeoutMS );
			//	}
			//}
			//return false;
		}
		void unlock() {
			LeaveCriticalSection(&crit);
		}
	};

*/
class Lock{
	protected:
		HANDLE mutex;
	public:
		Lock() {
			mutex = CreateMutex( NULL, false, NULL );
			}
		~Lock() {
			CloseHandle( mutex );
			}
		bool lock( int timeoutMS = INFINITE ) {
			if ( timeoutMS < 0 )
				timeoutMS = INFINITE;
			return ( WAIT_TIMEOUT != WaitForSingleObject( mutex, timeoutMS ) );
			}
		void unlock() {
			ReleaseMutex( mutex );
			}
	};

#else

#include <pthread.h>
#include <unistd.h>

class Lock{
 protected:
  pthread_mutex_t mutex;
 public:
  Lock() {
    pthread_mutexattr_t attr;
    pthread_mutexattr_init( &attr );
    pthread_mutexattr_settype( &attr, PTHREAD_MUTEX_RECURSIVE );
    pthread_mutex_init( &mutex, &attr );
    pthread_mutexattr_destroy( &attr );
  }
  ~Lock() {
    pthread_mutex_destroy( &mutex );
  }
  bool lock( int timeoutMS = -1 ) {
    // not sure if I can get real timeouts...
    if ( timeoutMS < 0 ) {
      pthread_mutex_lock( &mutex );
      return true;
    }
    else {
      while ( timeoutMS > 0 )
	if ( 0 == pthread_mutex_trylock( &mutex ) ) {
	  return true;
	} else {
	  sleep( 100 );     // of course it will sleep longer, so this is bad for serious timeout usage
	  timeoutMS -= 100; 
	}
      }
    return false;
  }
  void unlock() {
    pthread_mutex_unlock( &mutex );
  }
};

#endif


class AutoMutex{
	private:
		Lock mutex;
		volatile int refs;
		//std::vector<void*> whov;
		//void **who;
	public:
		AutoMutex() : refs(0) {}
		
		void incref() {
			lock(-1); // , 0);
			++refs;
			unlock(); //  0 );
			}
		
		void decref() {
			lock(-1); // , 0);
			--refs;
			unlock(); // 0 );
			
			if ( ! refs )
				delete this;
			
			}
		
		bool lock( int timeoutMS ) { // , void *caller ) {
			if ( mutex.lock( timeoutMS ) ) {
				//whov.push_back( caller );
				//who = &(whov[0]);
				return true;
				}
			return false;
			}
		
		void unlock() { //  void *caller ) {
			//std::vector<void*>::iterator whoi = whov.end();
			//do    {
			//	  if ( whoi == whov.begin() )
			//		  throw "unlock()er apparently didn't lock() first";
			//	  --whoi;
			//}
			//while ( *whoi != caller );
			//whov.erase( whoi );
			//who = whov.size() ? &(whov[0]) : 0;
			
			mutex.unlock();
			}
	};


class AutoMutexHolder
// can safely be switched from one to another or none
{
public:
	AutoMutex *mutex;
	volatile int count;
	//void *caller;

	AutoMutexHolder( AutoMutex *initMutex ) // , void *owner )
		: mutex( initMutex ), count( 0 ) // , caller( owner )
	{
		if ( mutex )
			mutex->incref();
	}
	~AutoMutexHolder()
	{
		if ( mutex )
			mutex->decref();
	}

	void switchMutex( AutoMutex *newMutex ) // must have lock, newMutex can be NULL
	{
		if ( newMutex ) {
			newMutex->lock( -1 ); // , caller );
			newMutex->incref();
		}

		AutoMutex *oldMutex = mutex;
		mutex = newMutex;

		if ( oldMutex ) {
			oldMutex->unlock(); //  caller );
			oldMutex->decref();
		}
	}

	bool lock( int timeoutMS = -1 ) {
		AutoMutex *myMutex = mutex;
		if ( ! myMutex ) {
			++count;
			return true;
		}

		if ( ! myMutex->lock( timeoutMS ) ) // , caller ) )
			return false;

		if ( mutex != myMutex ) { // switched out from under us
			myMutex->unlock(); //  caller );
			return lock( timeoutMS );
		}
		++count;
		return true;
	}

	void unlock() {
		--count;
		if ( mutex )
			mutex->unlock(); //  caller );
	}
};


#endif

