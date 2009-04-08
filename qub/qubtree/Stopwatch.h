#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <sys/types.h>
#include <sys/timeb.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <sys/time.h>
  #include <unistd.h>
#endif

class Stopwatch
{
 private:
  struct timeb init;
  struct timeb now;
 public:
  Stopwatch( bool started=true ) {
    ftime( &init );
  }
  int ms() {
    int elapsed = 0;
    ftime( &now );
    while ( now.time-- > init.time )
      elapsed += 1000;
    elapsed += now.millitm - init.millitm;
    return elapsed;
  }
  void waitMS( int msec ) {
#ifdef _WIN32
	  if ( msec > 0 )
		Sleep( msec );
#else
    struct timeval tv;
    tv.tv_sec = msec / 1000;
    tv.tv_usec = 1000 * (msec - tv.tv_sec * 1000);
    select( 0, 0, 0, 0, &tv );
#endif
  }
  void waitUntilMS( int msec ) {
    int rem = ms() - msec;
    if ( rem > 0 )
      waitMS( rem );
  }
};

#endif

  
