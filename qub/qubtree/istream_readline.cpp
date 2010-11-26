#include "istream_readline.h"
#include <cstring>

// no lame length restrictions here:
//   keeps doubling buflen until it fits; returns realloc'd buffer
char *istream_readline( istream& in, char *buf, int &buflen, int &linelen ){
	char newline;
	
	// bug in lib: empty line sets fail():
	if ( in.peek() == '\n' ) {
		in.get(newline);
		buf[0] = 0;
		return buf;
		}
	
	int at = 0;
	in.get( buf, buflen );
	
	while ( (! in.eof()) && in.fail() ) { // didn't find '\n'
		in.clear();
		at = buflen - 1;
		
		char *newbuf = new char[ 2 * buflen ];
		memcpy( newbuf, buf, at );
		delete [] buf;
		buf = newbuf;
		buflen = 2 * buflen;
		
		in.get( buf + at, buflen - at );
		}
	if ( in )
		in.get(newline);
	
	linelen = strlen( buf );
	if ( buf[ linelen - 1 ] == '\r' ) {
		buf[ linelen - 1 ] = 0;
		linelen--;
		}
	
	return buf;
	}

istream_linereader::istream_linereader( int initlen ){
	buflen = initlen;
	buf = new char[ initlen ];
	}

istream_linereader::~istream_linereader(){
	delete [] buf;
	}

char *istream_linereader::readLine( istream& in, int& len ){
	buf = istream_readline( in, buf, buflen, len );
	return buf;
	}
