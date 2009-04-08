#ifndef ISTREAM_READLINE_H
#define ISTREAM_READLINE_H

#include <iostream>
using namespace std;

char *istream_readline( istream& in, char *buf, int &buflen, int &linelen );

class istream_linereader{
	public:
		istream_linereader( int initlen = 512 );
		virtual ~istream_linereader();
		char *readLine( istream& in, int& len );

	private:
		char *buf;
		int buflen;
	};

#endif
