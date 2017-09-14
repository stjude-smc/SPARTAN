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
