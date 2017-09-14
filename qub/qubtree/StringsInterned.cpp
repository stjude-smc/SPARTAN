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

#include "StringsInterned.h"
#include <string.h>

StringsInterned::StringsInterned()
{
	emptyString = new char[1];
	emptyString[0] = 0;
}

StringsInterned::~StringsInterned()
{
	delete [] emptyString;
}

char *StringsInterned::get(const char *s)
{
  string ss(s);
  return get( ss );
}

char *StringsInterned::get(string &s)
{
	if ( s[0] == 0 )
		return emptyString;

	StringInterned *sRec = NULL;

	lock.lock();

	StringsInternedMap::iterator ii = strings.find(s);
	if ( ii == strings.end() ) {
		sRec = new StringInterned;
		sRec->count = 0;
		sRec->str = new char[ sizeof(void*) + s.size() + 1 ];
		* ((StringInterned **) (sRec->str)) = sRec;
		strcpy(sRec->str + sizeof(void*), &(s[0]));

		strings[s] = sRec;
	}
	else {
		sRec = ii->second;
	}

	sRec->count++;

	lock.unlock();

	return sRec->str + sizeof(void*);
}

char *StringsInterned::lookup(const char *s)
{
  string ss( s );
  return lookup( ss );
}

char *StringsInterned::lookup(string &s)
{
	if ( s[0] == 0 )
		return emptyString;

	StringInterned *sRec = NULL;

	lock.lock();

	StringsInternedMap::iterator ii = strings.find(s);
	if ( ii != strings.end() ) {
		sRec = ii->second;
	}

	lock.unlock();

	return sRec ? (sRec->str + sizeof(void*)) : NULL;
}

char *StringsInterned::dup(const char *s)
{
	if ( s == emptyString )
		return emptyString;

	StringInterned *sRec = * (StringInterned **) (s - sizeof(void*));
	lock.lock();
	sRec->count++;
	lock.unlock();

	return sRec ? (sRec->str + sizeof(void*)) : NULL;
}

void StringsInterned::release(const char *s)
{
	if ( s == emptyString )
		return;

	StringInterned *sRec = * (StringInterned **) (s - sizeof(void*));
	lock.lock();
	sRec->count--;
	if ( sRec->count == 0 ) {
		strings.erase( strings.find(s) );
		delete [] sRec->str;
		delete sRec;
	}
	lock.unlock();
}
