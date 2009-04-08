#include "StringsInterned.h"


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
