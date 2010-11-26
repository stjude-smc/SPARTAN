#ifndef STRINGS_INTERNED_H
#define STRINGS_INTERNED_H

#include <string>
#include <cstring>
#include <AutoMutex.h>

#if defined(_WIN32)
#include "notemplatewarning.h"
#include <hash_map>
using stdext::hash_map;
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;
namespace __gnu_cxx
{
        template<> struct hash< std::string >
        {
                size_t operator()( const std::string& x ) const
                {
                        return hash< const char* >()( x.c_str() );
                }
        };
}
#endif

using std::string;


/*

  Goals:
			keep one unique copy of each string in use
			be able to compare those strings for equality by pointer alone
			keep a reference count for each string so they can be freed
			strings transparently compatible with classic "char *" style
			quick lookup (via hash-table)

*/


typedef struct
{
	int count;
	char *str; // the first sizeof(void*) bytes are a pointer to this struct -- the string follows.
	           // the user's "char *" is 4 bytes past str.  Thus we can get this struct (e.g. to inc/dec count)
	           // directly from the char*, without hashmap lookup.
} StringInterned;

typedef hash_map<string, StringInterned *> StringsInternedMap;

class StringsInterned
{
private:
	Lock lock;
	StringsInternedMap strings;
	char *emptyString; // too common for the hash-map

public:
	StringsInterned();
	virtual ~StringsInterned(); // for now it is up to you to release all your strings first

	char *get(const char *s); // increments count
	char *lookup(const char *s); // doesn't increment count; may return NULL if none

	// the string& variants are faster unless you already have a char*
	char *get(string &s); // increments count
	char *lookup(string &s); // doesn't increment count; may return NULL if none

	// these operate on the interned char* returned by get or lookup
	char *dup(const char *s); // increments count
	void release(const char *s); // decrements count, possibly frees the string's memory
};


#endif
