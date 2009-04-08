#ifndef QUB_QTR_IMPL_H
#define QUB_QTR_IMPL_H

//
// Define API decoration for direct importing of DLL references.
//

#if defined(_WIN32)
#define QTR_DLLEXPORT __declspec(dllexport)
#define QTR_DLLIMPORT __declspec(dllimport)
#else
#define QTR_DLLEXPORT
#define QTR_DLLIMPORT
#endif

#if !defined(_QTR_)
#define QTR_API QTR_DLLIMPORT
#else
#define QTR_API QTR_DLLEXPORT
#endif

#include <stdio.h>
#include <string>

#if defined(_WIN32)
//#include <windows.h>
#else

#endif

//----- flags on disk
#define QTR_FLAG_PRELOAD             0x000002
#define QTR_FLAG_DATA_IN_NODE        0x000004
#define QTR_FLAG_MATRIX              0x000008
// #define NODE_FLAG_UNSIGNED            0x000001
// #define QTR_FLAG_NONWRITABLE         0x000010

//----- flags in memory only 
#define QTR_FLAG_CHANGED             0x001000
#define QTR_FLAG_CHILDCHANGED        0x002000

//----- Data types stored 
// For QTR_NumCol implementation, all types UCHAR forward are considered numeric.
enum QTR_DataType { 
			QTR_TYPE_EMPTY, QTR_TYPE_UNKNOWN, QTR_TYPE_POINTER, 
			QTR_TYPE_STRING,
		    QTR_TYPE_UCHAR, QTR_TYPE_CHAR, 
			QTR_TYPE_USHORT, QTR_TYPE_SHORT, 
			QTR_TYPE_UINT, QTR_TYPE_INT,
		    QTR_TYPE_ULONG, QTR_TYPE_LONG, 
			QTR_TYPE_FLOAT, QTR_TYPE_DOUBLE, QTR_TYPE_LDOUBLE }; 

#define QTR_LOAD_NONE 0xffffffff

extern int QTR_DataSize[];

//----- Placed at start of file as identifier 
#define QTR_MAGIC               "QUB_(;-)_QFS"
#define QTR_MAGIC_LEN           12

#define QTR_NRSV 7

#ifdef CAUTOMUTEX_NEEDED

//=============================================================================
// CAutoMutex and CAutoMutexHolder provide locking mechanism for trees.  Each 
// node has a CAutoMutexHolder pointer to a CAutoMutex object or NULL.  Usually
// the leaf nodes point to the root mutex.  The CAutoMutex is reference counted 
// for automatic deletion.
// The Mutex can safely be switched from one Holder to another or set to null.
class CAutoMutex {
	private:
		int refs;
		HANDLE h_mutex;
		int m_iGlobalIndex;

	public:
		CAutoMutex() ;
		~CAutoMutex();
		void incref();
		void decref();
		bool W_decref();
		bool lock( int timeoutMS=-1 );
		void unlock();
	};

class CAutoMutexHolder {
	public :
		CAutoMutex * m_pMutex;
		CAutoMutexHolder( bool lCreate ) ;
		~CAutoMutexHolder();
		void switchMutex( CAutoMutex * newMutex );
		bool lock( int timeoutMS = -1 );
		void unlock();
	};

#endif

//=============================================================================
//----- Header for each tree node on disk.
typedef struct {
	unsigned char flags[3]; //*
	unsigned char dataType; //*
	unsigned int dataSize; //* moved
	unsigned int dataCount; //* moved
	unsigned int dataPos; //* moved
	unsigned int childPos;
	unsigned int siblingPos;
	unsigned char reserved[ QTR_NRSV ]; //* consolidated, eroded
	unsigned char nameLen;
	// name goes here
	} QTR_OnDisk;


typedef struct {
	FILE * handle;
	char * path;
	unsigned int frontier;
	} QTR_File;


typedef struct QTR_Impl {
	QTR_OnDisk ondisk;
	
	int refs;
	
	char *name;
	
	struct QTR_Impl *parent, *child, *sibling;
	
	unsigned int  loadStart;
	unsigned int  loadEnd;
	void *data; // the loaded portion
	
	QTR_File *    file;
	unsigned int  fileStart;
	
	void *        m_mutex;	// ptr to object of type CAutoMutex
	int           readers; //* # of readers, or NODE_WRITING
	
	struct QTR_Impl *root; // 0 for root nodes
	
	} QTR_Impl;


extern "C" {
QTR_API char *         QTR_LookupName( const char *name ); // get the canonical interned char * for quick comparison instead of strcmp

	// Locked increment of refcount.  Returns new refcount.
QTR_API int            QTR_INCREF( QTR_Impl *impl );
	// locked decrement of refcount, and QTR_Free if 0.  Returns new refcount.
QTR_API int            QTR_DECREF( QTR_Impl *impl );

	// Lock using impl->m_mutex
QTR_API int            QTR_Lock( QTR_Impl *impl, int timeoutMS );
QTR_API void           QTR_Unlock( QTR_Impl *impl );

//TODOC
	//----- root nodes
QTR_API QTR_Impl *     QTR_Create( const char *name );
QTR_API QTR_Impl *     QTR_Clone( QTR_Impl *other, int deep );
QTR_API QTR_Impl *     QTR_Open( const char *path, int readOnly );
QTR_API void           QTR_Close( QTR_Impl *impl );
QTR_API int            QTR_Save( QTR_Impl *impl );
QTR_API int            QTR_SaveAs( QTR_Impl *impl, const char *path ); // if called on non-root, removes it to be a root first
QTR_API int            QTR_SaveAsTemp( QTR_Impl *impl );
// saveAs doesn't actually save; it just opens the file for writing

QTR_API QTR_Impl *     QTR_CreateChild( QTR_Impl *parent, QTR_Impl *presib, const char *name );
QTR_API QTR_Impl *     QTR_InsertClone( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *other, int deep );

QTR_API void           QTR_InsertChild( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child );
QTR_API void           QTR_RemoveChild( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child );

QTR_API int            QTR_Flag( QTR_Impl *impl, int mask );
QTR_API void           QTR_SetFlag( QTR_Impl *impl, int mask, int val );

	//----- set CHANGED flag and parent->CHILDCHANGED flag w. propagate up
QTR_API void           QTR_Changed( QTR_Impl *impl );
	//----- set CHILDCHANGED flag w. propagate up
QTR_API void           QTR_ChildChanged( QTR_Impl *impl );

QTR_API char *         QTR_Name( QTR_Impl *impl );
QTR_API void           QTR_SetName( QTR_Impl *impl, const char *name );

QTR_API int            QTR_IsInFile( QTR_Impl *impl );
QTR_API char *         QTR_FilePath( QTR_Impl *impl );

QTR_API int            QTR_SetupData( QTR_Impl *impl, QTR_DataType dtype, unsigned int size, unsigned int count );
QTR_API int            QTR_SetupNumData( QTR_Impl *impl, QTR_DataType type, unsigned int nr, unsigned int nc );
QTR_API int            QTR_SetupStringData( QTR_Impl *impl, int count );
QTR_API int            QTR_ResizeData( QTR_Impl *impl, unsigned int newnr );
QTR_API int            QTR_ClearData( QTR_Impl *impl );

QTR_API unsigned int   QTR_NumCol( QTR_Impl *impl );

QTR_API int            QTR_LoadRows( QTR_Impl *impl, unsigned int firstRow, unsigned int lastRow, int doRead );
QTR_API int            QTR_UnloadRows( QTR_Impl *impl, int doWrite );

QTR_API unsigned int   QTR_GetRows( QTR_Impl *impl, void *buf, unsigned int first, unsigned int last );
QTR_API unsigned int   QTR_SetRows( QTR_Impl *impl, const void *buf, unsigned int first, unsigned int last );


// conversion to/from a contiguous buffer in ondisk format
QTR_API unsigned int   QTR_ToBufferSize( QTR_Impl *impl );
QTR_API unsigned int   QTR_ToBuffer( QTR_Impl *impl, char *buffer );
QTR_API QTR_Impl *     QTR_FromBuffer( const char *buffer );


}

#endif
