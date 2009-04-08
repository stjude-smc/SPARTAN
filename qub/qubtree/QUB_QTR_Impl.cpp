#ifdef _USRDLL
  #define _QTR_
#endif

#include "notemplatewarning.h"

#include "QUB_QTR_Impl.h"
#include <set>
#include <AutoMutex.h>
#include "StringsInterned.h"

#include <stdexcept>
#include <iostream>

#define IMIN(a,b) ((a) < (b) ? (a) : (b))  
#define IMAX(a,b) ((a) > (b) ? (a) : (b))  

const long QTR_COPYING_MAXBYTES=1000000;	//TODO - test with 2^16 bytes - faster ? Smaller size should allow it to use preread buffer.

using namespace std;

int QTR_DataSize[] = { 0, 4, 4, 1, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 10 };

//=============================================================================// node name strings are interned
// TODO - why as a pointer, not an object ? 
StringsInterned * namesInterned = NULL;


char * QTR_LookupName( const char *name ){
	if ( namesInterned == NULL )
		namesInterned = new StringsInterned;
	return namesInterned->lookup( name );
	}

#if defined(_WIN32)

#ifdef CAUTOMUTEX_NEEDED

CAutoMutex :: CAutoMutex() {
	refs=0;
	h_mutex = CreateMutex( NULL, false, NULL );
	}

CAutoMutex :: ~CAutoMutex() {
	CloseHandle( h_mutex );
	}

void CAutoMutex :: incref() {
	lock();
	++refs;

	unlock();
	}

void CAutoMutex :: decref() {
	lock();
	--refs;
	unlock();

	if ( ! refs )
		delete this;
	}

bool CAutoMutex :: W_decref() {
	lock();
	--refs;
	unlock();

	if ( ! refs ) {
		//delete this;
		return true;
		}
	return false;
	}

bool CAutoMutex :: lock( int timeoutMS/*=-1*/ ) { 
	return ( WAIT_TIMEOUT != WaitForSingleObject( h_mutex, timeoutMS < 0 ? INFINITE : timeoutMS) );
	}

void CAutoMutex :: unlock() {
	ReleaseMutex( h_mutex );
	}


//=============================================================================
CAutoMutexHolder :: CAutoMutexHolder( bool lCreate )  {
	m_pMutex=NULL;
	if( lCreate ) {
		m_pMutex = new CAutoMutex;
		m_pMutex->incref();
		}
	}

CAutoMutexHolder :: ~CAutoMutexHolder() {
	if ( m_pMutex )
		m_pMutex->decref();
		/*
		if ( m_pMutex->W_decref() ) { 
			delete m_pMutex;
			m_pMutex=NULL;
			}
		*/
	}

void CAutoMutexHolder :: switchMutex( CAutoMutex * newMutex ) { // must have lock, newMutex can be NULL
	if ( newMutex ) {
		newMutex->lock( -1 ); 
		newMutex->incref();
		}

	CAutoMutex *oldMutex = m_pMutex;
	m_pMutex = newMutex;
	if ( oldMutex ) {
		oldMutex->unlock();
		
		oldMutex->decref();
/*		if ( oldMutex->W_decref() ) { 
			delete oldMutex;
			oldMutex=NULL;
			}
*/
		}
	}

bool CAutoMutexHolder :: lock( int timeoutMS /*= -1*/ ) {
	CAutoMutex * myMutex = m_pMutex;
	if ( ! myMutex )
		return true;
	
	if ( ! myMutex->lock( timeoutMS ) ) // , caller ) )
		return false;
	
	if ( m_pMutex != myMutex ) { // switched out from under us
		myMutex->unlock(); //  caller );
		return lock( timeoutMS );
		}
	return true;
	}

void CAutoMutexHolder :: unlock() {
	if ( m_pMutex )
		m_pMutex->unlock(); //  caller );
	}

#endif

#endif

//=============================================================================
// Local functions 

	//----- struct QTR_Impl default construction
	QTR_Impl * New_QTR_Impl(int iRefs);

	//----- Free a node and all children.
	// Calls RemoveChild for all children then ClearData (write) for parent
	void QTR_Free( QTR_Impl *impl );
	
	//----- Read node and children.  Data too, if PRELOAD flag is set.
	// Result is -1:error  0:success  1:SomeChildReadFailed.
	int QTR_Read( QTR_Impl *impl ); // using file/root and fileStart; requires something from create or createchild
	
	//----- Add child as the next sibling of presib under parent parent.
	void QTR_Link( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child );
	//----- Remove child from parent and previous sibling presib 
	int QTR_Unlink( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child );
	int QTR_Unlink( QTR_Impl * child );
	
	//----- Copy data <bytes> from file <from> to file <to>
	void QTR_CopyF2FinChunks( FILE * hSource, unsigned int iSourcePos, FILE * hDest, unsigned int iDestPos, unsigned int bytes );
	
	//todoc
	int QTR_MoveData( QTR_Impl *impl, QTR_File *from, QTR_File *to );

	//todoc
	int QTR_ImitateData( QTR_Impl *impl, QTR_Impl *other );
	
	//----- Set fileStart = 0 in impl and all its descendants.
	void QTR_ClearFileStarts( QTR_Impl *impl );
	
	//----- Set root in impl and all its descendants.
	void QTR_SetRoots( QTR_Impl *impl, QTR_Impl *root );
	
	//----- WriteData 
	int QTR_DumpLoadedToDisk( QTR_Impl *impl );
	
	//----- ReadData : Load rows of data in range impl->[loadStart,loadEnd] to impl->data
	int QTR_LoadFromDisk( QTR_Impl *impl ); // according to load*
	
	//----- Returns position in disk file, appending to file if necc.
	unsigned int QTR_FindFilePos( QTR_Impl *impl );
	
	//----- close, cleanup and delete QTR_FILE 
	void QTR_FreeFile( QTR_File * file );
	
	//----- Read <last> rows of data from file and send in blocks of 1MB to QTR_SetRows
	void QTR_SetRowsFromFile( QTR_Impl *impl, unsigned int loc, unsigned int last );


	
//-----------------------------------------------------------------------------
QTR_Impl * New_QTR_Impl(int iRefs) { 
	QTR_Impl * impl = new QTR_Impl;
	memset( impl, 0, sizeof(QTR_Impl) );
	*((int*) & impl->ondisk.flags) = (QTR_FLAG_PRELOAD | QTR_FLAG_CHANGED);	// default for new nodes
	impl->refs = iRefs;			// if 0, don't return a new ref, just leave one with parent
	impl->loadStart = impl->loadEnd = QTR_LOAD_NONE;
	return impl;
	}

//------------------------------------------------------------------------------
void QTR_Free( QTR_Impl *impl ) {
	if ( impl==NULL || impl->parent!=NULL )
		return;
	
	while ( impl->child ) {
		QTR_RemoveChild( impl, 0, impl->child );
	}

	QTR_ClearData( impl );
	
	if ( impl->file )
		QTR_Close( impl );

	delete (AutoMutexHolder *) impl->m_mutex;

	if ( impl->name )
		namesInterned->release( impl->name );

	delete impl;
}

//------------------------------------------------------------------------------
int QTR_Read( QTR_Impl *impl ) {	// using file/root and fileStart
	//----- Read 32 byte node header from disk.
	fseek( impl->root->file->handle, impl->fileStart, SEEK_SET );
	if ( 1 != fread( &(impl->ondisk), sizeof(QTR_OnDisk), 1, impl->root->file->handle ) ) {
		memset(&(impl->ondisk), 0, sizeof(QTR_OnDisk));
		return -1;
	}

	//----- Clear existing name 
	if ( impl->name ) {
		namesInterned->release( impl->name );
		impl->name = 0;
	}

	//----- Read node name to buffer then copy to node.
	if ( impl->ondisk.nameLen ) {
		char *tmpName = new char[ impl->ondisk.nameLen + 1 ];
		tmpName[ impl->ondisk.nameLen ] = 0;
		if ( 1 != fread(tmpName, impl->ondisk.nameLen, 1, impl->root->file->handle) ) {
			delete [] tmpName; 
			memset(&(impl->ondisk), 0, sizeof(QTR_OnDisk));
			return -1;
			} 
			impl->name = namesInterned->get( tmpName );
			delete [] tmpName;
		}
	else	// if name length is 0, set name to interned("")
		impl->name = namesInterned->get(""); 

	//----- if data exists and flag is preload - load data 
	if ( impl->ondisk.dataCount ) {
		if ( QTR_Flag( impl, QTR_FLAG_PRELOAD ) ) 
		  if ( QTR_LoadRows(impl, 0, impl->ondisk.dataCount - 1, 1) < 0 )
		    throw std::runtime_error("failed to load data");
	}
	else
		QTR_SetFlag( impl, QTR_FLAG_PRELOAD, 1 );

	//----- Clear modified flags 
	QTR_SetFlag( impl, QTR_FLAG_CHANGED, 0 );
	QTR_SetFlag( impl, QTR_FLAG_CHILDCHANGED, 0 );

	//----- Read all children 
	if ( impl->ondisk.childPos ) {
		QTR_Impl *lastChild = 0, *child;
		unsigned int nextChildPos = impl->ondisk.childPos;
		while ( nextChildPos ) {
			child = QTR_CreateChild( impl, lastChild, 0 );
			child->fileStart = nextChildPos;
			try {
			  QTR_Read(child);
			  nextChildPos = child->ondisk.siblingPos;
			  lastChild = child;
			} catch (std::exception& e) {
			  //QTR_RemoveChild(impl, lastChild, child);
			  throw e;
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
void QTR_Link( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child ) {
	if ( presib ) {
		child->sibling = presib->sibling;
		presib->sibling = child;
	}
	else {
		child->sibling = parent->child;
		parent->child = child;
	}
	child->parent = parent;
	QTR_INCREF( child );

	//TODO - set parent Changed and Parent->Parent ChildChanged, 
	// then sets presib->changed and presib->parent ChildChanged
	// Which both propagate up same line.  Also, should child->changed 
	// be set instead in first call ( since its siblings change ) 
	QTR_Changed( parent );
	QTR_ChildChanged( parent );
	if ( presib )
		QTR_Changed( presib );
}

//-----------------------------------------------------------------------------
int QTR_Unlink( QTR_Impl * child ){
	QTR_Impl * parent = child->parent;

	QTR_Impl * presib = NULL;
	if ( parent ) { 
		// scan linked list of <child> ->parent->child for prior sibling to <child>
		presib = parent->child;
		while ( presib && (presib->sibling != child) )
			presib = presib->sibling;
		}

	return QTR_Unlink( parent, presib, child );
	}

int QTR_Unlink( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child ){
	if ( presib ) {
		presib->sibling = child->sibling;
		}
	else {
		parent->child = child->sibling;
		}
	child->sibling = 0;
	child->parent = 0;

	QTR_Changed( parent );
	if ( presib )
		QTR_Changed( presib );
		//TODO : if QTR_Changed of presib then do not qtr_Changed of parent...?

	return QTR_DECREF( child );
	}

//-----------------------------------------------------------------------------
void QTR_CopyF2FinChunks( FILE * hSource, unsigned int iSourcePos, 
						 FILE * hDest, unsigned int iDestPos, unsigned int bytes ){

	unsigned int chunksize = IMIN( QTR_COPYING_MAXBYTES, bytes );
	char *buf = new char[ chunksize ];
	
	if ( 0 != fseek(hSource, iSourcePos, SEEK_SET)
		|| 0 != fseek(hDest, iDestPos, SEEK_SET) )
		return;

	while ( bytes ) {
		if ( bytes < chunksize )
			chunksize = bytes;
			fread(buf, chunksize, 1, hSource);
			fwrite(buf, chunksize, 1, hDest );
		bytes -= chunksize;
		}
	
	delete [] buf;
	}

//-----------------------------------------------------------------------------
int QTR_MoveData( QTR_Impl *impl, QTR_File *from, QTR_File *to ){
	int result = 0;

	if ( from == to )
		return 0;

	if ( ! impl->ondisk.dataCount ) {
	}
	else if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {
	}
	
	else if ( from && to ) { // one file to another
		unsigned int fromloc = impl->ondisk.dataPos;

		if ( ! fromloc ) { // nothing on disk

		}
		else {
			impl->ondisk.dataPos = to->frontier;
			to->frontier += impl->ondisk.dataSize * impl->ondisk.dataCount;
			
			if ( ! impl->data ) { // nothing loaded
				QTR_CopyF2FinChunks( from->handle, fromloc, to->handle, impl->ondisk.dataPos, impl->ondisk.dataSize * impl->ondisk.dataCount );
			}
			else {
				if ( impl->loadStart ) {
					QTR_CopyF2FinChunks( from->handle, fromloc, to->handle, impl->ondisk.dataPos, impl->ondisk.dataSize * impl->loadStart );
				}
				if ( impl->loadEnd + 1 < impl->ondisk.dataCount ) {
					QTR_CopyF2FinChunks( from->handle, fromloc + impl->ondisk.dataSize * (impl->loadEnd + 1),
										 to->handle, impl->ondisk.dataPos + impl->ondisk.dataSize * (impl->loadEnd + 1),
										 impl->ondisk.dataSize * (impl->ondisk.dataCount - impl->loadEnd) );
				}
			}
		}
	}

	else if ( from ) { // file to memory   
		unsigned int lf = impl->loadStart;
		unsigned int ll = impl->loadEnd;
		QTR_LoadRows( impl, 0, impl->ondisk.dataCount - 1, 1 );
		if ( impl->data ) { // successfully loaded -- enough memory
			impl->ondisk.dataPos = (unsigned long) impl->data;
			impl->data = ((char *) impl->data) + lf * impl->ondisk.dataSize;
			impl->loadStart = lf;
			impl->loadEnd = ll;
		}
		else {
			impl->ondisk.dataPos = 0;
			impl->ondisk.dataType = QTR_TYPE_EMPTY;
			impl->ondisk.dataSize = 0;
			impl->ondisk.dataCount = 0;
			result = -1;
		}
	}

	else { // memory to file
		QTR_LoadRows( impl, 0, impl->ondisk.dataCount - 1, 1 );
		impl->ondisk.dataPos = 0;

		// should really unload portions to disk, and restore prev.load.bounds, but...
		impl->loadStart = 0;
		impl->loadEnd = impl->ondisk.dataCount - 1;
	}

	QTR_Impl *child = impl->child;
	while ( child ) {
		result = result ? result : QTR_MoveData( child, from, to );
		child = child->sibling;
	}

	return result;
}

//-----------------------------------------------------------------------------
int QTR_ImitateData( QTR_Impl *impl, QTR_Impl *other ){
	QTR_SetFlag( impl, QTR_FLAG_PRELOAD, QTR_Flag( other, QTR_FLAG_PRELOAD ) );

	if ( ! other->ondisk.dataCount )
		return QTR_ClearData( impl );

	if ( ! QTR_SetupData( impl, (QTR_DataType) other->ondisk.dataType, other->ondisk.dataSize, other->ondisk.dataCount ) ) {
		unsigned int rowBytes = impl->ondisk.dataSize;
		unsigned int rowsAtOnce = QTR_COPYING_MAXBYTES / rowBytes;
		if ( ! rowsAtOnce )
			rowsAtOnce = 1;

		unsigned int r, remain;

		for ( r=0, remain=impl->ondisk.dataCount; r<impl->ondisk.dataCount; r += rowsAtOnce, remain -= rowsAtOnce ) {
			if ( remain < rowsAtOnce )
				rowsAtOnce = remain;

			if ( ! QTR_LoadRows( impl, r, r + rowsAtOnce - 1, 0 ) )
				QTR_GetRows( other, impl->data, r, r + rowsAtOnce - 1 );
			else
				return -1;
		}

		if ( other->loadStart != QTR_LOAD_NONE )
			QTR_LoadRows( impl, other->loadStart, other->loadEnd, true );
		else
			QTR_UnloadRows( impl, true );

		return 0;
	}
	else
		return -1;
}

//-----------------------------------------------------------------------------
void QTR_ClearFileStarts( QTR_Impl *impl ){
	impl->fileStart = 0;
	QTR_Changed( impl );
	
	QTR_Impl *child = impl->child;
	while ( child ) {
		QTR_ClearFileStarts( child );
		child = child->sibling;
	}
}

//-----------------------------------------------------------------------------
void QTR_SetRoots( QTR_Impl *impl, QTR_Impl *root ){
	((AutoMutexHolder *) impl->m_mutex)->switchMutex( ((AutoMutexHolder *) root->m_mutex)->mutex );
	impl->root = root;

	QTR_Impl *child = impl->child;
	while ( child ) {
		QTR_Lock( child, -1 );
		QTR_SetRoots( child, root );
		QTR_Unlock( child );
		child = child->sibling;
	}
}

//-----------------------------------------------------------------------------
unsigned int QTR_FindDataPos( QTR_Impl *impl ){
	if ( ! impl )
		return 0;

	if ( ! impl->root->file )
		return 0;

	if ( ! impl->ondisk.dataCount )
		return 0;

	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) )
		return impl->ondisk.dataPos;

	if ( ! impl->ondisk.dataPos ) {
		impl->ondisk.dataPos = impl->root->file->frontier;
		impl->root->file->frontier += impl->ondisk.dataSize * impl->ondisk.dataCount;
	}

	return impl->ondisk.dataPos;
}

//-----------------------------------------------------------------------------
int QTR_DumpLoadedToDisk( QTR_Impl *impl ){
	unsigned int rowBytes = impl->ondisk.dataSize;

	QTR_FindDataPos( impl );

	if ( (0==fseek(impl->root->file->handle, impl->ondisk.dataPos + rowBytes * impl->loadStart, SEEK_SET))
		&& (rowBytes == fwrite(impl->data, impl->loadEnd-impl->loadStart+1, rowBytes, impl->root->file->handle)) )
			return 0;

		return -1;
}

//-----------------------------------------------------------------------------
int QTR_LoadFromDisk( QTR_Impl *impl ){// according to load*
	if ( ! impl->ondisk.dataPos )
		return 0;

	unsigned int rowBytes = impl->ondisk.dataSize;

	if( (0 == fseek(impl->root->file->handle, impl->ondisk.dataPos + rowBytes * impl->loadStart, SEEK_SET) ) 
		&& (rowBytes == fread(impl->data, (impl->loadEnd-impl->loadStart+1), rowBytes, impl->root->file->handle)) )
		return 0;		// SUCCESS

	return -1;	// FAILURE
}

//-----------------------------------------------------------------------------
unsigned int QTR_FindFilePos( QTR_Impl *impl ){
	if ( ! impl )
		return 0;

	if ( ! impl->root->file )
		return 0;

	if ( ! impl->fileStart ) {
		impl->fileStart = impl->root->file->frontier;
		impl->root->file->frontier += sizeof(QTR_OnDisk) + impl->ondisk.nameLen;
	}

	return impl->fileStart;
}

//-----------------------------------------------------------------------------
void QTR_FreeFile( QTR_File *file ){
	fclose(file->handle);
	delete [] file->path;
	delete file;
}

//-----------------------------------------------------------------------------
void QTR_SetRowsFromFile( QTR_Impl *impl, unsigned int loc, unsigned int last ){
	if( impl->ondisk.dataCount <= 0 ) 
		return;

	QTR_Lock( impl, -1 );
	
	if ( last >= impl->ondisk.dataCount )		
		last = impl->ondisk.dataCount - 1;
	
	unsigned int rowBytes = impl->ondisk.dataSize;
	unsigned int rowsAtOnce = IMAX(1, QTR_COPYING_MAXBYTES / rowBytes) ;
	unsigned int iRow, iRows;

	if ( 0==fseek(impl->root->file->handle, loc, SEEK_SET ) ) {
		char * buf = new char[ rowsAtOnce * rowBytes ];
		for( iRow=0; iRow<=last; iRow += rowsAtOnce ) {
			iRows = IMIN( last-iRow+1, rowsAtOnce );
			if( iRows==fread(buf, rowBytes, iRows, impl->root->file->handle) )
				QTR_SetRows( impl, buf, iRow, iRow+iRows-1 );
			}
		delete [] buf;
		}
	
	QTR_Unlock( impl );
	}


void QTR_Stat(QTR_Impl *impl) {
  std::cerr << impl->name << "(" << impl->refs << "): " << impl->ondisk.dataCount << ", " << impl->child << ", " << impl->sibling << std::endl;
}


//===================================================================================
// Exported functions

//-----------------------------------------------------------------------------
QTR_API int QTR_INCREF( QTR_Impl *impl ){
	QTR_Lock( impl, -1 );
	int iRefs = ++impl->refs;
	QTR_Unlock( impl );
	return iRefs;		
}

//-----------------------------------------------------------------------------
QTR_API int QTR_DECREF( QTR_Impl *impl ){
	QTR_Lock( impl, -1 );
	int refsRemaining = --impl->refs;
	// QTR_Stat(impl);
	QTR_Unlock( impl );
	if ( ! refsRemaining )
		QTR_Free( impl );
	return refsRemaining;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_Lock( QTR_Impl *impl, int timeoutMS ){
	return ((AutoMutexHolder *) impl->m_mutex)->lock( timeoutMS );
}

//-----------------------------------------------------------------------------
QTR_API void QTR_Unlock( QTR_Impl *impl ){
	((AutoMutexHolder *) impl->m_mutex)->unlock();
}

//-----------------------------------------------------------------------------
QTR_API QTR_Impl * QTR_Create( const char *name ){
	if ( namesInterned == NULL )
		namesInterned = new StringsInterned;

	QTR_Impl * impl = New_QTR_Impl(1);
	impl->ondisk.nameLen = name ? strlen(name) : 0;
	impl->name           = namesInterned->get( name ? name : "" ); // new char[ impl->ondisk.nameLen + 1 ];
	//impl->m_mutex        = new CAutoMutexHolder( true );
	impl->m_mutex          = new AutoMutexHolder( new AutoMutex );
	impl->root           = impl;

	return impl;
}

//-----------------------------------------------------------------------------
QTR_API QTR_Impl * QTR_Clone( QTR_Impl *other, int deep ){
	if ( ! other )
		return 0;

	QTR_Lock( other, -1 );

	QTR_Impl * impl = New_QTR_Impl(1);
	impl->ondisk.nameLen = other->ondisk.nameLen;
	impl->name = namesInterned->dup( other->name );
	impl->m_mutex = new AutoMutexHolder( new AutoMutex );
	impl->root = impl;

	QTR_ImitateData( impl, other );

	if ( deep ) {
		QTR_Impl *lastChild = 0;
		QTR_Impl *otherChild = other->child;
		while ( otherChild ) {
			lastChild = QTR_InsertClone( impl, lastChild, otherChild, deep );
			otherChild = otherChild->sibling;
		}
	}

	QTR_Unlock( other );

	return impl;
}

//-----------------------------------------------------------------------------
QTR_API QTR_Impl * QTR_Open( const char *path, int readOnly ){
    QTR_Impl *impl = NULL;
    QTR_File *file = new QTR_File;
    try {
	const char *perm = readOnly ? "rb" : "r+b";
	file->handle = fopen(path, perm);
	if ( ! file->handle ) {
		delete file;
		return 0;
	}

	char magic[ QTR_MAGIC_LEN + 1 ];
	magic[ QTR_MAGIC_LEN ] = 0;
	int haveMagic = fread( magic, QTR_MAGIC_LEN, 1, file->handle );

	if ( (! haveMagic) || (strcmp(magic, QTR_MAGIC) != 0) ) {
		fclose( file->handle );
		delete file;
		return 0;
	}

	impl = QTR_Create( 0 );
	impl->file = file;
	impl->fileStart = QTR_MAGIC_LEN;
	
	QTR_Read( impl );

	int pathlen = strlen(path);
	file->path = new char[ pathlen+1 ];
	strcpy( file->path, path );

	fseek(file->handle, 0, SEEK_END);
	file->frontier = ftell(file->handle);

	return impl;
    } catch (...) {
        if ( impl )
	  QTR_DECREF( impl );
        return NULL;
    }
}

//-----------------------------------------------------------------------------
QTR_API void QTR_Close( QTR_Impl *impl ){
	if ( impl->file ) {
		QTR_Lock( impl, -1 );

		QTR_MoveData( impl, impl->file, 0 );
		QTR_FreeFile( impl->file );
		impl->file = 0;
		impl->fileStart = 0;

		QTR_Unlock( impl );
	}
}

//-----------------------------------------------------------------------------
QTR_API int QTR_Save( QTR_Impl *impl ){
	if ( ! impl->root->file )
		return 0;

	QTR_Lock( impl, -1 );

	FILE *f = impl->root->file->handle;

	int changed = QTR_Flag( impl, QTR_FLAG_CHANGED );
    int childChanged = QTR_Flag( impl, QTR_FLAG_CHILDCHANGED );
	QTR_SetFlag( impl, QTR_FLAG_CHANGED, 0 );
	QTR_SetFlag( impl, QTR_FLAG_CHILDCHANGED, 0 );
  
	if ( changed ) {
		// first make sure the child/sibling pos on disk is correct
		impl->ondisk.childPos = QTR_FindFilePos( impl->child );
		impl->ondisk.siblingPos = QTR_FindFilePos( impl->sibling );
		impl->ondisk.dataPos = QTR_FindDataPos( impl );

		fseek(f, impl->fileStart, SEEK_SET);
		fwrite(&(impl->ondisk), sizeof(QTR_OnDisk), 1, f);
		if ( impl->ondisk.nameLen )
			fwrite(impl->name, impl->ondisk.nameLen, 1, f);
		
		if ( impl->data && ! QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) )
			QTR_DumpLoadedToDisk( impl );

		// in case they set up non-preload data and never did anything with it;
		// make sure it's not pointing at file begin
		// (even though this might point past end of file)
		if ( impl->ondisk.dataCount && ! impl->ondisk.dataPos )
			QTR_FindDataPos( impl );
	}

	if ( childChanged ) {
		QTR_Impl *child = impl->child;
		while ( child ) {
			QTR_Save( child );
			child = child->sibling;
		}
	}

	QTR_Unlock( impl );

	if ( impl == impl->root )
		fflush( impl->file->handle );

	return 0;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_SaveAs( QTR_Impl *impl, const char *path ){
 // if called on non-root, removes it to be a root first
	QTR_File *file = new QTR_File;
	file->handle = fopen(path, "w+b");
	if ( ! file->handle ) {
		delete file;
		return -1;
	}

	fwrite(QTR_MAGIC, QTR_MAGIC_LEN, 1, file->handle);

	QTR_OnDisk decoyRoot;
	memset(&decoyRoot, 0, sizeof(QTR_OnDisk));
	decoyRoot.nameLen = impl->ondisk.nameLen;
	fwrite(&decoyRoot, sizeof(QTR_OnDisk), 1, file->handle);
	if ( decoyRoot.nameLen )
		fwrite( impl->name, decoyRoot.nameLen, 1, file->handle );

	int pathlen = strlen(path);
	file->path = new char[ pathlen + 1 ];
	strcpy( file->path, path );

	file->frontier = ftell(file->handle);

	QTR_Lock( impl, -1 );

	QTR_MoveData( impl, impl->root->file, file );
	QTR_ClearFileStarts( impl );
	
	if ( impl->root == impl ) {
		if ( impl->file ) {
			QTR_FreeFile( impl->file );
		}
	}
	else {
		QTR_Impl *oldParent = impl->parent;
		QTR_Lock( oldParent, -1 );
		QTR_Unlink( impl );
		((AutoMutexHolder *) impl->m_mutex)->switchMutex( new AutoMutex );
		QTR_SetRoots( impl, impl );
		QTR_Unlock( oldParent );
	}

	impl->file = file;
	impl->fileStart = QTR_MAGIC_LEN;

	QTR_Unlock( impl );

	return 0;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_SaveAsTemp( QTR_Impl *impl ){
 // if called on non-root, removes it to be a root first
	QTR_File *file = new QTR_File;
	file->handle = tmpfile();
	if ( ! file->handle ) {
		delete file;
		return -1;
	}

	fwrite(QTR_MAGIC, QTR_MAGIC_LEN, 1, file->handle);

	QTR_OnDisk decoyRoot;
	memset(&decoyRoot, 0, sizeof(QTR_OnDisk));
	decoyRoot.nameLen = impl->ondisk.nameLen;
	fwrite(&decoyRoot, sizeof(QTR_OnDisk), 1, file->handle);
	if ( decoyRoot.nameLen )
		fwrite( impl->name, decoyRoot.nameLen, 1, file->handle );

	file->path = new char[ 1 ];
	file->path[0] = 0;

	file->frontier = ftell(file->handle);

	QTR_Lock( impl, -1 );

	QTR_MoveData( impl, impl->root->file, file );
	QTR_ClearFileStarts( impl );
	
	if ( impl->root == impl ) {
		if ( impl->file ) {
			QTR_FreeFile( impl->file );
		}
	}
	else {
		QTR_Lock( impl->parent, -1 );
		QTR_Unlink( impl );
		((AutoMutexHolder *) impl->m_mutex)->switchMutex( new AutoMutex );
		QTR_SetRoots( impl, impl );
		QTR_Unlock( impl->parent );
	}

	impl->file = file;
	impl->fileStart = QTR_MAGIC_LEN;

	QTR_Unlock( impl );

	return 0;
}

//-----------------------------------------------------------------------------
QTR_API QTR_Impl *  QTR_CreateChild( QTR_Impl *parent, QTR_Impl *presib, const char *name ){
	QTR_Impl * impl = New_QTR_Impl(0);
	impl->ondisk.nameLen = name ? strlen(name) : 0;
	impl->m_mutex          = new AutoMutexHolder(0); // , impl);
	impl->name           = namesInterned->get( name ? name : "" ); // new char[ impl->ondisk.nameLen + 1 ];

	QTR_Lock( parent, -1 );
	QTR_Lock( impl, -1 );
	QTR_SetRoots( impl, parent->root );
	QTR_Link( parent, presib, impl );
	QTR_Unlock( impl );
	QTR_Unlock( parent );

	return impl;
}

//-----------------------------------------------------------------------------
QTR_API QTR_Impl * QTR_InsertClone( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *other, int deep ){
	if ( ! other )
		return 0;

	QTR_Lock( other, -1 );
	QTR_Lock( parent, -1 );

	// QTR_Impl *impl     = QTR_CreateChild( parent, presib, other->name );
	QTR_Impl * impl = New_QTR_Impl(0);
	impl->ondisk.nameLen = other->ondisk.nameLen;
	impl->m_mutex          = new AutoMutexHolder(0); // , impl);
	impl->name           = namesInterned->dup( other->name );

	QTR_Lock( parent, -1 );
	QTR_Lock( impl, -1 );
	QTR_SetRoots( impl, parent->root );
	QTR_Link( parent, presib, impl );
	QTR_Unlock( impl );
	QTR_Unlock( parent );

	QTR_ImitateData( impl, other );

	if ( deep ) {
		QTR_Impl *lastChild = 0;
		QTR_Impl *otherChild = other->child;
		while ( otherChild ) {
			lastChild = QTR_InsertClone( impl, lastChild, otherChild, deep );
			otherChild = otherChild->sibling;
		}
	}

	QTR_Unlock( parent );
	QTR_Unlock( other );

	return impl;
}

//-----------------------------------------------------------------------------
QTR_API void QTR_InsertChild( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child ){
	if ( ! child )
		return;

	QTR_Lock( parent, -1 );
	QTR_Lock( child, -1 );

	if ( parent->root->file != child->root->file ) {
		QTR_MoveData( child, child->root->file, parent->root->file );
		child->fileStart = 0;
	}
	
	if ( child->parent ) { // the inserted node already has a parent
		QTR_Unlink( child );
		QTR_Link( parent, presib, child );
		QTR_SetRoots( child, parent->root );
	}
	else { // child was a root node
		if ( child->file ) {
			QTR_FreeFile( child->file );
			child->file = 0;
		}
		QTR_SetRoots( child, parent->root );
		QTR_Link( parent, presib, child );
	}

	QTR_Unlock( child );
	QTR_Unlock( parent );
}

//-----------------------------------------------------------------------------
QTR_API void QTR_RemoveChild( QTR_Impl *parent, QTR_Impl *presib, QTR_Impl *child ){
	if ( (! parent) || (! child) || (parent != child->parent) )
		return;

	QTR_Lock( parent, -1 );
	QTR_Lock( child, -1 );
	
	if (   (presib == NULL && child != parent->child)
		|| (presib != NULL && presib->sibling != child) ) {
		presib = parent->child;
		while ( presib->sibling != child )
			presib = presib->sibling;
	}

	if ( child->refs > 1 ) {
		QTR_MoveData( child, parent->root->file, 0 );
		((AutoMutexHolder *) child->m_mutex)->switchMutex( new AutoMutex );
		QTR_SetRoots( child, child );
	}
	else {
		QTR_ClearData( child );
	}
	if ( child->refs > 1 ) {
		QTR_Unlink( parent, presib, child );
		child->fileStart = 0;
		QTR_Unlock( child );
	}
	else {
		QTR_Unlock( child );
		QTR_Unlink( parent, presib, child ); // will free -- do unlocked
	}

	QTR_Unlock( parent );
}

//-----------------------------------------------------------------------------
QTR_API int QTR_Flag( QTR_Impl *impl, int mask ) {
  int hasFlag = ( *((int *) & impl->ondisk.flags) & mask ) ? 1 : 0;
  return hasFlag;
}

//-----------------------------------------------------------------------------
QTR_API void QTR_SetFlag( QTR_Impl *impl, int mask, int val ) {
  QTR_Lock( impl, -1 );
  int& flags = *((int *) & impl->ondisk.flags);
  flags = (flags & ~ mask) | (val ? mask : 0);
  QTR_Unlock( impl );
  
  if ( ! (mask & (QTR_FLAG_CHANGED | QTR_FLAG_CHILDCHANGED) ) )
    QTR_Changed( impl );
}

//-----------------------------------------------------------------------------
QTR_API void QTR_Changed( QTR_Impl *impl ) {
  QTR_SetFlag( impl, QTR_FLAG_CHANGED, 1 );
  if ( impl->parent )
    QTR_ChildChanged( impl->parent );
}

//-----------------------------------------------------------------------------
QTR_API void QTR_ChildChanged( QTR_Impl *impl ) {
	//TODO - to avoid repeated walks up same chain, if ChildChanged 
	// is already set, assume all parents are also alredy set and exit.

  QTR_SetFlag( impl, QTR_FLAG_CHILDCHANGED, 1 );
  if ( impl->parent )
    QTR_ChildChanged( impl->parent );
}

//-----------------------------------------------------------------------------
QTR_API char * QTR_Name( QTR_Impl *impl ) {
	return impl->name;
}

//-----------------------------------------------------------------------------
QTR_API void QTR_SetName( QTR_Impl *impl, const char *name ) {
	if ( impl->file )
		return; // sorry, root nodes can't have name changes

	QTR_Lock( impl, -1 );

	int newNameLen = strlen(name);
	if ( impl->ondisk.nameLen < newNameLen )
		impl->fileStart = 0; // force re-location in file to a spot with enough space

	namesInterned->release( impl->name );
	impl->name = namesInterned->get( name );
	impl->ondisk.nameLen = newNameLen;

	QTR_Unlock( impl );
	QTR_Changed( impl );
}

//-----------------------------------------------------------------------------
QTR_API int QTR_IsInFile( QTR_Impl *impl ) {
	return impl->root->file ? 1 : 0;
}

//-----------------------------------------------------------------------------
QTR_API char * QTR_FilePath( QTR_Impl *impl ){
	return impl->root->file ? impl->root->file->path : 0;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_SetupData( QTR_Impl *impl, QTR_DataType dtype, unsigned int size, unsigned int count ){
	if ( ! (size && count) )
		return QTR_ClearData( impl );

	QTR_Lock( impl, -1 );
	QTR_Changed( impl );

	unsigned int curBytes = impl->ondisk.dataSize * impl->ondisk.dataCount;
	unsigned int reqBytes = size * count;
	bool curDataInNode = (0!=QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ));

	QTR_UnloadRows( impl, 0 );

	impl->ondisk.dataType = dtype;
	impl->ondisk.dataSize = size;
	impl->ondisk.dataCount = count;

	if ( curBytes && (curBytes < reqBytes) ) {
		if ( curDataInNode ) {
			QTR_SetFlag( impl, QTR_FLAG_DATA_IN_NODE, 0 );
			impl->ondisk.dataPos = 0;
		}
		else if ( impl->root->file ) {
			impl->ondisk.dataPos = 0;
		}
		else {
			delete [] (char *) impl->ondisk.dataPos;
		}
		curBytes = 0;
	}

	if ( ! curBytes ) {
		// if data fits in node...	Note : (size == 1) add space for a string's null terminator
		if ( size * count + (size==1?1:0) < sizeof(unsigned int) ) {	
			QTR_SetFlag( impl, QTR_FLAG_DATA_IN_NODE, 1 );
		}
		else if ( ! impl->root->file ) {
			char *data = new char[ reqBytes + 1 ];

			if ( data ) {
				data[ reqBytes ] = 0; // in case of string
				impl->ondisk.dataPos = (unsigned long) data;
			}
			else {
				impl->ondisk.dataType = QTR_TYPE_EMPTY;
				impl->ondisk.dataSize = 0;
				impl->ondisk.dataCount = 0;
				impl->ondisk.dataPos = 0;
			}		
		}
	}

	if ( QTR_Flag( impl, QTR_FLAG_PRELOAD ) )
		QTR_LoadRows( impl, 0, count-1, 0 );

	QTR_Unlock( impl );

	return (count == impl->ondisk.dataCount) ? 0 : -1;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_SetupNumData( QTR_Impl *impl, QTR_DataType type, unsigned int nr, unsigned int nc ){
	return QTR_SetupData( impl, type, nc * QTR_DataSize[ (int) type ], nr );
}

//-----------------------------------------------------------------------------
QTR_API int QTR_SetupStringData( QTR_Impl *impl, int count ) {
	return QTR_SetupData( impl, QTR_TYPE_STRING, 1, count );
}

//-----------------------------------------------------------------------------
QTR_API int QTR_ResizeData( QTR_Impl *impl, unsigned int newnr ){
	if ( ! newnr )
		return QTR_ClearData( impl );

	unsigned int prevnr = impl->ondisk.dataCount;
	if ( ! prevnr )
		return -1; // need valid type, size to resize with

	QTR_Changed( impl );

	if ( newnr <= prevnr ) {
		QTR_Lock( impl, -1 );
		impl->ondisk.dataCount = newnr;
		if ( impl->data && (impl->loadEnd >= newnr) )
			impl->loadEnd = newnr - 1;
		QTR_Unlock( impl );
		return 0;
	}

	int result = 0;
	QTR_Lock( impl, -1 );

	unsigned int rowBytes = impl->ondisk.dataSize;
	bool inNode = (0!=QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ));
	bool inMem = ! impl->root->file;
	unsigned int temp;
	char *prevData;
	unsigned int prevLoc;

	// unload:
	unsigned int restoreLoadFirst = impl->loadStart;
	unsigned int restoreLoadLast = impl->loadEnd;

	if ( QTR_Flag( impl, QTR_FLAG_PRELOAD ) && (restoreLoadFirst == 0) && (restoreLoadLast == prevnr - 1) ) {
		restoreLoadLast = newnr - 1;
	}
	QTR_UnloadRows( impl, 1 );

	// set aside old data:
	if ( inNode ) {
		prevData = (char *) &temp;
		temp = impl->ondisk.dataPos;
	}
	else if ( inMem ) {
		prevData = (char *) impl->ondisk.dataPos;
	}
	else {
		prevLoc = impl->ondisk.dataPos;
	}

	// set up new data
	impl->ondisk.dataCount = newnr;
	QTR_SetFlag( impl, QTR_FLAG_DATA_IN_NODE, 0 );
	if ( inMem ) {
		impl->ondisk.dataPos = (unsigned long) new char[ rowBytes * newnr + 1 ];
	}
	else {
		impl->ondisk.dataPos = 0;
	}

	if ( inMem && (! inNode) && (! impl->ondisk.dataPos) ) {
		impl->ondisk.dataPos = (unsigned long) prevData;
		impl->ondisk.dataCount = prevnr;
		if ( restoreLoadFirst != QTR_LOAD_NONE )
			QTR_LoadRows( impl, restoreLoadFirst, restoreLoadLast, 0 );

		result = -2;
	}
	else {
		// reload:
		if ( restoreLoadFirst != QTR_LOAD_NONE )
			QTR_LoadRows( impl, restoreLoadFirst, restoreLoadLast, 0 );

		// copy in old data:
		if ( inNode || inMem )
			QTR_SetRows( impl, prevData, 0, prevnr - 1 );
		else
			QTR_SetRowsFromFile( impl, prevLoc, prevnr-1 );

		// release old data:
		if ( inMem && ! inNode )
			delete [] prevData;
	}

	QTR_Unlock( impl );

	return result;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_ClearData( QTR_Impl *impl ){
	if ( ! impl->ondisk.dataCount )
		return 0;

	QTR_Lock( impl, -1 );
	QTR_UnloadRows( impl, 0 );
	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {
		QTR_SetFlag( impl, QTR_FLAG_DATA_IN_NODE, 0 );
	}
	else if ( ! impl->root->file ) {
		delete [] (char *) impl->ondisk.dataPos;
	}
	
	impl->ondisk.dataPos = 0;
	impl->ondisk.dataType = QTR_TYPE_EMPTY;
	impl->ondisk.dataSize = 0;
	impl->ondisk.dataCount = 0;

	QTR_Unlock( impl );
	QTR_Changed( impl );

	return 0;
}


//-----------------------------------------------------------------------------
QTR_API unsigned int QTR_NumCol( QTR_Impl *impl ){
	unsigned int numcol = 0;
	if ( ((int) impl->ondisk.dataType) > ((int) QTR_TYPE_UCHAR) )	// DataIsNumeric ? 
		numcol = impl->ondisk.dataSize / QTR_DataSize[ impl->ondisk.dataType ];
	if ( ! numcol )
		numcol = 1;
	return numcol;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_LoadRows( QTR_Impl *impl, unsigned int firstRow, unsigned int lastRow, int doRead ){
	QTR_Lock( impl, -1 );

	firstRow = IMAX( firstRow, 0 );
	firstRow = IMIN( firstRow, impl->ondisk.dataCount - 1 );
	lastRow = IMAX( firstRow, lastRow );
	lastRow = IMIN( lastRow, impl->ondisk.dataCount - 1 );

	unsigned int loadedFirst = impl->loadStart, loadedLast = impl->loadEnd;
	unsigned int rowCount = (lastRow - firstRow + 1);
	unsigned int byteCount = rowCount * impl->ondisk.dataSize;
	unsigned int byteStart = firstRow * impl->ondisk.dataSize;

	if ( (firstRow == loadedFirst) && (lastRow == loadedLast) ) {
		QTR_Unlock( impl );
		return 0;
	}

	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {
		impl->loadStart = firstRow;
		impl->loadEnd = lastRow;
		impl->data = (char *) & impl->ondisk.dataPos;
		QTR_Unlock( impl );
		return 0;
	}

	if ( ! impl->root->file ) {
		impl->loadStart = firstRow;
		impl->loadEnd = lastRow;
		impl->data = ((char *) impl->ondisk.dataPos) + byteStart;
		QTR_Unlock( impl );
		return 0;
	}

	if ( impl->data ) {
		if ( rowCount == (loadedLast - loadedFirst + 1) ) {
			QTR_DumpLoadedToDisk( impl );
		}
		else {
			QTR_UnloadRows( impl, 1 );
			delete [] (char*) impl->data;
			try {
			  impl->data = new char[ byteCount + 1 ];
			} catch (std::bad_alloc& ba) {
			  impl->data = NULL;
			}
		}
	}
	else {
	  try {
		impl->data = new char[ byteCount + 1 ];
	  } catch (std::bad_alloc& ba) {
	        impl->data = NULL;
	  }
	}

	if ( ! impl->data ) {
		impl->loadStart = impl->loadEnd = QTR_LOAD_NONE;
		QTR_Unlock( impl );
		return -2;
	}

	((char *) impl->data)[ byteCount ] = 0; // to make strings play nice
	impl->loadStart = firstRow;
	impl->loadEnd = lastRow;

	QTR_Unlock( impl );

	if ( doRead )
		return QTR_LoadFromDisk( impl );
	else
		return 0;
}

//-----------------------------------------------------------------------------
QTR_API int QTR_UnloadRows( QTR_Impl *impl, int doWrite ){
	if ( ! impl->data )
		return 0;

	QTR_Lock( impl, -1 );

	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {

	}
	else if ( impl->root->file ) {
		if ( doWrite ) {
			QTR_DumpLoadedToDisk( impl );
		}
		delete [] (char*) impl->data;
	}

	impl->loadStart = impl->loadEnd = QTR_LOAD_NONE;
	impl->data = 0;

	QTR_Unlock( impl );

	return 0;
}

//-----------------------------------------------------------------------------
QTR_API unsigned int QTR_GetRows( QTR_Impl *impl, void *buf, unsigned int first, unsigned int last ){
	QTR_Lock( impl, -1 );

	if ( first < 0 )							first = 0;
	if ( first >= impl->ondisk.dataCount )		first = impl->ondisk.dataCount - 1;
	if ( last < first )							last = first;
	if ( last >= impl->ondisk.dataCount )		last = impl->ondisk.dataCount - 1;

	if ( first < 0 ) { // (no data)
		QTR_Unlock( impl );
		return 0;
	}

	unsigned int rowBytes = impl->ondisk.dataSize;

	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {
		memcpy( buf, & impl->ondisk.dataPos, rowBytes * (last - first + 1) );
		QTR_Unlock( impl );
		return (last - first + 1);
	}
	else if ( impl->root->file ) {
		FILE *f = impl->root->file->handle;
		if ( ! impl->data ) {
			unsigned int rtn = 0;
			if ( ! fseek(f, impl->ondisk.dataPos + first * rowBytes, SEEK_SET) )
				rtn = fread(buf, rowBytes, last - first + 1, f);
			QTR_Unlock( impl );
			return rtn;
		}
		else {
			unsigned int bstart = IMAX( impl->loadStart, first );
			unsigned int bend = IMIN( impl->loadEnd, last );
			unsigned int d1end = IMIN( impl->loadStart, last + 1 ) - 1;
			unsigned int d2start = IMAX( impl->loadEnd + 1, first );
			char *ldata = (char *) impl->data;
			char *xbuf = (char *) buf;
			unsigned int rows = d1end - first + 1;
			unsigned int nread = 0;
			
			if ( rows > 0 ) {
				if ( ! fseek(f, impl->ondisk.dataPos + first * rowBytes, SEEK_SET) )
					nread += fread(xbuf, rowBytes, rows, f);
				xbuf += rowBytes * rows;
			}
  
			rows = bend - bstart + 1;
			if ( rows > 0 ) {
				memcpy( xbuf, ldata + rowBytes * (bstart - impl->loadStart), rowBytes * rows );
				xbuf +=  rowBytes * rows;
				nread += (bstart - impl->loadStart);
			}

			rows = last - d2start + 1;
			if ( rows > 0 ) {
				if ( ! fseek(f, impl->ondisk.dataPos + d2start * rowBytes, SEEK_SET) )
					nread += fread(xbuf, rowBytes, rows, f);
			}

			QTR_Unlock( impl );
			return nread;
		}
	}
	else {
		memcpy( buf, ((char *) impl->ondisk.dataPos) + first * rowBytes, (last - first + 1) * rowBytes );
		QTR_Unlock( impl );
		return ( last - first + 1 );
	}
}

//-----------------------------------------------------------------------------
QTR_API unsigned int QTR_SetRows( QTR_Impl *impl, const void *buf, unsigned int first, unsigned int last ){
	QTR_Lock( impl, -1 );

	if ( first < 0 )							first = 0;
	if ( first >= impl->ondisk.dataCount )		first = impl->ondisk.dataCount - 1;
	if ( last < first )							last = first;
	if ( last >= impl->ondisk.dataCount )		last = impl->ondisk.dataCount - 1;

	if ( first < 0 ) { // (no data)
		QTR_Unlock( impl );
		return 0;
	}

	QTR_Changed( impl );
	unsigned int rowBytes = impl->ondisk.dataSize;

	if ( QTR_Flag( impl, QTR_FLAG_DATA_IN_NODE ) ) {
		memcpy( & impl->ondisk.dataPos, buf, rowBytes * (last - first + 1) );
		QTR_Unlock( impl );
		return (last - first + 1);
	}
	else if ( impl->root->file ) {
		FILE *f = impl->root->file->handle;
		QTR_FindDataPos( impl );

		if ( ! impl->data ) {
			unsigned int rtn = 0;
			if ( ! fseek(f, impl->ondisk.dataPos + first * rowBytes, SEEK_SET) )
				rtn = fwrite(buf, rowBytes, last - first + 1, f);
			QTR_Unlock( impl );
			return rtn;
		}
		else {
			unsigned int bstart = IMAX( impl->loadStart, first );
			unsigned int bend = IMIN( impl->loadEnd, last );
			unsigned int d1end = IMIN( impl->loadStart, last + 1 ) - 1;
			unsigned int d2start = IMAX( impl->loadEnd + 1, first );
			char *ldata = (char *) impl->data;
			char *xbuf = (char *) buf;
			unsigned int rows;
			unsigned int nwritten = 0;
  
			rows = d1end - first + 1;
			if ( rows > 0 ) {
				if ( ! fseek(f, impl->ondisk.dataPos + first * rowBytes, SEEK_SET) )
					nwritten += fwrite(xbuf, rowBytes, rows, f);
				xbuf += rowBytes * rows;
			}
  
			rows = bend - bstart + 1;
			if ( rows > 0 ) {
				memcpy( ldata + rowBytes * (bstart - impl->loadStart), xbuf, rowBytes * rows );
				xbuf +=  rowBytes * rows;
				nwritten += (bstart - impl->loadStart);
			}

			rows = last - d2start + 1;
			if ( rows > 0 ) {
				if ( ! fseek(f, impl->ondisk.dataPos + d2start * rowBytes, SEEK_SET) )
					nwritten += fwrite(xbuf, rowBytes, rows, f);
			}

			QTR_Unlock( impl );
			return nwritten;
		}
	}
	else {
		memcpy( ((char *) impl->ondisk.dataPos) + first * rowBytes, buf, (last - first + 1) * rowBytes );
		QTR_Unlock( impl );
		return ( last - first + 1 );
	}
}



unsigned int QTR_ToBufferSizeRecur( QTR_Impl *impl )
{
	unsigned int size = 0;
	size += sizeof(QTR_OnDisk);
	size += impl->ondisk.nameLen;
	if ( impl->ondisk.dataCount && ! QTR_Flag(impl, QTR_FLAG_DATA_IN_NODE) )
		size += impl->ondisk.dataSize * impl->ondisk.dataCount;

	for (QTR_Impl *child = impl->child; child; child = child->sibling) {
		size += QTR_ToBufferSizeRecur( child );
	}

	return size;
}

QTR_API unsigned int   QTR_ToBufferSize( QTR_Impl *impl )
{
	if ( ! impl )
		return 0;

	unsigned int size = 0;
	size += QTR_MAGIC_LEN;
	size += QTR_ToBufferSizeRecur( impl );
	return size;
}

unsigned int QTR_ToBufferRecur( QTR_Impl *impl, char *buffer, int offset, bool root )
{
	memcpy(buffer + offset, &(impl->ondisk), sizeof(QTR_OnDisk));
	QTR_OnDisk& onbuf = * ((QTR_OnDisk*)(buffer + offset));
	offset += sizeof(QTR_OnDisk);

	memcpy(buffer + offset, impl->name, onbuf.nameLen);
	offset += onbuf.nameLen;

	if ( impl->ondisk.dataCount ) {
		if ( ! QTR_Flag(impl, QTR_FLAG_DATA_IN_NODE) )  {
			onbuf.dataPos = offset;
			QTR_GetRows(impl, buffer + offset, 0, onbuf.dataCount-1);
			offset += onbuf.dataSize * onbuf.dataCount;
		}
	}
	else {
		onbuf.dataPos = 0; // just to be safe
	}

	if ( impl->child ) {
		onbuf.childPos = offset;
		offset = QTR_ToBufferRecur( impl->child, buffer, offset, false );
	}
	else {
		onbuf.childPos = 0; // to be safe
	}
	if ( impl->sibling && ! root ) {
		// there can be no sibling for the root node
		onbuf.siblingPos = offset;
		offset = QTR_ToBufferRecur( impl->sibling, buffer, offset, false );
	}
	else {
		onbuf.siblingPos = 0;
	}

	return offset;
}

QTR_API unsigned int   QTR_ToBuffer( QTR_Impl *impl, char *buffer )
{
	if ( ! impl )
		return 0;

	memcpy(buffer, QTR_MAGIC, QTR_MAGIC_LEN);
	return QTR_ToBufferRecur( impl, buffer, QTR_MAGIC_LEN, true );
}

int QTR_FlagFromBuffer(QTR_OnDisk& onbuf, int mask)
{
	return ( ( *((int *) & onbuf.flags) & mask ) ? 1 : 0 );
}

// pre: impl is positioned in the tree and named
// offset is positioned just after node name (where data should be)
// op: copy reserved, data, flags; recur children and siblings
void QTR_FromBufferRecur( const char *buffer, QTR_Impl *impl, QTR_OnDisk& onbuf, unsigned int offset )
{
	for ( int i=0; i<QTR_NRSV; ++i )
		impl->ondisk.reserved[i] = onbuf.reserved[i]; //?

	if ( onbuf.dataCount ) {
		QTR_SetFlag(impl, QTR_FLAG_PRELOAD, QTR_FlagFromBuffer(onbuf, QTR_FLAG_PRELOAD));
		if ( QTR_FlagFromBuffer(onbuf, QTR_FLAG_DATA_IN_NODE) ) {
			QTR_SetFlag(impl, QTR_FLAG_DATA_IN_NODE, 1);
			impl->ondisk.dataCount = onbuf.dataCount;
			impl->ondisk.dataPos = onbuf.dataPos;
			impl->ondisk.dataSize = onbuf.dataSize;
			impl->ondisk.dataType = onbuf.dataType;
			QTR_LoadRows(impl, 0, onbuf.dataCount-1, 1/*doRead*/);
		}
		else {
			QTR_SetupData(impl, QTR_DataType(onbuf.dataType), onbuf.dataSize, onbuf.dataCount);
			QTR_SetRows(impl, buffer + offset, 0, onbuf.dataCount-1);
			offset += onbuf.dataSize * onbuf.dataCount;
		}
	}
	
	if ( onbuf.childPos ) {
		offset = onbuf.childPos;
		QTR_OnDisk& child_onbuf = * ((QTR_OnDisk*)(buffer + offset));
		offset += sizeof(QTR_OnDisk);

		char *child_name = new char[child_onbuf.nameLen+1];
		memcpy(child_name, buffer+offset, child_onbuf.nameLen);
		offset += child_onbuf.nameLen;
		child_name[child_onbuf.nameLen] = 0;
		QTR_Impl *child = QTR_CreateChild(impl, NULL, child_name);
		delete [] child_name;

		QTR_FromBufferRecur( buffer, child, child_onbuf, offset );
	}
	if ( onbuf.siblingPos ) {
		offset = onbuf.siblingPos;
		QTR_OnDisk& sibling_onbuf = * ((QTR_OnDisk*)(buffer + offset));
		offset += sizeof(QTR_OnDisk);

		char *sibling_name = new char[sibling_onbuf.nameLen+1];
		memcpy(sibling_name, buffer+offset, sibling_onbuf.nameLen);
		offset += sibling_onbuf.nameLen;
		sibling_name[sibling_onbuf.nameLen] = 0;
		QTR_Impl *sibling = QTR_CreateChild(impl->parent, impl, sibling_name);
		delete [] sibling_name;

		QTR_FromBufferRecur( buffer, sibling, sibling_onbuf, offset );
	}
}

QTR_API QTR_Impl *     QTR_FromBuffer( const char *buffer )
{
	unsigned int offset = 0;
	if (0 != strncmp(buffer, QTR_MAGIC, QTR_MAGIC_LEN))
		return NULL;
	offset += QTR_MAGIC_LEN;
	
	QTR_OnDisk& onbuf = * ((QTR_OnDisk*)(buffer + offset));
	offset += sizeof(QTR_OnDisk);

	char *name = new char[onbuf.nameLen+1];
	memcpy(name, buffer+offset, onbuf.nameLen);
	offset += onbuf.nameLen;
	name[onbuf.nameLen] = 0;
	QTR_Impl *root = QTR_Create( name );
	delete [] name;

	QTR_FromBufferRecur( buffer, root, onbuf, offset );

	return root;
}


