#include "notemplatewarning.h"

#ifdef _USRDLL
  #define _QTR_
#endif

#include "QUB_Tree.h"
#include <strstream>
#include <fstream>


void QTR_API trickMSVC6intoExportingTemplateMemberFunctions() {
	QUB_Tree crap;

	crap.setData(QTR_TYPE_UCHAR,     (unsigned char)     0  );
	crap.setData(QTR_TYPE_CHAR,      (char)              0  );
	crap.setData(QTR_TYPE_USHORT,    (unsigned short)    0  );
	crap.setData(QTR_TYPE_SHORT,     (short)             0  );
	crap.setData(QTR_TYPE_UINT,      (unsigned int)      0  );
	crap.setData(QTR_TYPE_INT,       (int)               0  );
	crap.setData(QTR_TYPE_ULONG,     (unsigned long long)0  );
	crap.setData(QTR_TYPE_LONG,      (long long)         0  );
	crap.setData(QTR_TYPE_FLOAT,     (float)             0.0);
	crap.setData(QTR_TYPE_DOUBLE,    (double)            0.0);
	crap.setData(QTR_TYPE_LDOUBLE,   (long double)       0.0);

	crap.dataAs( 0, 0, (unsigned char)   0   );
	crap.dataAs( 0, 0, (char)            0   );
	crap.dataAs( 0, 0, (unsigned short)  0   );
	crap.dataAs( 0, 0, (short)           0   );
	crap.dataAs( 0, 0, (unsigned int)    0   );
	crap.dataAs( 0, 0, (int)             0   );
	crap.dataAs( 0, 0, (unsigned long long)   0   );
	crap.dataAs( 0, 0, (long long)            0   );
	crap.dataAs( 0, 0, (float)           0.0 );
	crap.dataAs( 0, 0, (double)          0.0 );
	crap.dataAs( 0, 0, (long double)     0.0 );

	crap.dataAs( 0, (unsigned char)   0   );
	crap.dataAs( 0, (char)            0   );
	crap.dataAs( 0, (unsigned short)  0   );
	crap.dataAs( 0, (short)           0   );
	crap.dataAs( 0, (unsigned int)    0   );
	crap.dataAs( 0, (int)             0   );
	crap.dataAs( 0, (unsigned long long)   0   );
	crap.dataAs( 0, (long long)            0   );
	crap.dataAs( 0, (float)           0.0 );
	crap.dataAs( 0, (double)          0.0 );
	crap.dataAs( 0, (long double)     0.0 );
}

// static methods to create nodes and (pure) comments
QUB_Tree QUB_Tree::Create( string name ) {
	QUB_Tree node = QTR_Create( name.c_str() );
	QTR_DECREF( node.getImpl() );
	return node;
}

QUB_Tree QUB_Tree::CreateComment( string text ) {
  QUB_Tree com = QUB_Tree::Create( "QTR_COMMENT" );
  com.setData( text );
  return com;
}

QUB_Tree QUB_Tree::Open( string path, bool readOnly ) {
	QTR_Impl *impl = QTR_Open( path.c_str(), readOnly );
	QUB_Tree node( impl );
	if ( impl )
		QTR_DECREF( node.getImpl() );
	return node;
}

QUB_Tree QUB_Tree::ReadText( string path )
{
	std::ifstream in;
	in.open( path.c_str() );
	return QUB_Tree::CreateFromStream( in );
}

QUB_Tree QUB_Tree::CreateFromString( char *s )
{
	std::istrstream in( s );
	return QUB_Tree::CreateFromStream( in );
}

QUB_Tree QUB_Tree::CreateFromBuffer( const char *buffer )
{
	QTR_Impl *impl = QTR_FromBuffer( buffer );
	QUB_Tree node( impl );
	if ( impl )
		QTR_DECREF( node.getImpl() );
	return node;
}

// QUB_Tree_Text: QUB_Tree   QUB_Tree::CreateFromStream( std::istream& in )

//------------------- construct, copy, compare
QUB_Tree::QUB_Tree() // null
			: impl( 0 ) {
	}

QUB_Tree::QUB_Tree( QTR_Impl *n ) // refer to n
			: impl( n ) { 
	if ( impl )
		QTR_INCREF( impl ); 
	}

QUB_Tree::QUB_Tree( const QUB_Tree &n ) // refer to same as n
			: impl(n.getImpl()) { 
	if ( impl ) 
		QTR_INCREF( impl ); 
	}

QUB_Tree::~QUB_Tree() {
	if ( impl )
		QTR_DECREF( impl );
}


QUB_Tree& QUB_Tree::operator= (const QUB_Tree &n) { // refer to same as n
	if ( impl )                              
		QTR_DECREF( impl );
	impl = n.getImpl();         
	if ( impl )  
		QTR_INCREF( impl );
	return *this;
	}

bool QUB_Tree::equals (const QUB_Tree& n) const{
  return ( impl == n.impl );
}

bool QUB_Tree::save(){
	if ( impl )
		return (0 <= QTR_Save(impl));
	return false;
}

bool QUB_Tree::saveAs( string path ){
	if ( impl )
		if ( 0 <= QTR_SaveAs( impl, path.c_str() ) )
			return save();
	return false;
}

bool QUB_Tree::saveAsTemp(){
	if ( impl )
		if ( 0 <= QTR_SaveAsTemp( impl ) )
			return save();
	return false;
}

bool QUB_Tree::close(){
	if ( impl )
		QTR_Close( impl );
	return true;
}

bool QUB_Tree::inFile(){
	return (NULL!=impl) && (NULL!=impl->root->file);
	/*
	if ( impl )
		return impl->root->file;
	return false;
	*/
	}

string QUB_Tree::path() const
{
	string path;
	if ( impl && impl->root->file )
		path = impl->root->file->path;
	return path;
}

bool QUB_Tree::saveTextCopy( string path ) const
{
	if ( impl ) {
		std::ofstream out;
		out.open( path.c_str() );
		if ( out ) {
			toStream( out, "" );
			return true;
		}
	}
	return false;
}

string QUB_Tree::toString() const {
	string s;
	if ( impl ) {
		std::ostrstream out;
		toStream( out, "" );
		out << '\0';
		s = out.str();
		out.rdbuf()->freeze(0);
	}
	return s;
}

unsigned int QUB_Tree::toBufferSize() const {
	return QTR_ToBufferSize( impl );
}

unsigned int QUB_Tree::toBuffer(char *buffer) const {
	return QTR_ToBuffer( impl, buffer );
}

QUB_Tree QUB_Tree::clone( bool deep ) const {
	QUB_Tree cl;
	if ( impl ) {
		cl = QUB_Tree( QTR_Clone( impl, deep ) );
		QTR_DECREF( cl.getImpl() );
	}
  return cl;
}
  
QUB_Tree QUB_Tree::appendChild( string name ){
	QUB_Tree child;
	if ( impl ) {
		QUB_TreeIter lastChild = end();
		--lastChild;

		child = QTR_CreateChild( impl, lastChild->getImpl(), name.c_str() );
	}
	return child;
}

QUB_Tree QUB_Tree::appendChild( QUB_Tree newchild ){
	if ( impl )
		end().insert( newchild );
	return newchild;
}

QUB_Tree QUB_Tree::appendClone( QUB_Tree orig, bool deep ){
	QUB_Tree clone;
	if ( impl ) {
		QUB_TreeIter lastChild = end();
		--lastChild;

		clone = QTR_InsertClone( impl, lastChild->getImpl(), orig.getImpl(), deep );
	}
	return clone;
}

QUB_Tree QUB_Tree::insertChild( string name ){
	QUB_Tree child;
	if ( impl )
		child = QTR_CreateChild( impl, 0, name.c_str() );
	return child;
}

QUB_Tree QUB_Tree::insertChild( QUB_Tree newchild ){
	if ( impl )
		QTR_InsertChild( impl, 0, newchild.getImpl() );
	return newchild;
}

QUB_Tree QUB_Tree::insertClone( QUB_Tree orig, bool deep ){
	QUB_Tree clone;
	if ( impl )
		clone = QTR_InsertClone( impl, 0, orig.getImpl(), deep );
	return clone;
}

QUB_Tree QUB_Tree::insertChild( QUB_Tree after, string name ){
	QUB_Tree child;
	if ( impl )
		child = QTR_CreateChild( impl, after.getImpl(), name.c_str() );
	return child;
}

QUB_Tree QUB_Tree::insertChild( QUB_Tree after, QUB_Tree newchild ){
	if ( impl )
		QTR_InsertChild( impl, after.getImpl(), newchild.getImpl() );
	return newchild;
}

QUB_Tree QUB_Tree::insertClone( QUB_Tree after, QUB_Tree orig, bool deep )
{
	QUB_Tree clone;
	if ( impl )
		clone = QTR_InsertClone( impl, after.getImpl(), orig.getImpl(), deep );
	return clone;
}

void QUB_Tree::removeChild( QUB_Tree remchild )
{
	if ( impl && remchild.getImpl() ) {
		QUB_Tree prev = child();
		if ( remchild == prev ) // remchild is first child
			QTR_RemoveChild( impl, 0, remchild.getImpl() );
		else {
			while ( ! ( prev.sibling() == remchild ) ) {
				prev = prev.sibling();
			}
			QTR_RemoveChild( impl, prev.getImpl(), remchild.getImpl() );
		}
	}
}

// get iterators

QUB_Tree::iterator QUB_Tree::children() const {         return QUB_TreeIter( *this ); }
QUB_Tree::iterator QUB_Tree::end()      const {         return QUB_TreeIter( *this, QTR_ITER_END ); }

QUB_Tree::iterator QUB_Tree::find( string name ) const {
	iterator ch = children();
	if ( (*ch).name() != name )
		ch.next( name );
	return ch;
}

QUB_Tree::iterator QUB_Tree::rfind( string name ) const {
	iterator ch = end();
	ch.prev( name );
	return ch;
}

// get neighbors in tree

QUB_Tree QUB_Tree::parent() const {
  if ( impl && impl->parent )
    return QUB_Tree( impl->parent );
  return QUB_Tree();
}

QUB_Tree QUB_Tree::child() const {
  if ( impl && impl->child )
    return QUB_Tree( impl->child );
  return QUB_Tree();
}

QUB_Tree QUB_Tree::child( string name ) const {
	return * (find(name));
}

QUB_Tree QUB_Tree::sibling() const {
  if ( impl && impl->sibling )
    return QUB_Tree( impl->sibling );
  return QUB_Tree();
}

QUB_Tree QUB_Tree::sibling( string name ) const {
	QUB_Tree sib = sibling();
	while ( ! sib.isNull() )
		if ( sib.name() == name )
			return sib;
		else
			sib = sib.sibling();

	return QUB_Tree();
}

QUB_Tree QUB_Tree::siblingSameName() const {
	QUB_Tree sib = sibling();
	while ( ! sib.isNull() )
		if ( sib.impl->name == impl->name )
			return sib;
		else
			sib = sib.sibling();

	return QUB_Tree();
}

QUB_Tree QUB_Tree::operator[] (string childName) const
{
	iterator it = find( childName );
	if ( (*it).isNull() )
		it.insert( childName );
  
	return *it;
}

// get/set flags

unsigned int & QUB_Tree::flags() const {
  if ( ! impl ) return * ( (unsigned int *) impl );
  return * ( (unsigned int *) (& (impl->ondisk.flags)) );
}

bool QUB_Tree::isPreload() const {
	if ( impl )
		return (0!=QTR_Flag( impl, QTR_FLAG_PRELOAD ));
	return 0;
	}
 
void QUB_Tree::setPreload( bool pl ) {
  if ( impl )
    QTR_SetFlag( impl, QTR_FLAG_PRELOAD, pl );
}

void QUB_Tree::setChanged( bool ch ) {
	if ( impl )
		if ( ch )
			QTR_Changed( impl );
		else
		    QTR_SetFlag( impl, QTR_FLAG_CHANGED, 0 );
}

int QUB_Tree::lock( int timeoutMS )
{
  if ( impl )
    return QTR_Lock( impl, timeoutMS );
  return 0;
}

void QUB_Tree::unlock()
{
  if ( impl )
    QTR_Unlock( impl );
}

// get/set name

string QUB_Tree::name() const {
  string nm;
  if ( impl )
	  nm = QTR_Name( impl );
  return nm;
}

void QUB_Tree::setName( string name ) {
  if ( impl )
    QTR_SetName( impl, name.c_str() );
}

// raw data

bool QUB_Tree::dataIs( QTR_DataType type ) const {
  if ( (impl == 0) && (type == QTR_TYPE_EMPTY) )
    return true;
  
  return ( dataType() == type );
}

void *QUB_Tree::data( bool autoload ) {
  if ( impl ) {
    if ( autoload && (loadedCount() <= 0) )
      loadData();
    return impl->data;
  }
  return 0;
}

bool QUB_Tree::setData( string data )
{
	return setData( QTR_TYPE_STRING, 1, data.length(), (char*) data.c_str() );
}

bool QUB_Tree::setNumData( QTR_DataType type, unsigned int rows, unsigned int cols, void *inits )
{
	if ( impl )
		if (0 <= QTR_SetupNumData( impl, type, rows, cols )) {
			if ( inits )
				QTR_SetRows( impl, inits, 0, rows-1 );
			return true;
		}
	return false;
}

bool QUB_Tree::setNumData( QTR_DataType type, unsigned int rows, unsigned int cols, void **inits )
{
	if ( impl )
		if (0 <= QTR_SetupNumData( impl, type, rows, cols )) {
			if ( inits ) {
				for ( unsigned int r = 0; r < rows; ++r )
					QTR_SetRows( impl, inits[r], r, r );
			}
			return true;
		}
	return false;
}

bool QUB_Tree::setData( QTR_DataType type, unsigned int size, unsigned int count, void *inits )
{
	if ( impl )
		if (0 <= QTR_SetupData( impl, type, size, count )) {
			if ( inits )
				QTR_SetRows( impl, inits, 0, dataRows() - 1 );
			return true;
		}
	return false;
}

bool QUB_Tree::resizeData( unsigned int newrowcount ) {
  if ( impl )
    return (0 <= QTR_ResizeData( impl, newrowcount ));
  return false;
}

bool QUB_Tree::clearData() {
  if ( impl )
    return (0 <= QTR_ClearData( impl ));
  return false;
}

bool QUB_Tree::loadData( unsigned int firstRow, unsigned int lastRow, bool doRead ) {
  if ( impl ) {
    if ( lastRow == QTR_ALL )
      lastRow = dataRows() - 1;
    return (0 <= QTR_LoadRows( impl, firstRow, lastRow, doRead ));
  }
  return false;
}

bool QUB_Tree::unloadData( bool doWrite ) {
  if ( impl )
    return (0 <= QTR_UnloadRows( impl, doWrite ));
  return false;
}

unsigned int QUB_Tree::readDataInto( void *buf, unsigned int firstRow, unsigned int lastRow ) const
{
  if ( impl )
    return QTR_GetRows( impl, buf, firstRow, lastRow );
  return false;
}

unsigned int QUB_Tree::readDataInto( void **buf, unsigned int firstRow, unsigned int lastRow ) const
{
	unsigned int nread = 0;
	if ( impl ) {
		for ( unsigned int r = firstRow; r <= lastRow; ++r )
			nread += QTR_GetRows( impl, buf[r-firstRow], r, r );
	}
	return nread;
}

unsigned int QUB_Tree::writeDataFrom( void *buf, unsigned int firstRow, unsigned int lastRow )
{
  if ( impl )
    return QTR_SetRows( impl, buf, firstRow, lastRow );
  return false;
}

unsigned int QUB_Tree::writeDataFrom( void **buf, unsigned int firstRow, unsigned int lastRow )
{
	unsigned int nwrote = 0;
	if ( impl ) {
		for ( unsigned int r = firstRow; r <= lastRow; ++r )
			nwrote += QTR_SetRows( impl, buf[r-firstRow], r, r );
	}
	return nwrote;
}


// copy arrays of different types by casting elements.
// the arrays can be of different lengths
// if toLen > fromLen:
//   if there's a defaults array of toLen:
//     the remainder are filled from defs
//   else:
//     the remainder are uninitialized

// defLen == toLen unless defs == 0
template<class To, class From>
void CopyCastFill( To *into, To *defs, From *from,
		   unsigned int toLen, unsigned int fromLen )
{
  for ( unsigned int i=0; i<toLen; ++i ) {
    if ( i < fromLen )
      into[i] = (To) from[i];
    else if ( defs )
      into[i] = defs[i];
    else // no more elements in from, and no defaults
      break;
  }
}
/*
// copy a node's data with the appropriate casts
// intolen can be different from node's datalen

template<class T>
bool CopyCastFillFrom( T *into, unsigned int intolen, T *defs, QUB_Tree from,
		       int timeoutMS )
{
  unsigned int datalen = from.dataCount();
  
  if ( from.isNull() ) {
    // go away
  }
  else {
    char *fromdat = new char[ from.dataSize() * datalen ];
    bool success = false;
    if ( from.readDataInto( fromdat, 0, datalen - 1, timeoutMS ) ) {
      switch ( from.dataType() ) {
      case QTR_TYPE_EMPTY:
	// simple set?
	break;
      case QTR_TYPE_UNKNOWN:
	// treat like T* ?
	break;
      case QTR_TYPE_POINTER:
	CopyCastFill( into, defs, (unsigned int*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_STRING:
	CopyCastFill( into, defs, (char*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_UCHAR:
	CopyCastFill( into, defs, (unsigned char*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_CHAR:
	CopyCastFill( into, defs, (char*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_USHORT:
	CopyCastFill( into, defs, (unsigned short*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_SHORT:
	CopyCastFill( into, defs, (short*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_UINT:
	CopyCastFill( into, defs, (unsigned int*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_INT:
	CopyCastFill( into, defs, (int*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_ULONG:
	CopyCastFill( into, defs, (unsigned long long*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_LONG:
	CopyCastFill( into, defs, (long long*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_FLOAT:
	CopyCastFill( into, defs, (float*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_DOUBLE:
	CopyCastFill( into, defs, (double*)fromdat, intolen, datalen );
	break;
      case QTR_TYPE_LDOUBLE:
	CopyCastFill( into, defs, (long double*)fromdat, intolen, datalen );
	break;
      }
      success = true;
    }
    delete [] fromdat;
    return success;
  }
  return false;
}
*/
template<class T>
bool CopyCastFillRowsFrom( T *into, QUB_Tree from, unsigned int firstRow, unsigned int lastRow )
{
  unsigned int rowlen = from.dataCols();
  
  if ( from.isNull() ) {
    return false; // go away
  }
  else {
	  char *fromdat = new char[ from.dataSize() ];
	  for ( unsigned int r = firstRow; r <= lastRow; ++r, into += rowlen ) {
		  from.readDataInto( fromdat, r, r );
		  switch ( from.dataType() ) {
			  case QTR_TYPE_EMPTY:
				// simple set?
				break;
			  case QTR_TYPE_UNKNOWN:
				// treat like T* ?
				break;
			  case QTR_TYPE_POINTER:
				CopyCastFill( into, (T*) 0, (unsigned int*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_STRING:
				CopyCastFill( into, (T*) 0, (char*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_UCHAR:
				CopyCastFill( into, (T*) 0, (unsigned char*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_CHAR:
				CopyCastFill( into, (T*) 0, (char*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_USHORT:
				CopyCastFill( into, (T*) 0, (unsigned short*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_SHORT:
				CopyCastFill( into, (T*) 0, (short*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_UINT:
				CopyCastFill( into, (T*) 0, (unsigned int*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_INT:
				CopyCastFill( into, (T*) 0, (int*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_ULONG:
				CopyCastFill( into, (T*) 0, (unsigned long long*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_LONG:
				CopyCastFill( into, (T*) 0, (long long*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_FLOAT:
				CopyCastFill( into, (T*) 0, (float*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_DOUBLE:
				CopyCastFill( into, (T*) 0, (double*)fromdat, rowlen, rowlen );
				break;
			  case QTR_TYPE_LDOUBLE:
				CopyCastFill( into, (T*) 0, (long double*)fromdat, rowlen, rowlen );
				break;
		  }
	  }
	  delete [] fromdat;
  }
  return true;
}
/*
  // if is T, resize, apply defaults to growth zone.
  // else alloc T*, copy existing or defaults, setData
bool QUB_Tree::enforceInt( unsigned int count, int *defaults, int timeoutMS )
{
  if ( ! impl )
    return false;
  
  if ( dataType() == QTR_TYPE_INT ) {
    int oldlen = dataCount();
    if ( oldlen != count ) {
      if ( beginWrite( timeoutMS ) ) {
	resizeDataInWrite( count );
	loadDataInWrite();
	if ( defaults ) {
	  for ( unsigned int i=oldlen; i<count; i++ )
	    ( (int*) data() )[i] = defaults[i];
	}
	endWrite();
	return true;
      }
	  return false;
    }
	return true;
  } 
  else {
    int *tmpArr = new int[ count ];
    bool rslt = false;
    if ( CopyCastFillFrom( tmpArr, count, defaults, *this, timeoutMS ) )
      rslt = setData( QTR_TYPE_INT, count, tmpArr, timeoutMS );
    delete [] tmpArr;
    return rslt;
  }
  
  // return false;
}    

bool QUB_Tree::enforceFloat( unsigned int count, float *defaults, int timeoutMS )
{
  if ( ! impl )
    return false;
  
  if ( dataType() == QTR_TYPE_FLOAT ) {
    int oldlen = dataCount();
    if ( oldlen != count ) {
      if ( beginWrite( timeoutMS ) ) {
		resizeDataInWrite( count );
		loadDataInWrite();
		if ( defaults ) {
			for ( unsigned int i=oldlen; i<count; i++ )
				( (float*) data() )[i] = defaults[i];
		}
		endWrite();
		return true;
	  }
	  else {
		  return false;
	  }
    }
	else {
		return true;
	}
  } 
  else {
    float *tmpArr = new float[ count ];
    bool rslt = false;
    if ( CopyCastFillFrom( tmpArr, count, defaults, *this, timeoutMS ) ) {
      rslt = setData( QTR_TYPE_FLOAT, count, tmpArr, timeoutMS );
    }
    delete [] tmpArr;
    return rslt;
  }
  
  // return false;
}    

bool QUB_Tree::enforceDouble( unsigned int count, double *defaults, int timeoutMS )
{
  if ( ! impl )
    return false;
  
  if ( dataType() == QTR_TYPE_DOUBLE ) {
    int oldlen = dataCount();
    if ( oldlen != count ) {
      if ( beginWrite( timeoutMS ) ) {
	resizeDataInWrite( count );
	loadDataInWrite();
	if ( defaults ) {
	  for ( unsigned int i=oldlen; i<count; i++ )
	    ( (double*) data() )[i] = defaults[i];
	}
	endWrite();
	return true;
      }
	  return false;
    }
	return true;
  } 
  else {
    double *tmpArr = new double[ count ];
    bool rslt = false;
    if ( CopyCastFillFrom( tmpArr, count, defaults, *this, timeoutMS ) )
      rslt = setData( QTR_TYPE_DOUBLE, count, tmpArr, timeoutMS );
    delete [] tmpArr;
    return rslt;
  }
  
  // return false;
}

/*
bool QUB_Tree::enforceDouble( unsigned int count, double *defaults, int timeoutMS )
{
  if ( ! impl )
    return false;
  
  if ( ! beginWrite( timeoutMS ) )
    return false;
  
  if ( ! (loadDataInWrite( 0, QTR_ALL, true )
	  || (loadStart() == 0 && loadEnd() == (dataCount() - 1)) ) ) {
    endWrite();
    return false;
  }
  
  if ( dataIs( QTR_TYPE_DOUBLE ) ) {
    int oldlen = dataCount();
    if ( oldlen != count ) {
      resizeDataInWrite( count );
      
      if ( defaults ) {
	for ( unsigned int i=oldlen; i<count; i++ )
	  ( (double*) data() )[i] = defaults[i];
      }
    }
  }
  else {  
    double *tmpArr = new double[ count ];
    CopyCastFillFrom( tmpArr, count, defaults, *this );
    setDataInWrite( QTR_TYPE_DOUBLE, count, tmpArr );
    delete [] tmpArr;
  }
  
  endWrite();
  return true;
}
*/

string QUB_Tree::dataAsString( bool convert ) const
{
	string datstr;
	if ( ! dataCount() )
		return datstr;

	if ( convert ) {
		std::ostrstream out;
		dataToStream( out, "" );
		out << '\0';
		datstr = ( out.str() );
		out.rdbuf()->freeze(0);
		return datstr;
	}
	datstr.resize( dataSize() * dataCount() );
	readDataInto( &(datstr[0]), 0, dataRows() - 1 );
	return datstr;
}

void QUB_Tree::dataAsStrings( vector<string>& lines, bool convert ) const
{
  string data = dataAsString( convert );
  int start, next, nextCR, len;
  
  start = 0;
  next = 0;
  while ( start < int(data.length()) ) {
    next = data.find( "\n", start );
    if ( next == string::npos )
      next = data.length();
    
    nextCR = next;
    if ( data[ next - 1 ] == '\r' )
      nextCR--;
      
    len = nextCR - start;
    lines.push_back( string( data.substr(start, len) ) );
    
    start = next + 1;
  }
}

int QUB_Tree::dataAsInt( int def, int r, int c )
{
	int val = def;
	if ( (r < int(dataRows())) && (c < int(dataCols())) ) {
		if ( dataType() == QTR_TYPE_INT ) {
			if ( dataCols() == 1 )
				readDataInto( &val, r, r );
			else {
				int *row = new int[ dataCols() ];
				readDataInto( row, r, r );
				val = row[c];
				delete [] row;
			}
		}
		else {
			val = dataAs( r, c, val );
		}
	}
	return val;
}

double QUB_Tree::dataAsDouble( double def, int r, int c )
{
	double val = def;
	if ( (r < int(dataRows())) && (c < int(dataCols())) ) {
		if ( dataType() == QTR_TYPE_DOUBLE ) {
			if ( dataCols() == 1 )
				readDataInto( &val, r, r );
			else {
				double *row = new double[ dataCols() ];
				readDataInto( row, r, r );
				val = row[c];
				delete [] row;
			}
		}
		else {
			val = dataAs( r, c, val );
		}
	}
	return val;
}
	/*
  double val = def;
  if ( beginRead( timeoutMS ) ) {
    CountedArrayPtr<double> arr = dataAsDoubles( timeoutMS );
    if ( arr.size() )
      val = arr[0];
    endRead();
  }
  return val;
}*/

// fill your array, converting if needed/possible
unsigned int QUB_Tree::getDataAsInts( int *buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_INT )
		return readDataInto( buf, firstRow, lastRow );
	else if ( CopyCastFillRowsFrom( buf, *this, firstRow, lastRow ) )
		return (lastRow - firstRow + 1);
	else
		return 0;
}

unsigned int QUB_Tree::getDataAsInts( int **buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_INT )
		return readDataInto( buf, firstRow, lastRow );
	
	unsigned int nread = 0;
	for ( unsigned int r = firstRow; r <= lastRow; ++r )
		if ( CopyCastFillRowsFrom( buf[r - firstRow], *this, r, r ) )
			++nread;
	return nread;
}

unsigned int QUB_Tree::getDataAsFloats( float *buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_FLOAT )
		return readDataInto( buf, firstRow, lastRow );
	else if ( CopyCastFillRowsFrom( buf, *this, firstRow, lastRow ) )
		return (lastRow - firstRow + 1);
	else
		return 0;
}

unsigned int QUB_Tree::getDataAsFloats( float **buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_FLOAT )
		return readDataInto( buf, firstRow, lastRow );
	
	unsigned int nread = 0;
	for ( unsigned int r = firstRow; r <= lastRow; ++r )
		if ( CopyCastFillRowsFrom( buf[r - firstRow], *this, r, r ) )
			++nread;
	return nread;
}

unsigned int QUB_Tree::getDataAsDoubles( double *buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_DOUBLE )
		return readDataInto( buf, firstRow, lastRow );
	else if ( CopyCastFillRowsFrom( buf, *this, firstRow, lastRow ) )
		return (lastRow - firstRow + 1);
	else
		return 0;
}

unsigned int QUB_Tree::getDataAsDoubles( double **buf, unsigned int firstRow, unsigned int lastRow )
{
	if ( dataType() == QTR_TYPE_DOUBLE )
		return readDataInto( buf, firstRow, lastRow );
	
	unsigned int nread = 0;
	for ( unsigned int r = firstRow; r <= lastRow; ++r )
		if ( CopyCastFillRowsFrom( buf[r - firstRow], *this, r, r ) )
			++nread;
	return nread;
}



// ya know what?  forget the optimal (direct pointer if correct type)
// just give 'em a copy.
/* god this stuff is just so wrong for a dll
CountedArrayPtr<int> QUB_Tree::dataAsInts( int timeoutMS )
{
  unsigned int len = dataCount();
  CountedArrayPtr<int> intData( new int[len], len );
  
  if ( dataType() == QTR_TYPE_INT )
    readDataInto( intData.ptr(), 0, len - 1, timeoutMS );
  else
    CopyCastFillFrom( intData.ptr(), len, (int*) NULL, *this, timeoutMS );
  
  return intData;
}

CountedArrayPtr<double> QUB_Tree::dataAsDoubles( int timeoutMS )
{
  unsigned int len = dataCount();
  CountedArrayPtr<double> dblData( new double[len], len );
  
  if ( dataType() == QTR_TYPE_DOUBLE )
    readDataInto( dblData.ptr(), 0, len - 1, timeoutMS );
  else
    CopyCastFillFrom( dblData.ptr(), len, (double*) NULL, *this, timeoutMS );
  
  return dblData;
}

  int QUB_Tree::getNodeDataSize( QTR_DataType dtype )
{
	return QTR_DataSize[dtype];
}

CountedArrayPtr<double> QUB_Tree::dataAsDoubles()
{
  CountedArrayPtr<double> doubleData;
  unsigned int len = dataCount();
  loadData();
  
  if ( len == 0 ) {
    // leave doubleData empty/null
  }
  else if ( dataIs( QTR_TYPE_DOUBLE ) ) {
    doubleData = CountedArrayPtr<double>( (double*)data(), len, false/ *nodelete* / );
  }
  else {
    doubleData = CountedArrayPtr<double>( new double[len], len, true );
    startRead();
    CopyCastFillFrom( doubleData.ptr(), len, (double*)0, *this );
    endRead();
  }
  
  return doubleData;
}
*/


// same-line comment management

bool QUB_Tree::hasLineComment() const
{
  QUB_Tree com = *find("QTR_LINE_COMMENT");
  if ( ! com.isNull() )
    return true;
  return false;
}

string QUB_Tree::getLineComment() const
{
  QUB_Tree com = *find("QTR_LINE_COMMENT");
  if ( ! com.isNull() )
    return com.dataAsString();
  return "";
}

void QUB_Tree::setLineComment( string linecom )
{
	if ( 0 == linecom.length() )
		find("QTR_LINE_COMMENT").remove();
	else
		(*this)["QTR_LINE_COMMENT"].setData( linecom );
}

QTR_API bool operator== (const QUB_Tree& lhs, const QUB_Tree& rhs)
{
	return lhs.equals( rhs );
}

// ************************** QUB_TreeIter *******************************

QUB_TreeIter::QUB_TreeIter() // all null
  : parent( QUB_Tree() ), node( QUB_Tree() ), ix( 0 ), prevNodes( new vector<QUB_Tree> )
{}

QUB_TreeIter::QUB_TreeIter( const QUB_Tree &p )
  : parent( p ), node( p.child() ), ix( 0 ), prevNodes( new vector<QUB_Tree> )
{}

QUB_TreeIter::QUB_TreeIter( const QUB_Tree &p, int what )
  : parent( p ), prevNodes( new vector<QUB_Tree> )
{
  if ( what == QTR_ITER_END ) {
    ix = -1; // don't figure it out yet
  }
  else if ( ! p.isNull() ) {
    node = QUB_Tree( p.child() );
    ix = 0;
  }
}

QUB_TreeIter::QUB_TreeIter( const QUB_TreeIter& copy )
  : parent( copy.parent ), node( copy.node ), ix( copy.ix ),
    prevNodes( new vector<QUB_Tree> )
{
	* ((vector<QUB_Tree> *)prevNodes) = * ((vector<QUB_Tree> *)(copy.prevNodes));
}

QUB_TreeIter::~QUB_TreeIter()
{
	delete ((vector<QUB_Tree> *)prevNodes);
}

QUB_TreeIter& QUB_TreeIter::operator= (const QUB_TreeIter& copy) {
  parent = copy.parent;
  node = copy.node;
  ix = copy.ix;
  * ((vector<QUB_Tree> *)prevNodes) = * ((vector<QUB_Tree> *)(copy.prevNodes));
  return *this;
}

bool QUB_TreeIter::equals (const QUB_TreeIter &it) const {
  return ( parent == it.parent
	   && (    (node.isNull() && it.node.isNull()) // end/empty
		|| (ix == it.ix)) );
}

QUB_Tree QUB_TreeIter::operator* () {
  return node;
}

QUB_Tree* QUB_TreeIter::operator-> () {
  return &node;
}

QUB_TreeIter& QUB_TreeIter::operator++ () {
  if ( ! node.isNull() ) {
    ((vector<QUB_Tree> *)prevNodes)->push_back( node );
    node = node.sibling();
    ++ix;
  }
  return *this;
}

QUB_TreeIter& QUB_TreeIter::operator-- () {
  if ( ix == -1 ) {
    // constructed at end: go through whole list
    node = parent.child();
    ix = 0;
    while ( ! node.isNull() )
      ++(*this);
  }
  else if ( ix > 0 && ((vector<QUB_Tree> *)prevNodes)->empty() ) {
    // constructed from a mono and already -- as bi: refind
    int targix = ix;
    ix = 0;
    node = parent.child();
    while ( ix < targix )
      ++(*this);
  }
  
  if ( ! ((vector<QUB_Tree> *)prevNodes)->empty() ) {
    node = ((vector<QUB_Tree> *)prevNodes)->back();
    ((vector<QUB_Tree> *)prevNodes)->pop_back();
    --ix;
  }
  return *this;
}

QUB_TreeIter& QUB_TreeIter::next( string name ) {
  char *iName = QTR_LookupName( name.c_str() );

  while ( ! node.isNull() ) {
    ++(*this);
	QTR_Impl *chImpl = node.getImpl();
    if ( chImpl && (chImpl->name == iName) )
      break;
  }
  return *this;
}

QUB_TreeIter& QUB_TreeIter::nextSameName() {
  char *name = node.isNull() ? NULL : node.getImpl()->name;

  while ( ! node.isNull() ) {
    ++(*this);
	QTR_Impl *chImpl = node.getImpl();
    if ( chImpl && (chImpl->name == name) )
      break;
  }
  return *this;
}

QUB_TreeIter& QUB_TreeIter::prev( string name ) {
  char *iName = QTR_LookupName( name.c_str()  );

  while ( ix != 0 ) {
    --(*this);
	QTR_Impl *chImpl = node.getImpl();
    if ( chImpl && (chImpl->name == iName) )
      return *this; // skip GO_TO_END
  }
  node = QUB_Tree(0);
  ((vector<QUB_Tree> *)prevNodes)->clear();
  ix = -1;
  return *this;
}

QUB_TreeIter& QUB_TreeIter::prevSameName() {
  char *name = node.isNull() ? NULL : node.getImpl()->name;

  while ( ix != 0 ) {
    --(*this);
	QTR_Impl *chImpl = node.getImpl();
    if ( chImpl && (chImpl->name == name) )
      return *this; // skip GO_TO_END
  }
  node = QUB_Tree(0);
  ((vector<QUB_Tree> *)prevNodes)->clear();
  ix = -1;
  return *this;
}

int QUB_TreeIter::index() {
  if ( ix == -1 ) {
    --(*this);
    ++(*this);
  }
  return ix;
}

QUB_Tree QUB_TreeIter::getParent() {
  return parent;
}

void QUB_TreeIter::insert( QUB_Tree n ) {
	if ( ix == -1 ) {
		--(*this);
		++(*this);
	}
	
  if ( ! parent.isNull() ) {
	   QTR_InsertChild( parent.getImpl(), ix ? prevNode().getImpl() : 0, n.getImpl() );
       node = n;
  }
}

void QUB_TreeIter::insert( string name ) {
	if ( ix == -1 ) {
		--(*this);
		++(*this);
	}
	if ( ! parent.isNull() ) {
		node = QTR_CreateChild( parent.getImpl(), ix ? prevNode().getImpl() : 0, name.c_str() );
	}
}

QUB_Tree QUB_TreeIter::remove() {
  QUB_Tree removed = node;
  
	if ( ix == -1 ) {
		--(*this);
		++(*this);
	}

  if ( ! node.isNull() ) {
    if ( ix == 0 ) {
      QTR_RemoveChild( parent.getImpl(), 0, node.getImpl() );
      node = parent.child();
    }
    else {
      QUB_Tree prv( prevNode() );
      QTR_RemoveChild( parent.getImpl(), prv.getImpl(), node.getImpl() );
      node = prv.sibling();
    }
  }
  
  return removed;
}

QUB_Tree QUB_TreeIter::prevNode() {
  if ( ix == -1 ) {
    --(*this);
    ++(*this);
  }
  else if ( ix > 0 && ((vector<QUB_Tree> *)prevNodes)->empty() ) {
    // constructed from a mono then -- as bi: refind
    int targix = ix;
    ix = 0;
    node = parent.child();
    while ( ix < targix )
      ++(*this);
  }
  
  if ( ! ((vector<QUB_Tree> *)prevNodes)->empty() )
    return ((vector<QUB_Tree> *)prevNodes)->back();
  else
    return QUB_Tree(); // null
}

QTR_API bool operator== (const QUB_TreeIter& lhs, const QUB_TreeIter& rhs)
{
	return lhs.equals( rhs );
}

QTR_API bool operator!= (const QUB_TreeIter& lhs, const QUB_TreeIter& rhs)
{
	return ! ( lhs == rhs );
}



QUB_TreeMonoIter::QUB_TreeMonoIter() // all null
{}

QUB_TreeMonoIter::QUB_TreeMonoIter( const QUB_Tree &p )
  : QUB_TreeIter(p)
{}

QUB_TreeMonoIter::QUB_TreeMonoIter( const QUB_Tree &p, int what )
  : QUB_TreeIter(p, what)
{}

QUB_TreeMonoIter::QUB_TreeMonoIter( const QUB_TreeIter& copy )
  : QUB_TreeIter( copy )
{
	if ( ((vector<QUB_Tree> *)prevNodes)->size() ) {
		QUB_Tree prv = ((vector<QUB_Tree> *)prevNodes)->back();
		((vector<QUB_Tree> *)prevNodes)->clear();
		((vector<QUB_Tree> *)prevNodes)->push_back( prv );
	}
}

QUB_TreeMonoIter::~QUB_TreeMonoIter()
{}

QUB_TreeMonoIter& QUB_TreeMonoIter::operator= (const QUB_TreeIter& copy) {
  parent = copy.parent;
  node = copy.node;
  ix = copy.ix;
  if ( ((vector<QUB_Tree> *)(copy.prevNodes))->size() )
	((vector<QUB_Tree> *)prevNodes)->push_back( ((vector<QUB_Tree> *)(copy.prevNodes))->back() );
  return *this;
}

QUB_TreeIter& QUB_TreeMonoIter::operator++ () {
  if ( ! node.isNull() ) {
    ((vector<QUB_Tree> *)prevNodes)->clear();
	((vector<QUB_Tree> *)prevNodes)->push_back( node );
    node = node.sibling();
    ++ix;
  }
  return *this;
}
  
QUB_TreeIter& QUB_TreeMonoIter::operator-- () {
  throw "-- is for bidirectional iterators only";
  // return *this;
}

QUB_TreeIter& QUB_TreeMonoIter::prev( string name ) {
  throw "prev is for bidirectional iterators only";
  // return *this;
}

QUB_TreeIter& QUB_TreeMonoIter::prevSameName() {
  throw "prev is for bidirectional iterators only";
  // return *this;
}

int QUB_TreeMonoIter::index() {
  return ix;
}

QUB_Tree QUB_TreeMonoIter::prevNode() {
  if ( ! ((vector<QUB_Tree> *)prevNodes)->empty() )
    return ((vector<QUB_Tree> *)prevNodes)->back();
  else
    return QUB_Tree(); // null
}

