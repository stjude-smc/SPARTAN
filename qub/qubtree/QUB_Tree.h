#ifndef QUB_TREE_H
#define QUB_TREE_H

// classes which wrap the C API in QUB_QTR_Impl.h
//   they handle the ref.counting for you,
//   and add convenience methods.

// along the lines of NodeClasses.h but with
//   simpler access control and data access
//   missing enforce*
//   broken text features

#include "QUB_QTR_Impl.h"
#include <stdio.h>

#include <string>
#include <vector>
using namespace std;

#include "istream_readline.h"

#define QTR_ITER_END -1

#define QTR_ALL 0xffffffff

class QTR_API QUB_TreeIter;

class QTR_API QUB_Tree
{
 public:
	// never say "new QUB_Tree"; use QUB_Tree::Create instead
	static QUB_Tree Create( string name );
	static QUB_Tree CreateComment( string text );
	static QUB_Tree Open( string path, bool readOnly = false );
	static QUB_Tree ReadText( string path );
	static QUB_Tree CreateFromString( char *s );
	static QUB_Tree CreateFromStream( std::istream& in );
	static QUB_Tree CreateFromBuffer( const char *buffer );

	QUB_Tree(); // null
	QUB_Tree( QTR_Impl *n ); // refer to n
	QUB_Tree( const QUB_Tree &n ); // refer to same as n
	virtual ~QUB_Tree();

	QUB_Tree& operator= (const QUB_Tree &n); // refer to same as n

	bool equals (const QUB_Tree& n) const;

	bool isNull() const         { return (impl == 0); }

	bool save();
	bool saveAs( string path );
	bool saveAsTemp();
	bool close();
	bool inFile();
	string path() const;

	bool saveTextCopy( string path ) const;
	string toString() const;
	void toStream( std::ostream& out, string indent ) const;
	
	unsigned int toBufferSize() const;
	unsigned int toBuffer(char *buffer) const;

	QUB_Tree clone( bool deep=true ) const;

	QUB_Tree appendChild( string name );
	QUB_Tree appendChild( QUB_Tree newchild );
	QUB_Tree appendClone( QUB_Tree orig, bool deep=true );
	QUB_Tree insertChild( string name );
	QUB_Tree insertChild( QUB_Tree newchild );
	QUB_Tree insertClone( QUB_Tree orig, bool deep=true );
	QUB_Tree insertChild( QUB_Tree after, string name );
	QUB_Tree insertChild( QUB_Tree after, QUB_Tree newchild );
	QUB_Tree insertClone( QUB_Tree after, QUB_Tree orig, bool deep=true );

	void removeChild( QUB_Tree child ); // please use iterators instead for speed

	typedef QUB_TreeIter iterator;
	iterator children() const;
	iterator end()      const;
	iterator find( string name ) const;
	iterator rfind( string name ) const;

	QUB_Tree parent() const;
	QUB_Tree child() const;
	QUB_Tree child( string name ) const;
	QUB_Tree sibling() const;
	QUB_Tree sibling( string name ) const;
	QUB_Tree siblingSameName() const;
	QUB_Tree operator[] (string childName) const; // creates if none

	unsigned int & flags() const;
	bool isPreload() const;
	void setPreload( bool pl=true );
  
	// after editing data from a node in an open file,
	// call setChanged() to ensure the data's written on file.save()
	void setChanged( bool ch=true );

	int lock( int timeoutMS = -1 );
	void unlock();

	string name() const;
	void setName( string name );

	// information about the node's data
	bool dataIs( QTR_DataType type ) const;

	QTR_DataType dataType()  const {
		return (QTR_DataType) ( impl ? impl->ondisk.dataType : 0 );
	}
	unsigned int dataSize()  const { return ( impl ? impl->ondisk.dataSize : 0 ); }
	unsigned int dataCount() const { return dataRows() * dataCols(); }
	unsigned int dataCols()  const { return ( impl ? QTR_NumCol( impl ) : 0 ); }
	unsigned int dataRows()  const { return ( impl ? impl->ondisk.dataCount : 0 ); }
	unsigned int loadedFirst() const { return ( impl ? impl->loadStart : QTR_LOAD_NONE ); }
	unsigned int loadedLast() const   { return ( impl ? impl->loadEnd : QTR_LOAD_NONE ); }
	unsigned int loadedCount() const {
		if ( loadedFirst() == QTR_LOAD_NONE )
			return 0;
		else
			return ( loadedLast() - loadedFirst() + 1 );
	}
	unsigned int loadedFirstItem() const { return ( impl ? (impl->loadStart * dataCols()) : QTR_LOAD_NONE ); }
	unsigned int loadedLastItem() const   { return ( impl ? ((impl->loadEnd+1) * dataCols() - 1): QTR_LOAD_NONE ); }
	unsigned int loadedItemCount() const {
		if ( loadedFirst() == QTR_LOAD_NONE )
			return 0;
		else
			return ( loadedLastItem() - loadedFirstItem() + 1 );
	}

	// raw data access
	void *data( bool autoload = true );  // autoload: load it all if nothing's loaded.

	bool setData( string data );
	bool setNumData( QTR_DataType type, unsigned int rows, unsigned int cols, void *inits = 0 );
	bool setNumData( QTR_DataType type, unsigned int rows, unsigned int cols, void **inits );
	bool setData( QTR_DataType type, unsigned int size, unsigned int count, void *inits = 0 );
  
	template<class T>
	bool setData( QTR_DataType type, T init ) {
		if ( impl )   return setNumData( type, 1, 1, &init );
		else          return false;
	}

	// for reading text files -- don't even use it in a multithread / startWrite situation
	void setDataAsExpr( const char *expr, unsigned int flags, bool isUnsigned );

	bool resizeData( unsigned int newrowcount );
	bool clearData();
  
	bool loadData( unsigned int firstRow=0, unsigned int lastRow=QTR_ALL, bool doRead=true );
	bool unloadData( bool doWrite=true );

	unsigned int readDataInto( void *buf, unsigned int firstRow, unsigned int lastRow ) const;
	unsigned int readDataInto( void **buf, unsigned int firstRow, unsigned int lastRow ) const;
	unsigned int writeDataFrom( void *buf, unsigned int firstRow, unsigned int lastRow ); 
	unsigned int writeDataFrom( void **buf, unsigned int firstRow, unsigned int lastRow ); 
  
  
  // enforceT(count, defs):
  //   if the data isn't already of type T, length count:
  //     sets data to T, count
  //     copies and converts values from prev. data where possible,
  //     uses defaults where conversion is impossible
  // bool enforceInt( unsigned int count, int *defaults=0, int timeoutMS = -1 );
  // bool enforceFloat( unsigned int count, float *defaults=0, int timeoutMS = -1 );
  // bool enforceDouble( unsigned int count, double *defaults=0, int timeoutMS = -1 );
  
	// getting data as a specified type, either converting arrays or getting as char arrays
	void dataToStream( std::ostream& out, string indent, bool forFile = false ) const; // forFile includes line comments and '\'
	string dataAsString( bool convert=true ) const;
	void dataAsStrings( vector<string>& lines, bool convert=true ) const; // broken on CR and/or LF

	// these r and c should be absolute -- the same regardless of load bounds
	int dataAsInt( int def=0, int r=0, int c=0 );
	double dataAsDouble( double def=0.0, int r=0, int c=0 );

	template <class T>
	T& dataAs( int i, T dummy ) { // T must match dataType and i = r*ncol + c must be loaded
		T *dat = (T *) data();
		return dat[i - loadedFirst() * dataCols()];
	}

	template <class T>
	T& dataAs( int r, int c, T dummy ) { // T must match dataType and r,c must be loaded
		T *dat = (T *) data();
		return dat[ (r - loadedFirst()) * dataCols() + c ];
	}

	// fill your array, converting if needed/possible
	unsigned int getDataAsInts( int *buf, unsigned int firstRow, unsigned int lastRow );
	unsigned int getDataAsInts( int **buf, unsigned int firstRow, unsigned int lastRow );
	unsigned int getDataAsFloats( float *buf, unsigned int firstRow, unsigned int lastRow );
	unsigned int getDataAsFloats( float **buf, unsigned int firstRow, unsigned int lastRow );
	unsigned int getDataAsDoubles( double *buf, unsigned int firstRow, unsigned int lastRow );
	unsigned int getDataAsDoubles( double **buf, unsigned int firstRow, unsigned int lastRow );

	// manipulate same-line comments in the text-file
	bool hasLineComment() const;
	string getLineComment() const;
	void setLineComment( string linecom );
  
	// QUB_TreeImpl (the underlying representation) can be used with "Node.h"
	QTR_Impl *getImpl() const   { return impl; } // borrowed ref

protected:

	void readChildren( std::istream& in, istream_linereader& rdr );

	QTR_Impl *impl;
};

QTR_API bool operator== (const QUB_Tree& lhs, const QUB_Tree& rhs);

// a lot like the stl iterators; it iterates over the children of a node.

class QTR_API QUB_TreeIter
{
 public:
  QUB_TreeIter(); // of children of null
  QUB_TreeIter( const QUB_Tree &p );
  QUB_TreeIter( const QUB_Tree &p, int what );
  QUB_TreeIter( const QUB_TreeIter& copy );
  virtual ~QUB_TreeIter();

  QUB_TreeIter& operator= (const QUB_TreeIter& copy);

  bool equals (const QUB_TreeIter &it) const;

  QUB_Tree operator* (); // QUB_Tree n = *iter; (the current node)
  QUB_Tree* operator-> ();
  virtual QUB_TreeIter& operator++ (); // move iter to next node
  virtual QUB_TreeIter& operator-- ();
  
  QUB_TreeIter& next( string name ); // next with name (not including current)
  QUB_TreeIter& nextSameName();
  virtual QUB_TreeIter& prev( string name );
  virtual QUB_TreeIter& prevSameName();

  // find by bool functor( QUB_Tree & )
  template<class Test>
  QUB_TreeIter& nextByTest( Test& found ) {
    while ( ! node.isNull() ) {
      ++(*this);
      if ( found( node ) )
	break;
    }
    return *this;
  }
  
  virtual int index(); // correct if nobody else is changing the tree

  virtual void insert( QUB_Tree n );
  virtual void insert( string name );
  virtual QUB_Tree remove();

  QUB_Tree getParent();
  // better than (*iter).parent() since it still works if iter == parent.end()

// protected:
  virtual QUB_Tree prevNode();

  QUB_Tree parent; // valid or always end
  QUB_Tree node; //   null if constr.END, end, or empty
  int ix;       //    -1 : constr.END
  void *prevNodes;
};

class QTR_API QUB_TreeMonoIter : public QUB_TreeIter
{
public:
  QUB_TreeMonoIter(); // of children of null
  QUB_TreeMonoIter( const QUB_Tree &p );
  QUB_TreeMonoIter( const QUB_Tree &p, int what );
  QUB_TreeMonoIter( const QUB_TreeIter& copy );
  virtual ~QUB_TreeMonoIter();
  QUB_TreeMonoIter& operator= (const QUB_TreeIter& copy);

  QUB_TreeIter& operator++ (); // skip prev stuff
  QUB_TreeIter& operator-- (); // the rest throw exceptions
  QUB_TreeIter& prev( string name );
  QUB_TreeIter& prevSameName();
  
  int index(); // don't try to find ix of end, just return -1

  QUB_Tree prevNode();
};
  
QTR_API bool operator== (const QUB_TreeIter& lhs, const QUB_TreeIter& rhs);
QTR_API bool operator!= (const QUB_TreeIter& lhs, const QUB_TreeIter& rhs);


#ifdef __cplusplus    /*`__cplusplus' is #defined iff compiler is C++*/
extern "C" {
#endif

QTR_API char *QTR_DataToCString( QTR_Impl *impl );
QTR_API char *QTR_ToCString( QTR_Impl *impl );
QTR_API void QTR_FreeDataCString( char *cstr );
QTR_API QTR_Impl* QTR_FromString( char *cstr );

QTR_API QTR_Impl* QTR_FromTextFile( char *path );
QTR_API void QTR_SaveTextFile( QTR_Impl *impl, char *path );
QTR_API QTR_Impl* QTR_FromTBL( char *path );

#ifdef __cplusplus
}
#endif


#endif
