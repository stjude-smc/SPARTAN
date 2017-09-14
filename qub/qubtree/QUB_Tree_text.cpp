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

#include "notemplatewarning.h"

#ifdef _USRDLL
  #define _QTR_
#endif

#include "QUB_Tree.h"
#include <strstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "CountedPtr.h"
#include "istream_readline.h"

// reading text

class QTR_TextLine
{
	void readFlags( char* &line, int &len );
	void readMultiNum( std::istream& in, istream_linereader& rdr );
	void readMultiStr( std::istream& in, istream_linereader& rdr );
	void finishReadingData( std::istream& in, istream_linereader& rdr );
public:

	QTR_TextLine( std::istream& in, istream_linereader& rdr );

	bool isEmpty;

	string name;

	bool hasData;
	string data;
	bool nopreload, isUnsigned, isMatrix, isString;

	bool isComment, hasComment;
	string comment;

	bool hasOpen, hasClose;

	bool ok;
};

void QTR_TextLine::readMultiNum( std::istream& in, istream_linereader& rdr )
{
	int len;
	char *line;
	char *stopper = 0;
	
	do {
		line = rdr.readLine( in,  len );
		stopper = strchr(line, ')');
		data.append("\r\n");
		data.append( line, stopper ? stopper - line + 1
								   : len );
	}
	while ( ! stopper );
}

void QTR_TextLine::readMultiStr( std::istream& in, istream_linereader& rdr )
{
	int len;
	char *line;
	
	data.resize(data.length()-1);

	do {
		line = rdr.readLine( in,  len );
		data.append("\r\n");
		data.append( line, len );
	}
	while ( strrchr(line, '\\') == (line + len - 1) );
}

void QTR_TextLine::readFlags( char* &line, int &len )
{
	if ( strncmp(line, "NO_PRELOAD", 10) == 0 ) {
		nopreload = true;
		line += 10;
		len -= 10;
	}
	else if ( strncmp(line, "UNSIGNED", 8) == 0 ) {
		isUnsigned = true;
		line += 8;
		len -= 8;
	}
	else if ( strncmp(line, "MATRIX", 6) == 0 ) {
		isMatrix = true;
		line += 6;
		len -= 6;
	}
	else if ( strncmp(line, "STRING", 6) == 0 ) {
		isString = true;
		line += 6;
		len -= 6;
	}
	else {
		return;
	}

	while ( (line[0] == ' ') || (line[0] == '\t') ) {
		++line;
		--len;
	}

	readFlags( line, len );
}

void QTR_TextLine::finishReadingData( std::istream& in, istream_linereader& rdr )
{
	if ( data[ data.length() - 1 ] == '\\' )
		readMultiStr( in, rdr );
	//else if ( data =~ /^[ \t]*\(/ )
	else if ( data.find(')') == string::npos ) {
		int openpos = data.find("(");
		if ( openpos != string::npos ) {
			bool openIsFirst = true;
			for ( int i=0; i<openpos; ++i )
				if ( (data[i] != ' ') && (data[i] != '\t') )
					openIsFirst = false;
			if ( openIsFirst )
				readMultiNum( in, rdr );
		}
	}
}

QTR_TextLine::QTR_TextLine( std::istream& in, istream_linereader& rdr ) : 
		isEmpty( false ), 
		hasData( false ), 
		nopreload( false ), 
		isUnsigned( false ), 
		isMatrix( false ), 
		isString( false ), 
		isComment( false ), 
		hasComment( false ), 
		hasOpen( false ), 
		hasClose( false ), 
		ok( false ) {

	int len;
	char *line = rdr.readLine( in, len );
	char *pos;

	if ( len || in )
		ok = true;
	else
		return;

	while ( (line[0] == ' ') || (line[0] == '\t') ) {
		++line;
		--len;
	}

	if ( strchr(line, '{') ) {
		hasOpen = true;
		return;
	}
	if ( strchr(line, '}') ) {
		hasClose = true;
		return;
	}

	if ( pos = strstr(line, "#\\") ) {
		hasComment = true;
		comment = string(pos+2);
		pos[0] = 0;
		len = pos - line;
	}

	if ( ! len ) {
		if ( hasComment )
			isComment = true;
		else
			isEmpty = true;
	}
	else {
		readFlags( line, len );

		pos = strchr(line, '=');
		if ( pos ) {
			data = pos + 1;
			
			do {
				--pos;
			} while ( (pos >= line) && ((*pos == ' ') || (*pos == '\t')) );
			
			len = pos - line + 1;
			name = string(line, len);

			hasData = (0!=data.length());
			finishReadingData( in, rdr );
		}
		else if ( line[0] == '(' ) {
			data = line;
			hasData = true;
			finishReadingData( in, rdr );
		}
		else {
			name = line;
			isMatrix = isUnsigned = isString = false;
		}
	}

	//  if (/(.*) ?=(.*)/) then name, data, hasData = \1, \2, true
	//                     else name, hasData = line, false
}


QUB_Tree QUB_Tree::CreateFromStream( std::istream& in )
{
	QUB_Tree tree;
	istream_linereader rdr;

	while ( in ) {
		QTR_TextLine line( in, rdr );

		if ( line.name.length() )
			tree = QUB_Tree::Create( line.name );
		
		if ( line.hasOpen ) {
			if ( tree.isNull() )
				tree = QUB_Tree::Create("");
			tree.readChildren( in, rdr );
			break;
		}
	}

	return tree;

	// read up to and including line with { or eof
	// discern the name from that stretch
	// create a node with that name
	// read children until line with } -- tree.readChildren( in, istream_linereader );
}


void QUB_Tree::readChildren( std::istream& in, istream_linereader& rdr )
{
	QUB_Tree lastChild, actualLastChild; // actual may be a comment, but last is the last non-comment; the one that gets children

	while ( in ) {
		QTR_TextLine line( in, rdr );
		if ( ! line.ok )
			break;
		if ( line.isEmpty )
			continue;

		if ( line.hasClose )
			break;
		if ( line.hasOpen )
			lastChild.readChildren( in, rdr );
		else {
			if ( line.isComment ) {
				actualLastChild = insertChild( actualLastChild, "QTR_COMMENT" );
			}
			else {
				actualLastChild = lastChild = insertChild( actualLastChild, line.name );
				if ( line.hasComment )
					lastChild.setLineComment( line.comment );
				if ( line.nopreload )
					lastChild.setPreload( false );
				if ( line.isString ) {
					lastChild.setData( line.data );
				}
				else if ( line.hasData ) {
					lastChild.setDataAsExpr( line.data.c_str(), line.isMatrix ? QTR_FLAG_MATRIX : 0, line.isUnsigned );

					// if that expr was a few multiplexed nodes, then actualLastChild is no longer with us
					if ( actualLastChild.parent().isNull() ) {
						QUB_TreeIter it = end();
						--it;
						actualLastChild = lastChild = *it;
					}
				}
			}
		}
	}
}

// *************** read expressions from strings ********************

void QTR_SetDataAsExpr( QTR_Impl *impl, const char *expr, unsigned short flags, int isUnsigned );

void QUB_Tree::setDataAsExpr( const char *expr, unsigned int flags, bool isUnsigned )
{
	if ( impl )
		QTR_SetDataAsExpr( impl, expr, flags, isUnsigned );
}

bool QTR_SkipWhite( std::istream &rdr )
{
	int c;
	while ( ! rdr.eof() && ( (c = rdr.get()) != -1 ) ) {
		switch (c) {
			case ' ':
			case '\n':
			case '\r':
			case '\t':
				break;
			default:
				rdr.putback(c);
				return true;
		}
	}
	return false;
}

template<class T>
string QTR_MakeString(T& x)
{
	ostrstream os;
	os << x << '\0';
	string s( os.str() );
	os.rdbuf()->freeze(0);
	return s;
}

template<class ITYPE>
class QTR_DorIArr
{
public:
	vector<ITYPE> ivec;
	vector<double> dvec;
	bool isI, isD;
	QTR_DorIArr() { isI = true; isD = false; }
	void addDefault() {
		if ( isI )
			ivec.push_back( 0 );
		else if ( isD )
			dvec.push_back( 0.0 );
	}
	void toDouble() {
		if ( isI ) {
			vector<double> nd;
			int i;
			for ( i=0; i<int(ivec.size()); ++i )
				nd.push_back( ivec[i] );
			ivec.clear();
			dvec = nd;
			isI = false;
			isD = true;
		}
	}
	typedef CountedPtr< QTR_DorIArr<ITYPE> > Ptr;
};

template<class ITYPE>
istream& operator>> (istream& in, QTR_DorIArr<ITYPE>& arr)
{
	in >> ws;
	if ( ! in.good() )
		return in;
  
	int i;
	double d;

	in >> d;

	if ( in.fail() ) {
		arr.dvec.clear();
		arr.isD = false;
	}
	else if ( arr.isI ) {
		i = (int) d;
		if ( d - i == 0 ) {
			arr.ivec.push_back( i );
		}
		else {
			arr.toDouble();
			arr.dvec.push_back( d );
		}
	}
	else {
		arr.dvec.push_back( d );
	}
  
	if ( in.good() ) {
		in >> ws;
		if ( in.fail() )
			in.clear();
	}

	return in;
}

template<class T>
void QTR_SetDataT( QUB_Tree node, QTR_DataType type, vector< vector<T>* >& rows )
{
	int datlen = rows.size() * rows[0]->size();
	QTR_SetupData( node.getImpl(), type, QTR_DataSize[type], datlen );

	T *loadbuf = (T *) node.data();
	int dat_i, i, j;

	for ( i=0, dat_i=0; i<int(rows.size()); ++i ) {
		vector<T>& row = *rows[i];
		for ( j=0; j<int(row.size()); ++j, ++dat_i )
			loadbuf[dat_i] = row[j];
	}
}

void QTR_Matricize( QUB_Tree node, int ncol )
{
	if ( node.dataCount() == 0 )
		return;
	if ( ((int) node.dataType()) < ((int) QTR_TYPE_UCHAR) )
		return;
	if ( int(node.dataSize()) != QTR_DataSize[ node.dataType() ] )
		return;
	if ( node.dataCount() % ncol )
		return;

	QTR_Impl *impl = node.getImpl();
	if ( impl->loadStart != QTR_LOAD_NONE )
		if ( (impl->loadStart % ncol) || (impl->loadEnd+1 % ncol) )
			node.unloadData();

	impl->ondisk.dataCount /= ncol;
	impl->ondisk.dataSize *= ncol;

	if ( impl->loadStart != QTR_LOAD_NONE ) {
		impl->loadStart /= ncol;
		impl->loadEnd = (impl->loadEnd+1 / ncol) - 1;
	}
}

template<class INT> // int or unsigned int
bool QTR_SetDataAsMatrix( QUB_Tree node, const char *s, QTR_DataType itype ) {
	int failed = false;

	int buflen = 256;
	int linelen = 0;
	istrstream allIn( s );

	istream_linereader oLineReader;

	vector< typename QTR_DorIArr<INT>::Ptr > rows;
	// while get next line
	//   if line not empty
	//     push back row from text
	while ( allIn ) {
		typename QTR_DorIArr<INT>::Ptr row( new QTR_DorIArr<INT> );
		istrstream lineIn( oLineReader.readLine( allIn, linelen ) );
		while ( lineIn.good() )
			lineIn >> *row;
		if ( row->ivec.size() || row->dvec.size() ) // read only numbers, at least one
			rows.push_back( row );
		if ( lineIn.fail() ) {
			failed = true;
			break;
		}
	}
	if ( rows.size() == 0 ) {
		failed = true;
	}
	if ( ! failed ) {
		// make line arrays be uniform type and length
		bool all_i = true;
		int maxlen = 0;

		typename QTR_DorIArr<INT>::Ptr row;
		int i;

		for ( i=0; i<int(rows.size()); ++i ) {
			row = typename QTR_DorIArr<INT>::Ptr( rows[i] );
			int isize = row->ivec.size();
			int dsize = row->dvec.size();

			if ( all_i ) {
				if ( isize ) {
					if ( maxlen < isize )
						maxlen = isize;
				}
				else
					all_i = false;
			}
			else {
				if ( maxlen < dsize )
					maxlen = dsize;
			}
		}
    
		for ( i=0; i<int(rows.size()); ++i ) {
			row = typename QTR_DorIArr<INT>::Ptr( rows[i] );
			if ( ! all_i ) {
				row->toDouble();
				row->dvec.resize( maxlen );
			}
			else
				row->ivec.resize( maxlen );
		}
    
		// set data: concat(rows)
		if ( all_i ) {
			vector< vector<INT>* > rowvecs;
			for ( i=0; i<int(rows.size()); ++i )
				rowvecs.push_back( &(rows[i]->ivec) );
			QTR_SetDataT( node, itype, rowvecs );
		}
		else {
			vector< vector<double>* > rowvecs;
			for ( i=0; i<int(rows.size()); ++i )
				rowvecs.push_back( &(rows[i]->dvec) );
			QTR_SetDataT( node, QTR_TYPE_DOUBLE, rowvecs );
		}

		QTR_Matricize( node, node.dataCount() / rows.size() );
	}
  
	return ! failed;
}

template<class INT>
bool QTR_SetDataAsList( QUB_Tree node, const char *s, QTR_DataType itype )
{
	istrstream in( s );
	typename QTR_DorIArr<INT>::Ptr vals( new QTR_DorIArr<INT> );
	while ( in.good() )
		in >> *vals;
  
	if ( vals->ivec.size() )
		node.setNumData( itype, vals->ivec.size(), 1, &(vals->ivec[0]) );
	else if ( vals->dvec.size() )
		node.setNumData( QTR_TYPE_DOUBLE, vals->dvec.size(), 1, &(vals->dvec[0]) );
	else
		return false;
  
	return true;
}

// no empty strings or plain-old-lists
// if replacing (name empty, no children) must have parent
template<class INT>
bool QTR_SetDataAsLists( QUB_Tree node, const char *s, QTR_DataType itype )
{
	int failed = false;

	int buflen = 256;
	int linelen = 0;
	istrstream allIn( s );
	vector<string> hdrs;
	string hdr;

	istream_linereader oLineReader;

	istrstream hdrIn( oLineReader.readLine( allIn, linelen ) );
	while ( hdrIn.good() ) {
		hdrIn >> hdr;
		if ( ! hdrIn.fail() )
			hdrs.push_back( hdr );
	}
  
	int nhdr = hdrs.size();
	if ( nhdr == 0 )
		failed = true;
  
	vector< typename QTR_DorIArr<INT>::Ptr > cols;
	for ( int i=0; i<nhdr; ++i )
		cols.push_back( typename QTR_DorIArr<INT>::Ptr( new QTR_DorIArr<INT> ) );
  
	// read each row column by column,
	// failing on non-numbers,
	// adding defaults on premature row-end
	while ( allIn && ! failed ) {
		istrstream lineIn( oLineReader.readLine( allIn, linelen ) );
		for ( int i=0; i<nhdr; ++i ) {
			if ( lineIn.eof() ) {
				cols[i]->addDefault();
			}
			else {
				lineIn >> * cols[i];
				if ( lineIn.fail() ) {
					failed = true;
					break;
				}
			}
		}
	}
  
	if ( ! failed ) {
		QUB_TreeMonoIter inserter;

		// decide if these are my children or my replacements
		if ( node.name().length() || ! node.child().isNull() ) {
			inserter = node.children();
			while ( ! (*inserter).isNull() )
				++inserter;
		}
		else {
			inserter = node.parent().children();
			while ( (! (*inserter).isNull()) && ( (*inserter).getImpl() != node.getImpl() ) )
				++inserter;
			inserter.remove(); // self
		}
    
		// insert node per col
		for ( int i=0; i<nhdr; ++i ) {
			QUB_Tree cnode = QUB_Tree::Create(hdrs[i]);
			inserter.insert( cnode );
			++inserter;
      
			if ( cols[i]->isI )
				cnode.setNumData( itype, cols[i]->ivec.size(), 1, &(cols[i]->ivec[0]) );
			else // isD
				cnode.setNumData( QTR_TYPE_DOUBLE, cols[i]->dvec.size(), 1, &(cols[i]->dvec[0]) );
		}
	}
  
	return ! failed;
}

template<class INT>
void QTR_SetDataAsExprT( QUB_Tree node, const char *expr, unsigned short flags, QTR_DataType itype )
{
	if ( node.isNull() )
		return;
  
	int len = strlen(expr);
	if ( len == 0 ) {
		node.clearData();
		return;
	}
  
	bool done = false;

	int first = 0;
	int firstChar = expr[first];
	while ( first < len && ((firstChar == ' ') || (firstChar == '\t')) )
		firstChar = expr[ ++first ];
	int lastChar = expr[len - 1];
	while ( len && ((lastChar == ' ') || (lastChar == '\t')) )
		lastChar = expr[ (--len) - 1 ];

	string contents;
	bool hasParen = false;
	if ( firstChar == '(' && lastChar == ')' ) {
		contents = string( expr, first+1, len - 2 );
		hasParen = true;
	}
	else
		contents = expr;
	istrstream allIn( contents.c_str() );
  
	if ( flags & QTR_FLAG_MATRIX ) {
		done = QTR_SetDataAsMatrix<INT>( node, contents.c_str(), itype );
		if ( ! done )
			flags &= ~ QTR_FLAG_MATRIX;
	}
	if ( ! done ) {
		done = QTR_SetDataAsList<INT>( node, contents.c_str(), itype );

		if ( (! done) && hasParen ) {
			done = QTR_SetDataAsLists<INT>( node, contents.c_str(), itype );
		}
	}
  
	if ( ! done )
		node.setData( expr );
	else if ( node.dataCount() ) // no flags if it was replaced
		QTR_SetFlag( node.getImpl(), flags, 1 );

	if ( ! node.isPreload() )
		node.unloadData();
}

void QTR_SetDataAsExpr( QTR_Impl *impl, const char *expr, unsigned short flags, int isUnsigned )
{
	QUB_Tree node( impl );
	if ( isUnsigned )
		QTR_SetDataAsExprT<unsigned int>( node, expr, flags, QTR_TYPE_UINT );
	else
		QTR_SetDataAsExprT<int>( node, expr, flags, QTR_TYPE_INT );
}





// writing text

int QTR_WriteArrays( QTR_Impl *impl, ostream& out, string indent );

int QTR_WriteTextFlags( QUB_Tree node, std::ostream& out )
{
	int plusindent = 0;

	if ( node.dataCount() && ! node.isPreload() ) {
		out << "NO_PRELOAD ";
		plusindent += 11;
	}
	QTR_DataType type = node.dataType();
	if ( (type == QTR_TYPE_UCHAR) || (type == QTR_TYPE_USHORT)
		|| (type == QTR_TYPE_UINT) || (type == QTR_TYPE_ULONG) ) {
		out << "UNSIGNED ";
		plusindent += 9;
	}
	if ( type == QTR_TYPE_STRING ) {
		out << "STRING ";
		plusindent += 7;
	}
	else if ( node.dataCols() > 1 ) {
		out << "MATRIX ";
		plusindent += 7;
	}

	return plusindent;
}

void QUB_Tree::toStream( std::ostream& out, string indent ) const
{
	bool hasLinecom = hasLineComment();
	int plusindent = 0;

	out << indent;

	plusindent += QTR_WriteTextFlags( *this, out );

	out << name();
	plusindent += name().length();

	if ( dataCount() ) {
		out << " =";
		plusindent += 2;

		string datindent = indent;
		datindent.resize( datindent.length() + plusindent );
		for ( int i=int(indent.length()); i<int(datindent.length()); ++i )
			datindent[i] = ' ';

		dataToStream( out, datindent, true );
	}
	else {
		if ( hasLinecom )
			out << "#\\" << getLineComment();
	}
	out << "\r\n";

	QUB_Tree ch = child();
	if ( hasLinecom )
		ch = ch.sibling();
	if ( ! ch.isNull() ) { // ok there are more children than the line comment
		out << indent << "{\r\n";
		string chindent = indent + "\t";
		
		for ( QUB_TreeMonoIter ci = children(); ! ci->isNull(); ++ci ) {
			int nskip;
			do {
				nskip = QTR_WriteArrays( ci->getImpl(), out, chindent );
				for ( int i=0; i<nskip; ++i )
					++ci;
			} while ( nskip );

			if ( (! ci->isNull()) && (ci->name() != "QTR_LINE_COMMENT") ) {
				ci->toStream( out, chindent );
			}
		}

		out << indent << "}\r\n";
	}
}

QTR_API ostream& operator<< (ostream& out, const QUB_Tree node)
{
	node.toStream(out, "");
	return out;
}


// ********************* array and matrix output *******************

// load partially?

#define QTR_MPLEX_MIN   10
#define QTR_LIST_SPILL  10

class QTR_xWritableSeq
{
public:
  virtual ~QTR_xWritableSeq() {}
  virtual void writeNext( ostream& out ) {}
  virtual operator bool() { return false; }
};

typedef CountedPtr<QTR_xWritableSeq> QTR_xSeqPtr;

template<class T>
class QTR_xDataWriter : public QTR_xWritableSeq
{
public:
  T *data;
  int i, N;
  
  QTR_xDataWriter( QUB_Tree n )
    : i( 0 ), N( n.dataCount() )
  {
    n.loadData(); // make sure it's loaded fully
    data = (T*) n.data();
  }
  
  virtual void writeNext( ostream& out ) {
    if ( i < N )
      out << data[i++];
  }
  
  virtual operator bool() {
    return ( i < N );
  }
};

template<class T>
class QTR_xCharWriter : public QTR_xWritableSeq // ignores '\0'
{
public:
  T *data;
  int i, N;
  
  QTR_xCharWriter( QUB_Tree n )
    : i( 0 ), N( n.dataCount() )
  {
    n.loadData(); // make sure it's loaded fully
    data = (T*) n.data();
  }
  
  virtual void writeNext( ostream& out ) {
	  while ( (i < N) && ! data[i] )
		  ++i;
	  if ( i < N )
		  out << data[i++];
  }
  
  virtual operator bool() {
    return ( i < N );
  }
};

template<class T>
class QTR_xByteWriter : public QTR_xWritableSeq
{
public:
  T *data;
  int i, N;
  
  QTR_xByteWriter( QUB_Tree n )
    : i( 0 ), N( n.dataCount() )
  {
    n.loadData(); // make sure it's loaded fully
    data = (T*) n.data();
  }
  
  virtual void writeNext( ostream& out ) {
    if ( i < N )
      out << (int) data[i++];
  }
  
  virtual operator bool() {
    return ( i < N );
  }
};

QTR_xSeqPtr QTR_GetDataSeq( QUB_Tree node )
{
  if ( node.dataIs( QTR_TYPE_UNKNOWN ) || node.dataIs( QTR_TYPE_UCHAR ) )
    return QTR_xSeqPtr( new QTR_xByteWriter<unsigned char>( node ) ); // wrong oh well
  else if ( node.dataIs( QTR_TYPE_CHAR ) ) // no unicode
    return QTR_xSeqPtr( new QTR_xByteWriter<char>( node ) );
  else if ( node.dataIs( QTR_TYPE_STRING ) ) // no unicode
    return QTR_xSeqPtr( new QTR_xCharWriter<char>( node ) );
  else if ( node.dataIs( QTR_TYPE_POINTER ) || node.dataIs( QTR_TYPE_UINT ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<unsigned int>( node ) );
  else if ( node.dataIs( QTR_TYPE_USHORT ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<unsigned short>( node ) );
  else if ( node.dataIs( QTR_TYPE_SHORT ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<short>( node ) );
  else if ( node.dataIs( QTR_TYPE_INT ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<int>( node ) );
  else if ( node.dataIs( QTR_TYPE_ULONG ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<unsigned long>( node ) );
  else if ( node.dataIs( QTR_TYPE_LONG ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<long>( node ) );
  else if ( node.dataIs( QTR_TYPE_FLOAT ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<float>( node ) );
  else if ( node.dataIs( QTR_TYPE_DOUBLE ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<double>( node ) );
  else if ( node.dataIs( QTR_TYPE_LDOUBLE ) )
    return QTR_xSeqPtr( new QTR_xDataWriter<long double>( node ) );
  
  return QTR_xSeqPtr( new QTR_xWritableSeq );
}

int QTR_WriteArrays( QTR_Impl *impl, ostream& out, string indent )
{
	QUB_Tree start( impl ), node;
	bool hasSigned = false;
	bool hasUnsigned = false;
	bool hasPreload = false;
	bool hasNoPreload = false;
	int narr = 0;
	int arrlen = 0;
	int i, j;

	// count number of matching consecutive arrays
	// by stopping at first string, too-short, or mismatch (size, signedness of int)
	for ( node = start; ! node.isNull(); node = node.sibling(), ++narr ) {
		if ( node.name().length() == 0 )
			break;
    
		if ( node.name().find(" ") != string::npos )
			break;
    
		if ( ! node.child().isNull() )
			break;

		if ( node.dataCols() > 1 )
			break;
    
		QTR_DataType type = node.dataType();

		if ( type == QTR_TYPE_STRING || type == QTR_TYPE_POINTER || type == QTR_TYPE_UNKNOWN )
			break;

		if ( arrlen == 0 )
			arrlen = node.dataCount();
		else if ( arrlen != int(node.dataCount()) )
			break;
    
		if ( arrlen < QTR_MPLEX_MIN )
			break;

		if ( type == QTR_TYPE_UCHAR || type == QTR_TYPE_USHORT || type == QTR_TYPE_UINT || type == QTR_TYPE_ULONG ) {
			hasUnsigned = true;
			if ( hasSigned )
				break;
		}
		else if ( type == QTR_TYPE_CHAR || type == QTR_TYPE_SHORT || type == QTR_TYPE_INT || type == QTR_TYPE_LONG ) {
			hasSigned = true;
			if ( hasUnsigned )
				break;
		}
    
		if ( node.isPreload() ) {
			hasPreload = true;
			if ( hasNoPreload )
				break;
		}
		else {
			hasNoPreload = true;
			if ( hasPreload )
				break;
		}
	}
  
	if ( narr > 1 ) {
		// get an output object of the right data type for each array
		vector<QTR_xSeqPtr> seqs;
		seqs.resize( narr );
		for ( i=0, node=start; i<narr; ++i, node = node.sibling() )
			seqs[i] = QTR_GetDataSeq( node.getImpl() );
    
		// write flags
		if ( hasUnsigned )
			out << "UNSIGNED ";
		if ( hasNoPreload )
			out << "NO_PRELOAD ";
    
		// begin
		out << indent << '(';
    
		// write headers
		for ( i=0, node=start; i<narr; ++i, node = node.sibling() ) {
			out << '\t' << node.name();
		}
    
		// write rows
		for ( j=0; j<arrlen; ++j ) {
			out << "\r\n" << indent;
			for ( i=0, node=start; i<narr; ++i, node = node.sibling() ) {
				out << '\t';
				seqs[i]->writeNext( out );
			}
		}
	
		// end
		out << " )\r\n";

		return narr;
	}
	else
		return 0;
}

void QUB_Tree::dataToStream( std::ostream& out, string indent, bool forFile ) const
{
	bool hasCom = forFile ? hasLineComment() : false;

	QTR_DataType type = dataType();
	// if is string, write it
	if ( type == QTR_TYPE_STRING ) {
		vector<string> multis;
		dataAsStrings( multis, false );
		int linecount = multis.size();
		int i;
    
		if ( linecount ) {
			// strip trailing '\0' placed by misguided QFS
			string& lastLine = multis[linecount - 1];
			i = lastLine.size() - 1;
			if ( (i>=0) && (lastLine[i] == '\0') ) {
			  lastLine.resize( i );
			}

			out << multis[0];
			if ( linecount > 1 ) {
				if ( forFile )
					out << '\\';
				if ( hasCom )
					out << "#\\" << getLineComment();
				out << "\r\n";
	
				for ( i=1; i<linecount-1; ++i ) {
					out << multis[i];
					if ( forFile )
						out << "\\";
					out << "\r\n";
				}
				out << multis[linecount-1];
			}
			else {
				if ( hasCom )
					out << "#\\" << getLineComment();
			}
		}
	}
	else {
		// get array writer:
		QTR_xSeqPtr seq( QTR_GetDataSeq(*this) );

		// check if matrix, get dimensions
		unsigned int ncol = dataCols();
		// unsigned int *dim = n.getMatrixDim();
		if ( ncol > 1 ) {
			out << '(';
			while ( *seq ) {
				for ( unsigned int i=0; *seq && i<ncol; ++i ) {
					out << '\t';
					seq->writeNext( out );
				}

				if ( *seq ) { // has more
					if ( hasCom ) {
						hasCom = false; // print it on first line only
						out << "#\\" << getLineComment();
					}
					out << "\r\n" << indent;
				}
				else {
					out << " )";
					if ( hasCom ) {
						out << "#\\" << getLineComment();
					}
				}
			}
		}
		else { // it's a list
			bool multi = dataCount() > QTR_LIST_SPILL;
			if ( multi )
				out << '(';
      
			while ( *seq ) {
				for ( int i=0; *seq && i<QTR_LIST_SPILL; ++i ) {
					if ( i )
						out << '\t';
					seq->writeNext( out );
				}
				if ( hasCom ) {
					hasCom = false;
					out << "#\\" << getLineComment();
				}
				if ( *seq ) {
					out << "\r\n" << indent;
				}
			}
      
			if ( multi )
				out << ')';
			if ( hasCom )
				out << "#\\" << getLineComment();
		}
	}
}

QTR_API char *QTR_DataToCString( QTR_Impl *impl )
{
	QUB_Tree node( impl );
	string cpstr = node.dataAsString( true );
	char *cstr = new char[ cpstr.length() + 1 ];
	memcpy( cstr, cpstr.c_str(), cpstr.length() + 1 );
	return cstr;
}

QTR_API char *QTR_ToCString( QTR_Impl *impl )
{
	ostrstream ost;
	ost << QUB_Tree(impl) << '\0';
	char *cstr = new char[ ost.pcount() ];
	memcpy( cstr, ost.str(), ost.pcount() );
	ost.rdbuf()->freeze(0);

	return cstr;
}

QTR_API void QTR_FreeDataCString( char *cstr )
{
  if ( cstr )
    delete [] cstr;
}

QTR_API QTR_Impl* QTR_FromString( char *cstr )
{
	QUB_Tree node = QUB_Tree::CreateFromString( cstr );
	QTR_Impl *impl = node.getImpl();
	QTR_INCREF( impl );
	return impl;
}

QTR_API QTR_Impl* QTR_FromTextFile( char *path )
{
	QUB_Tree node = QUB_Tree::ReadText( path );
	QTR_Impl *impl = node.getImpl();
	if ( impl )
		QTR_INCREF( impl );
	return impl;
}

QTR_API void QTR_SaveTextFile( QTR_Impl *impl, char *path )
{
	QUB_Tree node( impl );
	node.saveTextCopy( path );
}

QTR_API QTR_Impl* QTR_FromTBL( char *path )
{
	string s;
	int i;
	double d;

	QUB_Tree tbl = QUB_Tree::Create(path);
	QUB_Tree ds = tbl["DataSet"];
	QUB_TreeIter segi = ds.children();

	ifstream in;
	in.open( path );

	while ( ! in.eof() ) {
		int len; // points
		double start, sampling; // ms
		double ll;
		int iter;
		vector<double> amp, sd, occ, ltime, amatrix;
		vector<int> nevent;

		// read one seg
		do {
			if ( ! in ) {
				in.clear();
				in >> s; // leading junk, then "Segment:"
			}
			in >> i;     // Segment index (throwaway)
		}
		while ( ! in && ! in.eof() );

		if ( in.eof() )
			break;

		in >> s;       // Length:
		in >> len;
		in >> s;       // Starting
		in >> s;       // time
		in >> s;       // (ms):
		in >> start;
		in >> s;       // Sampling
		in >> s;       // duration
		in >> s;       // (ms):
		in >> sampling;
		in >> s;       // Likelihood:
		in >> ll;
		in >> s;       // Iterations:
		in >> iter;
		for ( i=0; i<6; ++i )
			in >> s;     // class amp xms occupancy lifetime(ms) nevent
    
		while ( in ) {
			in >> i;     // class;
			if ( ! in ) break; // no more classes

			in >> d;     amp.push_back( d );
			in >> d;     sd.push_back( d );
			in >> d;     occ.push_back( d );
			in >> d;     ltime.push_back( d );
			in >> i;     nevent.push_back( i );
		}
		in.clear();
    
		in >> s;       // Transition
		in >> s;       // Probability
		in >> s;       // Matrix:
		while ( in ) {
			in >> d;
			if ( in )
				amatrix.push_back( d );
		}
		in.clear();

		if ( amatrix.size() ) {
			// store seg
			QUB_Tree seg = QUB_Tree::Create("Segment");
			segi.insert( seg );
			++segi;

			seg["length"].setData( QTR_TYPE_INT, len );
			seg["start"].setData( QTR_TYPE_DOUBLE, start );
			seg["sampling"].setData( QTR_TYPE_DOUBLE, sampling );
			seg["LL"].setData( QTR_TYPE_DOUBLE, ll );
			seg["iterations"].setData( QTR_TYPE_INT, iter );
			seg["amp"].setNumData( QTR_TYPE_DOUBLE, amp.size(), 1, &(amp[0]) );
			seg["sd"].setNumData( QTR_TYPE_DOUBLE, sd.size(), 1, &(sd[0]) );
			seg["occupancy"].setNumData( QTR_TYPE_DOUBLE, occ.size(), 1, &(occ[0]) );
			seg["lifetime"].setNumData( QTR_TYPE_DOUBLE, ltime.size(), 1, &(ltime[0]) );
			seg["nevent"].setNumData( QTR_TYPE_INT, nevent.size(), 1, &(nevent[0]) );
			int nstate = (int) sqrt( double(amatrix.size()) );
			seg["a matrix"].setNumData( QTR_TYPE_DOUBLE, nstate, nstate, &(amatrix[0]) );
		}
	}

	QTR_Impl *tbli = tbl.getImpl();
	QTR_INCREF( tbli );
	return tbli;
}



