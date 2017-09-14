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

#ifdef _USRDLL
  #define _QTR_
#endif
#include "QTR_NodeFunc.h"

QTR_NodeCFunc::QTR_NodeCFunc( QTR_NodeFunction fn )
  : func( fn )
{}

QUB_Tree QTR_NodeCFunc::NodeCall( QUB_Tree arg )
{
	return func( arg );
}

QTR_NodeNodeCFunc::QTR_NodeNodeCFunc( QTR_NodeNodeFunction fn )
  : func( fn )
{}

QUB_Tree QTR_NodeNodeCFunc::NodeNodeCall( QUB_Tree arg1, QUB_Tree arg2 )
{
	return func( arg1, arg2 );
}

#ifdef _WIN32

QTR_NodePasFunc::QTR_NodePasFunc( QTR_INodeCFunc* fn )
  : func( fn )
{
	func->AddRef();
}

QTR_NodePasFunc::~QTR_NodePasFunc()
{
	func->Release();
}

QUB_Tree QTR_NodePasFunc::NodeCall( QUB_Tree arg )
{
	QTR_Impl *res = func->NodeCall( arg.getImpl() );
	QUB_Tree result( res );
	if ( res )
		QTR_DECREF( res );
	return result;
}

QTR_NodeNodePasFunc::QTR_NodeNodePasFunc( QTR_INodeNodeCFunc* fn )
  : func( fn )
{
	func->AddRef();
}

QTR_NodeNodePasFunc::~QTR_NodeNodePasFunc()
{
	func->Release();
}

QUB_Tree QTR_NodeNodePasFunc::NodeNodeCall( QUB_Tree arg1, QUB_Tree arg2 )
{
	QTR_Impl *res = func->NodeNodeCall( arg1.getImpl(), arg2.getImpl() );
	QUB_Tree result( res );
	if ( res )
		QTR_DECREF( res );
	return result;
}

#endif
