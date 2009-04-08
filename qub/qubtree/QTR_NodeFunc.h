#ifndef QTR_NODE_FUNC_H
#define QTR_NODE_FUNC_H

#include "QUB_Tree.h"
#include "CountedPtr.h"

class QTR_API QTR_NodeFunc
{
public:
	virtual ~QTR_NodeFunc() {}
	virtual QUB_Tree NodeCall( QUB_Tree arg ) = 0;
};

class QTR_API QTR_NodeNodeFunc
{
public:
	virtual ~QTR_NodeNodeFunc() {}
	virtual QUB_Tree NodeNodeCall( QUB_Tree arg1, QUB_Tree arg2 ) = 0;
};

typedef CountedPtr<QTR_NodeFunc> QTR_NodeFptr;
typedef CountedPtr<QTR_NodeNodeFunc> QTR_NodeNodeFptr;

// wrapping c functions

typedef QUB_Tree (*QTR_NodeFunction)( QUB_Tree );
typedef QUB_Tree (*QTR_NodeNodeFunction)( QUB_Tree, QUB_Tree );

class QTR_API QTR_NodeCFunc : public QTR_NodeFunc
{
public:
	QTR_NodeCFunc( QTR_NodeFunction fn );
	virtual QUB_Tree NodeCall( QUB_Tree arg );
private:
	QTR_NodeFunction func;
};

class QTR_API QTR_NodeNodeCFunc : public QTR_NodeNodeFunc
{
public:
	QTR_NodeNodeCFunc( QTR_NodeNodeFunction fn );
	virtual QUB_Tree NodeNodeCall( QUB_Tree arg1, QUB_Tree arg2 );
private:
	QTR_NodeNodeFunction func;
};

// do-nothing stand-ins
class QTR_NoNodeFunc : public QTR_NodeFunc
{
public:
	QUB_Tree NodeCall( QUB_Tree x ) { return x; }
};

class QTR_NoNodeNodeFunc : public QTR_NodeNodeFunc
{
public:
	QUB_Tree NodeNodeCall( QUB_Tree x, QUB_Tree y ) { return x; }
};

#if defined(_WIN32)
#include <Unknwn.h>

// delphi callbacks: interface, then wrapper to use QUB_Tree

// struct because of default public inheritance and access
class QTR_API QTR_INodeCFunc : public IUnknown
{
public:
	virtual QTR_Impl* _stdcall NodeCall( QTR_Impl *arg ) = 0;
};

class QTR_API QTR_INodeNodeCFunc : public IUnknown
{
public:
	virtual QTR_Impl* _stdcall NodeNodeCall( QTR_Impl *arg1, QTR_Impl *arg2 ) = 0;
};

class QTR_API QTR_NodePasFunc : public QTR_NodeFunc
{
public:
	QTR_NodePasFunc( QTR_INodeCFunc* fn );
	virtual ~QTR_NodePasFunc();
	virtual QUB_Tree NodeCall( QUB_Tree arg );
private:
	QTR_INodeCFunc* func;
};

class QTR_API QTR_NodeNodePasFunc : public QTR_NodeNodeFunc
{
public:
	QTR_NodeNodePasFunc( QTR_INodeNodeCFunc* fn );
	virtual ~QTR_NodeNodePasFunc();
	virtual QUB_Tree NodeNodeCall( QUB_Tree arg1, QUB_Tree arg2 );
private:
	QTR_INodeNodeCFunc* func;
};

#endif


#endif
