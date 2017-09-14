/* Copyright 1998-2012 Research Foundation State University of New York */

/* This file is part of QUB Express.                                     */

/* QUB Express is free software; you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation, either version 3 of the License, or     */
/* (at your option) any later version.                                   */

/* QUB Express is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */

/* You should have received a copy of the GNU General Public License,    */
/* named LICENSE.txt, in the QUB Express program directory.  If not, see */
/* <http://www.gnu.org/licenses/>.                                       */

#ifndef LOADDATACB_PAGER_H
#define LOADDATACB_PAGER_H

#include <iostream>

#include <QUB_Tree.h>

class QUB_Pager {
public:
	int dataCount;
	int pageSize;
	int pageStart;
	double *data;

	QUB_Pager(int dataCount, int pageSize, int pageStart, double *data);
	virtual ~QUB_Pager();
	virtual void nextPage();
	virtual void loadPageContaining( int ix );

	int indexWithinPage( int true_ix ){
		return true_ix % pageSize;
	}		
};

class LoadDataCB_Pager : public QUB_Pager {
	public:
		int segFirst, segLast;
		QUB_Tree ds, seg;
		QTR_Callback LoadData;
		int *loadBounds;
		std::ostream& err;
		int offset;

		LoadDataCB_Pager( QUB_Tree segNode, QTR_Callback loadDataCB, int pagesize, std::ostream& milerr, int samplesOffset=0 );
		virtual ~LoadDataCB_Pager();

		void nextPage();
		void loadPageContaining( int ix );
};

#endif
