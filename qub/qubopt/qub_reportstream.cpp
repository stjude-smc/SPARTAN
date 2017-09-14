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

#include "qubopt.h"
#include "qub_reportstream.h"

static QUB_Tree QTR_DoCallback(QTR_Callback cb, QUB_Tree arg) {
  QTR_Impl *impl = NULL;
  QUB_Tree result;
  if ( cb ) {
    impl = cb(arg.getImpl());
    if ( impl ) {
      result = impl;
      QTR_DECREF(impl);
    }
  }
  return result;
}


int mil_reportbuf::overflow( int ch )
{
	if ( ch != EOF ) {
		if ( ch == '\n' ) {
			if ( callback ) {
				QUB_Tree msg = QUB_Tree::Create("");
				msg.setData( lineout.str() );

				QTR_DoCallback(callback, msg);
			}
			lineout.str("");
			lineout.clear();
		}
		else {
			lineout << (char) ch;
		}
	}
	return ch;
}

// *************************************************************************

