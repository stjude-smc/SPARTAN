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

#ifndef QUB_REPORTSTREAM_H
#define QUB_REPORTSTREAM_H

#include <sstream>

class mil_reportbuf : public std::streambuf
{
public:
	mil_reportbuf()
	{ setp(0,0), setg(0,0,0); }

	void setTarget( QTR_Callback reportCB )
	{ callback = reportCB; }

protected:
	virtual int overflow( int ch );
private:
	QTR_Callback callback;
	ostringstream lineout;
};

class mil_reportstream : public std::ostream
{
public:
	mil_reportstream( QTR_Callback reportCB )
		: std::basic_ostream<char,std::char_traits<char> >( new mil_reportbuf )
	{
		((mil_reportbuf *) rdbuf())->setTarget( reportCB );
	}
        virtual ~mil_reportstream()
        {
                delete rdbuf();
        }
};



#endif
