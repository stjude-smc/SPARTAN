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
#ifndef EXPREVAL_H
#define EXPREVAL_H

#include "qubopt.h"



class CExpression { 
	private : 
		int m_iTokens, m_iStackSize;
		int * m_aTokens;		// Expression parsed and reordered w/o parentheses
		int * m_aRaw;			// Expression parsed in original order 
		double * m_aValues;
		double * m_aStack;
		bool IsFunction( int iToken );
		bool IsValue( int iToken );
		bool IsBinary( int iToken );
		int ParseFind( int iToken1, int iToken2, int i0, int i1 ); 
		int ParseFindL( int iToken1, int iToken2, int i0, int i1 );
		void ParseShiftL( int i0, int i1, int iShift=1 );
		void ParseReduce( int i0, int i1 );
		bool NextToken( const char * pchExpr, int * piPos, int * piToken, int iPrevToken, double * pnValue );

	public :
		CExpression( const char * pchExpr );
		CExpression();
		~CExpression();
			
		bool Parse( const char * pchExpr ); 
		double Evaluate( double t );

		#ifdef _DEBUG
		//----- This function displays the parse results for debug purposes
		void Text( char * achBuffer, int iBufSize );
		#endif
	};


//----- 6/2005 - Expression Evaluation engine 
	/*
	The following two functions provide a simple interface to a global CExpression 
	object which can be used for c style dll access.
	*/

#ifdef __cplusplus
extern "C" {
#endif


	QUBOPT_API void LoadExpression( const char * pchExpr );
	QUBOPT_API double EvaluateExpression( double t );

#ifdef __cplusplus
}
#endif

	/* 
	Expression parser / evaluator.
	Evaluate f(t) for provided t.  
	f = Expression of + - * / ^ ln exp sin cos ( ) abs sgn
	'-' preceded by operator or '(' is negate, otherwise it is minus
	Ex : 
		CExpression o;
		o.Parse("t^3-60*t^2+1100*t-6000+abs(sin(t/5)*sin(t/4)*100+90*cos(t/6)*cos(t/3))");
		double t1=o.Evaluate(0.532);
		double t2=o.Evaluate(9876.5432);
	*/
#endif
