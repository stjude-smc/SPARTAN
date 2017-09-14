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

//#include <math.h>
#ifdef _DEBUG
#include <stdio.h>
#endif

#include <string.h>
#include "ExprEval.h"


//----- Special interface for DLL
CExpression g_oExpression;
extern "C" QUBOPT_API void LoadExpression( const char * pchExpr ) {
	g_oExpression.Parse( pchExpr );
	}

extern "C" QUBOPT_API double EvaluateExpression( double t ) { 
	return g_oExpression.Evaluate(t);
	}

//---------------------------------------------------------------------------------
const int TOKEN_t = -1;				// passed parameter t
const int TOKEN_ADD = -2;			// Binary...
const int TOKEN_SUB = -3;
const int TOKEN_MUL = -4; 
const int TOKEN_DIV = -5;
const int TOKEN_POW = -6;
const int TOKEN_SIN = -7;			// Function...
const int TOKEN_COS = -8;
const int TOKEN_LN  = -9;
const int TOKEN_EXP = -10;
const int TOKEN_SGN = -11;
const int TOKEN_ABS = -12;
const int TOKEN_OP  = -13;			// Other...
const int TOKEN_CP  = -14;
const int TOKEN_NULL = -9999;
const int TOKEN_VALUE = 0; // Token >=0 : Index to double array 

const int FCOUNT = 6;
char * g_apchFunctionNames[FCOUNT]= { "sgn", "abs", "sin", "cos", "ln", "exp" };
int g_aiFunctionTokens[FCOUNT]= { TOKEN_SGN, TOKEN_ABS, TOKEN_SIN, TOKEN_COS,
								  TOKEN_LN, TOKEN_EXP };

const int OPCOUNT = 8;		// note includes non operators t,(,) ? 
char g_achOperators[OPCOUNT] = {'+','-','*','/','^','(',')','t'};
int g_aiOperatorTokens[OPCOUNT] = { TOKEN_ADD, TOKEN_SUB, TOKEN_MUL, TOKEN_DIV, 
									TOKEN_POW, TOKEN_OP, TOKEN_CP, TOKEN_t };


//--- constructor / destructor 
CExpression :: CExpression() { 
	m_aTokens=NULL;
	m_aValues=NULL;
	m_aStack=NULL;
	}

CExpression :: CExpression( const char * pchExpr ) { 
	m_aTokens=NULL;
	m_aValues=NULL;
	m_aStack=NULL;
	Parse( pchExpr );
	}

CExpression :: ~CExpression() { 
	if( m_aTokens != NULL ) 
		delete [] m_aTokens;
	if( m_aValues != NULL ) 
		delete [] m_aValues;
	if( m_aStack != NULL ) 
		delete [] m_aStack;
	}

// Set piToken and pnValue for next token or return false.  pchExpr must be lowercase.
bool CExpression :: NextToken( const char * pchExpr, int * piPos, int * piToken, int iPrevToken, double * pnValue ) { 
	int ilen, isub, i, ii;
	char ch;

	// skip whitespace
	while( *(pchExpr+(*piPos))==' ' || *(pchExpr+(*piPos))=='\t' ) 
		++(*piPos);

	// Check for empty string
	ilen = (int) strlen(pchExpr+(*piPos));
	if( ilen == 0 ) 
		return false;

	// Check for unary funcs ( sgn, abs, sin, cos, ln, exp )
	for( i=0; i<FCOUNT; ++i ) {
	  ii= (int) strlen(g_apchFunctionNames[i]);
		if( ilen>=3 && memcmp(pchExpr+*piPos,g_apchFunctionNames[i],ii)==0 ) {
			*piToken = g_aiFunctionTokens[i];
			(*piPos)+=ii;
			return true;
			}
		}
			
	// Check for operators and other single char tokens
	ch = *(pchExpr+*piPos); 
	for( i=0; i<OPCOUNT; ++i ) {
		if( ch==g_achOperators[i] && ( ch!='-' || IsValue(iPrevToken) || iPrevToken==TOKEN_CP ) ) {
			*piToken = g_aiOperatorTokens[i];
			++(*piPos);
			return true;
			}
		}

	// check for all others 
	if( ch=='-' || ch=='.' || ( ch>='0' && ch<='9') ) { 
		bool lDecimal=false;
		*piToken=TOKEN_VALUE;
		isub=0; 
		while( strchr("0123456789.",*(pchExpr+(*piPos)+isub+1) )!=NULL ) {
			if( *(pchExpr+(*piPos)+isub+1)=='.' ) {
				if( lDecimal ) 
					return false;
				lDecimal=true;
				}
			++isub;
			}
			
		char * achBuf = new char[isub+2];
		memcpy( achBuf, pchExpr+(*piPos), isub+1 );
		achBuf[isub+1]=0;
		*pnValue = atof(pchExpr+(*piPos));
		*piPos+=isub+1;
		delete [] achBuf;
		return true;
		}
	return false;
	}

bool CExpression :: Parse( const char * pchExpr ) { 
	if( m_aTokens != NULL ) 
		delete [] m_aTokens;
	if( m_aValues != NULL ) 
		delete [] m_aValues;
	if( m_aStack != NULL ) 
		delete [] m_aStack;

	//----- count tokens to create array 
	int iTokens=0, iPos=0;
	int iValues=0;
	double nValue;
	int iToken=TOKEN_NULL;
	while( NextToken( pchExpr, &iPos, &iToken, iToken, &nValue ) ) {
		if( iToken==0 ) 
			++iValues;
		++iTokens;
		}

	//----- save tokens and values to arrays
	m_aRaw=new int[iTokens];
	m_aValues=new double[iValues];
	iPos=0;
	iValues=0;
	iToken=0;
	int iPrevToken=TOKEN_NULL;
	while( NextToken( pchExpr, &iPos, &m_aRaw[iToken], iPrevToken, &nValue ) ) {
		iPrevToken=m_aRaw[iToken];
		if( m_aRaw[iToken]==TOKEN_VALUE ) {
			m_aRaw[iToken]=iValues;
			m_aValues[iValues++]=nValue;
			}
		++iToken;
		}

	//----- reorder tokens into m_aTokens array 
	m_aTokens=new int[iTokens];
	for( iToken=0; iToken<iTokens; ++iToken ) 
		m_aTokens[iToken]=TOKEN_NULL;
	ParseReduce( 0, iTokens-1 );
	delete [] m_aRaw;

	//----- Finish Token Array processing 
	m_iStackSize=0;
	int iStack=0;
	m_iTokens=0;
	for( iToken=0; iToken<iTokens; ++iToken ) {
		//-- pack m_aTokens array 
		if( m_aTokens[iToken]!=TOKEN_NULL ) 
			m_aTokens[m_iTokens++]=m_aTokens[iToken];
		//-- Calculate stacksize and pre allocate 
		if ( IsBinary(m_aTokens[iToken]) ) 
			iStack--;
		else if( IsValue(m_aTokens[iToken]) ) 
			if( (++iStack) > m_iStackSize ) 
				m_iStackSize=iStack;
		}
	m_aStack=new double[m_iStackSize];
	
	return true;
	}

int CExpression :: ParseFind( int iToken1, int iToken2, int i0, int i1 ) { 
	for( int i=i0; i<=i1; ++i ) 
		if( m_aRaw[i]==iToken1 || m_aRaw[i]==iToken2 ) 
			return i;
	return -1;
	}

int CExpression :: ParseFindL( int iToken1, int iToken2, int i0, int i1 ) { 
	for( int i=i1; i>=i0; --i ) 
		if( m_aRaw[i]==iToken1 || m_aRaw[i]==iToken2 ) 
			return i;
	return -1;
	}

void CExpression :: ParseShiftL( int i0, int i1, int iShift ) { 
	for( int i=i0; i<=i1; ++i ) 
		m_aTokens[i-iShift]=m_aTokens[i];
	}

void CExpression :: ParseReduce( int i0, int i1 ) { 
	//----- Repeatedly Evaluate innermost parentheses
	int iCP, iOP;
	while( (iCP=ParseFind(TOKEN_CP, TOKEN_CP, i0, i1)) >= 0 )  { 
		iOP = ParseFindL( TOKEN_OP, TOKEN_OP,  i0, iCP-1 );
		//-- reduce expr in (...)
		ParseReduce( iOP+1, iCP-1 );
		//-- remove parentheses
		m_aRaw[iOP]=TOKEN_NULL;
		m_aRaw[iCP]=TOKEN_NULL;
		//-- if parentheses are a function eval, copy the function
		if( iOP>0 && IsFunction(m_aRaw[iOP-1]) ) { 
			m_aTokens[iCP]=m_aRaw[iOP-1];
			m_aRaw[iOP-1]=TOKEN_NULL;
			}
		}

	//-----Evaluate Binary :  +,- then *,/ then ^
	int iOpPos;
	if ( (iOpPos=ParseFindL( TOKEN_ADD, TOKEN_SUB, i0, i1))>0 
			|| (iOpPos=ParseFindL( TOKEN_MUL, TOKEN_DIV, i0, i1))>0 
			|| (iOpPos=ParseFindL( TOKEN_POW, TOKEN_POW, i0, i1))>0 ) { 
		ParseReduce( i0, iOpPos-1 );
		ParseReduce( iOpPos+1, i1 );
		ParseShiftL( iOpPos+1, i1 );
		m_aTokens[i1]=m_aRaw[iOpPos];
		m_aRaw[iOpPos]=TOKEN_NULL;
		}

	//----- Evaluate Values 
	else if( IsValue(m_aRaw[i0]) ) {
		m_aTokens[i0]=m_aRaw[i0];
		m_aRaw[i0]=TOKEN_NULL;
		}
	}

double CExpression :: Evaluate( double t ) { 
	int iStack=0;

	for( int i=0; i<m_iTokens; ++i ) {
		switch(	m_aTokens[i] ) { 
			case TOKEN_t :
				m_aStack[iStack++]=t;
				break;
			case TOKEN_ADD :
				m_aStack[iStack-2]+=m_aStack[iStack-1];
				iStack--;
				break;
			case TOKEN_SUB :
				m_aStack[iStack-2]-=m_aStack[iStack-1];
				iStack--;
				break;
			case TOKEN_MUL :
				m_aStack[iStack-2]*=m_aStack[iStack-1];
				iStack--;
				break;
			case TOKEN_DIV :
				m_aStack[iStack-2]/=m_aStack[iStack-1];
				iStack--;
				break;
			case TOKEN_POW :
				m_aStack[iStack-2]=pow(m_aStack[iStack-2],m_aStack[iStack-1]);
				iStack--;
				break;
			case TOKEN_EXP :
				m_aStack[iStack-1]=exp(m_aStack[iStack-1]);
				break;
			case TOKEN_SIN :
				m_aStack[iStack-1]=sin(m_aStack[iStack-1]);
				break;
			case TOKEN_COS :
				m_aStack[iStack-1]=cos(m_aStack[iStack-1]);
				break;
			case TOKEN_LN  :
				m_aStack[iStack-1]=log(m_aStack[iStack-1]);
				break;
			case TOKEN_ABS :
				m_aStack[iStack-1]=fabs(m_aStack[iStack-1]);
				break;
			case TOKEN_SGN :
				m_aStack[iStack-1]=(m_aStack[iStack-1]>=0?1.0:-1.0);
				break;
			default : 
				m_aStack[iStack++]=m_aValues[m_aTokens[i]];
			}
		}
	return m_aStack[0];
	}

bool CExpression :: IsFunction( int iToken ) { 
	return iToken==TOKEN_EXP || iToken==TOKEN_SIN || iToken==TOKEN_COS 
		|| iToken==TOKEN_LN || iToken==TOKEN_SGN || iToken==TOKEN_ABS;
	}
bool CExpression :: IsValue( int iToken ) {
	return iToken >= 0 || iToken==TOKEN_t ; 
	}
bool CExpression :: IsBinary( int iToken ) { 
	return iToken==TOKEN_ADD || iToken==TOKEN_SUB || iToken==TOKEN_MUL 
		|| iToken==TOKEN_DIV || iToken==TOKEN_POW;
	}


#ifdef _DEBUG
void CExpression :: Text( char * achBuffer, int iBufSize ) {
	int i=0, iPos=0;
	while( iPos+20<iBufSize && i<m_iTokens ){
		switch( m_aTokens[i] ) {
			case TOKEN_POW :							// EXPONENT power
				achBuffer[iPos++]='^';
				break;
			case TOKEN_ADD :								// +
				achBuffer[iPos++]='+';
				break;
			case TOKEN_SUB :								// -
				achBuffer[iPos++]='-';
				break;
			case TOKEN_MUL :								// *
				achBuffer[iPos++]='*';
				break;
			case TOKEN_DIV :								// /
				achBuffer[iPos++]='/';
				break;
			case TOKEN_t :								// t - time 
				achBuffer[iPos++]='t';
				break;
			case TOKEN_SIN  :
				memcpy(achBuffer+iPos,"sin",3);
				iPos+=3;
				break;
			case TOKEN_COS  :
				memcpy(achBuffer+iPos,"cos",3);
				iPos+=3;
				break;
			case TOKEN_LN   :
				memcpy(achBuffer+iPos,"ln",2);
				iPos+=2;
				break;
			case TOKEN_EXP  :
				memcpy(achBuffer+iPos,"exp",3);
				iPos+=3;
				break;
			case TOKEN_ABS  :
				memcpy(achBuffer+iPos,"abs",3);
				iPos+=3;
				break;
			case TOKEN_SGN  :
				memcpy(achBuffer+iPos,"sgn",3);
				iPos+=3;
				break;
			case TOKEN_NULL :
				memcpy(achBuffer+iPos,"null",4);
				iPos+=4;
				break;
			default : 
				if( m_aTokens[i] >= 0 ) { 
					sprintf(achBuffer+iPos, "%f", m_aValues[m_aTokens[i]]);
					iPos=strlen(achBuffer);
					}
				else 
					achBuffer[iPos++]='?';	
				break;
			}
		++i;
		achBuffer[iPos++]=' ';
		}
	achBuffer[iPos++]=0;
	}
#endif

