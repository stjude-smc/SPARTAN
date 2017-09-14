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

#include <string.h>
#include <iostream>

#include "dfp_optimize.h"

using namespace std;
using namespace fq;

#define  ITMAX 	   100
#define  ITPRT 	   1
#define  EPS 	      3.0e-8
#define	TOLF 	      1.0e-5
#define	TOLX 	      1.0e-15
#define	TOLG 	      5.0e-3
#define VERY_SMALL    1.0e-15

// Local Functions
bool dfp_lnsrch( fq::vector<double> &z, double ll, fq::vector<double> &g,
		 fq::vector<double> &xi, fq::vector<double> &znew,
		 double &llnew, fq::vector<double> &gnew,
		 double stpmax, int &check, int &nfcn, int &ndfcn,
		 int &flagRund, dfp_optimizable *optable );
	// returns true unless optable->optimize(...) returns false



dfp_result dfp_optimize( dfp_optimizable *optable,
		       int maxIter, double convLL, double convGrad, double stepMax ) {
  /*
  double dsqrarg;
  double dmaxarg1,dmaxarg2;
  double dminarg1,dminarg2;
  int iminarg1,iminarg2;
  int imaxarg1,imaxarg2;
  */
  dfp_result result;
  result.ll = 0.0;
  result.iter = 1; // or 0?
  result.nfcn = 0;
  result.ndfcn = 0;
  result.err = 0;
  if ( ! optable->getStartingZ( result.z ) ) {
    // Message( "blah blah blah invalid starting set\n" );
    result.hessian.resize( result.z.size(), result.z.size() );
	for (int hi=result.z.size()-1; hi>=0; --hi)
		for (int hj=result.z.size()-1; hj>=0; --hj)
			result.hessian[hi][hj] = 1.0;
    result.err = -666; // invalid starting point
    return result;
  }
  result.hessian.resize( result.z.size(), result.z.size() );
  matrix<double> &hessian = result.hessian;
  
  int      n = result.z.size();
  int      check,i,its,j,flagRund;
  double   fac,fad,fae,stpmax,sumdg,sumxi,sum=0.0;
  double   tolf,tolg, tmpLL;
  
  fq::vector<double> g( n, 0.0 );
  fq::vector<double> dg( n, 0.0 );
  fq::vector<double> hdg( n, 0.0 );
  fq::vector<double> pnew( n, 0.0 );
  fq::vector<double> gnew( n, 0.0 );
  fq::vector<double> xi( n, 0.0 );
  
lab1:	
   if ( ! optable->evaluate( result.z, g, result.ll,
			     result.nfcn, result.ndfcn ) ) {
     result.err = -1;  // killed; unsuccessful evaluate()
     return result;
   }
   
   //cerr << "reset hessian etc" << endl;
   for (i=0; i<n; i++) {
     for (j=0; j<n; j++)
       hessian[i][j]=0.0;
     hessian[i][i]=1.0;
     xi[i] = -g[i];
     sum += result.z[i]*result.z[i];
   }
   stpmax=stepMax*DMAX(sqrt(sum),(double)n);
   
   //cerr << "iterate" << endl;
   for (its=1; its<=maxIter; its++, result.iter++) {
     if ((its-1)/ITPRT*ITPRT == its-1) {
       //cerr << "check" << endl;
       if ( ! optable->checkContinue( result.iter, result.nfcn,
				      result.ndfcn, result.ll,
				      result.z, g ) ) {
	 result.err = -2;  // killed; checkContinue == false
	 return result;
       }
     }
     
     //cerr << "lnsrch" << endl;
     tmpLL = result.ll;
     if ( ! dfp_lnsrch( result.z, tmpLL, g, xi, pnew, result.ll,
			gnew, stpmax, check, result.nfcn, result.ndfcn,
			flagRund, optable ) ) {
       result.err = -1; // opt.evaluate returned false
       return result;
     }
     /*
       if ( posieigencount > 5 )
	   { 
		   FREEALL;
		   Message(pvoid,"Positive eigenvalue,modify dead time and run again!");
		   return 0; // killed
	   }

     */
     if (flagRund==1) {
       // Message("Warning: Roundoff error and DFPMIN restarts.");
       //cerr << "Warning: Roundoff error and DFPMIN restarts." << endl;
       goto lab1;
     }
     
	 if ( tmpLL != 0.0 ) {
		tolf=fabs(result.ll - tmpLL)/fabs(tmpLL);
	 }
	 else {
		 tolf = 1.0;
	 }
     tmpLL = result.ll;
     for (i=0; i<n; i++) {
       xi[i] = pnew[i]-result.z[i];
       result.z[i] = pnew[i];
       dg[i] = g[i];
       g[i]  = gnew[i];
     }
     for (tolg=0.0,i=0; i<n; i++) 
       if (tolg<fabs(g[i])) 
	 tolg=fabs(g[i]);
     if (tolf<convLL && tolg<convGrad && fabs(tmpLL)<1.0e10) {
         result.err = 0;   // successful
	 	 return result;
     }
     
     for (i=0; i<n; i++)
       dg[i]=g[i]-dg[i];
     for (i=0; i<n; i++) {
       hdg[i]=0.0;
       for (j=0; j<n; j++)
	 hdg[i] += hessian[i][j]*dg[j];
     }
     fac=fae=sumdg=sumxi=0.0;
     for (i=0; i<n; i++) {
       fac += dg[i]*xi[i];
       fae += dg[i]*hdg[i];
       sumdg += DSQR(dg[i]);
       sumxi += DSQR(xi[i]);
     }
     if (fac*fac > EPS*sumdg*sumxi) {
       fac=1.0/fac;
       fad=1.0/fae;
       for (i=0; i<n; i++)
	 dg[i]=fac*xi[i]-fad*hdg[i];
       for (i=0; i<n; i++) 
	 for (j=0; j<n; j++) 
	   hessian[i][j] += fac*xi[i]*xi[j] -
	                    fad*hdg[i]*hdg[j] +
	                    fae*dg[i]*dg[j];
     }
     for (i=0; i<n; i++) {
       xi[i]=0.0;
       for (j=0; j<n; j++)
	 xi[i] -= hessian[i][j]*g[j];
     }
   }
   
   // Message("Warning: Iteration in dfpmin exceeds the limit(%d)",maxIter);
   //cerr << "Warning: Iteration in dfpmin exceeds the limit (" << maxIter
   // << ")" << endl;
   result.err = -3;      // exceeds the maximum number of iterations
   return result;
}

#undef ITMAX
#undef ITPRT
#undef TOLF
#undef TOLX
#undef TOLG
#undef EPS


#define	ALF  1.0e-4
#define	TOLX 1.0e-7

bool dfp_lnsrch( fq::vector<double> &z, double ll, fq::vector<double> &g,
		 fq::vector<double> &xxi, fq::vector<double> &znew,
		 double &llnew, fq::vector<double> &gnew,
		 double stpmax, int &check, int &nfcn, int &ndfcn,
		 int &flagRund, dfp_optimizable *optable )
{
  int 	   i; // flagGrad = 1
  int      n = z.size();
  double   a,alam,alam2,alamin,b,disc,f2,fold2,
           rhs1,rhs2,slope,sum,temp,test,tmplam;
  //double dmaxarg1,dmaxarg2;

  fq::vector<double> xi = xxi;
  
  check = flagRund = 0;
  
  for (sum=0.0, i=0; i<n; i++)
    sum += xi[i]*xi[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=0; i<n; i++)
      xi[i] *= stpmax/sum;
  
  for (slope=0.0, i=0; i<n; i++)
    slope += g[i]*xi[i];
  test=VERY_SMALL; // was 0.0; 9.25.01 Chris
  for (i=0; i<n; i++) {
    temp=fabs(xi[i])/DMAX(fabs(z[i]),1.0);
    if (temp > test)
      test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  
  for (;;) {
    for (i=0; i<n; i++)
      znew[i] = z[i] + alam*xi[i];
    // if ( posieigencount > 5 ) return;
    optable->evaluate( znew, gnew, llnew, nfcn, ndfcn );
    // (n,p,f,gnew,nf,ndf,flagGrad,pvoid);
    if (alam < alamin) {
      for (i=0; i<n; i++)
	znew[i]=z[i];
      check=1;
      return true;
    }
    else if (llnew <= ll+ALF*alam*slope)
      return true;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(llnew-ll-slope));
      else {
	rhs1 = llnew-ll-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0)
	  tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) {
	    flagRund=1;
	    return true;
	  }
	  else
	    tmplam=(-b+sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam) tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = llnew;
    fold2=ll;
    alam=DMAX(tmplam,0.1*alam);
  }
  // return true;
}

#undef ALF
#undef TOLX

class dfp_optable_dll : public dfp_optimizable
{
public:
	dfp_optable_dll( void *call_obj, int npar, double *initPar,
					 dfp_dll_func xfunc, dfp_dll_check xcheck )
		: caller( call_obj ), initZ( npar ), func( xfunc ), check( xcheck )
	{
		memcpy( &(initZ[0]), initPar, npar * sizeof(double) );
	}

  virtual bool getStartingZ( fq::vector<double> &z )
  {
	  z = initZ;
	  return true;
  }

  virtual bool evaluate( fq::vector<double> &z, fq::vector<double> &gz,
			   double &ll, int &nf, int &ndf, bool flagGrad=true )
  {
	  return ( ! func( caller, &(z[0]), &(gz[0]), &ll, &nf, &ndf, (int)flagGrad ) );
  }

  virtual bool checkContinue( int iter, int nf, int ndf, double ll,
			      fq::vector<double> &z, fq::vector<double> &gz )
  {
	  if ( check )
		return ( ! check( caller, iter, nf, ndf, ll, &(z[0]), &(gz[0]) ) );
	  else return true;
  }

private:
    void *caller;
	fq::vector<double> initZ;
	dfp_dll_func func;
	dfp_dll_check check;
};

extern "C" QUBOPT_API
double dfp_opt_dll( void *call_obj, int npar, double *initPar, double **hessianOut,
				    dfp_dll_func func, dfp_dll_check check,
					int maxIter, double convLL, double convGrad, double stepMax,
					int *err )
{
	dfp_optable_dll optable( call_obj, npar, initPar, func, check );
	dfp_result res = dfp_optimize( &optable, maxIter, convLL, convGrad, stepMax );

	*err = res.err;
	
	for ( int i=0; i<npar; ++i )
		for ( int j=0; j<npar; ++j )
			hessianOut[i][j] = res.hessian[i][j];
	
	return res.ll;
}
