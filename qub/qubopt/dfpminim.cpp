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

#include "qublib.h"

//----- Local Functions 
void lnsearch0(int n, double *xold, double fold, double *g, double *p, 
    double *x, double *gnew, double *f, int *icase, 
    void (*func)(int, double*, double*, double*, int*, MODEL, MAKEMODEL, UPDATEDISPLAY), 
	MODEL, MAKEMODEL, UPDATEDISPLAY);

//------ 
void dfpminim(int n, double *p, double *fret, double *var, int *iter,
			void (*func)(int, double*, double*, double*, int*, MODEL model, MAKEMODEL, UPDATEDISPLAY),
			void (*oupt)(int, int, double*, double, double*, MODEL model, UPDATEMODEL, MSGOUTPUT),
			MODEL model, MAKEMODEL makemodel, UPDATEMODEL updatemodel, 
			UPDATEDISPLAY updatedisplay, MSGOUTPUT msgoutput, int *processflag) {
	
	const double EPS=3.0e-8;
	const double TOLG=5.0e-3;
	const double TOLF=1.0e-2;
	const int IT_PRT=1;
	const int IT_MAX=100;
    int		icase, idef, i, its, j;
    double  fac, fad, fae, fp, sumdg, sumxi;	//, test, tmp
    double  *dg, *g, *hdg, **hessin, *pnew, *gnew, *xi, testg;
    
    dg     = dalloc1(n);
    g      = dalloc1(n);
    hdg    = dalloc1(n); 
    hessin = dalloc2(n,n);
    pnew   = dalloc1(n);  
    gnew   = dalloc1(n);
    xi     = dalloc1(n);
    
    //First call - calculate starting function value and gradient
	(*func)(n, p, &fp, g, &idef, model, makemodel, updatedisplay);
	
    for( i=0; i<n; i++ ) {
		for( j=0; j<n; j++)
			hessin[i][j] = 0.0;
		hessin[i][i] = 1.0;
		xi[i] =- g[i];
		}
	
	for( its=1; its<=IT_MAX; its++) {			// main loop over the iterations
		if (*processflag != PROCESS_STARTED)
			break;
		
		*iter = its;
		
		for( i=0; i<n; i++)
			var[i] = hessin[i][i];
		
		// always TRUE  since IT_PRT=1.  Perhaps this is iterations per output, 
		// so equivalent to : (its-1)%IT_PRT==0 
		if( (its-1)/IT_PRT*IT_PRT == its-1)		
			(*oupt)(its, n, p, fp, g, model, updatemodel, msgoutput);
		
		lnsearch0(n, p, fp, g, xi, pnew, gnew, fret, &icase, func, model, makemodel, updatedisplay);
		
		// test for convergence on zero gradient
		testg=0.0;
		for( i=0; i<n; i++ )  
			testg = max( fabs(g[i]) , testg) ;
		if ( fabs(*fret-fp)<TOLF && testg < TOLG) {
			*fret = fp;
			break;
			}
		
		if (*processflag != PROCESS_STARTED) {
			break;
			}
		
		if (icase == 1)
			break;

		if (icase == 2) {
			for( i=0; i<n; i++ ) {
				for( j=0; j<n; j++)
					hessin[i][j] = 0.0;
				hessin[i][i] = 1.0;
				}
			} 
		else {
			for( i=0; i<n; i++) {
				xi[i] = pnew[i] - p[i];
				dg[i] = gnew[i] - g[i];
				}
			
			for( i=0; i<n; i++) {
				hdg[i] = 0.0;
				for( j=0; j<n; j++)
					hdg[i] += hessin[i][j]*dg[j];
				}
			
			fac = fae = sumdg = sumxi = 0.0;	// calculate dot products for the denominators
			for( i=0; i<n; i++) {
				fac   += dg[i]*xi[i];
				fae   += dg[i]*hdg[i];
				sumdg += sqr(dg[i]);
				sumxi += sqr(xi[i]);
				}
			
			if (fac*fac > EPS*sumdg*sumxi) {	// skip update if fac not sufficiently positive
				fac = 1.0/fac;
				fad = 1.0/fae;
				for( i=0; i<n; i++)				// The vector that makes BFGS different from DFP
					dg[i] = fac*xi[i] - fad*hdg[i];
				for( i=0; i<n; i++) {
					for( j=0; j<n; j++)
						hessin[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					}
				}
			}
		
		for( i=0; i<n; i++) {		// calculate the next direction to go
			xi[i]= 0.0;
			for( j=0; j<n; j++)
				xi[i] -= hessin[i][j]*gnew[j];
			p[i] = pnew[i];
			g[i] = gnew[i];
			}
		fp = *fret;
		}
	
	if( its>=IT_MAX ) 
		msgoutput("\r\tToo many iterations in dfpminim\r\n");

	free(dg); 
	free(g); 
	free(hdg); 
	free(pnew); 
	free(gnew); 
	free(xi); 
	free2((char**)hessin);   
	}

//----------------------------------------------------------------------------------
// given an ndimemsional point xold[1..n], the value and gradient there : fold, g[]
// and a direction p[], finds a new point x[1..n] along direction p.
void lnsearch0(int n, double *xold, double fold, double *g, double *p, 
			 double *x, double *gnew, double *f, int *icase, 
			 void (*func)(int, double*, double*, double*, int*, 
			 MODEL, MAKEMODEL, UPDATEDISPLAY), 
			 MODEL model, MAKEMODEL makemodel, UPDATEDISPLAY updatedisplay){
	const double ALF = 1.0e-4;	// ensures sufficient decrease in function value 
	const double TOLX = 1.0e-7;	// convergence criteion for delta x
	
    int		i, idef, ilam;
    double	a, alam, alam2, alamin, b, disc, f2, fold2, tmp;
    double	rhs1, rhs2, slope, test, tmplam;
	
	slope=0.0;
    for( i=0; i<n; i++) 
		slope += g[i]*p[i];
	
	test=0.0;					// compute lamda.min 
    for( i=0; i<n; i++) {
		tmp = fabs(p[i]) / max(fabs(xold[i]), 1.0);
		test = max( test, tmp );
		}
	
    alamin = TOLX / test;
    alam = 1.0;   
    ilam = 0;
    for (;;) {
		for( i=0; i<n; i++) 
			x[i] = xold[i] + alam*p[i];
		(*func)(n, x, f, gnew, &idef, model, makemodel, updatedisplay);
		if (idef == 0) 
			alam *= 0.68;
		else {
			if (alam < alamin) {		// convergence on delta x.  For zero finding,
				for (i = 0; i < n; i++) // the caller should verify the convergence.
					x[i] = xold[i];
				*icase = 1;
				return;
				} 

			if (*f <= fold + ALF*alam*slope) {	// sufficient function decrease.  Backtrack.
				*icase = 0;
				return;
				}

			if (ilam == 0) {			// first time
				tmplam = -slope / (2.0*(*f - fold - slope));
				ilam++;
				} 
			else {						// subsequent backtracks
				rhs1 = *f - fold - alam*slope;
				rhs2 = f2 - fold2 - alam2*slope;
				a = (rhs1 / (alam*alam) - rhs2 / (alam2*alam2)) / (alam - alam2);
				b = (-alam2*rhs1 / (alam*alam) + alam*rhs2 / (alam2*alam2)) / (alam - alam2);
				if (a == 0.0) 
					tmplam = -slope / (2.0*b);
				else {
					disc = b*b - 3.0*a*slope;
					if (disc < 0.0) {
						*icase = 2;
						return;
						} 
					tmplam = (- b + sqrt(disc)) / (3.0*a);
					}
				if (tmplam > 0.5*alam) 
					tmplam = 0.5*alam;
				}
			alam2 = alam;
			f2 = *f;
			fold2 = fold;
			alam = max(tmplam, 0.1*alam);
			}
		}
	}  
