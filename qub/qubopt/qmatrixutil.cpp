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

#include <limits.h>

#include "qmatrixutil.h"
#include "qublib.h"

using namespace std;
using namespace fq;

//----- Local Functions 
void qabc(int na, int nb, int nc, int *ida, int *idb, int *idc, 
		  double **qab, double **qbb, double **qbc, double td, 
		  int **path, int npath, double **qac, double ***dqac_q);
void mil_mexp(double **q, double ***dq_x, int n, int m, double td, double **a, double ***da_x, int ie);
int nterm(double **q, int n, double t);	   


/*
 * Convention for multidimentional array indexing: (xtoq, qtomutiq, mexp ) 
 *    a[r][c], r=0..nrow-1, c=0..ncol-1  
 *            ~ a[r*ncol + c]
 *    a[d][r][c], d=0..ndepth-1, r=0..nrow-1, c=0..ncol-1    
 *            ~ a[d*nrow*ncol + r*ncol + c]
 */


//---------------------------------------------------------------
/* 
 * calculate the initial probabilities of the meta-states 
 */
void CalcMultiPr(int nchannel, int nstate, int nmetastate, int **metastate, 
				 double *pr, double* mprx){
	int i, j;
	double tmp, dum;
	
	for (i=0; i<nmetastate; i++){
		tmp = lnfactrl(nchannel);
		for (j=0; j<nstate; j++) {
			dum = (pr[j] < 1.0e-5) ? 1.0e-10 : pr[j];
			tmp += metastate[i][j]*log(dum) - lnfactrl(metastate[i][j]);
			}
		mprx[i] = exp(tmp);
		}
	}

void CalcMultiPr_mpl(int nchannel, int nstate, int nmutistate, int **mutistate,  
					 double *pr, double* mprx) {
	int i, j;
	double tmp, dum;
	
	for (i=0; i<nmutistate; i++) {
		tmp = lnfactrl(nchannel);
		for (j=0; j<nstate; j++) {
			dum = (pr[j] < 1.0e-5) ? 1.0e-10 : pr[j];
			tmp += mutistate[i][j+1]*log(dum) - lnfactrl(mutistate[i][j+1]);
			}
		mprx[i] = exp(tmp);
		}
	}

double lnfactrl(int n){
	int i;
	double tmp = 0;
	for (i=n; i>=1; i--)
		tmp += log(double(i));
	return tmp;
	}
	


//--------------------------------------------------------------------------------
// The CChannelCombine class iterates through all combinations of states to make 
// a multi channel model.  Ex( Ch=2, States=3 ) : [2 0 0] [1 1 0] [0 2 0] [1 0 1]
// [0 1 1] [0 0 2].  A c function equivalent is also in qublib nexcom0()
//--------------------------------------------------------------------------------
CChannelCombine :: CChannelCombine(int iChannels, int iStates ) { 
	m_aStates = new int[iStates+1];
	m_iChannels=iChannels;
	m_iStates=iStates;
	m_aStates[1]=iChannels;
	m_iTemp=iChannels;
	m_iStateIndex=0;
	for( int iState=2; iState<=iStates; ++iState)
		m_aStates[iState]=0;
	m_lContinue=true;
	}

bool CChannelCombine :: GetNextCombine( ) { 
	if( m_aStates[m_iStates] == m_iChannels )
		return false;
	if(m_iTemp > 1) 
		m_iStateIndex=0;
	m_iStateIndex++;
	m_iTemp = m_aStates[m_iStateIndex];
	m_aStates[m_iStateIndex] = 0;
	m_aStates[1] = m_iTemp - 1;
	m_aStates[m_iStateIndex + 1] = m_aStates[m_iStateIndex + 1] + 1;
	return true;
	}

int CChannelCombine :: State0(int i){
	return m_aStates[i+1];
	}

CChannelCombine :: ~CChannelCombine() { 
	delete [] m_aStates;
	}




/*
 * Multi-channel Markov model 
 */

void metaMarkov(int nchannel, int nstate, int npath, int **path,
				int *clazz, float *ratio,
				int *nmetastate, int *nmetapath, int *nmetaclass, float *metaamp,
				int **metastate, int **metapath, int *metaclass, int **subpath){
	
	int	i,j,k,l,m,n,s,t;

	if ( nstate <= 0 )
		return;

	// iterate through combinations of states / channels 
	*nmetastate=0;		
	CChannelCombine oMeta(nchannel, nstate);
	do {
		for ( j=0; j<nstate; ++j )
			metastate[*nmetastate][j] = oMeta.State0(j);
		++(*nmetastate);
		} while( oMeta.GetNextCombine() );

	
    // metastate class -- borrowed from SKM by Chris 3/19/4
	fq::vector<float> ampOfMetastate(*nmetastate);
	for ( n=0; n<*nmetastate; ++n )
		for ( j=0; j<nstate; ++j )
		  ampOfMetastate[n] += ((float) metastate[n][j]) * ratio[ clazz[j] ];
		
	double openDirection;
	try {
		openDirection = ratio[1];
		} 
	catch (...) {
		openDirection = 1.0;
		}
	
	eclass(ampOfMetastate, *nmetastate, *nmetaclass, metaamp, metaclass, openDirection);
	
	/* 
	number of metapaths, the start and end nodes of each path, and the start 
	and end nodes of individual channel transition involved in each metapath 
	*/
	*nmetapath=0;
	for (m=0; m<*nmetastate; m++) {
		for (n=0; n<*nmetastate; n++) {
			s=0;
			for (j=0; j<nstate; j++) {
				s=s+abs(metastate[m][j]-metastate[n][j]);
				}
			if(s==2) { // if you can go from metastate m to n by changing the state of one channel
				for (j=0; j<nstate; j++) {
					t=metastate[m][j]-metastate[n][j];
					if (t == 1) k=j;    // the transition is from state k to state l
					if (t ==-1) l=j;
					}
				for (i=0; i<npath; i++) {
					if (path[i][0] == k && path[i][1] == l) {  	// if this connection exists:
						metapath[*nmetapath][0]=m;
						metapath[*nmetapath][1]=n;
						subpath[*nmetapath][0]=k;
						subpath[*nmetapath][1]=l;
						(*nmetapath)++;
						}
					}
				}
			}
		}
	
	}
	    

/*
 * Constraint routine
 */
int freePar(int mc, int nx, double *x, double **Mtx, double *vct, int *nz, double *z, ostream &msgOut){
	const double TOL = 1.0e-5;
	const double MAD_SMALL = 1.0e-9;
	int     i,j,k,t; 
	double  **a,*b,**u,**v,*w,*h,wmax,s;
	
	/* no constraint */
	if (mc==0) {
		*nz=nx;
		for (i=0; i<nx; i++) {
			for (j=0; j<*nz; j++) 
				Mtx[i][j]=0.0;
			Mtx[i][i]=1.0;
			vct[i]=0.0;
			z[i]=x[i];
			}
		return TRUE;
		}
	
	/* initialization */
	a = double_alloc2D(mc+1,nx+1);
	b = double_alloc1D(nx+1);
	u = double_alloc2D(mc+1,nx+1);
	v = double_alloc2D(nx+1,nx+1);
	w = double_alloc1D(nx+1);
	h = double_alloc1D(mc);
	
	for (i=0; i<mc; i++) {
		b[i]=vct[i];
		for (j=0; j<nx; j++)
			a[i][j]=Mtx[i][j];
		}
	
	/* check if initial guess meets constraints */
	for (i=0; i<mc; i++) {
		s=0.0;
		for (j=0; j<nx; j++) 
			s=s+a[i][j]*x[j];
		if(fabs(s-b[i])>1.0e-2) {
			msgOut << "Fatal err: Initial guess doesn't meet constraints: " << (s-b[i]) << endl;
			return FALSE;
			}
		}
	
	/* SVD */
	for (i=0; i<mc; i++) 
		for (j=0; j<nx; j++) 
			u[i][j]=a[i][j];
	
	//-----------------------------------------------------------------------
	//----- SVD 
	svdecomp0(u,mc,nx,w,v);

	//----- Set wmax = max(w[j=0..nx]).  Our version of svdcmp
	// performs a sort at the end so this will always be w[0].
	wmax=fabs(w[0]);
	for( j=1; j<nx; j++)
		if(fabs(w[j])>wmax) 
			wmax=fabs(w[j]);

	t=0; 
	while(fabs(w[t])>wmax*TOL) 
		t++;
	
	/* confirm SVD results */   
	if ( t != mc ) {
		msgOut << "Too many constraints.  " << (mc - t) << " constraints may be linear combinations of the others." << endl;
		// msgOut << "Fatal err: Bad SVD in freepar" << endl;
		return FALSE;
		}

	for (i=t; i<nx; i++) {
		if(fabs(w[i])>wmax*TOL){ // w[t] changed to w[i] Chris 6/2/2
			msgOut << "Fatal err: Bad SVD in freepar" << endl;
			return FALSE;
			}
		}

	for (i=0; i<mc; i++) {
		for (j=mc; j<nx; j++) {
			if(fabs(u[i][j])>TOL) {
				msgOut << "Fatal err: Bad SVD in freepar" << endl;
				return FALSE;
				}
			}
		}

	for (i=0; i<mc; i++) { // make sure that u's columns are orthogonal unit vectors
		for (j=0; j<mc; j++) { // changed from nx to mc Chris 6/2/2
			s=0.0;
			for (k=0; k<mc; k++) 
				s=s+u[k][i]*u[k][j];
			if( (i!=j&&fabs(s)>TOL) || (i==j&&fabs(s-1.0)>TOL) ) {
				msgOut << "Fatal err: Bad SVD in freepar" << endl;
				return FALSE;
				}
			}
		}
	for (i=0; i<nx; i++) { // make sure that v's columns are orthogonal unit vectors
		for (j=0; j<nx; j++) {
			s=0.0;
			for (k=0; k<nx; k++) 
				s=s+v[k][i]*v[k][j];
			if( (i!=j && fabs(s)>TOL) || (i==j && fabs(s-1.0)>TOL) ){
				msgOut << "Fatal err: Bad SVD in freepar" << endl;
				return FALSE;
				}
			}
		}
	for (i=0; i<nx; i++) {
		for (j=0; j<nx; j++) {
			s=0.0;
			for (k=0; k<nx; k++) 
				s=s+v[i][k]*v[j][k];
			if( (i!=j&&fabs(s)>TOL) || (i==j&&fabs(s-1.0)>TOL) ){
				msgOut << "Fatal err: Bad SVD in freepar" << endl;
				return FALSE;
				}
			}
		}  
	
	/* set Mtx and vct */
	for (i=0; i<mc; i++) {
		h[i]=0.0;
		for (j=0; j<mc; j++) 
			h[i]=h[i]+u[j][i]*b[j];
		if (i < t) 
			h[i]=h[i]/w[i]; // i always less than t anyway?
		if (i >= t && fabs(h[i])>TOL){
			msgOut << "Fatal err: Bad SVD in freepar" << endl;
			return FALSE;
			}
		}
	
	for (i=0; i<nx; i++) {
		vct[i]=0.0;
		for (j=0; j<t; j++) 
			vct[i]=vct[i]+v[i][j]*h[j];
		}
	
	for (i=0; i<nx; i++) {
		for (j=t; j<nx; j++) {
			if ( fabs(v[i][j]) < MAD_SMALL ) // 5/31/2 Chris:  some weird svd was strewing 1.0e-18 when using loop balance
				Mtx[i][j-t] = 0.0;
			else
				Mtx[i][j-t]=v[i][j];
			}
		}
		
	*nz=nx-mc;
	
	/* get z from x */
	for (i=0; i<*nz; i++) {
		z[i]=0.0;
		for (j=0; j<nx; j++)
			z[i]=z[i]+Mtx[j][i]*x[j];
		}

	// TESTING: scale z to start at 1.0 (unless it's zero, stay at zero)
	// effectively, define diagonal matrix S with entries from z (or 1 instead of 0)
	// so x = Mtx * S * 1
	// Then (Mtx*S) can replace Mtx, and the params start as an all 1 vector (except any zeros remain)
	for (i=0; i<*nz; ++i) {
		if ( fabs(z[i]) > MAD_SMALL ) {
			for (j=0; j<nx; ++j)
				Mtx[j][i] *= z[i];
			z[i] = 1.0;
		}
	}
	
	free_2D( (char **)a);
	free(b);
	free_2D( (char **)v);
	free_2D( (char **)u);
	free(w);
	free(h);
	return TRUE;
	}


/*
 * Conversion routines
 */
int mil_ztox(int nz, int nx, double **Mtx, double *vct, double *z, double *x,
			  double **dx_z, double *xmin, double *xmax){
	int i,j;
	
	for (i=0; i<nx; i++) 
		x[i]=0.0;
	
	for (i=0; i<nx; i++) {
		for (j=0; j<nz; j++) 
			x[i] += Mtx[i][j]*z[j];
		x[i] += vct[i];
		}
	
	for (i=0; i<nx; i++) 
		for (j=0; j<nz; j++) 
			dx_z[i][j]=Mtx[i][j];
		
	for (j=0; j<nx; j++) 
		if( fabs(xmax[j]-xmin[j])>1.0e-5 && (x[j]<xmin[j] || x[j]>xmax[j]) )
			return FALSE;
	return TRUE;
	}

void qtometaq(QUBOPT_VAR_NOT_USED int nstate, int npath, int **path, int nmetastate,
			  int nmetapath, int **metapath, int **metastate, int **subpath,
			  double **q, double **metaq, double **dmetaq_q) {	
	int	i,j,k,l,m,n,ii,jj;
	
	for (i=0; i<nmetastate; i++) 
		for (j=0; j<nmetastate; j++) 
			metaq[i][j]=0.0;
		
	for (n=0; n<nmetapath; n++) {
		k=metapath[n][0];
		l=metapath[n][1];
		i=subpath[n][0];
		j=subpath[n][1];
		metaq[k][k] -= ( metaq[k][l] = metastate[k][i]*q[i][j] );

		for (m=0; m<npath; m++) {
			ii=path[m][0];
			jj=path[m][1];
			dmetaq_q[n][m]=dirac(i,ii)*dirac(j,jj)*metastate[k][i];
			}
		} 
	}

void aggregate(int nstate, int nclass, int *clazz, int *ngroup, int **index){
	int	ic,j,m,n;
	
	for (ic=0; ic<nclass; ic++) {
	    for (m=n=0,j=0; j<nstate; j++) 
	        if (clazz[j]==ic) {
	            index[ic][n]=j;
	            n++;
	            }
	        else {
	            index[ic][nstate+m]=j;
	            m++;
				}
	    ngroup[ic]=n;
		}
	}
	           	            
/*********************************************************************
**  	Decompose Q into Q=SDS1 with S1=S^-1 
**********************************************************************/
int qspctrm(int n, double **q, double* wr, double **s, double **s1,
			 int &posieigencount, ostream &msgOut){
	int     i,j;
	double  wmax;
	double  * wi = double_alloc1D(n+1);

	/* compute eigens of Q; return error message if no convergence */ 
	//if ( n > 1 )
		eigen( q,n,wr,wi,s);
	//else {
	//    if ( ! eigen( q,n,wr-1,wi-1,s,msgOut) ) {
	//	    free(wi);
	//		msgOut << "Err: No convergence in eigen calculations" << endl;
	//		return FALSE;
	//		}
	//}
	/* return if any complex eigenvalues found */
	wmax=0.0;
    for(i=0; i<n; i++) 
		if (fabs(wi[i])>wmax) 
			wmax=fabs(wi[i]);
	free(wi);
    if ( wmax > 1.0e-5 ) {
        msgOut << "Err: Complex eigenvalues found: " << wmax << endl;
        return FALSE;
		}
	
	/* return if any positive eigenvalues found */
	wmax=0.0;
    for (i=0; i<n; i++) 
		if (wr[i]>wmax) 
			wmax=wr[i];
    if ( wmax > 1.0e-5 ) {
        msgOut << "Err: Positive eigenvalues found: " << wmax << endl;
		posieigencount++;
		return FALSE;
		}
	    
	/* compute the inverse matrix of S */
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
			s1[i][j]=s[i][j];

	if ( ! gaussj0(s1,n,NULL,0,msgOut))
		return FALSE;

	return TRUE;
	}


void qtoqe(int nstate, int npath, int **path, int *clazz, int nclass,
		   double td, double **q, double **qe, double ***dqe_q){
	int     i,j,k,l,ic,jc,si,sj,na,nb,nc; 
	double  tmp;
	
	l=0;
	for (ic=0; ic<nclass; ic++) {
		for (na=0,i=0; i<nstate; i++) 
			if (clazz[i]==ic) na++;
		na=IMAX(na,nstate-na);
		l=IMAX(l,na);
		}

	l=nstate;
	int * ida    = int_alloc1D(l);
	int * idb    = int_alloc1D(l);
	int * idc    = int_alloc1D(l);
	double ** qac    = double_alloc2D(l,l);
	double ** qcc    = double_alloc2D(l,l);
	double ** qca    = double_alloc2D(l,l);
	double ** qcb    = double_alloc2D(l,l);
	double ** qab    = double_alloc2D(l,l);
	double ** qbb    = double_alloc2D(l,l);
	double ** u      = double_alloc2D(l,l);
	double ** v      = double_alloc2D(l,l);
	double ** e      = double_alloc2D(l,l);
	double *** de_q   = double_alloc3D(l,l,npath);
	double *** dqbb_q = double_alloc3D(l,l,npath);
	double *** dqab_q = double_alloc3D(l,l,npath);
   
	/* correct diagonal blocks */
	
	for (ic=0; ic<nclass; ic++) {
		for (na=nc=0,i=0; i<nstate; i++){
			if (clazz[i]==ic)
				ida[na++]=i; 
			else
				idc[nc++]=i; 
			}

		for (i=0; i<na; i++) 
			for (j=0; j<nc; j++) 
				qac[i][j]=q[ida[i]][idc[j]];
		for (i=0; i<nc; i++) 
			for (j=0; j<nc; j++) 
				qcc[i][j]=q[idc[i]][idc[j]];
		for (i=0; i<nc; i++) 
			for (j=0; j<na; j++) 
				qca[i][j]=q[idc[i]][ida[j]];
					
		qabc(na,nc,na,ida,idc,ida,qac,qcc,qca,td,path,npath,e,de_q);
					
		for (i=0; i<na; i++) 
			for (j=0; j<na; j++) {
				si=ida[i];
				sj=ida[j];
				qe[si][sj]=q[si][sj]+e[i][j];
				for (k=0; k<npath; k++) {
					if( si==path[k][0] && sj==path[k][1] ) 
						tmp = 1.0;
					else if( si==path[k][0] && sj==si )    
						tmp =-1.0;
					else 				
						tmp = 0.0;
					dqe_q[si][sj][k]=tmp+de_q[i][j][k];
					}
				}
		}
	
	
	/* correct off-diagonal blocks */
	for (ic=0; ic<nclass; ic++) 
        for (jc=0; jc<nclass; jc++) 
			if (ic!=jc) {
				for (na=0,i=0; i<nstate; i++)
					if (clazz[i]==ic) {
						ida[na]=i; 
						na++;
						}
				for (nb=0,i=0; i<nstate; i++)
					if (clazz[i]==jc) {
						idb[nb]=i; 
						nb++;
						}   
				for (nc=0,i=0; i<nstate; i++)
					if (clazz[i]!=ic&&clazz[i]!=jc) {
						idc[nc]=i;
						nc++;
						}
				for (i=0; i<na; i++) 
					for (j=0; j<nc; j++) 
						qac[i][j]=q[ida[i]][idc[j]];
				for (i=0; i<nc; i++) 
					for (j=0; j<nc; j++) 
						qcc[i][j]=q[idc[i]][idc[j]];
				for (i=0; i<nc; i++) 
					for (j=0; j<nb; j++) 
						qcb[i][j]=q[idc[i]][idb[j]];
								
				qabc(na,nc,nb,ida,idc,idb,qac,qcc,qcb,td,path,npath,e,de_q);
				
				for (i=0; i<na; i++) 
					for (j=0; j<nb; j++) {
						si=ida[i];
						sj=idb[j];
						qab[i][j]=q[si][sj]+e[i][j];
						for (k=0; k<npath; k++) {
							if( si==path[k][0] && sj==path[k][1] ) 
								tmp= 1.0;
							else 				
								tmp= 0.0;
							dqab_q[i][j][k]=tmp+de_q[i][j][k];
							}
						}
					
				for (i=0; i<nb; i++) 
					for (j=0; j<nb; j++) {
						si=idb[i];
						sj=idb[j];
						qbb[i][j]=q[si][sj];
						for (k=0; k<npath; k++) { 
							if( si==path[k][0] && sj==path[k][1] ) 
								tmp= 1;
							else if( si==path[k][0] && sj==si )    
								tmp=-1;
							else 				
								tmp= 0;
							dqbb_q[i][j][k]=tmp;
							}
						}    

				mil_mexp(qbb,dqbb_q,nb,npath,td,e,de_q,0);
				mxm(na,nb,nb,qab,e,u);
				
				for (i=0; i<na; i++) 
					for (j=0; j<nb; j++) 
						qe[ida[i]][idb[j]]=u[i][j];
						
				for (k=0; k<npath; k++) {
					for (i=0; i<na; i++) 
						for (j=0; j<nb; j++) 
							u[i][j]=dqab_q[i][j][k];

					mxm(na,nb,nb,u,e,v);
					for (i=0; i<na; i++) 
						for (j=0; j<nb; j++) 
							dqe_q[ida[i]][idb[j]][k] =v[i][j];
							
					for (i=0; i<nb; i++) 
						for (j=0; j<nb; j++) 
							u[i][j]=de_q[i][j][k];
					mxm(na,nb,nb,qab,u,v);
					for (i=0; i<na; i++) 
						for (j=0; j<nb; j++) 
							dqe_q[ida[i]][idb[j]][k]+=v[i][j];
					}
				}
			
	free(ida);
	free(idb);
	free(idc);
	free_2D( (char **)qac);
	free_2D( (char **)qcc);
	free_2D( (char **)qca);
	free_2D( (char **)qcb);
	free_2D( (char **)qab);
	free_2D( (char **)qbb);
	free_2D( (char **)u);
	free_2D( (char **)v);
	free_2D( (char **)e);
	free_3D( (char ***)de_q);
	free_3D( (char ***)dqbb_q);
	free_3D( (char ***)dqab_q);
	}	        

void qabc(int na, int nb, int nc, int *ida, int *idb, int *idc, 
		  double **qab, double **qbb, double **qbc, double td, 
		  int **path, int npath, double **qac, double ***dqac_q) {
	int	i,j,k,l,si,sj;
	double tmp;
	
	if (nb==0) {
		for (i=0; i<na; i++)
			for (j=0; j<nc; j++) {
				qac[i][j]=0.0;
				for (k=0; k<npath; k++) 
					dqac_q[i][j][k]=0.0;
				}
		return;
		}

	l=max(max(na,nb),nc)+1;

	double ** a      = double_alloc2D(l,l);
	double ** u      = double_alloc2D(l,l);
	double ** v      = double_alloc2D(l,l);
	double *** dqbb_q = double_alloc3D(l,l,npath);
	double *** da_q   = double_alloc3D(l,l,npath);
	    
	for (k=0; k<npath; k++)  
		for (i=0; i<nb; i++) 
			for (j=0; j<nb; j++) {
				si=idb[i];
				sj=idb[j];
				if( si==path[k][0] && sj==path[k][1] ) 
					tmp = 1;
				else if( si==path[k][0] && sj==si )    
					tmp = -1;
				else 				
					tmp = 0;
				dqbb_q[i][j][k]=tmp;
				}    

    mil_mexp(qbb,dqbb_q,nb,npath,td,a,da_q,1);

	mxm(na,nb,nb,qab,a,u);
	mxm(na,nb,nc,u,qbc,qac);
	
	for (k=0; k<npath; k++) {
		for (i=0; i<na; i++) 
			for (j=0; j<nb; j++) 
				v[i][j] = ( ida[i]==path[k][0] && idb[j]==path[k][1] ) ? 1 : 0;

		mxm(na,nb,nb,v,a,u);
		mxm(na,nb,nc,u,qbc,v);

		for (i=0; i<na; i++) 
			for (j=0; j<nc; j++) 
				dqac_q[i][j][k]  = v[i][j];
			
		for (i=0; i<nb; i++) 
			for (j=0; j<nb; j++) 
				v[i][j]=da_q[i][j][k];

		mxm(na,nb,nb,qab,v,u);
		mxm(na,nb,nc,u,qbc,v);

		for (i=0; i<na; i++) 
			for (j=0; j<nc; j++) 
				dqac_q[i][j][k] += v[i][j];
			
		for (i=0; i<nb; i++) 
			for (j=0; j<nc; j++) 
				v[i][j] = ( idb[i]==path[k][0] && idc[j]==path[k][1] ) ? 1 : 0;

		mxm(nb,nb,nc,a,v,u);
		mxm(na,nb,nc,qab,u,v);

		for( i=0; i<na; i++ ) 
			for( j=0; j<nc; j++) 
				dqac_q[i][j][k] += v[i][j];
		}
	
	free_2D( (char **)a);
	free_2D( (char **)u);
	free_2D( (char **)v);
	free_3D( (char ***)dqbb_q);
	free_3D( (char ***)da_q);
	}
	    
void mil_mexp(double **q, double ***dq_x, int n, int m, double td, 
			  double **a, double ***da_x, int ie) {
	int	i,j,k,l,it,itmax;
	double s;
	double ** u   = double_alloc2D(n+1,n+1);
	double ** u1  = double_alloc2D(n+1,n+1);
	double *** v   = double_alloc3D(n+1,n+1,m);
	double *** v1  = double_alloc3D(n+1,n+1,m);
	
	for (i=0; i<n; i++)    
		for (j=0; j<n; j++) {
			a[i][j] = 0.0;
			u[i][j] = dirac(i,j)*pow(td,ie);
			for( k=0; k<m; k++ ) 
				da_x[i][j][k]=v[i][j][k]=0.0;
			}
	
	itmax = nterm(q,n,td);	
       
    for (it=0; it <= itmax; it++) {
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
				a[i][j] += u[i][j];
				for (k=0; k<m; k++) da_x[i][j][k] += v[i][j][k];
				}
			
		mxm(n,n,n,u,q,u1);
			
		for (k=0; k<m; k++)  
			for (i=0; i<n; i++)         
				for (j=0; j<n; j++) {
					s=0.0;
					for( l=0; l<n; l++) 
						s += q[i][l]*v[l][j][k] + dq_x[i][l][k]*u[l][j];
					v1[i][j][k]=s;
					}
				
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
				u[i][j] = u1[i][j]*td/(it+1+ie);
				for (k=0; k<m; k++) 
					v[i][j][k]=v1[i][j][k]*td/(it+1+ie);
				}
		}
	
	free_2D( (char **)u);
	free_2D( (char **)u1);
	free_3D( (char ***)v);
	free_3D( (char ***)v1);   
	}
 
//todo - RETURN AT LEAST 50 - WHY ?   s ERROR ON SUM(Q) > 1 ? 
int nterm(double **q, int n, double t) {	   
	const double EPS= 1.0E-50;
	int	i,j,m;
	double r=1.0,s=0.0;
		   
	for (i=0; i<n; i++) 
		for (j=0; j<n; j++) 
			s += t*fabs(q[i][j]);
	
	for( m=0; r/(1.0-r/(m+2))>=EPS; m++) // s=5 then r=5 then [expr]~= -7.5 : perhaps abs() would help.
		r*=s/(m+1);
	
	return max(50,m);
	}

/* -----------------------------------------------------------------------
	Calculate A = I*t^i0 + Qt^(i0+1)/(i0+1) + Q^2*t^(i0+2)/(i0+2) + ... 
	matrix exponential and its derivatives.
*/
void mexp(int n, float t, float* q, float* a, int nz, float* dq, float* da, int i0) {
	// nz could be zero, dq and da could be NULL
	int i,j,k,l;
	float s;
	float qsum=0.0;				// Sum of abs(t*sum(q)) = magnitude of Q
   
	// initialization
	float * u  = new float[n*n];
	float * v  = new float[n*n*nz];
	float * u1 = new float[n*n];
	float * v1 = new float[n*n*nz];
	
	for (i=0; i<n; i++)    
		for (j=0; j<n; j++) {
			a[i*n+j] = 0.0;
			u[i*n+j] = float((i==j)?pow(t,i0):0.0);
			for (k=0; k<nz; k++){
				da[i*n*nz+j*nz+k] = 0.0;
				v[i*n*nz+j*nz+k] = 0.0;
				}
			qsum+=float(fabs(q[i*n+j]));
			}
	qsum*=t;
		
	//----- series summation for minimum of 4 <iterm> terms until convergence 
	float r=1.0;
	for( int iterm=0; iterm<4 || r/(1.0-r/float(iterm+2))>=float(1.0e-10); iterm++) {
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
				a[i*n+j] += u[i*n+j];                        // a[i][j]
				for (k=0; k<nz; k++) 
					da[i*n*nz+j*nz+k] += v[i*n*nz+j*nz+k];    // da[i][j][k]
				}
		for (i=0; i<n; i++)
			for (j=0; j<n; j++) { 
				for (s=0.,l=0; l<n; l++) 
					s += u[i*n+l]*q[l*n+j];                   // u[i][l], q[l][j]
				u1[i*n+j]=s; 
				} 
		for (k=0; k<nz; k++)  
			for (i=0; i<n; i++)         
				for (j=0; j<n; j++) {
					for (s=0.,l=0; l<n; l++) 
						s += q[i*n+l]*v[l*n*nz+j*nz+k] + dq[i*n*nz+l*nz+k]*u[l*n+j];
					v1[i*n*nz+j*nz+k] = s;                    // v1[i][j][k]
					}
		for (i=0; i<n; i++)         
			for (j=0; j<n; j++) {
			  u[i*n+j] = u1[i*n+j]*t/(float(iterm+1+i0));
				for (k=0; k<nz; k++) 
				  v[i*n*nz+j*nz+k] = v1[i*n*nz+j*nz+k]*t/(float(iterm+1+i0)); // v[i][j][k]
				}
		
		// Convergence calculation pt. 2
		r *= qsum/float(iterm+1);
		}
	
	// free memory
	delete [] u;
	delete [] v;
	delete [] u1;
	delete [] v1;
	}

//-----------------------------------------------------------------------
float gaussian(float mean, float xms, float x){
	//const float DELT = float(0.0001);
	const float EPS = float(0.00001);
	if ( xms == 0 ) {
		if ( fabs(x-mean) < EPS )
			return 1.0 * SHRT_MAX; // warning this is a sub for infinity
		else
			return 0.0;
	}
	float tmp = (x-mean)/xms;
	tmp = float(exp(-tmp*tmp/2.0));  // this is where it would underflow
	tmp/= float(2.5066283*xms);
	return (tmp<float(1.0e-8)) ? float(1.0e-8) : tmp;
	}

//-----------------------------------------------------------------------
void qtomutiq(int nstate, QUBOPT_VAR_NOT_USED int npath, QUBOPT_VAR_NOT_USED int* path, int nmutistate, 
			  int* mutistate, int nmutipath, int* mutipath, float* q, 
			  float* mutiq, int nz, float* dq, float* dmutiq){   
	// nz could be zero 
	int  I,J,i,j,k,m;
	float  s;                
    
	// rate constants 
	for (I=0; I<nmutistate; I++) 
		for (J=0; J<nmutistate; J++) 
			mutiq[I*nmutistate+J] = 0.0;
		
	for (m=0; m<nmutipath; m++) {
		I = mutipath[m*4+0];
		J = mutipath[m*4+1];
		i = mutipath[m*4+2];
		j = mutipath[m*4+3];
		mutiq[I*nmutistate+J] = ((float) mutistate[I*(nstate+1)+i+1])*q[i*nstate+j];
		}
	
	for (I=0; I<nmutistate; I++) {
		for (s=0.,J=0; J<nmutistate; J++) 
			if (J != I) 
				s -= mutiq[I*nmutistate+J];
			mutiq[I*nmutistate+I] = s;
		}
    
	// derivatives to the z-variables 
	for (I=0; I<nmutistate; I++) 
		for (J=0; J<nmutistate; J++) 
			for (k=0; k<nz; k++) 
				dmutiq[I*nmutistate*nz+J*nz+k] = 0.0;
			
	for (m=0; m<nmutipath; m++) {
		I = mutipath[m*4+0];
		J = mutipath[m*4+1];
		i = mutipath[m*4+2];
		j = mutipath[m*4+3];
		for (k=0; k<nz; k++) 
			dmutiq[I*nmutistate*nz+J*nz+k] =                // dmutiq[I][J][k]
			  ((float) mutistate[I*(nstate+1)+i+1])*dq[i*nstate*nz+j*nz+k];  
		} 
	
	for (k=0; k<nz; k++) {
		for (I=0; I<nmutistate; I++) {
			for (s=0.,J=0; J<nmutistate; J++) 
				if (J != I) 
					s -= dmutiq[I*nmutistate*nz+J*nz+k];      // dmutiq[I][J][k]
				dmutiq[I*nmutistate*nz+I*nz+k]=s;               // dmutiq[I][I][k]
			}
		}
	}

//-----------------------------------------------------------------------
// q-matrix utility  - Fill in q and dq arrays.
// Note : drug[] could be NULL and nz could be zero
void xtoq(int nstate, int npath, int* path, float* drug, float volt, 
		  float* x, float* q, int nz, float* dx, float* dq) {
	int  i,j,k,m;
	float nQ;
	double c;
	
	//-- Zero the q and dq arrays 
	for (i=0; i<nstate; i++) {
		for (j=0; j<nstate; j++) {
			q[i*nstate+j] = 0.0;						// q[i][j]
			for (k=0; k<nz; k++) 
				dq[i*nstate*nz+j*nz+k] = 0.0;			// dq[i][j][k]         
			}
		}
		
	//-- Fill in rates and diagonal = -sum(rates) 
	for (m=0; m<npath; m++) {
		i = path[m*4+0];								// From State
		j = path[m*4+1];								// To State
		c = ((drug==NULL) ? 1.0 : drug[path[m*4+2]]);	// Drug index 
		// Fill in Q 
		nQ = float(c * exp(x[2*m]+volt*x[2*m+1]));
		q[i*nstate+j] = nQ;
		q[i*nstate+i] -= nQ;
		
		// Fill in dQ - derivatives to the z-variables  
		for (k=0; k<nz; k++)
			dq[i*nstate*nz+i*nz+k] -= 
				( dq[i*nstate*nz+j*nz+k] =  nQ * ( dx[2*m*nz+k] + volt*dx[(2*m+1)*nz+k] ));  
		}
	}

//-----------------------------------------------------------------------
// Slightly Similar to a numerical recipes algorithm eclazz, but this 
// returns a sorted aClasses array (though it is less efficient).
// aData - array of values.  iData - number of values.
// iClasses - output number of classes.  aClasses - output class values.
// ia[0..iData-1] - output  index to aClasses for each aData[0..iData-1]
// This function iteratively finds the minimum value in the list, then 
// considers all values within 1.0e-5 of that to be equivalent class and 
// removes them from the list until all data have been placed in classes.
void eclass(float * aData, int iData, int & iClasses, float * aClasses, 
			int * ia, double openDirection) {
	int  i,k;
	float amin;

	int iCopyPos;				// While scanning aClasses array, position being copied to ( = current - # of elements at current class already processed )
	bool lCompare = ( openDirection > 0 );	// whether to find min or max ( ascending or descending )

	int * aRemain = new int[iData];	// Array of indices of items to be classified
	for (i=0; i<iData; i++) 
		aRemain[i] = i;

	iClasses = 0;
	while( iData > 0 ) {
		// Find the minimum (or maximum ) value in the remaining aRemain list 
		k = aRemain[0];
		amin = aData[k];
		
		for( i=1; i<iData; i++ ) {
			k=aRemain[i];
			if( lCompare == (aData[k]<amin)) 
				amin=aData[k];
			}
		
		// add amin to aClasses[] array of classes and all aData[k]~= amin set ia[k]=iClasses
		aClasses[iClasses]=amin;
		iCopyPos=0;
		for (i=0; i<iData; i++) {
			k = aRemain[i]; 
			if (fabs(aData[k]-amin)<1.0e-5) 
				ia[k]=iClasses;					// include as result for this minimum ...
			else 
				aRemain[iCopyPos++]=k;				// ... or copy back to array for next iteration
			}
		iData=iCopyPos;

		iClasses++;
		}   
	
	delete [] aRemain;
	}


