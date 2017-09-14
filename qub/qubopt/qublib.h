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

// 10/2004 JB - combined all QubLib headers here;  Set up as precompiled header.

#include <float.h>	// For _clearfp and _fpreset
#include <math.h>
#include <stdio.h>
#include <stdlib.h> // for exit() in errexit
#include <string.h>
#include "qubopt.h"

//===== Process status 
	//#define  PROCESS_STOPPED    0
	//#define  PROCESS_STOPPING	1
	#define  PROCESS_STARTED	2
	//#define  PROCESS_STARTING	3

//===== Magic numbers 
	const double PI = 3.14159265358979;
	const double PI2 = 6.28318530717959;

//===== Formerly qubtypes.h
	typedef void *MODEL;

	typedef int  MAKEMODEL(MODEL model, int *segment, int *nchannel, int *nclass, int *nmetastate, int *nstate, int *npath, int *nar, int *state,  int **path, double *ratio, double *x, double *amp, double *xms, double **ar, double *pr);
	typedef int  MAKECONSTRAINTS(MODEL model, int **path, double *x, double **a, double *b, int *mc);
	typedef int  UPDATEMODEL(MODEL model, int *segment, int *nmetastate, double *amp, double *xms, double **ar, double *pr, double *x, double *dk, double *dd, double **a, int noupdatep);
	typedef int  UPDATEDISPLAY(int *firsttime);

	typedef void MSGSTATUS(char *msg);
	typedef void MSGOUTPUT(char *msg);
	typedef int STOPPROCESS(void);
	typedef int INPQUERY(char const *ACaption, char const *APrompt, double *Value, double low, double high); 

//===== MPL constants
	#define MAXSTATE	  100
	#define MAXPATH       50      
	#define MAXSEG        500
	#define MAXDRUG       5
	#define MAXAR         10   
	#define LOWLIMIT      1.0e-8
	#define UPLIMIT       1.0e8  
	#define MAXMUTICLASS  2000//20  
	//#define MAXMETASTATE  2000//1500
	//#define MAXMUTISTATE    30  
	//#define MAXMUTIPATH     20  

/*	------------------------------------------------------------------------- 
	convlv 
	Local functions : four1_0, twofft0, realft0
	-------------------------------------------------------------------------*/
	void convolve0(double *data, long n, double *respns, long m, int isign, double *ans);

/*	------------------------------------------------------------------------- 
	dfpmin 
	Local Functions : lnsrch
	-------------------------------------------------------------------------*/
	void dfpminim(int n, double *p, double *fret, double *var, int *iter, 
				void (*func)(int, double*, double*, double*, int*, MODEL, MAKEMODEL, UPDATEDISPLAY), 
				void (*oupt)(int, int, double*, double, double*, MODEL, UPDATEMODEL, MSGOUTPUT), 
				MODEL model, MAKEMODEL makemodel, UPDATEMODEL updatemodel, UPDATEDISPLAY updatedisplay, 
				MSGOUTPUT msgoutput, int *processflag);

/*	------------------------------------------------------------------------- 
	eigen 
	bypass or gsl
	-------------------------------------------------------------------------*/

	typedef int (*bypass_eigen_f)(int n, double *a, double *wr, double *wi, double *v, double *v_inv);
	extern "C" QUBOPT_API void set_bypass_eigen(bypass_eigen_f eigen_f);

	extern "C" QUBOPT_API void eigen(double **a, int n, double *wr, double *wi, double **v);

/*	------------------------------------------------------------------------- 
	func (MPL)
	-------------------------------------------------------------------------*/
	void dfpfunc(int nz, double *z, double *fz, double *gz, int *idef, MODEL model, MAKEMODEL makemodel, UPDATEDISPLAY updatedisplay);

/*	------------------------------------------------------------------------- 
	ioldt (MPL)
	i/o handlers for data, model and constraints 
	-------------------------------------------------------------------------*/
	int rdHeader(FILE *fp, int *nlead, double *dt, int *count);
	int rdSegHeader(FILE *fp, int *tstart, int *ndata);
	int rdSegData(FILE *fp, int count, int ndata, double *data);

/*	------------------------------------------------------------------------- 
	mbaum 
	-------------------------------------------------------------------------*/
        extern "C" QUBOPT_API void mbaum(int ndata, double *data, int nstate, double **a, int nmetastate, int **metastate, double *pr, int namp, double *amp, double *xms, int nar, double **ar, int nfir, double *fir, double *fret, int *state_out, MSGOUTPUT msgoutput, int* process);
	extern "C" QUBOPT_API void mgamma(int nar, int nfir, int nmetastate, int ndata, double **alpha, double **beta, double *scale);
		// mgamma s.b. in utility ? 
	extern "C" QUBOPT_API void mfward(int nstate, int nar, int nfir, int nmetastate, 
			int **metastate, int ndata, double *pr, double **a,
			double **b, double **alpha, double *scale);
	extern "C" QUBOPT_API void mbward(int nstate, int nar, int nfir, int nmetastate,
			int **metastate, int ndata, double **a, double **b,
			double *scale, double **beta);

/*	------------------------------------------------------------------------- 
	metamakv 
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API void  metamakv(int nstate, int namp, int **state, int nar, 
			int *nmetastate, int ***metastate);

/*	------------------------------------------------------------------------- 
	mexp 
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API void mexp(int n, double t, double **q, double **a, int nz, 
			double ***dq, double ***da, int i0);

/*	------------------------------------------------------------------------- 
	mutimakv 
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API void mutimakv(int nchannel, int nstate, int npath, 
			int *state, int **path, double *amp, 
			int *nmutistate, int *nmutipath, int *nmuticlass, 
			int ***mutistate, int ***mutipath);

/*	------------------------------------------------------------------------- 
	mvectb 
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API void mvectb(int ndata, double *data, double *amp, double *xms, 
			int nar, double **ar, int nfir, double *fir, 
			int nmetastate, int **metastate, double **b);

/*	------------------------------------------------------------------------- 
	qtomutiq
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API void qtomutiq(int nstate, int npath, int **path, int nmutistate, 
			int **mutistate, int nmutipath, int **mutipath, double **q, double **qq, int nz, double ***dq, double ***dqq);

/*	------------------------------------------------------------------------- 
	qublib : 
	-------------------------------------------------------------------------*/
	extern "C" QUBOPT_API int mpl(char **pmdlname, char **pldtname, char **pldtideal, 
			MODEL model, double voltage, double concentration, int segment, int currseg, 
			double frequency, int filtercoeffs, double datafreq, MAKEMODEL makemodel, 
			MAKECONSTRAINTS makeconstraints, UPDATEMODEL updatemodel, UPDATEDISPLAY updatedisplay, 
			MSGSTATUS msgstatus, MSGOUTPUT msgoutput, int *processflag, int errcode, double *lnlikelihood);
	extern "C" QUBOPT_API int amp(char **pmdlname, char **pldtname, char **pldtideal, 
			MODEL model, double voltage, double concentration, int segment, int currseg, 
			int noupdatep, double frequency, int filtercoeffs, double datafreq, MAKEMODEL makemodel, 
			UPDATEMODEL updatemodel, UPDATEDISPLAY updatedisplay, MSGSTATUS msgstatus, 
			MSGOUTPUT msgoutput, int *processflag, int *errcode, double *lnlikelihood);
	int filterintdata(int *source_data, int *dest_data, 
			int data_count, double dt, double freq/*Hz*/, int fft_count);
	int filterintdata_nfir(double freq, double dt);
	int filterdatad_nfir_notch(double transition_width, double dt);
	int filterdatad_notch(double *source_data, double *dest_data,
		    int data_count, double dt, double freq_lo, double freq_hi/*Hz*/, double transition_width, int fft_count);
	int filterdatad(double *data, double *fdata, int ndata, double dt, double freq);
	extern "C" QUBOPT_API double mlogl(int nar, int nfir, int ndata, double *scale);
	extern "C" QUBOPT_API void mdlogl(int nstate, int nar, int nfir, int nmetastate, 
			int **metastate, int ndata, double **alpha, double **beta, double **b, double **dll);

/*	------------------------------------------------------------------------- 
	svdecomp : Singular value decomposition for 0 based arrays 
	-------------------------------------------------------------------------*/
	void svdecomp0(double **a, int m, int n, double *w, double **v);
	void svdecomp0_mem(double **a, int m, int n, double *w, double **v, double *rv1/*n*/, double *aa/*m*/, double *vv/*n*/);
	void svdecomp0_x_mem(long double **a, int m, int n, long double *w, long double **v, long double *rv1/*n*/, long double *aa/*m*/, long double *vv/*n*/);

/*	------------------------------------------------------------------------- 
	utility : 
	Utility : Array allocation / free, some math funcs.
	-------------------------------------------------------------------------*/
	/* constants */
	#define SUCCESS 1
	#define FAILURE 0       
	#define TRUE    1
	#define FALSE   0   
         
	/* standard error handler */
	void errexit(char *msg);

	/* memory allocation */
	int *ialloc1(int size);
	double* dalloc1(int size);
	int **ialloc2(int row, int col);
	double** dalloc2(int row, int col);
	double*** dalloc3(int row, int col, int depth);
	long double ** xalloc2(int row, int col);
	void free2(char **ptr);
	void free3(char ***ptr);  

	/* math */
	#ifndef max
		#define max(a,b)    ((a)<(b) ? (b) : (a)) 
	#endif
	#ifndef min
		#define min(a,b)    ((a)>(b) ? (b) : (a))  
	#endif
	#define sqr(x)      ((x)*(x))
	#define dirac(a,b)  ((a)==(b) ? 1 : 0)      
	#define sign(a,b)   ((b)>=0 ? fabs(a) : -fabs(a))
	#define signl(a,b)   ((b)>=0 ? fabsl(a) : -fabsl(a))
	
	void swapn(double *a, double *b, int iItems);
	void swapnl(long double *a, long double *b, int iItems);
	void swap(double * p1, double *p2 );
	void swapl(long double * p1, long double *p2 );
	double pythag(double a, double b);
	long double pythagl(long double a, long double b);

/*	------------------------------------------------------------------------- 
	xtoq : Q matrix build
	Calculate Q matrix and dq matrix.
	-------------------------------------------------------------------------*/
	void xtoq(int nstate, int npath, int **path, double *C, double V, double *x, double **q, int nz, double **dx, double ***dq);

/*	------------------------------------------------------------------------- 
	ztox : (matrix multiply)
	Calculate Q matrix and dq matrix.
	-------------------------------------------------------------------------*/
	void ztox(int nz, int nx, double **a, double *b, double *z, double *x, double **dx);


