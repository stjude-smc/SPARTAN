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


//===================== AMP and MPL 

#include "qublib.h"

//----- LOCAL FUNCTIONS 
void ovrlap(double *data, long ndata, double *respns, long m, double *ans);
void ovrlapint(int *source_data, int *dest_data, int data_count, double *response, int response_count, int fft_count);
double *gsfir(double fc, int *nfir);
void wrap(int nh, double *h);   
void freepar0(int m, int n, double *x, double **a, double *b, int *nz, double *z);
void dfp_mpl_func(int nz, double *z, double *fz, double *gz, int *idef, 
			 MODEL model, MAKEMODEL makemodel, UPDATEDISPLAY updatedisplay);
void oupamp(int nclass, int nmutistate, double **qq, int namp, double *amp, double *xms, int nar, double **ar, 
			int nmetastate, double *pr, MSGOUTPUT msgoutput);
void oupdfp(int it, int nz, double *z, double fz, double *g, MODEL model, UPDATEMODEL updatemodel, MSGOUTPUT msgoutput);
void oupmdl(char *mdlname, int nz, double *z, double *vz, double fz, int it, double tt, MODEL model, UPDATEMODEL updatemodel, MSGOUTPUT msgoutput);

//----- Global variables 
int				mpl_nchannel, mpl_nclass, mpl_nstate, mpl_npath,
				*mpl_state = NULL, **mpl_path = NULL,
				mpl_nmutistate, mpl_nmuticlass, mpl_nmutipath,
				**mpl_mutistate = NULL, **mpl_mutipath = NULL, mpl_nfir,
				mpl_nseg, mpl_ndat[MAXSEG], mpl_nar,
				mpl_nfunc = 0, mpl_ndfunc = 0;

double			**mpl_xlimit = NULL, **mpl_Mtx = NULL, *mpl_vct = NULL, *mpl_ratio = NULL, *mpl_fir = NULL, mpl_fc,
				mpl_dt[MAXSEG], mpl_V[MAXSEG], *mpl_C[MAXSEG], *mpl_dat[MAXSEG],
				*mpl_pr[MAXSEG], *mpl_amp[MAXSEG], *mpl_xms[MAXSEG], **mpl_ar[MAXSEG];

char			mpl_ldtideal[512];
int				mpl_count;
int				mpl_tstart[MAXSEG];
int				mpl_leader;
int				mpl_firsttime;

//---- Error Codes
//#define  SIMULATION_ERROR	1
//#define  FILTER_ERROR		2
//#define  CANT_WRITE_FILE	3
#define  CANT_READ_FILE		4
#define  OUT_OF_MEMORY		5
#define  DONT_DISPLAY		6  
//#define  RESAMPLE_ERROR		7
//#define  USER_BREAK			8
#define  LOOP_UNBALANCE		9
#define  LOOP_ADJUSTED		10
#define  ERROR_MBAUM		11

//--------------------------------------------------------------------------
//-----All the memory should be freed, including the global blocks !!!
//--------------------------------------------------------------------------
#define AMP_FREE_ALL	free(state); free2((char**)path); free(ratio); free(x);	free(fir);\
						free2((char**)a); free2((char**)q); free2((char**)qq); free2((char**)dx);\
						free3((char***)dq); free3((char***)dqq); free3((char***)da);\
						free(data); free(dataf);\
						if (fpldt != NULL) fclose(fpldt); if (fpfir != NULL) fclose(fpfir);
						
#define MPL_FREE_ALL	free(mpl_state); free2((char**)mpl_path); free(mpl_ratio); free(x);	free(mpl_fir);\
						if (fpldt != NULL) fclose(fpldt); if (fpfir != NULL) fclose(fpfir);


//----------------------------------------------------------------------------
void oupdfp(int it, int nz, double *z, double fz, double *g, MODEL model, UPDATEMODEL updatemodel, MSGOUTPUT msgoutput){       
    int i, m; 
    double *x, **dx;   
	char buff[256];
	char buffer[4096];

	//-------------------------------------------------
	x  = dalloc1(2*mpl_npath);
    dx = dalloc2(2*mpl_npath, nz);
    
    ztox(nz, 2*mpl_npath, mpl_Mtx, mpl_vct, z, x, dx);  

	updatemodel(model, NULL, NULL, NULL, NULL, NULL, NULL, x, NULL, NULL, NULL, 0);
	//UPDATEMODEL(MODEL model, int *segment, int *nmetastate, double *amp, double *xms, double **ar, double *pr, double *x, double *dk, double *dd, double **a);
   
    sprintf(buffer, "\r\n"
					"    Iter: %d   Func: %d   Dfunc: %d   Optimizer: DFP\r\n"
					"    --------------------------------------------\r\n"
					"      Logll:    %.2f\r\r\n", 
					it, mpl_nfunc, mpl_ndfunc, -fz);

    for( m=1; m<mpl_npath; m++ ) {
        sprintf(buff, "      %s     %d->%d\t%9.2f%9.3f\r\n", 
					m==0?"      ":"Rates:",
					mpl_path[m][0] + 1, 
					mpl_path[m][1] + 1, 
					exp(x[2*m]), x[2*m + 1] ); 
        strcat(buffer, buff);
		}

	strcat(buffer, "      Grad:  ");
    for( i=0; i<nz; i++) {
		sprintf(buff, "%9.4f", g[i]);
		strcat(buffer, buff);
		}
	strcat(buffer, "\r\n");
    
    msgoutput(buffer);
     
    free(x);   
    free2((char**)dx);
	}    

//--------------------------------------------------------------------------
void oupmdl(QUBOPT_VAR_NOT_USED char *mdlname, int nz, double *z, double *vz, double fz, int it, double tt, MODEL model, UPDATEMODEL updatemodel, MSGOUTPUT msgoutput){   
    int i, j, m; 
    double tmp, *x, **dx, *a, *va, *b, *vb;   

	//-----------------------------------------
	x  = dalloc1(2*mpl_npath);
    dx = dalloc2(2*mpl_npath, nz);
    a  = dalloc1(mpl_npath);
    va = dalloc1(mpl_npath);
    b  = dalloc1(mpl_npath);
    vb = dalloc1(mpl_npath);
    
    ztox(nz, 2*mpl_npath, mpl_Mtx, mpl_vct, z, x, dx);  
    
    for (m = 0; m < mpl_npath; m++) {
		a[m] = exp(x[2*m]);

		tmp=0.0;
        for( i=0; i<nz; i++ )
			tmp += sqr(dx[2*m][i])*vz[i];
        va[m] = sqr(a[m])*tmp;
        b[m]  = x[2*m + 1];

		tmp=0.0;
        for( i=0; i<nz; i++ )
			tmp += sqr(dx[2*m + 1][i])*vz[i];
        vb[m] = tmp; 
		}
    
    for( m=0; m<mpl_npath; m++) {
        va[m] = sqrt(va[m]);
        vb[m] = sqrt(vb[m]); 
		}

	updatemodel(model, NULL, NULL, NULL, NULL, NULL, NULL, NULL, va, vb, NULL, 0);

    //----- build output string 
    char buff[256];
	char buffer[4096];
	strcpy(buffer,"\r\n    ============ OPTIMIZATION RESULTS ============\r\n"
					"      start  end        a       da      b      db\r\n");

    for( m=0; m<mpl_npath; m++ ) {
        i = mpl_path[m][0] + 1;
		j = mpl_path[m][1] + 1;
		sprintf(buff, "        %d     %d%11.2f%8.2f%8.3f%8.3f\r\n", i, j, a[m], va[m], b[m], vb[m]);
		strcat(buffer, buff);
		}

	sprintf(buff,	"\r\n"
					"      likelihood:  %8.2f\r\n"
					"      iterations:  %d\r\n"	
					"      func call:   %d\r\n"	
					"      grad call:   %d\r\n"	
					"      core time:   %d\r\n\r\n\r\n",
					-fz, it, mpl_nfunc, mpl_ndfunc, (int)tt);
	strcat(buffer, buff);
    msgoutput(buffer);

    free(x);
    free2((char**)dx);
    free(a); 
    free(va);
    free(b);
    free(vb);   
	}       
       
//	-----------------------------------------------------------------------------------
void oupamp(int nclass, int nmutistate, double **qq, int namp, double *amp, double *xms, int nar, double **ar, int nmetastate, double *pr, MSGOUTPUT msgoutput){
	int i, j;
	char buff[256];
	char buffer[4096];

	//----- Title
	strcat(buffer, "\r\n\tOptimal values (I-sd-ar):\r\n");

	//----- for each amp[] : amp, xms, ar[1..nar] 
    for( i=0; i<namp; i++ ) {
		sprintf(buff, "\t\tmulticlass %d: %.3f %.3f", i, amp[i], xms[i]);
	    strcat(buffer, buff);
        for (j=1; j<=nar; j++) {
			sprintf(buff, " %.3f",ar[i][j]);
			strcat(buffer, buff);
			}
	    strcat(buffer, "\r\n");
    	}
	
    //----- output qq[i][j]
    for (i=0; i<nmutistate; i++) {
		strcat(buffer, (i==0) ? "\tmulti-A: " : "\t         ");
		for (j=0; j<nmutistate; j++) {
			sprintf(buff, "%.3f ", qq[i][j]);
			strcat(buffer, buff);
			}
		strcat(buffer, "\r\n");
		}
	strcat(buffer, "\r\n");

	//----- for each class : amp, std, ar[1..nar]
    for (i=0; i<nclass; i++) {
		sprintf(buff, "class: %d  amp: %.3f  std: %.3f", i, amp[i], xms[i]);
		strcat(buffer, buff);
		if (nar > 0) {
			strcat(buffer, "\tar: ");
			for (j=1; j<=nar; j++) {
				sprintf(buff, "%.3f ", ar[i][j]);
				strcat(buffer, buff);
				}
			}
		strcat(buffer, "\r\n");
    	}  

	//----- initial state probabilities ?
    strcat(buffer, "pr:");
    for (i=0; i<nmetastate; i++) {
		sprintf(buff, "  %.3f", pr[i]);
		strcat(buffer, buff);
		}
    strcat(buffer, "\r\n\r\n"); 
    msgoutput(buffer);
	}

//	-----------------------------------------------------------------------------------
void ovrlapint(int *source_data, int *dest_data, int data_count, double *response, int response_count, int fft_count){
    long	i, nd, np; 
	double	*s, *r, *y;
	int		rc1 = response_count - 1, rc12 = rc1 / 2;

	//Watch for filter shift !!!
	nd = fft_count - rc1;
	
	//Do we really need double precision here ?
	s = dalloc1(fft_count);
    r = dalloc1(fft_count);
    y = dalloc1(2 * fft_count);

    try {
		for (np = 0; np < data_count; np += nd) {

			//Clear s[]
			for (i = 0; i < fft_count; i++)
				s[i] = 0.0;

			for (i = 0; i < nd; i++)
				if (np + i < data_count)
					s[rc12 + i] = source_data[np + i];

			for (i = 0; i < response_count; i++)
				r[i] = response[i];
				//Is the response modified ?! Why is it copied ?!

			try {
				convolve0(s, fft_count, r, response_count, 1, y);
				} 
			catch (...) {
				}

			for (i = 0; i < rc12; i++)
				if (np - i >= 1)
					dest_data[np - i - 1] += (int)y[rc12 - i - 1];
			
			for (i = 0; i < nd + rc12; i++)
				if (np + i < data_count)
					dest_data[np + i] += (int)y[rc12 + i];
			}
		} 
    catch (...) {
    }
}


int filterintdata(int *source_data, int *dest_data, int data_count, double dt, double freq/*Hz*/, int fft_count) {
    int	nh;
	double * h = NULL;
	int rtnVal = 0;
	
	try {
		h = gsfir(freq*dt, &nh);
		wrap(nh, h);
		while ( fft_count < nh )
			fft_count *= 2;
		ovrlapint(source_data, dest_data, data_count, h, nh, fft_count);
		} 
	catch (...) {
#ifdef _WIN32
		unsigned int fstat = _clearfp(); // just in case
		_fpreset();
#endif
		rtnVal = 1;
		}
	
	free(h);
	return rtnVal;
	}

	
extern "C" QUBOPT_API double mlogl(int nar, int nfir, int ndata, double *scale){
    double s=0.0;
    for( int t=nar+nfir; t<ndata; t++ )
		if ( scale[t] > 0.0 )
			s -= log(scale[t]);
    return s;
	}

extern "C" QUBOPT_API void mdlogl(int nstate, int nar, int nfir, int nmetastate,
										   int **metastate, int ndata, double **alpha, 
										   double **beta, double **b, double **dll) {
    int	I, J, i0, j0, k, n, t, t0=nar+nfir;
    double tmp;
    
    for( i0=0; i0<nstate; i0++)
		for( j0=0; j0<nstate; j0++)
			dll[i0][j0] = 0.0;

    for( t=t0; t<ndata-1; t++) {
		for (J = 0; J < nmetastate; J++) {
			n = metastate[J][t0 + 3 + nstate];
			j0 = metastate[J][0];
			tmp = beta[J][t + 1]*b[J][t + 1];
			for (k = 0; k < n; k++) {
				I = metastate[J][t0 + 4 + nstate + k];
				i0 = metastate[I][0];
				dll[i0][j0] += alpha[I][t]*tmp;
				}
			}
		}
	}

//	-----------------------------------------------------------------------------------
void ovrlap(double *data, long ndata, double *respns, long m, double *ans, int nfft){
    long i, nd, np; 
	double *s, *r, *y;

	nd = nfft - (m - 1);

	s = dalloc1(nfft);
    r = dalloc1(nfft);
    y = dalloc1(2*nfft);
    
	memset( ans, 0, ndata*sizeof(double));	// for( i=0; i<ndata; i++) ans[i]=0.0;

    for( np=0; np<ndata; np+=nd) {
		for( i=1; i<=nfft; i++) 
			s[i-1] = 0.0;
		for( i=0; i<nd; i++)
			if( np+i < ndata )
				s[(m-1)/2+i] = data[np+i];
		for( i=0; i<m; i++) 
			r[i]=respns[i];

		try {
			convolve0(s, nfft, r, m, 1, y);
			} 
		catch (...) {
			}

		for( i=0; i<(m-1)/2; i++)
			if (np-i>=1) 
				ans[np-i-1] += y[(m-1)/2 - i - 1];

		for( i=1; i<=nd+(m-1)/2; i++)
			if (np+i<=ndata) 
				ans[np+i-1] += y[(m-1)/2 + i - 1];
		}

	free(s);
    free(r);
    free(y);
	}

//	-----------------------------------------------------------------------------------
int filterdatad(double *data, double *fdata, int ndata, double dt, double freq) {
    int	nh;
	double * h = NULL;
	int rtnVal = 0;
	int nfft = 1024;
	
	try {
		h = gsfir(freq*1e3*dt, &nh);
		wrap(nh, h);
		while ( nfft < nh )
			nfft *= 2;
		ovrlap(data, ndata, h, nh, fdata, nfft);
		} 
	catch (...) {
#ifdef _WIN32
		unsigned int fstat = _clearfp(); // just in case
		_fpreset();
#endif
		rtnVal = 1;
		}

	free(h);
	return rtnVal;
	}

//	-----------------------------------------------------------------------------------
int filterintdata_nfir(double freq, double dt) {
    double fc = freq * dt;
    double sigma = 0.1325 / fc;
    int n = (int)(4 * sigma);               
    return 2*n + 1;
}

//	-----------------------------------------------------------------------------------
double *gsfir(double fc, int *nfir){
    int	k, n;
    double tmp, sigma, *h;
    
    sigma = 0.1325/fc;
    n = (int)(4*sigma);               
    *nfir = 2*n + 1;

    h = dalloc1(*nfir + 3);
	
    if (sigma < 0.6) {
		*nfir = 3;
		h[0] = h[2] = sqr(sigma)/2.0;
		h[1] = 1 - h[0] - h[2];
		} 
	else {
		h[n] = 1.0/(sqrt(2.0*PI)*sigma);
		for (k = 1; k <= n; k++) {
			tmp = k/sigma;
			tmp = sqr(tmp);
			h[n+k] = h[n-k] = h[n]*exp(-tmp/2.0);
			}
		}
    return h;
	}

//	-----------------------------------------------------------------------------------
//	nh must be an odd integer
//	This is a strange function.   Used in filtering ?
//	Example output  : wrap( [1,2,3,4,5,6,7] ) -> [4,5,6,7,6,5,4]
void wrap(int nh, double *h)  {
    int i, m = (nh-1)/2;
    
    for( i=0; i<=m; i++) 
		h[i] = h[m+i];

    for( i=1; i<=m; i++)
		h[m+i] = h[m-i+1];
	}


int filterdatad_nfir_notch(double transition_width, double dt)
{
	double delta_f = transition_width * dt;
    const double width_factor = 3.3;
       /* 'hamming': 3.3,
          'hann': 3.1,
          'blackman': 5.5,
          'rectangular': 2.0, */
    int ntaps = int(width_factor/delta_f + 0.5);
    return (ntaps & ~0x1) + 1;   // ensure it's odd
}

int filterdatad_notch(double *source_data, double *dest_data,
	int data_count, double dt, double freq_lo, double freq_hi/*Hz*/, double transition_width, int fft_count)
{
	int rtnVal = 0;
	double wc_lo = 2*PI*freq_lo*dt;
	double wc_hi = 2*PI*freq_hi*dt;
	int nh = filterdatad_nfir_notch(transition_width, dt);
	int middle = (nh - 1) / 2;
	int order = nh - 1;
	double fmax = 0.0;
	double *h = new double[nh];
	for (int i=0; i<nh; ++i) {
		double f = (0.53836 - 0.46164*cos((PI*i)/order)) /*Hamming window*/
					* ((i == middle)
						? (1.0 - ((wc_hi - wc_lo)/PI))
						: ((sin((i-middle)*wc_lo) - sin((i-middle)*wc_hi)) / ((i-middle)*PI)));
		h[i] = f;
		fmax += f;
	}
	for (int i=0; i<nh; ++i)
		h[i] /= fmax;

	try {
		while ( fft_count < nh )
			fft_count *= 2;
		ovrlap(source_data, data_count, h, nh, dest_data, fft_count);
		} 
	catch (...) {
#ifdef _WIN32
		unsigned int fstat = _clearfp(); // just in case
		_fpreset();
#endif
		rtnVal = 1;
		}
	
	delete [] h;
	return rtnVal;
	}







/*	------------------------------------------------------------------------
	Formulate Ax=b into x=Az+b. n must be greater than m, and 
	constraints should be independent of each other.
*/
void freepar0(int m, int n, double *x, double **a, double *b, int *nz, double *z) {
    int i, j; 
    double tmp, eps = 1.0e-8, **v, *w, *h, sum;
    
	//no constraint
    if (m == 0) {
		for( i=0; i<n; i++) {
			for( j=0; j < n; j++)
				a[i][j] = 0.0;

			a[i][i] = 1.0;
			b[i] = 0.0;
			z[i] = x[i];
			}
		*nz = n;
		return;
		}

	//check if initial guess meets constraints
    for( i=0; i<m; i++) {
		tmp=0.0;
		for( j=0; j<n; j++)
			tmp += a[i][j]*x[j];
		tmp -= b[i];
		
		//My checkpoint - can be division by zero
		if( b[i]!=0 )
			tmp /= b[i];
		if( fabs(tmp) > eps )
			errexit("\tInitial guesses don't meet constraints");
		}

	//SVD    
    v = dalloc2(n, n);
    w = dalloc1(n);
    h = dalloc1(m);
    
    svdecomp0(a, m, n, w, v);
    *nz = n - m;

    for( i=0; i<m; i++) {
		sum = 0.0;
		for( j=0; j<m; j++)
			sum += a[j][i]*b[j];
		h[i] = sum / w[i];
		}

    for( i=0; i<n; i++ ) {
		sum=0.0;
		for( j=0; j<m; j++)
			sum += v[i][j]*h[j];
		b[i]=sum;
		for( j=0; j<*nz; j++)
			a[i][j] = v[i][j + m];
		}
    
	//initial z
    for( i=0; i<*nz; i++) {
		sum = 0.0;
		for( j=0; j<n; j++)
			sum += a[j][i]*x[j];
		z[i]=sum;
		}
    
    free(w);
    free(h);
    free2((char**)v);
	}

//---------------------------------------------------------------------------
void dfp_mpl_func(int nz, double *z, double *fz, double *gz, int *idef, 
			 MODEL model, MAKEMODEL makemodel, UPDATEDISPLAY updatedisplay){       
    int		i, j, k, n, nmetastate, **metastate;
    double  **alpha, **beta, *scale, **b;


	//****************************************
	short *idl_buffer;
	unsigned short  idl_sampling, idl_count;
	double	idl_max;
	int		idl_imax;
	int		idl_ndata;
	int		idl_leader;
	FILE	*fpidl = NULL;

	//****************************************
    *idef = 1;

    double  *x   = dalloc1(2*mpl_npath); 
    double  **dx  = dalloc2(2*mpl_npath, nz);  
	ztox(nz, 2*mpl_npath, mpl_Mtx, mpl_vct, z, x, dx);
    
	//----- Check if any of the rate constants is outside limits.  Exit if so.
	for( i=0; i<2*mpl_npath; i++) {
		if (x[i] < mpl_xlimit[i][0] || x[i] > mpl_xlimit[i][1]) {
			*idef = 0;
			free((char*)x); 
			free2((char**)dx); 
			return;
			}
		}
     
    double  **q   = dalloc2(mpl_nstate, mpl_nstate);
    double  **qq  = dalloc2(mpl_nmutistate, mpl_nmutistate);
    double  **a   = dalloc2(mpl_nmutistate, mpl_nmutistate);
    double  ***dq  = dalloc3(mpl_nstate, mpl_nstate, nz);
    double  ***dqq = dalloc3(mpl_nmutistate, mpl_nmutistate, nz);
    double  ***da  = dalloc3(mpl_nmutistate, mpl_nmutistate, nz);
    double  **dfa = dalloc2(mpl_nmutistate, mpl_nmutistate); 

	*fz = 0.0;
    for( i=0; i<nz; i++ )
		gz[i] = 0.0;


	//****************************************************
	//Open the idealized ldt file
	try {
		fpidl = fopen(mpl_ldtideal, "wb");
		} 
	catch (...) {
		}

	//Write the file header
	idl_sampling = (short)(mpl_dt[0]/1.0e-6);
	idl_count    = (short) mpl_count;
	idl_leader   = mpl_leader*idl_sampling;
	fwrite(&idl_leader, sizeof(idl_leader), 1, fpidl);
	fwrite(&idl_sampling, sizeof(idl_sampling), 1, fpidl);
	fwrite(&idl_count, sizeof(idl_count), 1, fpidl);

	//****************************************************
	try {
		//-----for each segment
		for( n=0; n<mpl_nseg; n++ ) {
			//Concentration and Voltage are a function of each segment - optimization accross C and V
			//Rate constants to Q matrix
			xtoq(mpl_nstate, mpl_npath, mpl_path, mpl_C[n], mpl_V[n], x, q, nz, dx, dq); 

			//Q matrix to multi-Q matrix
			qtomutiq(mpl_nstate, mpl_npath, mpl_path, mpl_nmutistate, mpl_mutistate, mpl_nmutipath, mpl_mutipath, q, qq, nz, dq, dqq);

			//Calculate the exponentials (???)
			mexp(mpl_nmutistate, mpl_dt[n], qq, a, nz, dqq, da, 0); 

			//Check if a[i][j] is between zero and one
			for( i=0; i<mpl_nmutistate; i++) {
				for (j=0; j<mpl_nmutistate; j++) {
					if (a[i][j] >= 1.0 || a[i][j] < 0.0) {
						*idef = 0;  //say some error message
						return;
						}
					}
				}
		
			//Basically, calculate meta-states (???)
			metamakv(mpl_nmutistate, mpl_nmuticlass, mpl_mutistate, mpl_nar + mpl_nfir, &nmetastate, &metastate);
		
			if ((mpl_nar + mpl_nfir) != 0)
				makemodel(model, NULL, NULL, NULL, &nmetastate, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, mpl_pr[n]);
			else
				makemodel(model, NULL, NULL, NULL, NULL, &mpl_nstate, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, mpl_pr[n]);

			alpha = dalloc2(nmetastate, mpl_ndat[n]);
			beta  = dalloc2(nmetastate, mpl_ndat[n]);
			b     = dalloc2(nmetastate, mpl_ndat[n]);
			scale = dalloc1(mpl_ndat[n]);
		
			//Make b vector		
			mvectb(mpl_ndat[n], mpl_dat[n], mpl_amp[n], mpl_xms[n], mpl_nar, mpl_ar[n], mpl_nfir, mpl_fir, nmetastate, metastate, b);
		

			//Modified for MChannels; for one is nstate
			mfward(mpl_nmutistate, mpl_nar, mpl_nfir, nmetastate, metastate, mpl_ndat[n], mpl_pr[n], a, b, alpha, scale);
			mbward(mpl_nmutistate, mpl_nar, mpl_nfir, nmetastate, metastate, mpl_ndat[n], a, b, scale, beta);
		
			//???		
			*fz -= mlogl(mpl_nar, mpl_nfir, mpl_ndat[n], scale);
		
			mdlogl(mpl_nmutistate, mpl_nar, mpl_nfir, nmetastate, metastate, mpl_ndat[n], alpha, beta, b, dfa);

			//Modified for MChannels; for one is nstate
			for (k = 0; k < nz; k++) {
				for (i = 0; i < mpl_nmutistate; i++) {
					for (j = 0; j < mpl_nmutistate; j++) {
						gz[k] -= dfa[i][j]*da[i][j][k];
						}
					}
				}
		
			//----------------------------------------------------------------------
			//-----  Idealized  LDT
			mgamma(mpl_nar, mpl_nfir, nmetastate, mpl_ndat[n], alpha, beta, scale);

			//Allocate memory
			idl_buffer = new short[mpl_ndat[n] + 2*mpl_leader];

			//For each point, find the state where the probability, as given by beta, is maximum
			//Firstly, do the points for which there is no probability
			for (i = 0; i < mpl_nar + mpl_nfir; i++)
				idl_buffer[i + mpl_leader] = (short)(mpl_amp[n][0]*(mpl_count/100.0));

			for (i = mpl_nar + mpl_nfir; i < mpl_ndat[n]; i++) {
				idl_max  = beta[0][i];
				idl_imax = 0;
				for (j = 0; j < nmetastate; j++) {
					if (beta[j][i] > idl_max) {
						idl_max  = beta[j][i];
						idl_imax = j;
						}
					}

				//idl_imax represents a metastate
				idl_buffer[i + mpl_leader] = (short)(mpl_amp[n][metastate[idl_imax][1]]*(mpl_count/100.0));
				}

			for (i = 0; i < mpl_leader; i++) {
				idl_buffer[i] = idl_buffer[i + mpl_leader];
				idl_buffer[i + mpl_leader + mpl_ndat[n]] = idl_buffer[i + mpl_leader];
				}

			idl_ndata = mpl_ndat[n] + 2*mpl_leader;

			//Write segment header
			fwrite(&mpl_tstart[n], sizeof(mpl_tstart[n]), 1, fpidl);
			fwrite(&idl_ndata, sizeof(idl_ndata), 1, fpidl);
			//Write data points		  
			fwrite(idl_buffer, idl_ndata, sizeof(short), fpidl);

			delete [] idl_buffer;
			
			//------------------------------------------------------------------------
			free2((char**)metastate);
			free((char**)scale);
			free2((char**)alpha);
			free2((char**)beta);
			free2((char**)b);
			}
		} 
	catch (...) {
	}
	//************
	fclose(fpidl);
	
	//************
	mpl_nfunc++; 
	mpl_ndfunc++;     
	updatedisplay(&mpl_firsttime);
	
	free((char*)x); 
	free2((char**)dx); 
	free2((char**)q);
	free2((char**)qq); 
	free2((char**)a); 
	free3((char***)dq);
	free3((char***)dqq); 
	free3((char***)da); 
	free2((char**)dfa); 
}


//**************************************************************************************************
extern "C" int QUBOPT_API mpl(char **pmdlname, char **pldtname, char **pldtideal, MODEL model, double voltage, double concentration, int segment, QUBOPT_VAR_NOT_USED int currseg, double frequency, int filtercoeffs, double datafreq, MAKEMODEL makemodel, MAKECONSTRAINTS makeconstraints, UPDATEMODEL updatemodel, UPDATEDISPLAY updatedisplay, MSGSTATUS msgstatus, MSGOUTPUT msgoutput, int *processflag, QUBOPT_VAR_NOT_USED int errcode, double *lnlikelihood){
	int		n, nn, nnn, mc, nz, it, i, j, k;

    double  *x = NULL, *z = NULL, *vz = NULL, fz, tt = 0;
	double	*dataf = NULL;

	double	*ta = NULL, *tx = NULL, **tr = NULL;

	double	ts;

    char	msg[256];
    char	*buffer;
    
    FILE	*fpfir = NULL, *fpldt = NULL;   

    //-------------------------------------------------------
	*processflag = PROCESS_STARTED;

	sprintf(msg, "Allocating memory ...");
	msgstatus(msg);
	
	mpl_firsttime = 1;

	//global name of the ideal ldt file
	sprintf(mpl_ldtideal, "%s", *pldtideal);
	
	//----- Alocate memory 1 
	try {
		mpl_state = ialloc1(MAXSTATE);
		mpl_path  = ialloc2(MAXPATH,3); 
		mpl_ratio = dalloc1(MAXSTATE);
		x		  = dalloc1(2*MAXPATH);     
		} 
	catch (...) {
		}		

	if (mpl_state==NULL || mpl_path==NULL || mpl_ratio==NULL || x==NULL) {
		sprintf(msg, "Memory allocation error !");
		msgstatus(msg);
		MPL_FREE_ALL
		return OUT_OF_MEMORY;
		}

    //*********************************************************
	sprintf(msg, "Reading model ...");
	msgstatus(msg);

	//Read model ************************************************************************************
	if (makemodel(model, NULL, &mpl_nchannel, &mpl_nclass, NULL, &mpl_nstate, &mpl_npath, &mpl_nar, mpl_state, mpl_path, mpl_ratio, x, NULL, NULL, NULL, NULL)) {
		MPL_FREE_ALL
		return DONT_DISPLAY;
		}

	//***********************************************************************************************
	//Multimodel
	//Modified for MChannels; for one is disabled
    mutimakv(mpl_nchannel, mpl_nstate, mpl_npath, mpl_state, mpl_path, mpl_ratio, &mpl_nmutistate, &mpl_nmutipath, &mpl_nmuticlass, &mpl_mutistate, &mpl_mutipath);

	//Constraints
    mpl_Mtx = dalloc2(2*mpl_npath, 2*mpl_npath);
    mpl_vct = dalloc1(2*mpl_npath);
    sprintf(msg, "Reading constraints ...");
	msgstatus(msg);
	//input constraints
	//inpcns(mpl_nstate, mpl_npath, mpl_path, x, mpl_Mtx, mpl_vct, &mc);
	i = 0;
	do {
		for (i = 0; i < 2*mpl_npath; i++) {
			for (j = 0; j < 2*mpl_npath; j++) {
				mpl_Mtx[i][j] = 0;
				}
			mpl_vct[i] = 0;
			}

		i = makeconstraints(model, mpl_path, x, mpl_Mtx, mpl_vct, &mc);
		if (i == LOOP_UNBALANCE) {
			MPL_FREE_ALL
			return LOOP_UNBALANCE;
			} 
		else if (i == LOOP_ADJUSTED) {
			makemodel(model, NULL, &mpl_nchannel, &mpl_nclass, NULL, &mpl_nstate, &mpl_npath, &mpl_nar, mpl_state, mpl_path, mpl_ratio, x, NULL, NULL, NULL, NULL);
			}
		} while (i == LOOP_ADJUSTED);

	//**********************************************************************************************
	//Read data, dt, voltage and concentration #####################################################
	//
	//                 ONLY FOR ONE LDT FILE !!!
	//
	//                 ONLY ONE DRUG !!!		
	//
	//
	//Open ldt file
        try {
		fpldt = fopen(*pldtname, "rb");
		} 
	catch (...) {
		}

	if (fpldt == NULL) {
		sprintf(msg, "Cannot open data file !");
		msgstatus(msg);
		MPL_FREE_ALL
		return CANT_READ_FILE;
		}
	
	//Read ldt file header
	rdHeader(fpldt, &mpl_leader, &ts, &mpl_count);

	//Read ldt file segments
	n = 0; nn = 0;

	while(rdSegHeader(fpldt, &mpl_tstart[n], &mpl_ndat[n])) {
		if (((segment > -1) && (nn == segment)) || (segment == -1)) {
			//Read segment
			
			sprintf(msg, "Reading segment %d ...", nn);
			msgstatus(msg);
			
			//Sampling rate
			mpl_dt[n] = ts;
			
			//Voltage
			mpl_V[n]  = voltage;
			
			//Concentration
			mpl_C[n]  = dalloc1(MAXDRUG + 1);
			mpl_C[n][0] = concentration;
			
			//Data
			mpl_dat[n] = dalloc1(mpl_ndat[n]);
			
			rdSegData(fpldt, mpl_count, mpl_ndat[n], mpl_dat[n]);
			mpl_ndat[n] = mpl_ndat[n] - 2*mpl_leader;
			for (i = 0; i < mpl_ndat[n]; i++) {
				mpl_dat[n][i] = mpl_dat[n][mpl_leader + i];
				}
			
			//Filter data, if necessary
			if (datafreq > 0) {
				try {
					dataf = (double*)realloc(dataf, (mpl_ndat[n] + 1)*sizeof(double));
					if (filterdatad(mpl_dat[n], dataf, mpl_ndat[n], mpl_dt[n], datafreq) == 0) {
						for (i = 0; i < mpl_ndat[n]; i++) {
							mpl_dat[n][i] = dataf[i];
							}
						}
					} 
				catch (...) {
					}
				}
			
			//Allocate memory for Amp, Std, Ar, Pr
			//Modified for MChannels; for one is disabled
			
			//mpl_amp[n] = dalloc1(mpl_nmutistate + 1/*mpl_nmuticlass*/);
			//mpl_xms[n] = dalloc1(mpl_nmutistate + 1/*mpl_nmuticlass*/);
			//mpl_ar[n]  = dalloc2(mpl_nmutistate + 1/*mpl_nmuticlass*/, mpl_nar + 1);
			//mpl_pr[n]  = dalloc1(mpl_nmuticlass + 1);
			mpl_amp[n] = dalloc1(MAXMUTICLASS);
			mpl_xms[n] = dalloc1(MAXMUTICLASS);
			mpl_ar[n]  = dalloc2(MAXMUTICLASS, MAXAR + 1);
			mpl_pr[n]  = dalloc1(MAXMUTICLASS);
			
			
			for (i = 0; i < mpl_nmuticlass; i++) {
				mpl_amp[n][i] = 0;
				mpl_xms[n][i] = 0;
				mpl_ar[n][i][0] = 1;
				for (j = 1; j < (mpl_nar + 1); j++) {
					mpl_ar[n][i][j] = 0;
					}
				}
			for (i = 0; i < mpl_nmutistate; i++) {
				mpl_pr[n][i] = 0;
				}
			
			
			
			if (mpl_nchannel == 1) {		//One Channel
				if (segment > -1)
					//Curr segment
					makemodel(model, NULL, NULL, NULL, NULL, &mpl_nstate, NULL, &mpl_nar, NULL, NULL, NULL, NULL, mpl_amp[n], mpl_xms[n], mpl_ar[n], mpl_pr[n]);
				else
					//This segment
					makemodel(model, &n,   NULL, NULL, NULL, &mpl_nstate, NULL, &mpl_nar, NULL, NULL, NULL, NULL, mpl_amp[n], mpl_xms[n], mpl_ar[n], mpl_pr[n]);
				} 
			else {						//Multiple Channels
				ta  = dalloc1(mpl_nclass);
				tx  = dalloc1(mpl_nclass);
				tr  = dalloc2(mpl_nclass, mpl_nar + 1);
				
				if (segment > -1)	
					//Curr segment
					makemodel(model, NULL, NULL, NULL, NULL, &mpl_nstate, NULL, &mpl_nar, NULL, NULL, NULL, NULL, ta, tx, tr, mpl_pr[n]);
				else				
					//This segment
					makemodel(model, &n,   NULL, NULL, NULL, &mpl_nstate, NULL, &mpl_nar, NULL, NULL, NULL, NULL, ta, tx, tr, mpl_pr[n]);
				
				for( nnn=0; nnn<mpl_nmutistate; nnn++ ) {
					for (mpl_amp[n][nnn] = 0.0, mpl_xms[n][nnn] = 0, i = 0; i < mpl_nstate; i++) {
						k = mpl_state[i];
						mpl_amp[n][nnn] += (mpl_mutistate)[nnn][i + 1]*ta[k];
						mpl_xms[n][nnn] += (mpl_mutistate)[nnn][i + 1]*tx[k];
						for (j = 1; j < (mpl_nar + 1); j++) {
							mpl_ar[n][nnn][j] += (mpl_mutistate)[nnn][i + 1]*tr[k][j];
							}
						}
					for (j = 1; j < (mpl_nar + 1); j++) {
						mpl_ar[n][nnn][j] -= tr[0][j];
						}
					mpl_amp[n][nnn] -= ta[0];
					mpl_xms[n][nnn] -= tx[0];
					}
				
					/*
					for (nnn = 0; nnn < mpl_nmuticlass; nnn++) {
					for (mpl_amp[n][nnn] = 0.0, mpl_xms[n][nnn] = 0, i = 0; i < mpl_nstate; i++) {
					k = mpl_state[i];
					mpl_amp[n][nnn] += (mpl_mutistate)[nnn][i + 1]*ta[k];
					mpl_xms[n][nnn] += (mpl_mutistate)[nnn][i + 1]*tx[k];
					for (j = 1; j < (mpl_nar + 1); j++) {
					mpl_ar[n][nnn][j] += (mpl_mutistate)[nnn][i + 1]*tr[k][j];
					}
					}
					for (j = 1; j < (mpl_nar + 1); j++) {
					mpl_ar[n][nnn][j] -= tr[0][j];
					}
					mpl_amp[n][nnn] -= ta[0];
					mpl_xms[n][nnn] -= tx[0];
					}
				*/
				
				free(ta);
				free(tx);
				free((char**)tr);
				}
			
			//-----
			sprintf(msg, "\tsegment %d\r\n", n + 1);
			msgoutput(msg);
			buffer = (char*)calloc(4096, sizeof(char));
			for (i = 0; i < mpl_nclass; i++) {
				sprintf(msg, "\t\tmulticlass: %d  amp: %.3f  std: %.3f", i, mpl_amp[n][i], mpl_xms[n][i]);
				buffer = strcat(buffer, msg);
				if (mpl_nar > 0) {
					buffer = strcat(buffer, "  ar: ");
					for (j = 1; j <= mpl_nar; j++) {
						sprintf(msg, "%.3f ", mpl_ar[n][i][j]);
						buffer = strcat(buffer, msg);
						}
					}
				buffer = strcat(buffer, "\r\n");
				}
			msgoutput(buffer);
			free(buffer);
			
			n++;
			
			if (segment > -1)
				break;
			} 
		else {
			//Skip segment
			fseek(fpldt, mpl_ndat[n]*sizeof(short), SEEK_CUR);
			}
		
		nn++;
		}

    if (segment != -1)
		mpl_nseg = 1;
	else
		mpl_nseg = n;
		
	if (datafreq > 0)
		free(dataf);

	fclose(fpldt);
	fpldt = NULL;
	
	//##############################################################################################
	//*********************************** F I L T E R **********************************************
	//
    //nfir refers to the order here
    if (filtercoeffs == 0) {
		mpl_nfir = 0;
		try {
			mpl_fir  = dalloc1(mpl_nfir + 1);
			mpl_fir[0] = 1.0;
		} catch (...) {
		};
		if (mpl_fir == NULL) {
			sprintf(msg, "Memory allocation error !");
			msgstatus(msg);
			MPL_FREE_ALL
			return OUT_OF_MEMORY;
			}
		}	 
	else {
		double	*h = NULL, *r = NULL, fc;
		int	nh, ki;
		fc = frequency*1e3;

		h = gsfir(fc*ts, &nh);
		r = new double[nh + 1];
		r[0]=h[0];
		for( i=1; i<nh; i++)
			r[i] = r[i - 1] + h[i];
		free(h);

		ki = (nh - filtercoeffs) / 2;
		for (i = 0; i < filtercoeffs - 1; i++)
			r[i] = r[ki + i];
		r[filtercoeffs - 1] = 1.0;
		mpl_fir = dalloc1(filtercoeffs);
		mpl_fir[0] = r[0];
		for (i = 1; i < filtercoeffs; i++)
			mpl_fir[i] = r[i] - r[i - 1];
		mpl_nfir = filtercoeffs - 1;
		delete [] r;
		}
	
    mpl_xlimit = dalloc2(2*mpl_npath, 2); 
    
    for (i = 0; i < 2*mpl_npath; i++) {
		mpl_xlimit[i][0] = log(LOWLIMIT) + x[i];   
        mpl_xlimit[i][1] = log(UPLIMIT)  + x[i];   
		}
    
    z = dalloc1(2*mpl_npath);
    
    freepar0(mc, 2*mpl_npath, x, mpl_Mtx, mpl_vct, &nz, z);
    
    //Modified for MChannels; for one is active
	//mutimakv(mpl_nchannel, mpl_nstate, mpl_npath, mpl_state, mpl_path, mpl_ratio, &mpl_nmutistate, &mpl_nmutipath, &mpl_nmuticlass, &mpl_mutistate, &mpl_mutipath);
    
    vz = dalloc1(2*mpl_npath);
    
	sprintf(msg, "Optimizing rates ...");
	msgstatus(msg);

	fz = 0;
    dfpminim(nz, z, &fz, vz, &it, dfp_mpl_func, oupdfp, model, makemodel, updatemodel, updatedisplay, msgoutput, processflag);
           
    sprintf(msg, "Writing output ...");
	msgstatus(msg);

	oupmdl(*pmdlname, nz, z, vz, fz, it, tt, model, updatemodel, msgoutput);

	*lnlikelihood = -fz;

	
	MPL_FREE_ALL

	return 0;
	}



extern "C" int QUBOPT_API amp(QUBOPT_VAR_NOT_USED char **pmdlname, char **pldtname, QUBOPT_VAR_NOT_USED char **pldtideal, MODEL model, double voltage, double concentration, int segment, int currseg, int noupdatep, double frequency, int filtercoeffs, double datafreq, MAKEMODEL makemodel, UPDATEMODEL updatemodel, QUBOPT_VAR_NOT_USED UPDATEDISPLAY updatedisplay, MSGSTATUS msgstatus, MSGOUTPUT msgoutput, int *processflag, int *errcode, double *lnlikelihood)
{
    //**********************************************************************************************
	int     nchannel, nstate, npath, nclass,
	        nmutistate, nmutipath, nmuticlass,
			nar, nfir, nmetastate,
            nseg, nlead, ndata, tstart, count,
            i, j, k, n, s;
    
    int		*state=NULL, **path=NULL, **mutistate, **mutipath, **metastate;
    
	double  *ratio = NULL, *x = NULL, **q = NULL, **qq = NULL, **a = NULL,
            **dx = NULL, ***dq = NULL, ***dqq = NULL, ***da = NULL,
            *data = NULL, *dataf = NULL, *amp = NULL, *xms = NULL, **ar = NULL, *fir = NULL, *pr = NULL,
			*ta = NULL, *tx = NULL, **tr = NULL;

	double	dt, V, C[1], ll, lltot = 0;
    
    char	msg[256];
    
    FILE    *fpldt = NULL, *fpfir = NULL;
    //**********************************************************************************************           
    
    
	*errcode = 1;

	*processflag = PROCESS_STARTED;

	sprintf(msg, "Allocating memory ...");
	msgstatus(msg);
	
	
	
	//Alocate memory 1 ****************************************
	try {
		state = ialloc1(MAXSTATE);
		path  = ialloc2(MAXPATH,3); 
		ratio = dalloc1(MAXSTATE);
		x     = dalloc1(2*MAXPATH);     
		} 
	catch (...) {
		};		
	if (state==NULL || path==NULL || ratio==NULL || x==NULL) {
		sprintf(msg, "Memory allocation error !");
		msgstatus(msg);
		AMP_FREE_ALL
		return OUT_OF_MEMORY;
	};
    //*********************************************************
	


	sprintf(msg, "Reading model...");
	msgstatus(msg);
	
	
	
	//Read model ************************************************************************************

	if (makemodel(model, NULL, &nchannel, NULL, NULL, &nstate, &npath, &nar, state, path, ratio, x, NULL, NULL, NULL, NULL)) {
		AMP_FREE_ALL
		return DONT_DISPLAY;
	};
	//***********************************************************************************************

   

	//Open ldt file ****************************
	sprintf(msg, "Opening data file...");
	msgstatus(msg);
	try {
		fpldt = fopen(*pldtname, "rb");
		} 
	catch (...) {
		};
	if (fpldt == NULL) {
		sprintf(msg, "Cannot open data file !");
		msgstatus(msg);
		AMP_FREE_ALL
		return CANT_READ_FILE;
	};
	//******************************************

    

	//Read ldt file header **************
   	sprintf(msg, "Reading file header...");
	msgstatus(msg);
	rdHeader(fpldt, &nlead, &dt, &count); 
	//***********************************


    
    //nfir refers to the order here
    if (filtercoeffs == 0) {
		nfir = 0;
		try {
			fir  = dalloc1(nfir + 1);
			fir[0] = 1.0;
			} 
		catch (...) {
			};
		if (fir == NULL) {
			sprintf(msg, "Memory allocation error !");
			msgstatus(msg);
			AMP_FREE_ALL
			return OUT_OF_MEMORY;
		};
	} 

	else {
		double	*h = NULL, *r = NULL, fc;
		int		nh, ki;

		fc = frequency*1e3;

		h = gsfir(fc*dt, &nh);
		r = new double[nh + 1];
		r[0]=h[0];
		for( i=1; i<nh; i++)
			r[i] = r[i-1] + h[i];
		free(h);
		ki = (nh - filtercoeffs) / 2;
		for( i=0; i<filtercoeffs-1; i++)
			r[i] = r[ki + i];
		r[filtercoeffs - 1] = 1.0;
		fir = dalloc1(filtercoeffs);
		fir[0] = r[0];
		for ( i=1; i<filtercoeffs; i++ )
			fir[i] = r[i] - r[i-1];
		nfir = filtercoeffs - 1;
		delete [] r;
		}

   	sprintf(msg, "Applying mutimakv...");
	msgstatus(msg);
	mutimakv(nchannel, nstate, npath, state, path, ratio, &nmutistate, &nmutipath, &nmuticlass, &mutistate, &mutipath); 
   	sprintf(msg, "Applying metamakv...");
	msgstatus(msg);
    metamakv(nmutistate, nmuticlass, mutistate, nar + nfir, &nmetastate, &metastate);
	
	nclass = 0;
	for( i=0; i<nstate; i++)
		if (state[i] > nclass)
			nclass = state[i];
	nclass++;

	//Alocate memory 2 ****************************************
   	sprintf(msg, "Allocating memory...");
	msgstatus(msg);
	try {
		q   = dalloc2(nstate, nstate);
		qq  = dalloc2(nmutistate, nmutistate);
		dx  = dalloc2(2*npath, 1);
		dq  = dalloc3(nstate, nstate, 1);
		dqq = dalloc3(nmutistate, nmutistate, 1);
		a   = dalloc2(nmutistate, nmutistate);
		da  = dalloc3(nmutistate, nmutistate, 1);
		} 
	catch (...) {
		}

	if (q == NULL || qq == NULL || a == NULL || dx == NULL || q == NULL || dqq == NULL || da == NULL) {
		sprintf(msg, "Memory allocation error !");
		msgstatus(msg);
		AMP_FREE_ALL
		return OUT_OF_MEMORY;
		}
	
	//------------------------------------------------------------
    C[0] = concentration;
    V = voltage; 

    //Rate constants and derivatives *****************
   	sprintf(msg, "Applying xtoq...");
	msgstatus(msg);
	xtoq(nstate, npath, path, C, V, x, q, 1, dx, dq);  

    //Rate constants and derivatives multi-channel **********************************************
   	sprintf(msg, "Applying qtomutiq...");
	msgstatus(msg);
	qtomutiq(nstate, npath, path, nmutistate, mutistate, nmutipath, mutipath, q, qq, 1, dq, dqq); 
	
   	sprintf(msg, "Applying mexp...");
	msgstatus(msg);
	mexp(nmutistate, dt, qq, a, 1, dqq, da, 0);

	//Allocate memory 3 ******************************
   	sprintf(msg, "Allocating memory...");
	msgstatus(msg);
    try {
		amp = dalloc1(nmutistate + 1/*nmuticlass*/);  //see here how many ?!
		xms = dalloc1(nmutistate + 1/*nmuticlass*/);
		ar  = dalloc2(nmutistate + 1/*nmuticlass*/, nar + 1); 
		pr  = dalloc1(nmetastate + 1);
		} 
	catch (...) {
		}

	if (amp == NULL || xms == NULL || ar == NULL || pr == NULL) {
		sprintf(msg, "Memory allocation error !");
		msgstatus(msg);
		AMP_FREE_ALL
		return OUT_OF_MEMORY;
		}

	//************************************************
	for( i=0; i<nmuticlass; i++) {
		amp[i] = 0;
		xms[i] = 0;
		ar[i][0] = 1;
		for( j=1; j<(nar+1); j++ )
			ar[i][j] = 0;
		}
	for( i=0; i<nmetastate; i++ )
		pr[i] = 0;

   	sprintf(msg, "Reading amp, std and ar...");
	msgstatus(msg);
	
	if (nchannel == 1) {				//One Channel
		makemodel(model, NULL, NULL, NULL, &nmetastate, &nstate, NULL, NULL, NULL, NULL, NULL, NULL, amp, xms, ar, pr);
		} 
	else {						//Multiple Channels
		ta  = dalloc1(nclass);
		tx  = dalloc1(nclass);
		tr  = dalloc2(nclass, nar + 1); 
		//Initially, the amps for the current segment are read, so segment is NULL
		makemodel(model, NULL, NULL, NULL, &nmetastate, &nstate, NULL, NULL, NULL, NULL, NULL, NULL, ta, tx, tr, pr);
		for (n = 0; n < nmutistate; n++) {
			for (amp[n] = 0.0, xms[n] = 0, i = 0; i < nstate; i++) {
				k = state[i];
				amp[n] += (mutistate)[n][i + 1]*ta[k];
				xms[n] += (mutistate)[n][i + 1]*tx[k];
				for (j = 1; j < (nar + 1); j++) {
					ar[n][j] += (mutistate)[n][i + 1]*tr[k][j];
				};
			};
			for (j = 1; j < (nar + 1); j++) {
				ar[n][j] -= tr[0][j];
			};
			amp[n] -= ta[0];
			xms[n] -= tx[0];
		};
		free(ta);
		free(tx);
		free((char**)tr);
	};
	//****************************************************************************************

	data = NULL;//just in case, to be sure is NULL
	dataf = NULL;

	nseg = 1;
	while(rdSegHeader(fpldt, &tstart, &ndata)) {
		//Read segment
		if ( ((segment > -1) && (nseg == (segment + 1))) || (segment == -1)) {
			//Alocate memory for one segment of data
			sprintf(msg, "Allocating memory...");
			msgstatus(msg);

			try {
				data = (double*)realloc(data, (ndata + 1)*sizeof(double));
				} 
			catch (...) {
				}
			
			if (data == NULL) {
				sprintf(msg, "Memory allocation error !");
				msgstatus(msg);
				AMP_FREE_ALL
				return OUT_OF_MEMORY;
			};
	

			//Read one segment
   			sprintf(msg, "Processing segment %i...\r\n", nseg);
			msgoutput(msg);
			rdSegData(fpldt, count, ndata, data);
			ndata -= 2*nlead;
			for (i = 0; i < ndata; i++)
				data[i] = data[nlead + i];

			//Filter data, if necessary
			if (datafreq > 0) {
				try {
					dataf = (double*)realloc(dataf, (ndata + 1)*sizeof(double));
					if (filterdatad(data, dataf, ndata, dt, datafreq) ==0)
						for (i = 0; i < ndata; i++)
							data[i] = dataf[i];
					} 
				catch (...) {
					}
				}

			sprintf(msg, "Applying mbaum...");
			msgstatus(msg);
			try { 
				mbaum(ndata, data, nmutistate, a, nmetastate, metastate, pr, nmuticlass, amp, xms, nar, ar, nfir, fir, &ll, NULL, msgoutput, processflag);
				} 
			catch (...) {
				sprintf(msg, "Error in MBaum !");
				msgstatus(msg);
				AMP_FREE_ALL
				return ERROR_MBAUM;
				}

			lltot += ll;
   			//sprintf(msg, "Converting A matrix to Q matrix...");
			//msgstatus(msg);
			//atoq(a, nmutistate, dt, qq);
       			
			sprintf(msg, "Writing output...");
			msgstatus(msg);
			//oupamp(nmutistate, qq, nmuticlass, amp, xms, nar, ar, nmetastate, pr, msgoutput);
			oupamp(nclass, nmutistate,  a, nmuticlass, amp, xms, nar, ar, nmetastate, pr, msgoutput);

			s = nseg - 1;

			if ((nar != 0) || (nfir !=0)) {
				//Has metastates
				//only if all segments
				updatemodel(model, &s, &nmetastate, amp, xms, ar, pr, NULL, NULL, NULL, a, noupdatep);
				if (s == currseg)
					updatemodel(model, NULL, &nmetastate, amp, xms, ar, pr, NULL, NULL, NULL, NULL, noupdatep);
				} 
			else {
				//No metastates
				updatemodel(model, &s, NULL, amp, xms, ar, pr, NULL, NULL, NULL, a, noupdatep);
				if (s == currseg)
					updatemodel(model, NULL, NULL, amp, xms, ar, pr, NULL, NULL, NULL, NULL, noupdatep);
				}
		
			if (segment > -1)
				break;
			} 
		else {
			//Skip segment
			fseek(fpldt, ndata*sizeof(short), SEEK_CUR);
			}
		nseg++;
		}
		
	free(data);
	if (datafreq > 0) {
		free(dataf);
	};
	
	if (fpldt != NULL) fclose(fpldt);
	if (fpfir != NULL) fclose(fpfir);

	*errcode = 0;

	*lnlikelihood = lltot;

	free(state);
	free2((char**)path);
	free(ratio);
	free(x);
	free(fir);
	
	free2((char**)q);
	free2((char**)qq);
	free2((char**)dx);
	free3((char***)dq);
	free3((char***)dqq);
	free3((char***)da);
	
	try {
		free2((char**)a);//something is damaging a. Why ?!
		} 
	catch (...) {
		}
	
	return 0;
	}
