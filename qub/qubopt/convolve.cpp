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

// utilities for Amp, Mpl 

//----- local functions 
void fourier1_0(double *data, long nn, int isign);
void two_fft0(double *data1, double *data2, double *fft1, double *fft2, long n);
void real_fft0(double *data, long n, int isign);

//	-----------------------------------------------------------------------------------
//	Fast fourier transform : 0 based arrays.
//	Replaces data [0..2*nn-1] with ...
//		isign=1 : its discrete fourier transform
//		isign=-1: nn * its inverse discrete fourier transform
//	Data is a complex array of length nn or a real array of length nn
//	nn must be an integer power of 2 !!!
#define SWAP(a, b) tempr = (a); (a) = (b); (b) = tempr
void fourier1_0(double *data, long nn, int isign){
	long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;
	for( i=1; i<n; i+=2 ) {		// This is the bit reversal section of the routine
		if (j > i) {			// Exchange the two complex numbers
			SWAP( data[j-1], data[i-1]);
			SWAP( data[j], data[i]);
			}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
			}
		j += m;
		}

	//----- Danielson Lanczos 
	mmax = 2;
	while( n>mmax ) {					// outer loop executed log2(nn) times
		istep = mmax << 1;
		theta = ((double)isign) * ( PI2 / (double) mmax);	// initialize the trigonometric occurrence
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for( m=1; m<mmax; m+=2) {
			for( i=m; i<=n; i+=istep ) {
				j = i + mmax;			// Danielson Lanczos formula
				tempr = wr*data[j-1] - wi*data[j];	
				tempi = wr*data[j] + wi*data[j-1];
				data[j-1] = data[i-1] - tempr;
				data[j] = data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
				}
			wr = (wtemp=wr)*wpr - wi*wpi + wr;	// Trigonometric recurrence
			wi = wi*wpr + wtemp*wpi + wi;
			}
		mmax = istep;
		}
	}
#undef SWAP

//	-----------------------------------------------------------------------------------
//	two ffts - 0 based arrays.
//	Simultaneous fourier1() of two real arrays using real and imag of complex numbers alg.
//	Given data1[0..n-1] data2[0..n-1], call four1 and return fft1[0..2n-1],fft2[0..2n-1]
//	complex arrays of FFT's
//	n must be an integer power of 2 !! 
void two_fft0(double *data1, double *data2, double *fft1, double *fft2, long n){
	double rep, rem, aip, aim;
	long j, nn=n+n;

	// pack the two real arrays into one complex array
	for( j=0; j<n; ++j) {
		fft1[j+j] = data1[j];
		fft1[j+j+1] = data2[j];
		}

	// transform the complex array 
	fourier1_0(fft1, n, 1);

	fft2[0] = fft1[1];
	fft1[1] = fft2[1] = 0.0;

	for( j=2; j<=n; j+=2 ) {	// Use symmetries to separate the two transforms
		rep = 0.5 * (fft1[j] + fft1[nn-j]);
		rem = 0.5 * (fft1[j] - fft1[nn-j]);
		aip = 0.5 * (fft1[j+1] + fft1[nn+1-j]);
		aim = 0.5 * (fft1[j+1] - fft1[nn+1-j]);
		fft1[j] = rep;			// Ship them out in two complex arrays
		fft1[j+1] = aim;
		fft1[nn-j] = rep;
		fft1[nn+1-j] = -aim;
		fft2[j] = aip;
		fft2[j+1] = -rem;
		fft2[nn-j] = aip;
		fft2[nn+1-j] = rem;
		}
	}

//	-----------------------------------------------------------------------------------
//	realft : Numerical recipes in c - 1992, modifed for 0 based arrays.
//	FFT a single real data set:  data1(0:n-1)
void real_fft0(double *data, long n, int isign){
	long i, ii;
	double c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;

	theta = PI / (double) (n >> 1);

	if (isign == 1) {
		c2 = -0.5;
		fourier1_0(data, n >> 1, 1);
		} 
	else {
		c2 = 0.5;
		theta =- theta;
		}

	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;

	for( i=2; i<=(n>>2); i++) {
		ii=i+i;
		h1r = c1*( data[ii-2] + data[n+2-ii] );
		h1i = c1*( data[ii-1] - data[n+3-ii] );
		h2r =-c2*( data[ii-1] + data[n+3-ii] );
		h2i = c2*( data[ii-2] - data[n+2-ii] );
		data[ii-2] = h1r + wr*h2r - wi*h2i;
		data[ii-1] = h1i + wr*h2i + wi*h2r;
		data[n+2-ii] = h1r - wr*h2r + wi*h2i;
		data[n+3-ii] =-h1i + wr*h2i + wi*h2r;
		wr = (wtemp = wr)*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
		}

    if (isign == 1) {
		data[0] = (h1r = data[0]) + data[1];
		data[1] = h1r - data[1];
		} 
	else {
		data[0] = c1*( (h1r=data[0]) + data[1] );
		data[1] = c1*( h1r-data[1] );
		fourier1_0(data, n >> 1, -1);
		}
	}

/*
 * Convolve or deconvolve a real data set data[0..n-1] with a response 
 * function respns[0..n-1].
 *  1. n must be an integer power of two.
 *  2. data[] must be padded with zeros (at the end) in the calling 
 *     program, and the number of zeros should not be less than the 
 *     positive duration of the response function (not including t=0).
 *  3. m must be an odd integer.
 *  4. The response function must be stored in wrap-around order in  
 *     the first m elements of respons[].
 *  5. The result is returned in the first n components of ans[].
 *     However, ans[] must be supplied in the calling program with 
 *     dimensions [0..2*n-1].
 */
void convolve0(double *data, long n, double *respns, long m, int isign, double *ans){
	long i, no2;
	double dum, mag2, *fft;

	fft = dalloc1(2*(int)n);
   
    for( i=1; i<=(m-1)/2; i++)
		respns[n-i] = respns[m-i];

    for( i=(m+3)/2; i<=n-(m-1)/2; i++)
		respns[i-1]=0.0;

    two_fft0(data, respns, fft, ans, n);

    no2 = n >> 1;

    for( i=0; i<=n; i+=2) {
		if (isign == 1) {
			dum = ans[i];
			ans[i]	= (fft[i] * dum - fft[i+1] * ans[i+1]) / (double) no2;
			ans[i+1]= (fft[i+1] * dum + fft[i+1] * ans[i+1]) / (double) no2;
			} 
		else if (isign == - 1) {
			if ((mag2 = sqr(ans[i]) + sqr(ans[i+1])) == 0.0) 
				errexit("Deconvolving at response zero in convolve0.");
			dum = ans[i];
			ans[i]	= (fft[i] * dum + fft[i+1] * ans[i+1]) / mag2 / (double) no2;
			ans[i+1]= (fft[i+1] * dum - fft[i] * ans[i+1]) / mag2 / (double) no2;
			}
		}

	ans[1] = ans[n];
	real_fft0(ans, n, -1);
    
	free(fft);
	}

