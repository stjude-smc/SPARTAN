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

#include "histpdf.h"
using namespace std;

// reworked to handle multiple segments

// now the dwells can have the dead time already subtracted.
// i.e. they can be mil's working arrays


const float ZERO_SINGLE=float(1.0e-8);

void binning(float ts, int iSegs, int* ndwell, float** tdwell, int* nbin, float* bin, float tdead, double& rOut){
	int i, j;//, n;
	float tmin, tmax, dx, r, tdw;
	
	tmin=tmax=tdwell[0][0];
	for (i=0; i<iSegs; i++){
		for ( j=0; j<ndwell[i]; j++ ) {
			tdw=tdwell[i][j];
			tmin = (tmin < tdw) ? tmin : tdw;
			tmax = (tmax > tdw) ? tmax : tdw;
			}
		}

	tmin += tdead;
	tmax += tdead;
	
	dx = float(log10(tmax/tmin)/(float)(*nbin-1));   
	r = float(pow(10.0f, dx));
	r = (r>1+ZERO_SINGLE) ? r : float(1.01);
	rOut = r;

	bin[0] = tmin;
	for (i=1; i<*nbin; i++) 
		bin[i] = bin[i-1]*r;
	
	if (ts > ZERO_SINGLE) 
		for (i=0; i<*nbin; i++)
		  bin[i] = (float) (ts*floor(bin[i]/ts+0.5));
		
	for (i=1; i<*nbin; i++) 
		while( fabs(bin[i]-bin[i-1]) < ZERO_SINGLE ) {
			for (j=i; j<*nbin-1; j++) 
				bin[j] = bin[j+1];
			(*nbin)--;
			}
	}


void calhist(float ts, int ic, int iSegs, int* ndwell, int **idwell, float** tdwell,
			 int iBins, float* bin, float* hst, float tdead) {
	int iDwell, iBin, iSeg;
	int iPrevBin, iCurrBin, iValue; 
	float nPrevBin, nCurrBin, nValue;
	float halfts = float(0.5 * ts);
      
	for( iBin=0; iBin<iBins; iBin++ ) 
		hst[iBin] = 0.0;
	
	if ( ts > ZERO_SINGLE ) 
		for ( iSeg=0; iSeg<iSegs; iSeg++ ){
			for (iDwell=0; iDwell<ndwell[iSeg]; iDwell++) {
				iValue = int((tdwell[iSeg][iDwell] + tdead + halfts)/ts);
				iPrevBin = 0;
				if ( idwell[iSeg][iDwell] == ic ){
					for (iBin=0; iBin<iBins; iBin++) {
						iCurrBin = int(bin[iBin]/ts);
						if (iValue>iPrevBin && iValue<=iCurrBin) 
							hst[iBin]++;
						iPrevBin = iCurrBin;
						}
					}
				}
			}
	else 
		for ( iSeg=0; iSeg<iSegs; iSeg++ ){
			for (iDwell=0; iDwell<ndwell[iSeg]; iDwell++) {
				if ( idwell[iSeg][iDwell] == ic ){
  					nValue = tdwell[iSeg][iDwell] + tdead;
  					nPrevBin = 0.0;
					for (iBin=0; iBin<iBins; iBin++) { 
						nCurrBin = bin[iBin];
						if ( nValue>nPrevBin && nValue<=nCurrBin ) 
							hst[iBin]++; 
						nPrevBin = nCurrBin;
						} 
					}
				}
			}
		/* Better code if we can assume bin[n]<bin[n+1]
		for ( iSeg=0; iSeg<iSegs; iSeg++ ){
			for (iDwell=0; iDwell<ndwell[iSeg]; iDwell++) {
				if ( idwell[iSeg][iDwell] == ic ){
	  				nValue = tdwell[iSeg][iDwell] + tdead;
					for (iBin=0; iBin<iBins; iBin++) { 
						if( nValue<=bin[iBin]  ) {	//    || iBin==iBins-1
							hst[iBin]++; 
							break;
							}
						} 
					}
				}
			}
		*/
	}

void intersect(double a0, double a1, double b0, double b1, double& c0, double& c1)
{
	c0 = max(a0, b0);
	c1 = min(a1, b1);
}
/* oops this is only valid for 0<x<1 and maybe not precise enough

altho...the biggest factor appears to be "floating point consistency" in the optimizations.
apparently when you operate on different length float types the least significant bits
get carried along with incorrect results.  this code (even the NR replacement) really brings it out.

int factorial(int x)
{
	if ( x <= 0 )
		return 1;
	return x * factorial(x - 1);
}

double erf(double x)
{
	const double eps = 1e-12;
	long double term;
	long double sum = 0.0;
	int i = 0;
	do {
		term = pow(-1, i) * pow(x, 2*i+1) / ((2*i+1) * factorial(i));
		sum += term;
		++i;
	}
	while (term > eps);
	return (sum * 2.0 / sqrt(PI));
}
*/

/* 
Numerical Recipes p221:

... you might wish to consider the following
routine, based on Chebyshev fitting to an inspired guess as to the functional form:

double erfcc(double x)
// Returns the complementary error function erfc(x) with fractional error everywhere less than
// 1.2 × 10-7.
{
double t,z,ans;
z=fabs(x);
t=1.0/(1.0+0.5*z);
ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
t*(-0.82215223+t*0.17087277)))))))));
return x >= 0.0 ? ans : 2.0-ans;
}

double erf(double x)
{
	return (1.0f - erfcc(x));
}

double gauss_area(double mean, double sd, double from, double to)
{
	return (0.5 * (erf((to-mean)/(sd*sqrt(2.0))) - erf((from-mean)/(sd*sqrt(2.0)))));
}
*/
// referring to a triangular distribution, from (peak-1, 0) to (peak, 1) to (peak+1, 0)

double tri_rising_y(double peak, double x)
{
	return (x - peak + 1);
}

double tri_falling_y(double peak, double x)
{
	return - (x - peak - 1);
}

double tri_rising_area(double peak, double from, double to)
{
	double x0, x1;
	intersect(peak-1, peak, from, to, x0, x1);
	if ( x0 < x1 ) {
		double y0 = tri_rising_y(peak, x0);
		double y1 = tri_rising_y(peak, x1);
		return (0.5*(x1-x0)*(y0+y1));
	}
	else
		return 0.0;
}

double tri_falling_area(double peak, double from, double to)
{
	double x0, x1;
	intersect(peak, peak+1, from, to, x0, x1);
	if ( x0 < x1 ) {
		double y0 = tri_falling_y(peak, x0);
		double y1 = tri_falling_y(peak, x1);
		return (0.5*(x1-x0)*(y0+y1));
	}
	else
		return 0.0;
}

double tri_area(double peak, double from, double to)
{
	return tri_rising_area(peak, from, to) + tri_falling_area(peak, from, to);
}

#define WAY_TOO_MANY_SAMPLES 256000

void spread_events(float *bin, float *hst, int iBins, int samples, int count, float sampling, float r) {
	double a = bin[0] / r;
	double b; // , aa, bb;
	for ( int bi=0; bi<iBins; ++bi ) {
		b = bin[bi];
		//intersect(sam-1.0, sam+1.0, a/sampling, b/sampling, aa, bb);
		//if ( aa < bb ) {
			double area = tri_area(samples, a/sampling, b/sampling); // tri_area(sam, aa, bb); // gauss_area(sam, 0.5, aa, bb);
			hst[bi] += float(count * area);
		//}
		a = b;
	}
}

void calhist_smooth(int ic, int iSegs, int* ndwell, int **idwell, float** tdwell,
			 int iBins, float* bin, float* hst, float tdead, float sampling, float r) {
	vector<int> sampleCounts;
	vector<int> veryLongInSamples;

	// first bin the events by how many samples they were recorded as
	// fixed: treat extraordinarily long events separately, so we don't need 200 million empty bins
	for ( int iSeg=0; iSeg<iSegs; iSeg++ ) {
		for ( int iDwell=0; iDwell<ndwell[iSeg]; iDwell++ ) {
			if ( idwell[iSeg][iDwell] == ic ) {
				int samples = int(floor((tdwell[iSeg][iDwell]+tdead)/sampling));
				if ( samples > WAY_TOO_MANY_SAMPLES ) {
					veryLongInSamples.push_back(samples);
					continue;
				}
				if ( ((int) sampleCounts.size()) <= samples )
					sampleCounts.resize( samples+1 );
				sampleCounts[samples]++;
			}
		}
	}
	// then distribute those binned events among the log bins
	// using the empirical "fact" that events recorded as k samples long are actually between (k-1)dt and (k+1)dt,
	// distributed normally with mean of k*dt and std dev of dt/2
	// 
	// oops i think it is better approximated by a simple triangle
	for ( int bi=0; bi<iBins; ++bi )
		hst[bi] = 0.0;
	for ( size_t sam=1; sam<sampleCounts.size(); ++sam )
		spread_events(bin, hst, iBins, (int) sam, sampleCounts[sam], sampling, r);
	// now handle those set-aside way-too-long events
	for ( vector<int>::iterator si = veryLongInSamples.begin(); si != veryLongInSamples.end(); ++si )
		spread_events(bin, hst, iBins, *si, 1, sampling, r);
}
