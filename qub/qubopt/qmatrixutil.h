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

#ifndef QMATRIXUTIL_H
#define QMATRIXUTIL_H

#include "qubopt.h"

	class CChannelCombine { 
		private :
			int m_iChannels;
			int m_iStates;
			int m_iTemp;
			int m_iStateIndex;
			bool m_lContinue;
			int * m_aStates;
		public : 
			// Initialize and create the first combination
			CChannelCombine(int iChannels, int iStates);
			// Check if another combination exists 
			bool GetNextCombine();
			// Access the state elements with a 0 based index ( internally stored as 1 based array ) 
			int State0(int i);
			~CChannelCombine();
		};

	void CalcMultiPr(int nchannel, int nstate, int nmetastate, int **metastate, double *pr, double* mprx);
	void CalcMultiPr_mpl(int nchannel, int nstate, int nmutistate, int **mutistate, double *pr, double* mprx);
	double lnfactrl(int n);
	void metaMarkov(int nchannel, int nstate, int npath, int **path, int *class1, 
					float *ratio, int *nmetastate, int *nmetapath, int *nmetaclass, 
					float *metaamp, int **metastate, int **metapath, int *metaclass, 
					int **subpath);
	int freePar(int mc, int nx, double *x, double **Mtx, double *vct, int* nz, double *z, ostream &msgOut);
	int mil_ztox(int nz, int nx, double **Mtx, double *vct, double *z, double *x, double **dx_z, double *xmin, double *xmax);
	void qtometaq(int nstate, int npath, int **path, int nmetastate,
				  int nmetapath, int **metapath, int **metastate, int **subpath,
				  double **q, double **metaq, double **dmetaq_q);
	void aggregate(int nstate, int nclass, int *class1, int *ngroup, int **index);
	int qspctrm(int n, double **q, double* wr, double **s, double **s1, int &posieigencount, ostream &msgOut);
	void qtoqe(int nstate, int npath, int **path, int *class1, int nclass, double td, double **q, double **qe, double ***dqe_q);
	void mexp(int n, float t, float* q, float* a, int nz, float* dq, float* da, int i0) ;
	float gaussian(float mean, float xms, float x);
	void qtomutiq(int nstate, int npath, int* path, int nmutistate, 
				  int* mutistate, int nmutipath, int* mutipath, float* q, 
				  float* mutiq, int nz, float* dq, float* dmutiq);
	void xtoq(int nstate, int npath, int* path, float* drug, float volt, 
			  float* x, float* q, int nz, float* dx, float* dq);
	void eclass(float * aData, int iData, int & iClasses, float * aClasses, int * ia, double openDirection);


#endif
