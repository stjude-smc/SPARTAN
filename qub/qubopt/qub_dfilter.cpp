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

//Digital filtering routines

#ifdef _WIN32
  #include <windows.h>
  #include <time.h>
  #include <process.h>
  #define sleep Sleep
#else
  #include <stdlib.h>
  #include <unistd.h>
#endif

#include "qub_dfilter.h"
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;

//#define Debug

template<class Precision>
struct Filter_t {
	int Order;
	Precision * An;
	Precision * Bn;
	Precision * Dn;
	Precision * Di;
};

//This is Form II Transform - Potentially better computational behavior.

template<class Precision>
Precision DFIIT(Precision X, Filter_t<Precision> * F){

/*Arrays are dimensioned from 0 to order. The loop runs from 1 to order and
uses the [0] index last.*/

    int Order = F->Order;                       //Referenced every time through loop
    Precision *Dn = F->Dn;                      //Ditto
    Precision *Ai, *Bi, *Dk=F->Di;              //Define Result and pointers to indexed An, Bn & Dk

//First calculate the output value since this is used to update the delay line. This uses only 
//An[0], Bn[0] and the oldest Delay line entry (Pointed to by Dk).

    Precision Xf = (*F->An * X + *Dk) / *F->Bn;
    *Dk = 0.0;                                  //Zero for reuse by newest entry

//Loop to update the delays. Note that we calculate the value to replace the oldest (nth) delay from the
//sum of n+1th delay, X*An[n+1] and -Xf*B[n+1]. For the purposes of the loop, we replace the n+1th delay
//with this result. When the loop finishes, we "rotate" the delay ring buffer by one step to put all the
//delays in their proper places. There is a **KLUDGE** here: Note that the circularization of the delay
//line detects the pointer being less than the start of the array. You might think that we would want to
//add back "1" then add "Order". However, the Dn array is one shorter than the An and Bn arrays, so adding
//back order is just perfect.

    for (Ai = F->An + 1, Bi = F->Bn + 1; Ai <= F->An + Order; ++Ai, ++Bi){
        if (--Dk < Dn) Dk = Dk + Order;	              //Decrement Dk and circularize the delay line
        *Dk = *Dk + *Ai * X - *Bi * Xf;	              //Compute new delay value
    }

//When we exit the loop, we will have replaced the Dk that we used to calculate Xf with the newest
//(the last computation in the loop). Since we want the saved index (F->Di) to point to the oldest
//entry, we need to decrement Dk one last time to roll tht ring buffer by that last step. 

    if (--Dk < Dn) Dk = Dk + Order;	                //Decrement Dk and circularize the delay line
    F->Di = Dk;                                     //Then update F->Di for next time.
    return Xf;                                      //Return Filtered value.
}

//This is Form II.

template<class Precision>
Precision DFII(Precision X, Filter_t<Precision> * F){

/*Arrays are dimensioned from 0 to order. The loop runs from 1 to order and
uses the [0] index last.*/

    Precision Xf=0.0, *Ai, *Bi, *Dk=F->Di;      //Define Result and pointers to indexed An, Bn & Dk
    Precision *Dkp1;                            //For post-incrementing the Delay Line Pointer
    Precision Di=X;                             //Set Local Di to incoming X value
    int Order = F->Order;                       //Referenced every time through loop
    Precision *Dn = F->Dn;                      //Ditto

//Note the use of pointer arithmetic here. Ai=An+1 sets Ai to point at An(1). Ditto for Bi...
//Note then that ++Ai and ++Bi add enough to the pointer to get it pointing to the next array
//element. (We don't actually know how much this is - it depends upon the definition of "Precision")

    for (Ai = F->An + 1, Bi = F->Bn + 1; Ai <= F->An + Order; ++Ai, ++Bi){
        Di = Di - *Bi * *Dk;                    //Accumulate the delay... and
        Xf = Xf + *Ai * *Dk;                    //the output terms
        Dkp1=Dk;                                //Save the pointer before decrementing
        if (--Dk < Dn) Dk = Dk + Order;         //Decrement Dk and circularize the delay line
    }

//When we exit the loop, Dk will be the same as it was when we entered, which is pointing to one past
//the oldest entry, which we want to replace. So we undo the last decrement before storing the new delay.
//(Really we just use Dkp1 instead of Dk.)

    Di = Di / *F->Bn;                           //Scale the newest delay by Bn(0)
    Xf = Xf + *F->An * Di;                      //Then finish output sum by adding Di*An(0)
    *Dkp1 = Di;                                 //And finally set Dn[Index] to the Local Di
    F->Di = Dkp1;                               //Then update F->Di for next time.
    return Xf;                                  //Return Filtered value.
}

/*
Analyze a filter configuration file for the coefficients. File are in MatLab format. Strategy is to
build a new filter data structure (array of An, Bn, and Dn) using data from the input file.
Coefficients are initially stored in data vectors so they may be sorted and corrected for missing
poles or zeros. Once the coefficient data has been been fully constructed, then copy into conventional
arrays and assign the array pointers into the passed in "F" parameter. In version 3.0.188, the ability
to initialize the delayline (read in eigenvalues) was added. This should allow the program to start out
the filter with some initial condition ("IC"). This feature is experimental and assumes that the initial delay
line values are based on the steady state response of the filter to a unity input.
*/

template<class Precision>
bool DFII_Init(string FileName, Filter_t<Precision> * F, Precision IC){

//	cerr << "Entering DFII_Init with Filename of " << FileName << endl;

//First see if the file name is empty or --None--. This is a simple error check for some lame
//programmer who might just call this routine with a bad (or no) filter filename.

    if(FileName.empty())return false;

//So, we may actually have a valid file. Filter files must all be located in the common "Filters"
//directory, so file names were passed in "undecorated"and the path had to be added. As of version
//3.0.0.91, the file name is passed from client with the full path. There may be an issue here if
//we want true remote operation, since then the path to the filters must match exactly on both client
//and server. One solution would be to send entire contents of filter file in the message from the
//For now, we assume that the path structure will always match.client.

//Open the file, read only.

    ifstream FilterFile(FileName.c_str(),ios::in);
    if (FilterFile.fail()) return false;                      //No such file in Filters
    char TxtLine[128];                                        //Someplace to put input line
    vector<pair<int, double> > An;                            //Temp storage vector for An
    vector<pair<int, double> > Bn;                            //Temp storage vector for Bn
    vector<pair<int, double> > Dn;                            //Optional EigenValues for delay line
    int AMax=-1, BMax=-1, DMax=-1;                            //Input item Max indexes (for format error check)
    double Kn=1.0, Kd=1.0;                                    //N and D of gain term. Default is 1.0 if none found

    while (FilterFile.getline(TxtLine,128)){                  //Read an input line and check EOF
        if (FilterFile.gcount()<6) continue;                  //Ignore short lines: they may kill following logic
        char iABC=(char) (toupper(TxtLine[0])-'A');                    //Get 1st character as 'A'=0, 'B'=1, etc.
//      if (iABC<0 || iABC>3 || iABC==2) continue;            //Ignore lines not starting with A, B or D (>3.0.187)
   	    if (iABC<0 || iABC>3) continue;                       //Ignore if not starting with A, B, C or D (>3.0.191)

//Probably we should ignore all spaces to make the file format more robust.

        if (TxtLine[1]!='(') continue;                        //Ignore lines without '(' after A, B, D
        int j=atoi(&TxtLine[2])-1;                            //Get (n) & fixup "1 based" array index
        if (j<0) continue;                                    //No number there? Then ignore line.
        char * nTxt=strchr(TxtLine,'=');                      //Get pointer to '=' character
        if(nTxt==NULL) continue;                              //There is no equals character: Ignore line
        double C=atof(&nTxt[1]);                              //Get An[j], Bn[j], Dn[j] or C[0]
        switch (iABC) {                                       //Decide which and then store appropriately
        case 0:
            An.push_back(pair<int,double>(j,C));              //Store An[j]
            AMax=max(AMax,j);                                 //Keep track of largest index
            break;
        case 1:
            Bn.push_back(pair<int,double>(j,C));              //Store Bn[j]
            BMax=max(BMax,j);                                 //Keep track of largest index
            break;
        case 2:
            if (j==0){
                Kn=C;                                         //Save gain if found. Only accept C(1) as
                cerr << "Setting Gain Numerator to "<< Kn <<endl;
            }
            if (j==1){
                Kd=C;                                         //Numerator and C(2) as denominator
                cerr << "Setting Gain Denominator to "<< Kd <<endl;
            }
            break;
        case 3:
            Dn.push_back(pair<int,double>(j,IC*C));           //Store IC*Dn[j] (>3.0.187)
            DMax=max(DMax,j);                                 //Keep track of largest index
            break;
        }
    } //end while (FilterFile.getline(TxtLine,128))

    FilterFile.close();                                       //Close the file now
    sort(An.begin(), An.end());                               //put An's in order
    sort(Bn.begin(), Bn.end());                               //put Bn's in order
    sort(Dn.begin(), Dn.end());                               //put Dn's in order (>3.0.187)

/*Note that the filter's order is one less than the number of coefficients. An[0] and Bn[0] are actually
the overal filter gain terms (used at the end of the filter loop), i.e.,t he loop runs from 1 to Order - see
DF2 code. For this reason, we subtract 1 from the vector size (which is the number of items) to get the order.*/

//A bug was discovered in the 3.0.187 and previous implementations when filters did not have an equal number of
//poles and zero's. In this case, the sizes of An and Bn will not be equal. Since the assumption was made that
//they would always be equal, the old code set the "missing" coefficients to some random value (whatever was in
//the memory allocated to the vector). In some cases, this could have produced a seg-fault if the extra space
//allocated for the vector was not enough to hold the missing elements. We must test for this and set
//the array elements for the missing coefficients to zero.

    int ASize=(int)An.size();                        //Get size of the An vector (>3.0.187)
    int BSize=(int)Bn.size();                        //And the size of the Bn vector (>3.0.187)
    int DSize=(int)Dn.size();                        //And the size of the Eigenvalue array (>3.0.187)
    if (AMax+1!=ASize || BMax+1!=BSize || DMax+1!=DSize) {
        cerr << "Error in " << FileName << ": Non-contiguous or non-zero based parameter array definition found for a(n), b(n) or d(n)!\n";
        return false;
    }

//Size needed for Coef. arrays is the greater of ASize or BSize. The filter order is exactly one less than this.

    int TS=max(ASize,BSize);                    //This will be the size of each coef array (>3.0.187)
    int Order=TS-1;	                            //Define the order
    if (DSize>Order) cerr << "Warning: "<< FileName << " contains too many elements in the delay line initializer.\nExtra elements will be ignored.\n";

//The copy operation puts the data into one appropriate array and then passes the pointers to each part back.
//Note that if either i, j, or k exceeds <A|B|D>Size, then we need to assign "0.0" since that element does not exist.

    Precision *tAn=new Precision [TS];                  //Define new An
    Precision *tBn=new Precision [TS];                  //New Bn

//The delay line is complicated by the possibility that only some delay values may be given, in which case
//we want to propagate the last one to all the remaining delays. The easiest way to do this is to initialize
//the entire array to the last Dn value. Then if any others are present, they will get replaced in the loop
//below. Also note that the delay line length is equal to the filter order, not the number of coefficients, so
//we have to prevent assignments of Dn past the filter order.

    Precision *tDn=new Precision [Order];
    Precision DFault=(Precision) ((DSize>0) ? Dn[DSize-1].second:0.0);

/*In version 192, Net filter gain was added to allow precise correction of filter's DC gain to be exactly
 0 dB. The gain factor appears in the filter file as C(0) and C(1), where the gain, K, is C(1)/C(2). We
 include this in the filter's response by dividing all the numerator coefficients by the K. If either C(1)
 or C(2) is not found in the file, then they are defaulted to 1.0 so the value of K is will be
 defined as C(1) or 1/C(2), respectively. If both are missing K will be 1.0, which has no effect.*/

    double iGain=Kd/Kn;
//    cerr << "Gain Correction is " << iGain << endl;
    for (int i=0; i<TS; ++i){
        tAn[i]=(Precision) ((i<ASize) ? iGain*An[i].second:0.0);            //Copy ith An. If i=>ASize, set to 0.0 (>3.0.187)
        tBn[i]=(Precision) ((i<BSize) ? Bn[i].second:0.0);                  //Copy ith Bn. If i=>BSize, set to 0.0 (>3.0.187)
        if (i<Order) {                                        //Only do Dn[i] for i<Order!!!
            tDn[i]=(Precision) ((i<DSize) ? Dn[i].second:DFault);           //Copy ith Dn. If i=>Order, Use DFault
        }
    }

/* At long last, assign the data into the filter data structure with a pointer set. We need to acquire
a mutex on the filter to make sure that we don't change the coefficients while the filter is running.
It is then safe to delete the old arrays and assign the new ones. Always set the Delay Line index to
Dn[0] */

    delete [] F->An; F->An=tAn;	                              //Set up pointers
    delete [] F->Bn; F->Bn=tBn;
    delete [] F->Dn; F->Dn=tDn;
    F->Di=F->Dn;
    F->Order=Order;                                           //And also the filter order.
    return true;
}


// extern "C" wrapper, exported in qub_dfilter.h

extern "C" QUBOPT_API Filter32_p DFII32_New()
{
  Filter_t<float> *F = new Filter_t<float>;
  F->An = F->Bn = F->Dn = F->Di = NULL;
  F->Order = 0;
  return (Filter32_p) F;
}

extern "C" QUBOPT_API void DFII32_Free(Filter32_p F)
{
  if ( ! F ) return;
  Filter_t<float> *FF = (Filter_t<float>*) F;
  delete [] FF->An;
  delete [] FF->Bn;
  delete [] FF->Dn;
  delete FF;  
}

extern "C" QUBOPT_API int DFII32_Init(const char * FileName, Filter32_p F, float IC) //  IC = 0.0
{
  return DFII_Init(FileName, (Filter_t<float>*)F, IC);
}

extern "C" QUBOPT_API float DFII32(float X, Filter32_p F)
{
  return DFII(X, (Filter_t<float>*)F);
}

extern "C" QUBOPT_API float DFIIT32(float X, Filter32_p F)
{
  return DFIIT(X, (Filter_t<float>*)F);  
}


extern "C" QUBOPT_API Filter32_p DFII64_New()
{
  Filter_t<double> *F = new Filter_t<double>;
  F->An = F->Bn = F->Dn = F->Di = NULL;
  F->Order = 0;
  return (Filter64_p) F;
}

extern "C" QUBOPT_API void DFII64_Free(Filter32_p F)
{
  if ( ! F ) return;
  Filter_t<double> *FF = (Filter_t<double>*) F;
  delete [] FF->An;
  delete [] FF->Bn;
  delete [] FF->Dn;
  delete FF;
}

extern "C" QUBOPT_API int DFII64_Init(const char * FileName, Filter64_p F, double IC) // IC = 0.0
{
  return DFII_Init(FileName, (Filter_t<double>*)F, IC);
}

extern "C" QUBOPT_API double DFII64(double X, Filter64_p F)
{
  return DFII(X, (Filter_t<double>*)F);
}

extern "C" QUBOPT_API double DFIIT64(double X, Filter64_p F)
{
  return DFIIT(X, (Filter_t<double>*)F);
}

