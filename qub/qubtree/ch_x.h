#ifndef _CH_X_
#define _CH_X_


#include <stdio.h>
#include <list>
#include <matrix.h>
#include <string>


//
// type definitions
// 

typedef struct tag_mdlinf
{
   std::string mdlfile;

   fqvector<int> clazz; // class of state i		//**USED
   matrix<int> path;    // [from, to, drug] of path i		//**USED
   fqvector<double> x ; // [k0, k1] for each path		//**USED
   fqvector<double> pr;  // start prob of state i		//**USED

   int nstate;  // derived		//**USED
   int npath;   // derived // one-directional		//**USED
   //int nclazz;

   fqvector<double> xlimit[2];		//**USED
   matrix<double> mtx;		//**USED
   fqvector<double> vct;		//**USED
   fqvector<double> z;
   int nz;
   fqvector<double> gz;
   fqvector<double> var_z;
   fqvector<double> var_qa;
   fqvector<double> var_qb;
} mdlinf;

typedef struct
{
   fqvector<int> clazz; // class of state i	(formerly ics)	//**USED
   matrix<int> path;    // [from, to, fwd drug index, bkwd drug index] of path i //**USED
   matrix<double> x ; // [fwd u, fwd v, bkwd u, bkwd v] for each path	//**USED
   fqvector<double> pr;  // start prob of state i		//**USED

   int nclass;
   int nstate;  // derived		//**USED
   int npath;   // derived // one-directional		//**USED
   //int nclazz;
   fqvector<double> ptx;  // x start position of state i		//**USED
   fqvector<double> pty;  // y start position of state i		//**USED

   int nchannel;		//**USED

   fqvector<double> amp;       // single amplitudes		//**USED
   fqvector<double> xms;       // single standard deviatons 	//**USED
   fqvector<int> nar;       //**USED

   // int iconstraint;
   // int rconstraint;
   int nconstr;		//**USED
   fqvector<int> itype;	// constraint type 0:fix a  1:scale a  2:loop  3:fix b  4:scale b //**USED
   matrix<int> iconstr; // constraint state index //**USED
   fqvector<int> nloop;	// loop constraint size //**USED

} mdlfinf;


typedef struct tag_multmdl
{
   fqvector<int> clazz;
   matrix<int> path;
   matrix<int> state;
   matrix<int> cinf;
   int nclazz;
   int nstate;  // derived
} multmdl;


typedef struct tag_metamdl
{
   matrix<int> metastate;
   matrix<int> src;
   matrix<int> dst;
   
   int nmetastate;  // derived
   int nstate;      // temporary
   int nar;
   int nfir;
   int nclazz;
} metamdl;


class seginf 
{
public:
   fqvector<double> data;
   seginf() { }
};


class datinf 
{
public:
   std::string ldtfile;
   std::string parfile;

   fqvector<double> fir;
   int nfir;
   int nchannel;		//**USED
   double dt;
   
   fqvector<double> chold;
   double vhold;
   fqvector<double> c;
   double v;
   
   std::list<seginf*> lssi;

   multmdl mm;
   metamdl tm;

   fqvector<double> i;       // single amplitudes		//**USED
   matrix<double> r;       // single correlations		//**USED
   int ni;		//**USED
   int nr;		//**USED

   int iconstraint;
   int rconstraint;
   tensor<double> mtxr;
   matrix<double> vctr;
   matrix<double> mtxi;
   fqvector<double> vcti;
   fqvector<double> zr;
   fqvector<double> zi;
   int nzr;
   int nzi;
   fqvector<double> gzr;
   fqvector<double> gzi;
   fqvector<double> var_zr;
   fqvector<double> var_zi;
   fqvector<double> var_i;
   matrix<double> var_r;

   matrix<double> zai;     
   fqvector<double> tpi;            
   tensor<double> phi; 
   matrix<double> rho;      
   fqvector<double> sig;            
   double ll;
   datinf() { }
};


#ifdef CH_CODE_ENABLED


//
// Define API decoration for direct importing of DLL references.
//

#if defined(_WIN32)
#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)
#else
#define DLLEXPORT extern
#define DLLIMPORT
#endif

#if !defined(_CH_XLIB_)
#define CH_XAPI DLLIMPORT
#else
#define CH_XAPI DLLEXPORT
#endif
 
            // atoq.cpp

CH_XAPI void atoq (
            matrix<double>& a, 
            double dt, 
            matrix<double>& q);

            // freepar.cpp

CH_XAPI void freepar(
            matrix<double>& a, 
            fqvector<double>& b, 
            fqvector<double>& x, 
            fqvector<double>& z) ;

            // hmm.cpp

CH_XAPI double logl(
            metamdl& tm, 
            fqvector<double>& scale);

CH_XAPI void forward (
            metamdl& tm, 
            fqvector<double>& pr, 
            matrix<double>& a, 
            matrix<double>& b, 
            matrix<double>& alpha, 
            fqvector<double>& scale);

CH_XAPI void backward(
            metamdl& tm, 
            matrix<double>& a, 
            matrix<double>& b, 
            fqvector<double>& scale, 
            matrix<double>& beta);

CH_XAPI void gamma(
            metamdl& tm, 
            matrix<double>& alpha, 
            matrix<double>& beta, 
            fqvector<double>& scale);

CH_XAPI void vectb(
            metamdl& tm, 
            fqvector<double>& amp, 
            fqvector<double>& sd, 
            matrix<double>& ar, 
            fqvector<double>& fir, 
            fqvector<double>& data, 
            matrix<double>& b )  ;

CH_XAPI void dllda(
            metamdl& tm, 
            matrix<double>& alpha, 
            matrix<double>& beta, 
            matrix<double>& b, 
            matrix<double>& ga) ;

CH_XAPI void dlldp(
            metamdl& tm, 
            matrix<double>& beta, 
            matrix<double>& b, 
            fqvector<double>& gp) ;

CH_XAPI void dlldi(
            metamdl& tm, 
            fqvector<double>& amp, 
            fqvector<double>& sd, 
            matrix<double>& ar, 
            fqvector<double>& fir, 
            fqvector<double>& data, 
            matrix<double>& gama, 
            fqvector<double>& gi) ;

CH_XAPI void dlldn(
            metamdl& tm, 
            fqvector<double>& amp, 
            fqvector<double>& sd, 
            matrix<double>& ar, 
            fqvector<double>& fir, 
            fqvector<double>& data, 
            matrix<double>& gama, 
            fqvector<double>& gsd, 
            matrix<double>& gar) ;


#ifdef __cplusplus
extern "C" {
#endif
            // inpcns.cpp

//CH_XAPI void inpcns (mdlinf& mi);

            // inpldt.cpp

CH_XAPI void inpldt( std::list<std::string>& lsldt,
                    std::list<datinf*>& lsdi ) ;

CH_XAPI void rd_ldt_file ( std::string& filename, datinf& di);

            // inpmdl.cpp

//CH_XAPI void inpmdl (std::string& mdlfile, mdlinf& mi);

            // inpar.cpp

CH_XAPI void inppar (std::list<datinf*>& lsdi) ;
CH_XAPI void rd_par_file (std::string& filename, datinf& di);

            // ldtutil.cpp

CH_XAPI void rd_file_header (
            FILE* fp, 
            int& nlead, 
            double& dt, 
            int& count);

CH_XAPI bool rd_seg_info(
            FILE *fp, 
            int& tstart, 
            int& ndata) ;

CH_XAPI void rd_seg_data(
            FILE *fp, 
            int count, 
            int ndata, 
            double* data) ;

CH_XAPI void wt_file_header(
            FILE* fp, 
            int nlead, 
            double dt, 
            int count) ;
 
CH_XAPI void wt_segment(
            FILE *fp, 
            int tstart, 
            int ndata, 
            fqvector<double>& data, 
            int count) ;

#ifdef __cplusplus
}
#endif

            // metamakv,cpp

CH_XAPI void metamakv (
            fqvector<int>& clazz, 
            int nar, 
            int nfir, 
            metamdl& tm);


            // eclass.cpp

CH_XAPI void eclass(
            fqvector<double>& a, 
            fqvector<double>& c, 
            fqvector<int>& ia); 

CH_XAPI void eclass(
            fqvector<float>& a, 
            fqvector<float>& c, 
            fqvector<int>& ia);

CH_XAPI void eclass(
            fqvector<int>& a, 
            fqvector<int>& c, 
            fqvector<int>& ia) ; 

            // initp.cpp

CH_XAPI void ptomp(
            multmdl& mm, 
            fqvector<double>& pr, 
            fqvector<double>& mpr);

CH_XAPI void qtope(
            matrix<double>& q, 
            fqvector<double>& pe, 
            tensor<double>& dq, 
            matrix<double>& dpe);

CH_XAPI void qtope (matrix<double>& q, fqvector<double>& pe);

CH_XAPI void inimpr(
            mdlinf& mi, 
            multmdl& mm, 
            fqvector<double>& c, 
            double v, 
            fqvector<double>& x, 
            fqvector<double>& mpr, 
            matrix<double>& dx, 
            matrix<double>& dmpr);

CH_XAPI void inimpr(
            mdlinf& mi, 
            multmdl& mm, 
            fqvector<double>& c, 
            double v, 
            fqvector<double>& x, 
            fqvector<double>& mpr);

CH_XAPI void ptotp(
            metamdl& tm, 
            fqvector<double>& pr, 
            matrix<double>& a, 
            fqvector<int>& clazz,
            fqvector<double>& tpr, 
            matrix<double>& dpr, 
            tensor<double>& da, 
            matrix<double>& dtpr);

CH_XAPI void ptotp(
            metamdl& tm, 
            fqvector<double>& pr, 
            matrix<double>& a, 
            fqvector<int>& clazz,
            fqvector<double>& tpr);

            // multmakv.cpp

CH_XAPI void mstate (int nstate, int nchannel, multmdl& mm);
CH_XAPI void mclass (fqvector<int>& clazz, multmdl& mm);
CH_XAPI bool mpath (multmdl& mm, int I, int J, int& i, int& j);
CH_XAPI void mpath (matrix<int>& path, multmdl& mm);

            // ztoar.cpp

CH_XAPI void ztoi(
            matrix<double>& mtx, 
            fqvector<double>& vct, 
            fqvector<double>& z, 
            fqvector<double>& i, 
            matrix<double>& di);

CH_XAPI void ztoi(
            matrix<double>& mtx, 
            fqvector<double>& vct, 
            fqvector<double>& z, 
            fqvector<double>& i);

CH_XAPI void itomi(
            multmdl& mm, 
            fqvector<double>& u, 
            fqvector<double>& amp, 
            matrix<double>& du, 
            matrix<double>& damp);

CH_XAPI void itomi(
            multmdl& mm, 
            fqvector<double>& u, 
            fqvector<double>& amp);

CH_XAPI void ztor(
            int rconstraint, 
            fqvector<double>& z, 
            matrix<double>& r, 
            tensor<double>& dr);

CH_XAPI void ztor(
            int rconstraint, 
            fqvector<double>& z, 
            matrix<double>& r);

CH_XAPI void rtomr(
            multmdl& mm, 
            matrix<double>& r, 
            matrix<double>& mr, 
            tensor<double>& dr, 
            tensor<double>& dmr);

CH_XAPI void rtomr (
            multmdl& mm, 
            matrix<double>& r, 
            matrix<double>& mr);

CH_XAPI int rtoar (
            matrix<double>& r, 
            fqvector<double>& sd, 
            matrix<double>& ar);

CH_XAPI int rtoar (
            matrix<double>& r, 
            fqvector<double>& sd, 
            matrix<double>& ar, 
            tensor<double>& dr, 
            matrix<double>& dsd, 
            tensor<double>& dar);

CH_XAPI int rtoar (
            fqvector<double>& r, 
            double& sd, 
            fqvector<double>& ar, 
            matrix<double>& dr,
            fqvector<double>& dsd, 
            matrix<double>& dar);

            // ztoq.cpp

CH_XAPI void ztox (
            matrix<double>& a, 
            fqvector<double>& b, 
            fqvector<double>& z, 
            fqvector<double>& x, 
            matrix<double>& dx);

CH_XAPI void xtoq (
            matrix<int>& path, 
            fqvector<double>& c, 
            double v, 
            fqvector<double>& x, 
            matrix<double>& q, 
            matrix<double>& dx, 
            tensor<double>& dq) ;

CH_XAPI void qtomq (
            multmdl& mm, 
            matrix<double>& q, 
            matrix<double>& qq, 
            tensor<double>& dq, 
            tensor<double>& dqq);

CH_XAPI void ztox (
            matrix<double>& a, 
            fqvector<double>& b, 
            fqvector<double>& z, 
            fqvector<double>& x);

CH_XAPI void xtoq (
            matrix<int>& path, 
            fqvector<double>& c, 
            double v, 
            fqvector<double>& x, 
            matrix<double>& q) ;

CH_XAPI void qtomq (
            multmdl& mm, 
            matrix<double>& q, 
            matrix<double>& qq);

            // rtoma.cpp

CH_XAPI int rtoma (fqvector<double>& r, fqvector<double>& a); 
CH_XAPI int rtoma (matrix<double>& r, matrix<double>& ma);
 


#ifdef __cplusplus
extern "C" {
#endif
            // ranpack

CH_XAPI double rand1 (long& idum) ; 

            // minpack

CH_XAPI void frprmn ( 
            fqvector<double>& p, 
            double ftol, 
            double stpmax, 
            int& iter, 
            double& fret, 
            void* pvoid,
            double (*func)(fqvector<double>&, void*), 
            void (*dfunc)(fqvector<double>&, double&, fqvector<double>&, void*));

CH_XAPI void linmin (
            fqvector<double>& p, 
            fqvector<double>& xi, 
            int n, 
            double& fret, 
            double (*func)(fqvector<double>&, void*), 
            void* pvoid);

CH_XAPI void mnbrak (
            double& ax, double& bx, double& cx, 
            double& fa, double& fb, double& fc, 
            double (*func)(double, void*), 
            void* pvoid);

CH_XAPI void brent (
            double ax, 
            double bx, 
            double cx, 
            double tol, 
            double& xmin, 
            double& fmin, 
            double (*f)(double, void*), 
            void* pvoid);

CH_XAPI int dfpmin (
            fqvector<double>& p, 
            double& fret, 
            matrix<double>& hessin, 
            int& iter, 
            double stpmax, 
            double tolf, 
            double tolg, 
            int itmax,
            int (*func)(fqvector<double>&, double&, fqvector<double>&, void*), 
            void (*oupt)(int, fqvector<double>&, double, fqvector<double>&, void*), 
            void* pvoid) ;

CH_XAPI int lnsrch(
            fqvector<double>& pold, 
            double fold, 
            fqvector<double>& gold, 
            fqvector<double>& xi, 
            fqvector<double>& p, 
            double& f, 
            fqvector<double>& g, 
            int (*func)(fqvector<double>&, double&, fqvector<double>&, void*), 
            void* pvoid);

            // linpack

CH_XAPI void svdc( 
            matrix<double>& a, 
            fqvector<double>& w, 
            matrix<double>& v);
             
CH_XAPI void gaussj(
            matrix<double>& a, 
            fqvector<double>& b);
             
CH_XAPI void invert(
            matrix<double>& a);
 
CH_XAPI void balanc(
            matrix<double>& a, 
            int& low, 
            int& igh, 
            fqvector<double>& scale);

CH_XAPI void balbak ( 
            int n, int m, 
            int low, int igh, 
            fqvector<double>& scale, 
            matrix<double>& z);

CH_XAPI void elmhes(
            matrix<double>& a, 
            int low, int igh, 
            fqvector<int>& inf);

CH_XAPI void eltran (
            matrix<double>& a, 
            int low, 
            int igh, 
            fqvector<int>& inf, 
            matrix<double>& z);

CH_XAPI void mspec(
            matrix<double>& q, 
            fqvector<double>& wr, 
            fqvector<double>& wi, 
            matrix<double>& z, 
            matrix<double>& v);
           
            // others

CH_XAPI int zrhqr(fqvector<double>& a, fqvector<double>& rtr, fqvector<double>& rti);
CH_XAPI int rt2cof(fqvector<double>& rtr, fqvector<double>& rti, fqvector<double>& a);

#ifdef __cplusplus
}
#endif

             // Fourier transform, convolution

                   
            // filter design

//extern "C" void gsfir (double fc, int& nh, fqvector<double>& h); 
//extern "C" void wrap (int nh, fqvector<double>& h);
            
            // linear programming
/*
void simplx (
     double** a, int m1, int m2, int m3, int n, 
     int iposv[], int izrov[], std::pair<int,int>& iter,
     const double eps//=1.0e-8
	 );
*/
 
            // the following require c++ compilation
 
#ifndef __cplusplus
#error ERROR: C++ compiling needed
#endif

CH_XAPI int hqr (
            matrix<double>& h, 
            int low, 
            int igh, 
            fqvector<double>& lamdar, 
            fqvector<double>& lamdai);

CH_XAPI int hqr (
            matrix<double>& h, 
            int low, 
            int igh, 
            fqvector<double>& lamdar, 
            fqvector<double>& lamdai, 
            matrix<double>& z);

CH_XAPI void eigen( 
            matrix<double>& a, 
            fqvector<double>& wr, 
            fqvector<double>& wi);

CH_XAPI void eigen(
            matrix<double>& a, 
            fqvector<double>& wr, 
            fqvector<double>& wi, 
            matrix<double>& z);

CH_XAPI void mexp(
            matrix<double>& q, 
            double t, 
            matrix<double>& a, 
            int d);

CH_XAPI void mexp(
            matrix<double>& q, 
            double t, 
            matrix<double>& a, 
            int d, 
            tensor<double>& dq, 
            tensor<double>& da);

#endif /* CH_CODE_ENABLED */

#endif /*_CH_X_*/















