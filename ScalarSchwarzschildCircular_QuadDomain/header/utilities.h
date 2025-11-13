#include <stdio.h>
#include <stdlib.h>
//  #include <stdio_ext.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>
#include <dirent.h>
#include <errno.h>
#include <sys/utsname.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define Pi       3.14159265358979323846264338328
#define Pih      1.57079632679489661923132169164 // Pi/2
#define Piq      0.78539816339744830961566084582 // Pi/4
#define Third    0.33333333333333333333333333333 // 1/3
#define TwoThird 0.66666666666666666666666666667 // 2/3
#define Sqrt_2   1.41421356237309504880168872421 // Sqrt[2]

#define ERRTOL 0.05
#define TINY 1.0e-20
#define BIG  1.0e+10
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define NR_END 1
#define FREE_ARG char*
#define NMAX 200

typedef struct DCOMPLEX {double r,i;} dcomplex;

void nrerror(char error_text[]);
int      *ivector( int nl,  int nh);
double   *dvector( int nl,  int nh);
double  **dpvector(int nl, int nh);
double complex *cvector(int nl, int nh);
int     **imatrix( int nrl, int nrh, int ncl, int nch);
double  **dmatrix( int nrl, int nrh, int ncl, int nch);
double complex **cmatrix(int nrl, int nrh, int ncl, int nch);
double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_ivector(int       *v, int nl,  int nh);
void free_dvector(double    *v, int nl,  int nh);
void free_dpvector(double   **v, int nl,  int nh);
void free_cvector(double complex *v, int nl, int nh);
void free_imatrix(int      **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double   **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix(double complex **m, int nrl, int nrh, int ncl, int nch);
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
								int ndl, int ndh);
void fill0_dvector(double *X, int n0, int ntotal);
void fill0_ivector(int *X, int n0, int ntotal);
void fill0_cvector(double complex *X, int n0, int n1);
void fill0_dmatrix(double **X, int m0, int mtotal, int n0, int ntotal);
void fill0_imatrix(int **X, int m0, int mtotal, int n0, int ntotal);
void copy_dvector(double *aout, double *ain, int n0, int ntotal);
double norm1(double *v, int n);
double norm2(double *v, int n);
double scalarproduct(double *v, double *w, int n);

int minimum2(int i,int j);
int minimum3(int i,int j,int k);
int maximum2(int i,int j);
double dmaximum2(double a, double b);
double dminimum2(double a, double b);
int maximum3(int i,int j,int k);
int pow_int(int mantisse,int exponent);
double sinch(double x);
double Sqrt(double x);
double sqr(double x);



dcomplex Cadd(dcomplex a, dcomplex b);
dcomplex Csub(dcomplex a, dcomplex b);
dcomplex Cmul(dcomplex a, dcomplex b);
dcomplex RCmul(double x, dcomplex a);
dcomplex Cdiv(dcomplex a, dcomplex b);
dcomplex Complex(double re, double im);
dcomplex Conjg(dcomplex z);
double Cabs(dcomplex z);

dcomplex Csqrt(dcomplex z);
dcomplex Cexp(dcomplex z);
dcomplex Clog(dcomplex z);
dcomplex Csin(dcomplex z);
dcomplex Ccos(dcomplex z);
dcomplex Ctan(dcomplex z);
dcomplex Ccot(dcomplex z);
dcomplex Csinh(dcomplex z);
dcomplex Ccosh(dcomplex z);
dcomplex Ctanh(dcomplex z);
dcomplex Ccoth(dcomplex z);

void tridag(double a[], double b[], double c[], double r[], double u[],
   unsigned long n);
void ludcmp(double **a, int N, int *indx, double *d, int FLAG);
void lubksb(double **a, int N, int *indx, double *b, int FLAG);

double pythag(double a, double b);
void svdcmp(double **a, int m, int n, double w[], double **v);
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

void bandec(double **a, int n, int m1, int m2, double **al, int indx[], double *d);
void banbks(double **a, int n, int m1, int m2, double **al, int indx[], double b[]);

double plgndr(int l, int m, double x);
double plgndr_axis_normalised(int l, int m, double x);
double gammln(double xx);
double factrl(int n);
double dfactrl(double n);

double rd(double x, double y, double z);
double rf(double x, double y, double z);
double ellf(double phi, double ak);
double elle(double phi, double ak);

double binomial(double n, int k);
