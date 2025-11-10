#include <stdio.h>
// #include <stdio_ext.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

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
#define TINY 1.5e-38
#define BIG 3.0e37
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*
#define NMAX 2000

typedef struct DCOMPLEX {double r,i;} dcomplex;


void nrerror(char error_text[]);
int      *ivector( int nl,  int nh);
double   *dvector( int nl,  int nh);
double complex 	  *cvector(int nl, int nh);
double  **dpvector(int nl, int nh);
int     **imatrix( int nrl, int nrh, int ncl, int nch);
double  **dmatrix( int nrl, int nrh, int ncl, int nch);
double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_ivector(int       *v, int nl,  int nh);
void free_dvector(double    *v, int nl,  int nh);
void free_cvector(double complex *v, int nl, int nh);
void free_dpvector(double   **v, int nl,  int nh);
void free_imatrix(int      **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double   **m, int nrl, int nrh, int ncl, int nch);
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
								int ndl, int ndh);
void fill0_dvector(double *X, int n0, int ntotal);
void fill0_ivector(int *X, int n0, int ntotal);
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

void tridag(double a[], double b[], double c[], double r[], double u[],
   unsigned long n);
void ludcmp(double **a, int N, int *indx, double *d, int FLAG);
void lubksb(double **a, int N, int *indx, double *b, int FLAG);

double plgndr(int l, int m, double x);
double gammln(double xx);
double factrl(int n);

double rd(double x, double y, double z);
double rf(double x, double y, double z);
double ellf(double phi, double ak);
double elle(double phi, double ak);


