// source /home_local/ansorg/src/intel/bin/compilervars.sh intel64
// source /home_local/nb-oma/Downloads/src_Intel/bin/compilervars.sh intel64
// clear; icc -O2 *.c -o Solve_ODE_Lobatto

#include "utilities.h"
#include "spectral_utilities.h"

#define nFields 2
#define nPar 0
#define nDom 2

#define Newton_itmin  2      
#define Newton_itmax  10      
#define Newton_tol    5.e-9
#define Newton_verb   1

#define flag_NewCoord 0




typedef struct DERIVS_Space{ 
	double *d0, *d1, *d11;
	
} derivs;

typedef struct COMPLEX_DERIVS_Space{ 
	double complex d0, d1, d11;
	
} complex_derivs;

typedef struct COMPLEX_SIGMA_DERIVS_Space{ 
	double complex d0, dsigma, d2sigma;
	
} complex_sigma_derivs;

typedef struct PARAMETERS{ 
	int N[nDom], ntotal, Ntotal, mA[nDom+1];
	
	char SimName[200], grid[nDom][50], Equation[50];	

	double AnMR_kappa[nDom], AnMR_x_boundary[nDom];
	
	int ell, spin, m, effective_source_FLAG, 
	Nread_Seff_rPlus, Nread_Seff_rMinus,
	Nread_phiP_rPlus, Nread_phiP_rMinus;
	
	double q, r0_over_M, rho_0, rho_1, lambda_over_M, r_plus_over_M, r_minus_over_M, 
		 sigma[nDom+1], sigma0, f0, E0, M_Omega0, L0_over_M, eta,
		 sigma_minus, sigma_plus, rh_over_M,
		 jump_field[nFields][nDom-1], jump_derivative[nFields][nDom-1];
	
	double complex s, ds_dr0,  bar_kappa, dbar_kappa_dr0;

	FILE *fout; 
	int i_omp;
	
} parameters;


  
// Routines in "Allocate.c"
void allocate_derivs(derivs *v, int n);
void fill0_derivs(derivs v, int n);
void free_derivs(derivs *v, int n);

// Routines in "Derivatives.c"
void Spectral_derivative(char *grid, double *f, double *df, double *d2f, int N);
void Get_Derivatives(parameters par, derivs W);


// Routines in "newton.c"
int newton(parameters par, double *X);

// Routines in "IndexRoutines.c"
int Index_i(parameters par, int iDom, int i);
int Index(parameters par, int iDom, int iField, int i);
void Get_Index_Arrays(parameters *par);

// Routines in "parameteres.c"
void set_parameters(parameters *par, double r0, int l, int m, double eta, int N);

// Routines in "InitialGuess.c"
void Initial_Guess(parameters par, double *X);

// Routines in "debug.c"
void PrintVector(double *J, int na, int nb);
void PrintMatrix(double **J, int la, int lb, int ca, int cb);
void pause();
void create_directory(char *dir_name);


//Routines in "hyperboloidal_functions.c"
void func_rho(parameters par, double sigma, double *rho, double *drho_dsigma, double *d2rho_dsigma2);
void func_beta(parameters par, double sigma, double *beta, double *dbeta_dsigma);
void func_r_of_sigma(parameters par, double sigma, double *r, double *dr_dsigma);
void func_sigma_of_r(parameters par, double r, double *sigma, double *dsigma_dr);
void func_f(double r, double *f, double *df_dr);
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma);
void func_height(parameters par, double sigma, double *h, double *dh_dsigma);
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma);
void func_Z(parameters par, double sigma, double complex *Z, double complex *dlnZ_dsigma);
double complex Diff_Operator(parameters par, complex_derivs W, int iDom, int i, int FLAG_dr0);
void func_alpha2(parameters par, double sig, double complex *alpha2, double complex *dalpha2_dsig, double complex *dalpha2_dr0);
void func_alpha1(parameters par, double sig, double complex *alpha1, double complex *dalpha1_dr0);
void func_alpha0(parameters par, double sig, double complex *alpha0, double complex *dalpha0_dr0);
double complex Get_ModeSum_Source(parameters par, int iDom, int i);
double complex Get_Extended_Source(parameters par, complex_derivs W, int iDom, int i);


// Routines in "SelfForce_Functions.c"
double get_clm(int l, int m);
double complex Get_phi_From_v(parameters par, double sigma, double complex V);
double complex Get_psi_From_w(parameters par, double sigma, double complex V, double complex W);
double complex Get_psi_From_w_at_particle(parameters par, double complex V, double complex W);
double get_Ft_PN(parameters par);
// double get_dr0_Ft_PN(parameters par);
void get_Ft_SSF_lm(parameters par, double *X, double complex *Ft_SSF_lm_dom1, double complex *Ft_SSF_lm_dom2);
// void get_Ft_dr0_SSF_lm(parameters par, double *X, double complex *Ft_dr0_SSF_lm_dom1, double complex *Ft_dr0_SSF_lm_dom2);
void get_Flux_lm(parameters par, double *X, double *Flux_lm_hrzn, double *Flux_lm_scri, double *Flux_lm);
// void get_Flux_dr0_lm(parameters par, double *X, double *Flux_dr0_lm_hrzn, double *Flux_dr0_lm_scri, double *Flux_dr0_lm);
double get_Ft_FLux_lm(parameters par, double *X);
// double get_Ft_dr0_FLux_lm(parameters par, double *X);
void Sum_over_m(parameters par, double *X, double *sum_ell_Dom1, double *sum_ell_Dom2);
void Get_Bfield(parameters par, double *Bfield, double complex *Hyp_Bfield);
double get_Flux_PN(parameters par);
// double get_dr0_Flux_PN(parameters par);

// Routines in "FuncAndJacobian.c"
void copy_X_to_W(parameters par, double *X, derivs W);
void F_of_X(parameters par, double *X, double *F);
void Jacobian(parameters par, double *X, double **J);
void Jacobian_FD(parameters par, double *X, double **J);

// Routines in "Field_Eqs.c"
void get_ComplexField(parameters par, int indx_Re, int indx_Im, derivs W , complex_derivs *z);
void Field_Eqs(parameters par, int iDom, int i, derivs W, double *F);
void Boundary_Data_Left(parameters par, int iDom, int i, derivs W, double *F);
void Boundary_Data_Right(parameters par, int iDom, int i, derivs W, double *F);
void Transition_Condition_Left(parameters par, int iDom, int i, derivs W, double *F);
void Transition_Condition_Right(parameters par, int iDom, int i, derivs W, double *F);

//Routines in "output.c"
void output_cheb(parameters par, double *X);
void output_Solution(parameters par, double *X);
void output_Puncture(parameters par);
void output_Ret_Field_Boundary_data(parameters par, double *X);
void output_EffectiveSource(parameters par);
void output_phi_ell(parameters par, int ell_min, int ell_max, double *Sum_DomLeft, double *Sum_DomRight);
void output_phi(parameters par, int ell_min, int ell_max, double *Sum_Dom);
void output_Convergence(parameters par, double *X);


//Routines in "Coordinates.c"
double get_ChebyshevGrid_x(int i, int N, char *grid);
double get_x_from_chi(double x_boundary, double kappa, double chi);
double get_chi_from_x(double x_boundary, double kappa, double x);
void get_sigma(parameters par, int iDom, int i, double *sigma, double *dx_dsigma);
double get_x_from_sigma(parameters par, int iDom, double sigma);
void get_SigmaDerv_from_SpecDerv_complex(parameters par, int iDom, int i, complex_derivs W, complex_sigma_derivs *U);
double func_ZlnZ(double z);
double z_of_zeta(double zeta);
double zeta_of_z(double z);


// Routines in "Solve_ODE_Lobatto.c"
int main() ;
