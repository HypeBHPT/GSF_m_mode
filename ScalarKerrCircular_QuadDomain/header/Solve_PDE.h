// source /home_local/ansorg/src/intel/bin/compilervars.sh intel64
// source /home_local/nb-oma/Downloads/src_Intel/bin/compilervars.sh intel64
// clear; icc -O2 *.c -o Solve_PDEs

#include "utilities.h"
#include "spectral_utilities.h"
#include "effsource.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define nDom 4
#define nFields 2
#define nPar 0


#define Newton_itmin  2  
#define Newton_itmax  5  
#define Newton_tol    1.e-12
#define Newton_verb   1 

#define bicgstab_itmax 200
#define bicgstab_decr 1.e-3
#define bicgstab_verb 1


#define FD_ORDER 4
#if FD_ORDER == 2
 #define STENCILSIZE 9 
#elif FD_ORDER == 4
 #define STENCILSIZE 25
#elif FD_ORDER == 6
 #define STENCILSIZE 49
#endif

typedef struct PARAMETERS{ 
	int N1[nDom], N2[nDom], n1[nDom], n2[nDom], n1n2[nDom+1], n_2D, ntotal, Ntotal, SOLVER_METHOD;//, N1_max[nDom], N2_max[nDom];
	double bicgstab_tol;
	
	
	char grid_1[nDom][50], grid_2[nDom][50], SimName[800];

	double 	*grid_chi_1[nDom], *grid_chi_2[nDom],
			AnMR_kappa_1[nDom], AnMR_x_boundary_1[nDom], AnMR_kappa_2[nDom], AnMR_x_boundary_2[nDom];
	
	//---- HYPERBOLOIDAL SELF FORCE PARAMETERS
	
	int spin, m, nbar, FLAG_Trajec;
	int Dom_scri, Dom_bulk, Dom_ptcl, Dom_hrzn;
	double q, a_over_M, kappa, r0_over_M, r_plus_over_M, r_minus_over_M, lambda_over_M, rh_over_M, rC_over_M, rho_0, rho_1, f0, E0, M_Omega0, L0_over_M, eta,
	sigma0, sigma_minus, sigma_plus, sigma_hrz, sigma_Cauchy;


	int N1_LoadPunc, N1_LoadSeff, N2_LoadSeff; //Resolution for loading external data
	int N1_PuncSeff, N2_PuncSeff; //Resolution for inbuilt Puncture and Effective Source

	char grid_1_PuncSeff[50], grid_2_PuncSeff[50]; //grid for inbuilt Puncture and Effective Source
	
	double 	**Re_cheb_Seff, **Im_cheb_Seff, 
			**Re_cheb_PuncField, **Im_cheb_PuncField,
			*Re_cheb_phi_Punc, *Re_cheb_phi_Punc_sigma, *Re_cheb_phi_Punc_y,
			*Im_cheb_phi_Punc, *Im_cheb_phi_Punc_sigma, *Im_cheb_phi_Punc_y;

	double complex s;

	double rho_min;

	double prec;
	FILE *fout; 
	int i_omp,  CoordMap_FLAG, TEST_Func_FLAG;
	
} parameters;

typedef struct DERIVS_Space_VECTOR{ 
	double *d0, 
		   *d1, *d2,
		   *d11,  *d22, *d12,
		   *d111, *d112, *d122, *d222,
		   *d1111, *d1112, *d1122, *d1222, *d2222,
		   *d11222;

} derivs_2D;

typedef struct DERIVS_Space_{ 
	double d0, d1, d2,
		   d11,  d22, d12,
		   d111, d112, d122, d222,
		   d1111, d1112, d1122, d1222, d2222,
		   d11222;

} func_derivs_2D;

typedef struct SIGMA_DERIVS_Space_{ 
	double d0, dsigma, d2sigma, dy, d2y, d2sigmay;

} func_sigma_derivs_2D;

typedef struct COMPLEX_DERIVS_Space{ 
	double complex d0, d1, d2,
		   		   d11,  d22, d12,
		           d111, d112, d122, d222,
		           d1111, d1112, d1122, d1222, d2222,
				   d11222;
	
} complex_derivs_2D;

typedef struct COMPLEX_SIGMA_DERIVS_Space{ 
	double complex d0, dsigma, d2sigma, dy, d2y, d2sigmay;
	
} complex_sigma_derivs;



typedef struct JFD_COMPONENTS{ 
	double	**J_UL, **K_UL, **Kl_UL;
	int	*ncols_J_UL, **cols_J_UL, *iK_UL,
		m1, m2;
		
	double	**J_UR, **J_LL, **J_LR, **K_UR, **K_LR;
	int    *indx_LR;
	
	
	
} JFD_Components;

// Routines in "test_routines.c"
func_sigma_derivs_2D Test_Func(parameters par, double sigma, double y);
void output_test_function(parameters par);
double complex Test_Effective_Source(parameters par, int iDom, int i1, int i2);
void output_test_Effective_Source(parameters par);
void output_test_function_error(parameters par, double *X);
void Check_Derivative(parameters par, char *fn, int j1, int j2, derivs_2D W, int FLAG);

// Routines in "parameters.c"
void set_parameters(parameters *par, int N, int nbar, int m, double M_Omega0, double a_over_M);

// Routines in "solve_equations.c"
void get_Solution(parameters par, double *X, derivs_2D Sol);
int solve_equations(parameters par, double *X);

// Routines in "load_data.c"
void read_header(parameters par, FILE *fr, int *N1, int *N2);
void load_Puncture_at_Boundary(parameters *par);
void load_EffectiveSource(parameters *par);
void free_external_data(parameters *par);
void load_PunctureField(parameters *par);


// Routines in "SelfForce_functions.c"
double get_clm(int l, int m);
complex_derivs_2D get_complex_dervs_conjugate(complex_derivs_2D z);
double complex Get_phi_From_hyp_phi(parameters par, double sigma, double y, int FLAG_NS, double complex V);
void Get_SelfForce_Ft_m(parameters par, double *X, double complex *Ft_m_bulk, double complex *Ft_m_equator);
void Get_SelfForce_Fr_m(parameters par, double *X, double complex *Fr_m_bulk, double complex *Fr_m_equator);
void get_Puncture_EffectiveSource(parameters *par);
double complex Get_HypFunc_From_BLFunc(parameters par, double sigma, double y, int FLAG_NS, double complex BL_Func);
void get_HypDerv_from_BLDerv(parameters par, double sigma, double y, int FLAG_NS,  double complex d_phi_dr, double complex d_phi_dtheta, double complex *d_phi_dsigma, double complex *d_phi_dy);
void get_Derivatives_HypFunc_From_BLFunc(parameters par, double sigma, double y, int FLAG_NS, double complex BL_Func, double complex dBL_Func_dsigma, double complex dBL_Func_dy, double complex *dHypFunc_dsigma, double complex *dHypFunc_dy);
double complex Rescale_Source(parameters par, double sigma, double y,  double complex SourceIn);

 // Routines in "hyperboloidal_functions.c"
void func_rho(parameters par, double sigma, double *rho, double *drho_dsigma, double *d2rho_dsigma2);
void func_beta(parameters par, double sigma, double *beta, double *dbeta_dsigma);
void func_r_of_sigma(parameters par, double sigma, double *r, double *dr_dsigma);
void func_sigma_of_r(parameters par, double r, double *sigma, double *dsigma_dr);
void func_f(parameters par,double r_over_M, double *f, double *df_dr);
void func_Delta(parameters par, double r_over_M, double *Delta, double *d_Delta_dr);
void func_Sigma(parameters par, double r_over_M, double y, double *Sigma);
void get_x_infty(parameters par, double sigma, double *x_infty, double *dx_infty_dsigma);
void get_x_horizon(parameters par, double sigma, double *x_hrz, double *dx_hrz_hrz_dsigma);
void get_x_Cauchy(parameters par, double sigma, double *x_C, double *dx_C_dsigma);
void get_x_Regular(parameters par, double sigma, double *x_Reg, double *dx_Reg_dsigma);
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma);
void get_Chi_horizon(parameters par, double sigma, double *Chi_hrz, double *dChi_hrz_hrz_dsigma);
void get_Chi_Cauchy(parameters par, double sigma, double *Chi_C, double *dChi_C_dsigma);
void get_Chi_Regular(parameters par, double sigma, double *Chi_Reg, double *dChi_Reg_dsigma);
void func_tortoise_Chi(parameters par, double sigma, double *Chi, double *dChi_dsigma);
void func_Height(parameters par, double sigma, double *H, double *dH_dsigma);
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma);
void func_Z(parameters par, double sigma, double y, int FLAG_NS, double complex *Z, double complex *dlnZ_dsigma, double complex *dlnZ_dy);

void func_alpha2(parameters par, double sig, double complex *alpha2);
void func_alpha1(parameters par, double sig, double complex *alpha1);
void func_alpha0(parameters par, double sig, double y, double complex *alpha0);
void func_gamma2(parameters par, double y, double complex *gamma2);
void func_gamma1(parameters par, double y, double complex *gamma1);
double complex HyperboloidalEffectiveSource(parameters par, int iDom, int j1, int j2);
double complex Diff_Operator(parameters par, complex_derivs_2D W, int iDom, int j1, int j2, int FLAG_dr0);

// Routines in "Allocate.c"
void allocate_derivs_2D(derivs_2D *v, int n);
void fill0_derivs_2D(derivs_2D v, int n);
void free_derivs_2D(derivs_2D *v, int n);


// Routines in "newton_direct.c"
void Jacobian_FD(parameters par, double *X, double **J);
void Jacobian(parameters par, double *X, double **J);
void JacobianBand(parameters par, double *X, double **J);
int newton_direct(parameters par, double *X);
int newton_SVN(parameters par, double *X);

// Routines in "newton.c"
void get_BandMatrix(parameters par,double *X, JFD_Components *JFD);
void free_bandMatrix(parameters par,JFD_Components *JFD);
void Get_boundary_i1(parameters par, int jdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1, int *mdom);
void Get_boundary_i2(parameters par, int jdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1, int *mdom);
void Get_gridpoints(parameters par, int n, int *mdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1);
void Get_JFD_Matrix(parameters par, double *X, JFD_Components *JFD);
void Get_JFD_Components(parameters par, double *X, JFD_Components *JFD);
void PreCond(parameters par, JFD_Components JFD, double *X, double *b, double *Y, double *res);
void PreCondLU(parameters par, double **J, int *indx, double *X, double *b, double *Y, double *res);
int bicgstab(parameters par, double *X, double *DX, double *F, double *normres);
int newton(parameters par, double *X);
int newton_error(parameters par, double *X);


// Routines in "FieldEqs.c"
void get_ComplexField(parameters par, int indx_Re, int indx_Im, derivs_2D W , complex_derivs_2D *z);
void FieldEquations(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F);
void LinearFieldEquations(parameters par, derivs_2D w, derivs_2D Dw, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J);
void BoundaryCondition(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F);
void LinearBoundaryCondition(parameters par, derivs_2D w, derivs_2D Dw, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J);
void ParameterEquations(parameters par, derivs_2D w, double *Xpar, double *F);
void LinearParameterEquations(parameters par, derivs_2D w, derivs_2D Dw, double *Xpar, double *DXpar, double *J);
void FieldJumpEquations_Boundary_chi1(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *F);
void LinearFieldJumpEquations_Boundary_chi1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *J);
void FieldDerivativeJumpEquations_Boundary_chi1(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *F);
void LinearFieldDerivativeJumpEquations_Boundary_chi1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *J);
void FieldJumpEquations_Boundary_chi2(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *F);
void LinearFieldJumpEquations_Boundary_chi2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *J);
void FieldDerivativeJumpEquations_Boundary_chi2(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *F);
void LinearFieldDerivativeJumpEquations_Boundary_chi2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *J);


// Routines in "derivatives.c"
void Spectral_derivative(char *grid, double *f, double *df, double *d2f, int N);
void Get_Derivatives(parameters par, derivs_2D W);
void Get_DerivativesFinDif(parameters par, derivs_2D v);
void Get_DerivativesFinDif_grid(parameters par, int i, int j, int i_field, int idom, derivs_2D v);

// Routines in "FuncAndJacobian.c"
void copy_X_to_W(parameters par, double *X, derivs_2D W);
void get_W_at_Grid(parameters par, int iDom, int j1, int j2, derivs_2D W, derivs_2D W_atGrid);
void get_ExtraPar_from_X(parameters par, double *X, double *ExtraParameters);
void F_of_X(parameters par, double *X, double *F);
void J_times_DX(parameters par, double *X, double *DX, double *F);
void J_times_DX_FD(parameters par, double *X,double *DX, double *JDX);
void Jacobian(parameters par, double *X, double **J);


//Routine in get_InitialGuess.c
void get_InitialGuess_test(parameters par, double *X);
void get_InitialGuess(parameters par, double *X);
void input_Solution(parameters par, double *X);


//Routines in IndexRoutines.c
int Index_Boundary(char grid[], int i, int n);
void Get_Index_Arrays(parameters *par);
int Index(parameters par, int iDom, int iField, int j1, int j2);
void Get_Indices_From_Index(parameters par, int indx, int *idom, int *If, int *i, int *j);

// Routines in "output.c"
void output_Solution(parameters par, double *X);
void output_SpecCoef(parameters par, double *X);
void output_SolutionChebyshevPostProc(parameters par, double *X);
void output_SelfForce(parameters par, double *X);

void output_F_of_X(parameters par, double *X);
void output_EffectiveSource(parameters par);
void output_Seff_SpecCoef(parameters par, double **Seff, char *Seff_sector);
void output_DifferentialOperator(parameters par, double *X);
void A_of_X(parameters par, double *X, double *F);
void FieldEquations_AX(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F);
void output_SpecCoef_AX(parameters par, double *X);
void output_Solution_AX(parameters par, double *X);




void output_FieldParticle(parameters par, double *X);
void output_SolutionDerivatives(parameters par, double *X);


void output_PunctureField(parameters par);
void output_Puncture_at_Boundary(parameters par);


// Routines in "Coordinates.c"
int total_grid_points(char *grid, int N);
double get_grid(char *grid, int i, int N);
void Construct_grid(parameters *par);
void free_grid(parameters *par);
double get_Jacobian(parameters par, int iDom, int j1, int j2);
double get_x_from_chi(double x_boundary, double kappa, double chi);
double get_chi_from_x(double x_boundary, double kappa, double x);
void get_sigma(parameters par, int iDom, double chi_1, double chi_2, func_derivs_2D *sigma);
void get_y(parameters par, int iDom, double chi_1, double chi_2, func_derivs_2D *y);
void MapFirstDerivative_Dom_a_to_Dom_b(parameters par, 
	int iDom_a, double chi_a_1, double chi_a_2, complex_derivs_2D W_a, 
	int iDom_b, double chi_b_1, double chi_b_2, complex_derivs_2D *W_b);
void get_chi1_chi2_From_sigma_y(parameters par, int iDom, func_derivs_2D sigma, func_derivs_2D y, double *chi_1, double *chi_2);
void get_PhysDerv_from_SpecDerv_complex(parameters par, int iDom, int j1, int j2, complex_derivs_2D W, complex_sigma_derivs *U);

double get_sigma_from_x1(double kappa, double x1);
double get_tau_from_x2(parameters par, double x2);
double get_x1_from_sigma(double kappa, double z);
double get_x2_from_tau(parameters par, double x);
void get_SpecDerv_from_x1x2(parameters par, derivs_2D W, derivs_2D U, double x1, int n);
double get_AngleIncrement(char *grid, int N);

// Routine in debug.c
void pause();
void PrintVector(parameters par, double *J, int na, int nb);
void PrintMatrix(parameters par, char *fp, double **J, int la, int lb, int ca, int cb);
void create_directory(char *dir_name);
void Check_OS();

// Routines in "Solve_PDEs.c"
int main() ;


