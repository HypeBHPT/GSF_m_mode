#include "Solve_PDE.h"

//--------------------------------------------------------------------------------------
void func_rho(parameters par, double sigma, double *rho, double *drho_dsigma, double *d2rho_dsigma2){
  double rho_0 = par.rho_0, rho_1 = par.rho_1;


  *rho = rho_0 + rho_1*sigma;
  *drho_dsigma = rho_1;
  *d2rho_dsigma2 = 0.;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_beta(parameters par, double sigma, double *beta, double *dbeta_dsigma){
  double rho, drho, d2rho;
  func_rho(par, sigma, &rho, &drho, &d2rho);
  
  *beta = rho - sigma*drho;
  *dbeta_dsigma = - sigma*d2rho;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_r_of_sigma(parameters par, double sigma, double *r_over_M, double *dr_dsigma){
  double rho, drho, d2rho, beta, dbeta;
  
  func_rho(par, sigma, &rho, &drho, &d2rho);
  func_beta(par, sigma, &beta, &dbeta);
  
  
  *r_over_M = par.lambda_over_M*rho/sigma;
  *dr_dsigma = -par.lambda_over_M*beta/sqr(sigma);
  
  return;
}
//--------------------------------------------------------------------------------------
void func_sigma_of_r(parameters par, double r_over_M, double *sigma, double *dsigma_dr){
  
  *sigma = par.lambda_over_M*par.rho_0/(r_over_M - par.lambda_over_M*par.rho_1);
  *dsigma_dr = - par.lambda_over_M*par.rho_0/sqr(r_over_M - par.lambda_over_M*par.rho_1);
  
  return;
}
//--------------------------------------------------------------------------------------
void func_f(parameters par,double r_over_M, double *f, double *df_dr){
  double rh = par.rh_over_M, rc = par.rC_over_M, r = r_over_M, r2=sqr(r), r3 = r*r2;

  *f=(1.-rh/r)*(1.- rc/r);
  *df_dr= rh/r2 + rc/r2 - 2.*rc*rh/r3;

  return;
}
//--------------------------------------------------------------------------------------
void func_Delta(parameters par, double r_over_M, double *Delta, double *d_Delta_dr){
  double rh = par.rh_over_M, rc = par.rC_over_M, r = r_over_M;

  *Delta = (r- rh)*(r-rc);
  *d_Delta_dr = (r- rh) + (r-rc);

  return;
}
//--------------------------------------------------------------------------------------
void func_Sigma(parameters par, double r_over_M, double y, double *Sigma){
  double r2 = sqr(r_over_M), cos_theta_2=y, a = par.a_over_M, a2 = sqr(a);

  *Sigma = r2 + a2*cos_theta_2;
  
  return;
}
//--------------------------------------------------------------------------------------
void get_x_infty(parameters par, double sigma, double *x_infty, double *dx_infty_dsigma){
  double rh=par.rh_over_M, lambda = par.lambda_over_M, rh_over_lambda = rh/lambda, kappa = par.kappa, kappa2=sqr(kappa), rho0=par.rho_0;

   *x_infty =   rho0/sigma - rh_over_lambda*(1+kappa2)*log(sigma);
   *dx_infty_dsigma =  -rho0/sqr(sigma) - rh_over_lambda*(1+kappa2)/sigma;

   return;
}
//--------------------------------------------------------------------------------------
void get_x_horizon(parameters par, double sigma, double *x_hrz, double *dx_hrz_hrz_dsigma){
  double rh=par.rh_over_M, lambda = par.lambda_over_M, rh_over_lambda = rh/lambda, kappa = par.kappa, kappa2=sqr(kappa), Kh;

  Kh = 0.5*(1.-kappa2)/(1.+kappa2); //BH Horizon Surface Gravity
  *x_hrz = 0.5*rh_over_lambda/Kh*log(1-sigma);
  *dx_hrz_hrz_dsigma = -0.5*rh_over_lambda/Kh/(1-sigma);

  return;
}
//--------------------------------------------------------------------------------------
void get_x_Cauchy(parameters par, double sigma, double *x_C, double *dx_C_dsigma){
  double rh=par.rh_over_M, lambda = par.lambda_over_M, rh_over_lambda = rh/lambda, kappa = par.kappa, kappa2=sqr(kappa), Kc;

  Kc = -0.5*(1.-kappa2)/(kappa2*(1.+kappa2)); //Cauchy Horizon Surface Gravity
 
  if(kappa==0. || par.rho_1 != 0.){
    *x_C = 0.;
    *dx_C_dsigma = 0.;
  }
  else{
    *x_C =         0.5*rh_over_lambda/Kc*log(1-kappa2*sigma);
    *dx_C_dsigma = -kappa2*0.5*rh_over_lambda/Kc/(1-kappa2*sigma);
  }
  

  return;
}
//--------------------------------------------------------------------------------------
void get_x_Regular(parameters par, double sigma, double *x_Reg, double *dx_Reg_dsigma){

  *x_Reg = 0.;
  *dx_Reg_dsigma = 0.;

  return;
}
//--------------------------------------------------------------------------------------
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma){
  double x_infty, dx_infty_dsigma, x_hrz, dx_hrz_hrz_dsigma, x_C, dx_C_dsigma, x_Reg, dx_Reg_dsigma;

  get_x_infty(par, sigma, &x_infty, &dx_infty_dsigma);
  get_x_horizon(par, sigma, &x_hrz, &dx_hrz_hrz_dsigma);
  get_x_Cauchy(par, sigma, &x_C, &dx_C_dsigma);
  get_x_Regular(par, sigma, &x_Reg, &dx_Reg_dsigma);

  *x = x_Reg + x_infty + x_hrz + x_C;
  *dx_dsigma = dx_Reg_dsigma + dx_infty_dsigma + dx_hrz_hrz_dsigma + dx_C_dsigma;

  return;

}
//--------------------------------------------------------------------------------------
void get_Chi_horizon(parameters par, double sigma, double *Chi_hrz, double *dChi_hrz_hrz_dsigma){
  double kappa = par.kappa, kappa2=sqr(kappa);

  *Chi_hrz = kappa/(1.-kappa2)*log(1-sigma);
  *dChi_hrz_hrz_dsigma = -kappa/(1.-kappa2)/(1-sigma);

  return;
}
//--------------------------------------------------------------------------------------
void get_Chi_Cauchy(parameters par, double sigma, double *Chi_C, double *dChi_C_dsigma){
  double kappa = par.kappa, kappa2=sqr(kappa);

  if(par.rho_1 != 0.){
    *Chi_C = 0.;
    *dChi_C_dsigma = 0.;
  }
  else{
    *Chi_C = -kappa/(1.-kappa2)*log(1-kappa2*sigma);
    *dChi_C_dsigma = kappa2 * kappa/(1.-kappa2)/(1-kappa2*sigma);
  }

  return;
}
//--------------------------------------------------------------------------------------
void get_Chi_Regular(parameters par, double sigma, double *Chi_Reg, double *dChi_Reg_dsigma){
  *Chi_Reg = 0.;
  *dChi_Reg_dsigma = 0.;

  return;
}
//--------------------------------------------------------------------------------------
void func_tortoise_Chi(parameters par, double sigma, double *Chi, double *dChi_dsigma){
  double Chi_hrz, dChi_hrz_hrz_dsigma, Chi_C, dChi_C_dsigma, Chi_Reg, dChi_Reg_dsigma;

  get_Chi_horizon(par, sigma, &Chi_hrz, &dChi_hrz_hrz_dsigma);
  get_Chi_Cauchy(par, sigma, &Chi_C, &dChi_C_dsigma);
  get_Chi_Regular(par, sigma, &Chi_Reg, &dChi_Reg_dsigma);

  *Chi = Chi_hrz + Chi_C + Chi_Reg;
  *dChi_dsigma = dChi_hrz_hrz_dsigma + dChi_C_dsigma + dChi_Reg_dsigma;

  return;
}
//--------------------------------------------------------------------------------------
void func_Height(parameters par, double sigma, double *H, double *dH_dsigma){
 //height function w.r.t to BL coordinate t: Out-In Minimal Gauge

  double x_infty, dx_infty_dsigma, x_hrz, dx_hrz_hrz_dsigma, x_C, dx_C_dsigma, x_Reg, dx_Reg_dsigma;

  get_x_infty(par, sigma, &x_infty, &dx_infty_dsigma);
  get_x_horizon(par, sigma, &x_hrz, &dx_hrz_hrz_dsigma);
  get_x_Cauchy(par, sigma, &x_C, &dx_C_dsigma);
  get_x_Regular(par, sigma, &x_Reg, &dx_Reg_dsigma);

  *H = -x_Reg - x_infty + x_hrz + x_C;
  *dH_dsigma = -dx_Reg_dsigma - dx_infty_dsigma + dx_hrz_hrz_dsigma + dx_C_dsigma;

  return;
}
//--------------------------------------------------------------------------------------
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma){

  *Omega = sigma/par.lambda_over_M;
  *dOmega_dsigma = 1./par.lambda_over_M;
  return;
}
//--------------------------------------------------------------------------------------
void func_Z(parameters par, double sigma, double y, int FLAG_NS, double complex *Z, double complex *dlnZ_dsigma, double complex *dlnZ_dy){
  // NOTE: Z Factor with respect to Teukoslky Master Function in BL coordinates, defined in terms of Kinnersley tetrad
  //INPUT: angular coordinate y = cos(theta)^2. 
  //FLAG_NS Specifies if field is calcualte at north (FLAG_NS = 1): 0 <= theta <= pi/2, 0<=x <=1, x = sqrt(y), or south (FLAG_NS = -1) pi/2 <= theta <= pi, -1<=x <=0, x = - sqrt(y)
  //
  
  int ss = par.spin, m = par.m, d1 = abs(m), d2 = abs(m);
  double complex s=par.s;
  double scale_cte,
  cos_thetaOvertwo, sin_thetaOvertwo, dcos_thetaOvertwo_dy, dsin_thetaOvertwo_dy, 
  x = FLAG_NS*sqrt(y),
  r, dr_dsigma, Delta, d_Delta_dr, d_Delta_dsigma,
  Omega, dOmega,
  H, dH, Chi, dChi;

  

  func_r_of_sigma(par, sigma, &r, &dr_dsigma);
  func_Delta(par, r, &Delta, &d_Delta_dr);
  d_Delta_dsigma = dr_dsigma*d_Delta_dr;


  scale_cte = 1;//pow(rh, 1 + 2*ss)/pow(lambda, ss);
  cos_thetaOvertwo = sqrt(0.5*( 1. + x));
  dcos_thetaOvertwo_dy = 1./(8. * cos_thetaOvertwo * x);

  sin_thetaOvertwo = sqrt(0.5*( 1. - x));
  dsin_thetaOvertwo_dy = -1./(8. * sin_thetaOvertwo * x);
      
  func_Omega(par, sigma, &Omega, &dOmega);
  func_Height(par, sigma, &H, &dH);
  func_tortoise_Chi(par, sigma, &Chi, &dChi);
  
    
  *Z = scale_cte * Omega * pow(Delta, -ss) * cexp( s*H ) * pow(cos_thetaOvertwo, -d1) * pow(sin_thetaOvertwo, -d2) * cexp( I*m*Chi );
    
  *dlnZ_dsigma = dOmega/Omega - ss*d_Delta_dsigma/Delta + s*dH + I*m*dChi ;
  *dlnZ_dy =  -d1 *  dcos_thetaOvertwo_dy/cos_thetaOvertwo - d2 * dsin_thetaOvertwo_dy/sin_thetaOvertwo;

  // double complex out = *Z;
  // printf("%3.15e %3.15e, %3.15e %3.15e\n", sigma, y, creal(out), cimag(out));
  
  return;
}
//------------------------------------------------------------
void func_alpha2(parameters par, double sig, double complex *alpha2){
  double sig2 = sqr(sig);
    double kappa=par.kappa;
    double kappa2 = sqr(kappa);

    *alpha2 = -((-1+sig)*sig2*(-1+kappa2*sig));
    *alpha2 *=-1.; //NOTE: Mathematica Notebook had the opposite sign definition for Teukolskly operator


  return;
}
//------------------------------------------------------------
void func_alpha1(parameters par, double sig, double complex *alpha1){
  double sig2 = sqr(sig);
  double kappa=par.kappa;
  double kappa2 = sqr(kappa);
  int m=par.m;
  double lambda=par.lambda_over_M;
  double rh=par.rh_over_M;

  double complex s=par.s;
    
    *alpha1 = (2*rh*s*(-1 + (2 + kappa2 *(3 - 2*sig) -2*kappa2*kappa2*(-1 +sig))* sig2))/lambda + sig*(-2 + sig*(3 +2* I*m*kappa  +kappa2*(3-4*sig)));
    *alpha1 *=-1.; //NOTE: Mathematica Notebook had the opposite sign definition for Teukolskly operator 
 
  return;
}
//------------------------------------------------------------
void func_alpha0(parameters par, double sig, double y, double complex *alpha0){
  double  m=1.*par.m;
  double delta1=fabs(m);
  double delta2=fabs(m);
  double sig2 = sqr(sig);
  double kappa=par.kappa;
  double kappa2 = sqr(kappa);
  double lambda=par.lambda_over_M;
  double rh=par.rh_over_M;    
  double complex s=par.s;
  

    *alpha0 = (sqr(rh)*s*s*(4*(1 + sig) + kappa2*(7 + y + 4*kappa2 +4*(2 + 2*kappa2 + sqr(kappa2))*sig -4 *sqr(1 +kappa2)*sig2)))/sqr(lambda) + (2*rh*s*((2 + sqr(kappa2)* (2 - 3*sig) -3*kappa2*(-1 + sig))*sig +I*m*(kappa +2*kappa*(1 + kappa2)*sig) ))/lambda +(sqr(m) + delta1*(-1 + delta2) - delta2 +4*I*m *kappa *sig - 4*kappa2*sig2 +2*(1 + kappa2) *sig)/2;
    *alpha0 *=-1.; //NOTE: Mathematica Notebook had the opposite sign definition for Teukolskly operator

  return;
}
//------------------------------------------------------------
void func_gamma2(parameters par, double y, double complex *gamma2){
  
    *gamma2 = 4*y*(y-1);
    *gamma2 *=-1.; //NOTE: Mathematica Notebook had the opposite sign definition for Teukolskly operator

  return;
}
//------------------------------------------------------------
void func_gamma1(parameters par, double y, double complex *gamma1){
 double m =1.*par.m;
 double delta1=fabs(m);
 double delta2=fabs(m);
 *gamma1 =  -2.*(1.+y*(-3.+delta1+delta2));
 *gamma1 *=-1.; //NOTE: Mathematica Notebook had the opposite sign definition for Teukolskly operator
    
  return;
}

//------------------------------------------------------------
double complex Diff_Operator(parameters par, complex_derivs_2D W, int iDom, int j1, int j2, int FLAG_dr0){
  
  double chi_1, chi_2;
  func_derivs_2D sigma, y;
  chi_1 = par.grid_chi_1[iDom][j1];
  chi_2 = par.grid_chi_2[iDom][j2];
  get_sigma(par, iDom, chi_1, chi_2, &sigma);
  get_y(par, iDom, chi_1, chi_2, &y);

  
  double complex a2, a1, a0, g2, g1, A_f=0.;
  complex_sigma_derivs f;

  get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, W, &f);
    
  func_alpha0(par, sigma.d0, y.d0, &a0);
  func_alpha1(par, sigma.d0, &a1);
  func_alpha2(par, sigma.d0, &a2);
  func_gamma2(par, y.d0, &g2);
  func_gamma1(par, y.d0, &g1);

  A_f = a2*f.d2sigma + a1*f.dsigma + a0*f.d0 + g2*f.d2y + g1*f.dy;
    
  return A_f;
}
//--------------------------------------------------------------------------------------
double complex HyperboloidalEffectiveSource(parameters par, int iDom, int j1, int j2){
  double complex hyp_S_eff=0;
  
  if(iDom == par.Dom_ptcl){
    int N1_Seff=par.N1_PuncSeff, N2_Seff=par.N2_PuncSeff;
    
    double real_hS_eff, imag_hS_eff, chi_1=par.grid_chi_1[iDom][j1], chi_2=par.grid_chi_2[iDom][j2];

    real_hS_eff = Clenshaw_Chebyshev_2D(par.Re_cheb_Seff, N1_Seff, N2_Seff, chi_1, chi_2);
    imag_hS_eff = Clenshaw_Chebyshev_2D(par.Im_cheb_Seff, N1_Seff, N2_Seff, chi_1, chi_2);

    hyp_S_eff = real_hS_eff + I*imag_hS_eff;
  }  

  return hyp_S_eff;
}