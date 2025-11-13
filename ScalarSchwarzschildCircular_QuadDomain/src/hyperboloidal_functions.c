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
  
  *sigma = par.lambda_over_M*par.rho_0/r_over_M;
  *dsigma_dr = - par.lambda_over_M*par.rho_0/sqr(r_over_M);
  
  return;
}
//--------------------------------------------------------------------------------------
void func_f(double r_over_M, double *f, double *df_dr){     
  
  *f=(1.-2./r_over_M);
  *df_dr=2/sqr(r_over_M);

  return;
}
//--------------------------------------------------------------------------------------
void func_tortoise_x(parameters par, double sigma, double *x, double *dx_dsigma){
  

  *x= 2*( 1./sigma + log(1.-sigma) - log(sigma) )/par.lambda_over_M;

  double r,dr_dsigma;  
  func_r_of_sigma(par, sigma, &r, &dr_dsigma);    
    
  double f, df_dr;
  func_f(r, &f, &df_dr);

  double beta, dbeta;
  func_beta(par, sigma, &beta, &dbeta);
  
  double drstar_dr = 1./f, dx_drstar = 1./par.lambda_over_M;
  
  *dx_dsigma = dx_drstar*drstar_dr*dr_dsigma;
  
  return;
}
//--------------------------------------------------------------------------------------
void func_height(parameters par, double sigma, double *h, double *dh_dsigma){
 //height function w.r.t to Eddington Finkelstein coordinate v
  *h = -2*par.rho_0*(1./sigma - 2*log(sigma)/(par.lambda_over_M*par.rho_0)  );
  *dh_dsigma = 2*par.rho_0*(1./sqr(sigma) + 2./(sigma*par.lambda_over_M*par.rho_0)  );
  
  return ;
}
//--------------------------------------------------------------------------------------
void func_Height(parameters par, double sigma, double *H, double *dH_dsigma){
 //height function w.r.t to Eddington Finkelstein coordinate v
  double lambda=par.lambda_over_M, rh = par.rh_over_M;
  *H = rh/lambda * ( log(1-sigma) - 1./sigma + log(sigma));
  *dH_dsigma = rh/lambda * ( ( 1. -2*sqr(sigma) ) /(sqr(sigma)*(1.-sigma)));
  return ;
}
//--------------------------------------------------------------------------------------
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma){

  *Omega = sigma/par.lambda_over_M;
  *dOmega_dsigma = 1./par.lambda_over_M;
  return;
}
//--------------------------------------------------------------------------------------
void func_Z(parameters par, double sigma, double y, double complex *Z, double complex *dlnZ_dsigma){
  double complex s=par.s;
  double Omega, dOmega, x, h, dx, dh, m = 1.*par.m;
  
  func_Omega(par, sigma, &Omega, &dOmega);
  func_tortoise_x(par, sigma, &x, &dx);
  func_height(par, sigma, &h, &dh);
  
  *Z = Omega * cexp( s*(x+h) ) * pow(1.-y, -0.5*m);
  
  *dlnZ_dsigma = (dOmega/Omega + s*(dx+dh));
  
  
  return;
}
//------------------------------------------------------------
void func_alpha2(parameters par, double sig, double complex *alpha2, double complex *dalpha2_dsig, double complex *dalpha2_dr0){
  double sig2 = sqr(sig);

    *alpha2 = sig2*(1-sig);
    
    *dalpha2_dsig = 2*sig-3*sig2;
    *dalpha2_dr0 = 0.;

  return;
}
//------------------------------------------------------------
void func_alpha1(parameters par, double sig, double complex *alpha1, double complex *dalpha1_dr0){
  double complex s, ds_dr0;
  double sig2 = sqr(sig);

  s=par.s;
  ds_dr0 = par.ds_dr0;
  *alpha1 =  sig*(2-3*sig) + (1-2*sig2)*s;  
  

  *dalpha1_dr0 = ds_dr0*(1-2*sig2);  
 
  return;
}
//------------------------------------------------------------
void func_alpha0(parameters par, double sig, double complex *alpha0, double complex *dalpha0_dr0){
  double complex s, s2, ds_dr0;
  double spin=par.spin, spin2=sqr(spin), m=1.*par.m;

  s=par.s;
  ds_dr0 = par.ds_dr0;
  s2=s*s;

  *alpha0 = -(1+sig)*s2 - 2*sig*s - ( m*(m-1.) + sig*(1-spin2) );
  


  *dalpha0_dr0 = - 2*s*ds_dr0*(1+sig) - 2*ds_dr0*sig;

  return;
}
//------------------------------------------------------------
void func_gamma2(parameters par, double y, double complex *gamma2){
  
    *gamma2 = 4*y*(1-y);

  return;
}
//------------------------------------------------------------
void func_gamma1(parameters par, double y, double complex *gamma1){
  double m=1.*par.m;
    *gamma1 = 2. - 2.*y*(3. - 2.*m);
    
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

  
  double complex a2, a1, a0, g2, g1, da2, da2_dsig, da1, da0, A_f=0.;
  complex_sigma_derivs f;

  get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, W, &f);
    
  func_alpha0(par, sigma.d0, &a0, &da0);
  func_alpha1(par, sigma.d0, &a1, &da1);
  func_alpha2(par, sigma.d0, &a2, &da2_dsig, &da2);
  func_gamma2(par, y.d0, &g2);
  func_gamma1(par, y.d0, &g1);

  A_f = a2*f.d2sigma + a1*f.dsigma + a0*f.d0 + g2*f.d2y + g1*f.dy;
    
  return A_f;
}
//--------------------------------------------------------------------------------------
double complex HyperboloidalEffectiveSource(parameters par, int iDom, int j1, int j2){
  double complex hyp_S_eff=0;
  
  if(iDom == par.Dom_ptcl){
    int N1_Seff=par.N1_LoadSeff, N2_Seff=par.N2_LoadSeff;
    
    double real_hS_eff, imag_hS_eff, chi_1=par.grid_chi_1[iDom][j1], chi_2=par.grid_chi_2[iDom][j2];

    real_hS_eff = Clenshaw_Chebyshev_2D(par.Re_cheb_Seff, N1_Seff, N2_Seff, chi_1, chi_2);
    imag_hS_eff = Clenshaw_Chebyshev_2D(par.Im_cheb_Seff, N1_Seff, N2_Seff, chi_1, chi_2);

    hyp_S_eff = real_hS_eff + I*imag_hS_eff;
  }  

  return hyp_S_eff;
}