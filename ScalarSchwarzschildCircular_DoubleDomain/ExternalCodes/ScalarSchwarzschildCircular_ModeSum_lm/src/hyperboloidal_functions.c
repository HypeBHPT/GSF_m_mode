#include "Solve_ODE.h"
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
 
  *h = -2*par.rho_0*(1./sigma - 2*log(sigma)/(par.lambda_over_M*par.rho_0)  );
  *dh_dsigma = 2*par.rho_0*(1./sqr(sigma) + 2./(sigma*par.lambda_over_M*par.rho_0)  );
  
  return ;
}
//--------------------------------------------------------------------------------------
void func_Omega(parameters par, double sigma, double *Omega, double *dOmega_dsigma){

  *Omega = sigma/par.lambda_over_M;
  *dOmega_dsigma = 1./par.lambda_over_M;
  return;
}
//--------------------------------------------------------------------------------------
void func_Z(parameters par, double sigma, double complex *Z, double complex *dlnZ_dsigma){
  double complex s=par.s;
  double Omega, dOmega, x, h, dx, dh;
  
  func_Omega(par, sigma, &Omega, &dOmega);
  func_tortoise_x(par, sigma, &x, &dx);
  func_height(par, sigma, &h, &dh);
  
  *Z = Omega*cexp( s*(x+h) );
  
  *dlnZ_dsigma = dOmega/Omega + s*(dx+dh);
  
  
  return;
}
//------------------------------------------------------------
void func_alpha2(parameters par, double sig, double complex *alpha2, double complex *dalpha2_dsig, double complex *dalpha2_dr0){
  double sig2 = sqr(sig);

  if(par.spin == 0 && par.ell==0 && cabs(par.s)==0. && sig==0.){
    *alpha2 = sig*(1-sig);
    *dalpha2_dsig = 2-3*sig;
    *dalpha2_dr0 = 0.;
  }
  else{
    *alpha2 = sig2*(1-sig);
    *dalpha2_dsig = 2*sig-3*sig2;
    *dalpha2_dr0 = 0.;
  }

  return;
}
//------------------------------------------------------------
void func_alpha1(parameters par, double sig, double complex *alpha1, double complex *dalpha1_dr0){
  double complex s, ds_dr0;
  double sig2 = sqr(sig), spin=par.spin;

  s=par.s;
  ds_dr0 = par.ds_dr0;
  if(spin ==0 && par.ell==0 && cabs(s)==0. && sig==0.){
    *alpha1 =  (2-3*sig);  
    *dalpha1_dr0 = 0.;
  }
  else{
    if(strcmp( par.Equation,"ReggeWheeler")==0){
      *alpha1 =  sig*(2-3*sig) + (1-2*sig2)*s;  
      *dalpha1_dr0 = ds_dr0*(1-2*sig2);  
    }
    else if(strcmp( par.Equation,"Zerilli") ==0){
      *alpha1 =  sig*(2-3*sig) + (1-2*sig2)*s;  
      *dalpha1_dr0 = ds_dr0*(1-2*sig2); 
    }
    else if(strcmp( par.Equation,"BardeenPress") ==0){
      *alpha1 = + sig*(2-3*sig + spin*(2.-sig) ) + (1-2*sig2)*s;
      *dalpha1_dr0 = ds_dr0*(1-2*sig2);
    }
    else{
      printf("Error in func_alpha1: par.Equation has to be: ReggeWheeler / Zerilli / BardeenPress\n Equation was: %s\n", par.Equation);
      exit(1);
    }
  }

  return;
}
//------------------------------------------------------------
void func_alpha0(parameters par, double sig, double complex *alpha0, double complex *dalpha0_dr0){
  double complex s, s2, ds_dr0;
  double l=1.*par.ell, spin=par.spin, spin2=sqr(spin);

  s=par.s;
  ds_dr0 = par.ds_dr0;
  s2=s*s;
  if(spin ==0 && par.ell==0 && cabs(s)==0. && sig==0.){
    *alpha0 = -1.;
    *dalpha0_dr0 = 0.;
  }
  else{
    if(strcmp( par.Equation,"ReggeWheeler") ==0){
      *alpha0 = -(1+sig)*s2 - 2*sig*s - ( l*(l+1) + sig*(1-spin2));
      *dalpha0_dr0 = - 2*s*ds_dr0*(1+sig) - 2*ds_dr0*sig;
    }
    else if(strcmp( par.Equation,"Zerilli") ==0){
      double n = (l-1.)*(l+2.)/2.;
      *alpha0 = -(1+sig)*s2 - 2*sig*s - ( sig + 2.*n*( 1. + 4.*n*(3.+2.*n)/sqr(2.*n+3.*sig) )/3. );
      *dalpha0_dr0 = - 2*s*ds_dr0*(1+sig) - 2*ds_dr0*sig;
    }
    else if(strcmp( par.Equation,"BardeenPress") ==0){
      *alpha0 =  -(1+sig)*s2 - (2*sig - spin*(1.-sig) )*s - ( l*(l+1) + (sig-spin)*(1.+ sig) );
      *dalpha0_dr0 = - 2*s*ds_dr0*(1+sig) - ds_dr0*(2*sig - spin*(1.-sig) );
  }
  else{
      printf("Error in operator_A: par.Equation has to be: ReggeWheeler / Zerilli / BardeenPress\n Equation was: %s\n", par.Equation);
      exit(1);
    }
}

  return;
}
//------------------------------------------------------------
double complex Diff_Operator(parameters par, complex_derivs W, int iDom, int i, int FLAG_dr0){
  
  double sigma, dx_dsigma ;
  get_sigma(par, iDom, i, &sigma, &dx_dsigma);
  

  
  double complex a2, a1, a0, da2, da2_dsig, da1, da0, A_f=0.;
  complex_sigma_derivs f;

  get_SigmaDerv_from_SpecDerv_complex(par, iDom, i, W, &f);
    
  func_alpha0(par, sigma, &a0, &da0);
  func_alpha1(par, sigma, &a1, &da1);
  func_alpha2(par, sigma, &a2, &da2_dsig, &da2);

  if(FLAG_dr0==0) A_f = a2*f.d2sigma + a1*f.dsigma + a0*f.d0;  
  else if(FLAG_dr0==1) A_f = da2*f.d2sigma + da1*f.dsigma + da0*f.d0;  
  else {fprintf(stderr, "Error in function Diff_Operator\nFLAG_dr0 = %d not implementend\n", FLAG_dr0); exit(-1);}
  
  return A_f;
}
// -------------------------------------------------------------------------------
double complex Get_ModeSum_Source(parameters par, int iDom, int i){
  
  double complex HypSource = 0.;
  return HypSource;
}
// -------------------------------------------------------------------------------
double complex Get_Extended_Source(parameters par, complex_derivs W, int iDom, int i){
  
  double complex HypSource;
    
  HypSource = -Diff_Operator(par, W, iDom, i, 1);


  return HypSource;
}