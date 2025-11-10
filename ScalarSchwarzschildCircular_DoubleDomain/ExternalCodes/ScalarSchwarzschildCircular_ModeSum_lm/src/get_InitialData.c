#include "Solve_ODE.h"

// -------------------------------------------------------------------------------
void get_Initial_Data(parameters par, double z, double *V0, double *V0_sigma, double *W0)
{ 

  *V0   = 0;//z*(1.-z);
  *V0_sigma = 0;//1.-2*z;
  *W0   = 0;//0.;
  
  return;
	
}
// -------------------------------------------------------------------------------
double complex get_funcB(parameters par, double sigma, double complex s)
{ 
  double sigma2=sqr(sigma), V0, V0_sig, W0;  
  double complex B;
  
  get_Initial_Data(par, sigma, &V0, &V0_sig, &W0);
    
  B = (1.- 2.*sigma2*(1.+kappa*(1+kappa)*(1.-sigma)) )*V0_sig    
    -(1.+sigma*(1+kappa))*(1+kappa*(1+kappa)*(1-sigma))*(s*V0+W0)    
    -sigma*(2+kappa*(2-3*sigma)*(1+kappa))*V0;

  
  return B;
	
}
