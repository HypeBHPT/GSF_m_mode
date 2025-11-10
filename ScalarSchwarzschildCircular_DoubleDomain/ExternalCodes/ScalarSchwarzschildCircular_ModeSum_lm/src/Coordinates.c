#include "Solve_ODE.h"

double get_ChebyshevGrid_x(int i, int N, char *grid){

  double h;
         if(strcmp( grid,"Radau_RHS") ==0) 
      h=2.*Pi/(2*N+1);      
    else if(strcmp( grid,"Radau_LHS")==0)
      h=2.*Pi/(2*N+1);
    else if(strcmp( grid,"Gauss")==0)
      h=Pi/(N+1);
    else if(strcmp( grid,"Lobatto")==0)
      h=Pi/(N);
    else{
      fprintf(stderr, "Error in get_ChebyshevGrid_x: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto\n grid was: %s\n", grid);
      exit(1);
    }
  
  double x_ChebGrid = cos(i*h);

  return x_ChebGrid;
}
//----------------------------------------------------------------
double get_x_from_chi(double x_boundary, double kappa, double chi){
  //Input: x_boundary = -1 (Analytical Mesh Refinement at domain's Left Boundary)
  //       x_boundary =  1 (Analytical Mesh Refinement at domain's Right Boundary)

  double x;

  if(kappa==0.){
    x = chi;
  }
  else{
    if(1.-sqr(x_boundary)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundary was: %lf\n", x_boundary);
      exit(1);
    }
    x = x_boundary*( 1. - 2.* sinh(kappa* (1-chi*x_boundary))/sinh(2*kappa) );
}
  return x;
}
//----------------------------------------------------------------
double get_chi_from_x(double x_boundary, double kappa, double x){
  //Input: x_boundary = -1 (Analytical Mesh Refinement at domain's Left Boundary)
  //       x_boundary =  1 (Analytical Mesh Refinement at domain's Right Boundary)

  double chi;

  if(kappa==0.){
    chi = x;
  }
  else{
    if(1.-sqr(x_boundary)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundary was: %lf\n", x_boundary);
      exit(1);
    }
     chi = (1 - asinh(0.5*sinh(2*kappa)*(1. - x/x_boundary ))/kappa)/x_boundary ;
}
  return chi;
}
//----------------------------------------------------------------
void get_sigma(parameters par, int iDom, int i, double *sigma, double *dx_dsigma){
	
	int N=par.N[iDom-1];
  double chi, x, x_B = par.AnMR_x_boundary[iDom -1], kappa = par.AnMR_kappa[iDom-1];


  chi = get_ChebyshevGrid_x(i, N, par.grid[iDom-1]);
  x = get_x_from_chi(x_B, kappa, chi);

	*sigma = 0.5*(par.sigma[iDom]*(1.+ x) + par.sigma[iDom-1]*(1.- x));
	*dx_dsigma=2./(par.sigma[iDom]-par.sigma[iDom-1]);
	
	return;
}
//---------------------------------------------------------------------
double get_x_from_sigma(parameters par, int iDom, double sigma){
  
  double x, 
        sigma_sum = par.sigma[iDom] + par.sigma[iDom-1] , 
        sigma_dif = par.sigma[iDom] - par.sigma[iDom-1] ; 

  x =  (2*sigma - sigma_sum)/sigma_dif;
  return x;
}
//---------------------------------------------------------------------
void get_SigmaDerv_from_SpecDerv_complex(parameters par, int iDom, int i, complex_derivs W, complex_sigma_derivs *U){
  
  int N=par.N[iDom-1];
  double x_B = par.AnMR_x_boundary[iDom -1], kappa = par.AnMR_kappa[iDom-1], chi, sigma, dx_dsigma, dx_dsigma_2;

  chi = get_ChebyshevGrid_x(i, N, par.grid[iDom-1]);
  get_sigma(par, iDom, i, &sigma, &dx_dsigma);
  dx_dsigma_2 = sqr(dx_dsigma);
  
  //Map derivatives w.r.t to chi [spectral] into derivatives w.r.t to x [AnMr]
  if(kappa==0.){
    (*U).dsigma = W.d1;
    (*U).d2sigma = W.d11;
  }
  else{
    if(1.-sqr(x_B)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundary was: %lf\n", x_B);
      exit(1);
    }
    double J, J2;

    J= 2.*kappa* cosh(kappa* (1-chi*x_B))/sinh(2*kappa);
    J2 = -2.*sqr(kappa)*x_B* sinh(kappa* (1-chi*x_B))/sinh(2*kappa);

    (*U).dsigma = W.d1/J;
    (*U).d2sigma = (W.d11 - J2* (*U).dsigma )/sqr(J);
  }

  //Map derivatives w.r.t to x [AnMr] into derivatives w.r.t to sigma [Hyperboloidal]
  (*U).d0 = W.d0;
  (*U).dsigma *= dx_dsigma;
  (*U).d2sigma *= dx_dsigma_2;

  return;
}



//-------------------------------------------------------------
// COORDINATE TRANSFORMATION TO ENHANCE ACCURACY
// OF TERMS z*ln(z). MAP INTO NEW COORDINATE zeta
//-------------------------------------------------------------
double func_ZlnZ(double z){
	double ZlnZ = z==0.? 0.: z*log(z);
	return ZlnZ;
} 
//----------------------------------------------------
double z_of_zeta(double zeta){
  
  double z;
  
  if(flag_NewCoord==0)
    z= zeta;
  else
    z = zeta==0.? 0: exp(1.-1./zeta);
  
  return z;  
}
//----------------------------------------------------
double zeta_of_z(double z){
  
  double zeta;
  
  
  if(flag_NewCoord==0)
    zeta = z;
  else
    zeta = z==0? 0.: 1./(1.-log(z));
  
  return zeta;  
}
