#include "Solve_PDE.h"
//   ---------------------------------------------------------------
double get_grid(char *grid, int i, int N){
 	
	double arg, x;
	if(strcmp( grid,"Radau_RHS" ) ==0){
	  arg = Pi*i/(2*N+1);
	  x = cos(2.*arg);
	}
	else if(strcmp( grid,"Radau_LHS")==0){
	  arg = Pih - Pi*i/(2*N+1);
	  x = cos(2.*arg);
	}
	else if(strcmp( grid,"Gauss")==0){
	  arg = Pih*(i+0.5)/(N+1);
	  x = cos(2.*arg);
	}
	else if(strcmp( grid,"Lobatto")==0){
	  arg = Pih*i/N;
	  x = cos(2.*arg);
	}
	else if(strcmp( grid,"Fourier")==0){
	  x = 2.*Pi*i/(2*N+1);
	}
	else{ 
	  printf("Error in get_grid: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n",  grid);
	  exit(1);
	}	
	return x;
}
//----------------------------------------------------------------
void Construct_grid(parameters *par){
	int iDom, j1, j2, N1, N2;

	for(iDom=0; iDom < nDom; iDom++){
		N1 = (*par).N1[iDom];	(*par).grid_chi_1[iDom] = dvector(0,N1);
		N2 = (*par).N2[iDom]; (*par).grid_chi_2[iDom] = dvector(0,N2);		
		
		for(j1=0; j1<=N1; j1++) (*par).grid_chi_1[iDom][j1] = get_grid((*par).grid_1[iDom] , j1, N1);
		for(j2=0; j2<=N2; j2++) (*par).grid_chi_2[iDom][j2] = get_grid((*par).grid_2[iDom] , j2, N2);
	}

	return;
}
//----------------------------------------------------------------
void free_grid(parameters *par){
	int iDom, N1, N2;

	for(iDom=0; iDom < nDom; iDom++){
		N1 = (*par).N1[iDom];	free_dvector( (*par).grid_chi_1[iDom], 0, N1);
		N2 = (*par).N2[iDom]; free_dvector( (*par).grid_chi_2[iDom], 0, N2);		
		

	}

	return;
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
    chi = x_boundary*(1. - 0.5*sinh(2*kappa)*asinh( (1 - x*x_boundary) )/kappa);
  }
  return chi;
}
//----------------------------------------------------------------
void get_sigma(parameters par, int iDom, double chi_1, double chi_2, func_derivs_2D *sigma){
	double x1, x1_B = par.AnMR_x_boundary_1[iDom], kappa1 = par.AnMR_kappa_1[iDom],
  			 x2, x2_B = par.AnMR_x_boundary_2[iDom], kappa2 = par.AnMR_kappa_2[iDom];

  if(kappa1!=0 || kappa2!=0){
  	fprintf(stderr, "WARNING in get_sigma: AnMR not yet double-checked (iDom = %d, kappa1=%lf, kappa2=%lf)\n", iDom, kappa1, kappa2);
  }

  x1 = get_x_from_chi(x1_B, kappa1, chi_1);
  x2 = get_x_from_chi(x2_B, kappa2, chi_2);

  double eta = par.eta, eta_2 = sqr(eta), 
         sqrt_f0 = sqrt(par.f0), f0 = par.f0,
         sig0 = par.sigma0, sig0_2 = sqr(sig0), sig0_3 = sig0_2*sig0,
         rho, rho_2, drho_dx2, d2rho_dx2x2;

  switch(iDom){
  
  	case 0: //Near-Particle Vacuum 
		//Map sigma(x1,x2)
		(*sigma).d0  = sig0/(1. - eta * sig0 * sqrt_f0 * x1);
		//derivative x1
		(*sigma).d1  = sig0_2 * eta * sqrt_f0/sqr(1. - eta * sig0 * sqrt_f0 *x1);
	  //(drho_B_dTheta * cos_Theta - rho_B * sin_Theta)*dTheta_dx1;	
		//derivative x1 x1
		(*sigma).d11 = 2 * sig0_3 * eta_2 * f0/pow(1. - eta * sig0 * sqrt_f0 *x1, 3);
		//derivative x1 x2 
		(*sigma).d12 = 0.; 
		//derivative x2
		(*sigma).d2 = 0.; 
		//derivative x2 x2
		(*sigma).d22 = 0.; 
  	break;

  	case 1:	//Near-Particle Excision 
  	rho = (1.+x2)/2.; 
  	rho_2 = sqr(rho);
  	drho_dx2 = 1./2;
  	d2rho_dx2x2 = 0.;

  	//Map sigma(x1,x2)
  	(*sigma).d0  = sig0/(1. - eta * sig0 * sqrt_f0 * rho * x1);
		//derivative x1
		(*sigma).d1  = sig0_2 * eta * sqrt_f0 * rho/sqr(1. - eta * sig0 * sqrt_f0 * rho * x1);
		//derivative x1 x1
		(*sigma).d11  = 2 * sig0_3 * eta_2 * f0 * rho_2/pow(1. - eta * sig0 * sqrt_f0 * rho * x1, 3);
		//derivative x1 x2
		(*sigma).d12  =  sig0_2 * eta * sqrt_f0 * drho_dx2/sqr(1. - eta * sig0 * sqrt_f0 * rho * x1) 
                  + 2 * sig0_3 * eta_2 * f0 * rho * drho_dx2 * x1/pow(1. - eta * sig0 * sqrt_f0 * rho * x1, 3);
		//derivative x2
		(*sigma).d2  = sig0_2 * eta * sqrt_f0 * drho_dx2 * x1/sqr(1. - eta * sig0 * sqrt_f0 * rho * x1);
		//derivative x2 x2
		(*sigma).d22 =  sig0_2 * eta * sqrt_f0 * d2rho_dx2x2 * x1/sqr(1. - eta * sig0 * sqrt_f0 * rho * x1)
                 + 2 * sig0_3 * eta_2 * f0 * sqr(drho_dx2 * x1)/pow(1. - eta * sig0 * sqrt_f0 * rho * x1, 3);
	

		break;


  	default:
  		fprintf(stderr,"Error in get_sigma: iDom = %d does not exist", iDom);
  		exit(-1);

  }


	return;
}
//----------------------------------------------------------------
void get_y(parameters par, int iDom, double chi_1, double chi_2, func_derivs_2D *y){
	
	double x1, x1_B = par.AnMR_x_boundary_1[iDom], kappa1 = par.AnMR_kappa_1[iDom],
  			 x2, x2_B = par.AnMR_x_boundary_2[iDom], kappa2 = par.AnMR_kappa_2[iDom];

  if(kappa1!=0 || kappa2!=0){
  	fprintf(stderr, "Warning in get_y: AnMR not yet checked (iDom = %d, kappa1=%lf, kappa2=%lf)\n", iDom, kappa1, kappa2);
  }

  x1 = get_x_from_chi(x1_B, kappa1, chi_1);
  x2 = get_x_from_chi(x2_B, kappa2, chi_2);

  double eta = par.eta, eta_2 = sqr(eta), x1_2 = sqr(x1),
         sig0 = par.sigma0, sig0_2 = sqr(sig0),
         rho, drho_dx2, d2rho_dx2x2, rho_2, drho_2_dx2, d2rho_2_dx2x2;

   switch(iDom){
	  
	  case 0://Near-Particle Vacuum 
		//Map y(x1,x2)
		(*y).d0 = eta_2 * sig0_2 * 0.5*(1-x2)*( 1 - x1_2 ) +  0.5*(1+x2);
		//derivative x1
		(*y).d1 = - eta_2 * sig0_2 * (1-x2) * x1;
		//derivative x1x1
		(*y).d11 = - eta_2 * sig0_2 * (1-x2);
		//derivative x1 x2
		(*y).d12 = eta_2 * sig0_2 * x1;
		//derivative x2
		(*y).d2 = - eta_2 * sig0_2 * 0.5*( 1 - x1_2 ) +  0.5;
		//derivative x2 x2
		(*y).d22 = 0;
		break;

		case 1: //Near-Particle Excision 
			rho = (1+x2)/2.;
			drho_dx2 = 1./2;
			d2rho_dx2x2 = 0.;

			switch(par.CoordMap_FLAG){
				case 1:
					//Map y(x1,x2)
					(*y).d0 = eta_2 * sig0_2 * rho*(1 - x1_2);
					//derivative x1
					(*y).d1 = -2 * eta_2 * sig0_2 * rho * x1;
					//derivative x1x1
					(*y).d11 = -2 * eta_2 * sig0_2 * rho;
					//derivative x1 x2
					(*y).d12 = -2 * eta_2 * sig0_2 * drho_dx2 * x1;
					//derivative x2
					(*y).d2 = eta_2 * sig0_2 * drho_dx2 * (1 - x1_2);
					//derivative x2 x2
					(*y).d22 = eta_2 * sig0_2 * d2rho_dx2x2 * (1 - x1_2);
				break;

				case 2:
					rho_2 = sqr(rho);
					drho_2_dx2 = 2*rho*drho_dx2;
					d2rho_2_dx2x2 = 2*sqr(drho_dx2) + 2*rho*d2rho_dx2x2;

					//Map y(x1,x2)
					(*y).d0 = eta_2 * sig0_2 * rho_2*(1 - x1_2);
					//derivative x1
					(*y).d1 = -2 * eta_2 * sig0_2 * rho_2 * x1;
					//derivative x1x1
					(*y).d11 = -2 * eta_2 * sig0_2 * rho_2;
					//derivative x1 x2
					(*y).d12 = -2 * eta_2 * sig0_2 * drho_2_dx2 * x1;
					//derivative x2
					(*y).d2 = eta_2 * sig0_2 * drho_2_dx2 * (1 - x1_2);
					//derivative x2 x2
					(*y).d22 = eta_2 * sig0_2 * d2rho_2_dx2x2 * (1 - x1_2);
				break;

				default:
				fprintf(stderr,"Error in get_y: CoordMap_FLAG = %d does not exist", par.CoordMap_FLAG);
				exit(-1);
			}
		break;

	  default:
		fprintf(stderr,"Error in get_y: iDom = %d does not exist", iDom);
		exit(-1);

	}

	return;
}
//----------------------------------------------------
void MapFirstDerivative_Dom_a_to_Dom_b(parameters par, 
	int iDom_a, double chi_a_1, double chi_a_2, complex_derivs_2D W_a, 
	int iDom_b, double chi_b_1, double chi_b_2, complex_derivs_2D *W_b){
	
	func_derivs_2D sigma_a, sigma_b, y_a, y_b;
	double J_a;

	get_sigma(par, iDom_a, chi_a_1, chi_a_2, &sigma_a);
	get_y(par, iDom_a, chi_a_1, chi_a_2, &y_a);

	get_sigma(par, iDom_b, chi_b_1, chi_b_2, &sigma_b);	
	get_y(par, iDom_b, chi_b_1, chi_b_2, &y_b);

	J_a = sigma_a.d1*y_a.d2 - sigma_a.d2*y_a.d1;

	
	(*W_b).d1 =(  ( sigma_b.d1*y_a.d2 - sigma_a.d2*y_b.d1 )*W_a.d1  + ( sigma_a.d1*y_b.d1 - sigma_b.d1*y_a.d1 )*W_a.d2      )/J_a;
	
	(*W_b).d2 =(  ( sigma_b.d2*y_a.d2 - sigma_a.d2*y_b.d2 )*W_a.d1  + ( sigma_a.d1*y_b.d2 - sigma_b.d2*y_a.d1 )*W_a.d2      )/J_a;

	return;

}
//---------------------------------------------------------------
void get_chi1_chi2_From_sigma_y(parameters par, int iDom, func_derivs_2D sigma, func_derivs_2D y, double *chi_1, double *chi_2){
  double x1, x1_B = par.AnMR_x_boundary_1[iDom], kappa1 = par.AnMR_kappa_1[iDom],
         x2, x2_B = par.AnMR_x_boundary_2[iDom], kappa2 = par.AnMR_kappa_2[iDom];
         
  double sigma0 = par.sigma0, rho, dsigma = sigma.d0-sigma0,
         eta = par.eta, eta_2 = sqr(eta), sigma0_2 = sqr(sigma0), yB,
         sqrt_f0 = sqrt(par.f0), f0 = par.f0;

    switch(iDom){
	  
	  case 0:
  	x1 = (1./sigma0 - 1./sigma.d0)/(eta*sqrt_f0);

  	yB = eta_2 * sigma0_2 * ( 1. - sqr(x1) );
  	x2 = (2*y.d0 - (1.+yB))/(1.-yB);	
  	break;

  	case 1:
		switch(par.CoordMap_FLAG){
			case 1:
				rho = (sqrt_f0*y.d0*sigma.d0 + sqrt( f0*sqr(y.d0*sigma.d0) + sqr(2 *eta*(sigma.d0-sigma0)*sigma0)    ))/(2*sqrt_f0*eta_2*sigma0_2*sigma.d0);
			break;

			case 2:
				rho = sqrt(y.d0 + sqr(1. - sigma0/sigma.d0)/f0)/(eta*sigma0);
			break;

			default:
				fprintf(stderr,"Error in get_y: CoordMap_FLAG = %d does not exist\n", par.CoordMap_FLAG);
				exit(-1);
		}
  	
  	x1 =  (1. - sigma0/sigma.d0)/(eta*sqrt_f0*rho*sigma0);
  	x2 = 2*rho - 1.;

  	if(y.d0==0 && dsigma==0){
    	fprintf(stderr,"WARNING in get_chi1_chi2_From_sigma_y: Coordinate trasformation degenerates at sigma = sigma0 and y=0\n");
    	fprintf(stderr,"Interpolation maybe wrong as limit to singular point depends on direction\n");
    	x1 = 0.;
  	}

  	break;

  	default:
		fprintf(stderr,"Error get_chi1_chi2_From_sigma_y: iDom = %d does not exist", iDom);
		exit(-1);

  }

  *chi_1 = get_chi_from_x(x1_B, kappa1, x1);
  *chi_2 = get_chi_from_x(x2_B, kappa2, x2);
  return;
}
//---------------------------------------------------------------------
double get_Jacobian(parameters par, int iDom, int j1, int j2){
	int N1=par.N1[iDom], N2=par.N2[iDom];

	double chi_1 = get_grid(par.grid_1[iDom], j1, N1), chi_2 = get_grid(par.grid_2[iDom], j2, N2); 
	
	func_derivs_2D sigma, y;
	
	get_sigma(par, iDom, chi_1, chi_2, &sigma);  
  	get_y(par, iDom, chi_1, chi_2, &y);
	
	double det = sigma.d1*y.d2 - sigma.d2*y.d1;

	return det;
}
//---------------------------------------------------------------------
void get_PhysDerv_from_SpecDerv_complex(parameters par, int iDom, int j1, int j2, complex_derivs_2D W, complex_sigma_derivs *U){
  
  int N1=par.N1[iDom], N2=par.N2[iDom];
  double x_B_1 = par.AnMR_x_boundary_1[iDom], kappa_1 = par.AnMR_kappa_1[iDom], chi_1,
  		   x_B_2 = par.AnMR_x_boundary_2[iDom], kappa_2 = par.AnMR_kappa_2[iDom], chi_2,
  		   det, det_2;
  func_derivs_2D sigma, y;


  chi_1 = get_grid(par.grid_1[iDom], j1, N1);
  chi_2 = get_grid(par.grid_2[iDom], j2, N2); 

  get_sigma(par, iDom, chi_1, chi_2, &sigma);  
  get_y(par, iDom, chi_1, chi_2, &y);
  det = sigma.d1*y.d2 - sigma.d2*y.d1;
  det_2 = sqr(det);

  // if(fabs(det)<1.e-15){
  // 	fprintf(stderr, "Error in get_PhysDerv_from_SpecDerv_complex. Singular coordinate mapping: det = %3.15e\n", det);
  // 	exit(-1);
  // }
  
  	if(kappa_1!=0 && kappa_2==0){
    double J, J2, kappa=kappa_1, chi=chi_1, x_B=x_B_1;
    if(1.-sqr(x_B)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary_1 has to be -1 or 1.\n x_boundary_1 was: %lf\n", x_B);
      exit(1);
    }

    J= 2.*kappa* cosh(kappa* (1-chi*x_B))/sinh(2*kappa);
    J2 = -2.*sqr(kappa)*x_B* sinh(kappa* (1-chi*x_B))/sinh(2*kappa);

    W.d1 = W.d1/J;
    W.d11 = (W.d11 - J2* W.d1 )/sqr(J); //REMINDER NOTE: W.d1 WAS MODIFIED IN THE PREVIOUS LINE. EXPRESSION IS CORRECT AS IT TAKES THE NEW VALUE 
    W.d12 = W.d12/J;
  }
  else if(kappa_1==0 && kappa_2!=0){
  	double J, J2, kappa=kappa_2, chi=chi_2, x_B=x_B_2;
  	if(1.-sqr(x_B)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary_2 has to be -1 or 1.\n x_boundary_2 was: %lf\n", x_B);
      exit(1);
    }    

    J= 2.*kappa* cosh(kappa* (1-chi*x_B))/sinh(2*kappa);
    J2 = -2.*sqr(kappa)*x_B* sinh(kappa* (1-chi*x_B))/sinh(2*kappa);

    W.d2 = W.d2/J;
    W.d22 = (W.d22 - J2* W.d2 )/sqr(J); //REMINDER NOTE: W.d2 WAS MODIFIED IN THE PREVIOUS LINE. EXPRESSION IS CORRECT AS IT TAKES THE NEW VALUE 
    W.d12 = W.d12/J;
  }
  else if(kappa_1!=0 && kappa_2!=0){
    if(1.-sqr(x_B_1)!= 0. || 1.-sqr(x_B_2)!= 0. ){
      fprintf(stderr, "Error in get_x_from_chi: \n x_boundary has to be -1 or 1.\n x_boundaries were: %lf and  %lf\n", x_B_1, x_B_2);
      exit(1);
    }
    double J_1, J2_1, J_2, J2_2;

    J_1= 2.*kappa_1* cosh(kappa_1* (1-chi_1*x_B_1))/sinh(2*kappa_1);
    J2_1 = -2.*sqr(kappa_1)*x_B_1* sinh(kappa_1* (1-chi_1*x_B_1))/sinh(2*kappa_1);

    J_2= 2.*kappa_2* cosh(kappa_2* (1-chi_2*x_B_2))/sinh(2*kappa_2);
    J2_2 = -2.*sqr(kappa_2)*x_B_2* sinh(kappa_2* (1-chi_2*x_B_2))/sinh(2*kappa_2);

    W.d1 = W.d1/J_1;
    W.d11 = (W.d11 - J2_1* W.d1 )/sqr(J_1); //REMINDER NOTE: W.d1 WAS MODIFIED IN THE PREVIOUS LINE. EXPRESSION IS CORRECT AS IT TAKES THE NEW VALUE 
    W.d2 = W.d2/J_2;
    W.d22 = (W.d22 - J2_2* W.d2 )/sqr(J_2); //REMINDER NOTE: W.d2 WAS MODIFIED IN THE PREVIOUS LINE. EXPRESSION IS CORRECT AS IT TAKES THE NEW VALUE 
    W.d12 = W.d12/(J_1*J_2);
  }

  //Map derivatives w.r.t to (x1,x2) [AnMr] into derivatives w.r.t to (sigma,y) [Hyperboloidal]
  (*U).d0 = W.d0;
  
  (*U).dsigma = (y.d2 * W.d1 - y.d1 * W.d2)/det;

  (*U).dy = (-sigma.d2*W.d1 + sigma.d1*W.d2)/det;
  
  (*U).d2sigma = (  
  	sqr(y.d2)*W.d11 - 2*y.d1*y.d2*W.d12 + sqr(y.d1)*W.d22
  	- ( 
  			( sqr(y.d1)*y.d22 - 2*y.d1*y.d2*y.d12 + sqr(y.d2)*y.d11 ) * (*U).dy
  		+ ( sqr(y.d1)*sigma.d22 - 2*y.d1*y.d2*sigma.d12 + sqr(y.d2)*sigma.d11 ) * (*U).dsigma
  		)  
  	)/det_2;
  
  (*U).dy = (-sigma.d2*W.d1 + sigma.d1*W.d2)/det;
  
  (*U).d2y = (  
  	sqr(sigma.d2)*W.d11 - 2*sigma.d1*sigma.d2*W.d12 + sqr(sigma.d1)*W.d22
  	- ( 
  			( sqr(sigma.d1)*y.d22 - 2*sigma.d1*sigma.d2*y.d12 + sqr(sigma.d2)*y.d11 ) * (*U).dy
  		+ ( sqr(sigma.d1)*sigma.d22 - 2*sigma.d1*sigma.d2*sigma.d12 + sqr(sigma.d2)*sigma.d11 ) * (*U).dsigma
  		)  
  	)/det_2;
  
  (*U).d2sigmay = (  
  	-(sigma.d2*y.d2*W.d11 - (sigma.d1*y.d2 + sigma.d2*y.d1)*W.d12 + sigma.d1*y.d1*W.d22)
  	+( 
  			( sigma.d1*y.d1*y.d22 - (sigma.d1*y.d2 + sigma.d2*y.d1)*y.d12 + sigma.d2*y.d2*y.d11 ) * (*U).dy
  		+ ( sigma.d1*y.d1*sigma.d22 - (sigma.d1*y.d2 + sigma.d2*y.d1)*sigma.d12 + sigma.d2*y.d2*sigma.d11 ) * (*U).dsigma
  		)  
  	)/det_2;

  return;
}
//---------------------------------------------------
int total_grid_points(char *grid, int N){
  int n;
  
  	  if(strcmp( grid,"Fourier" ) ==0)
	    n=2*N+1;
	  else
	    n=N+1;
	  
  
  return n;
}

//   ---------------------------------------------------------------
double get_AngleIncrement(char *grid, int N){
       
 	
	double arg;
	if(strcmp( grid,"Radau_RHS" ) ==0){
	  arg = 2*Pi/(2*N+1);	  
	}
	else if(strcmp( grid,"Radau_LHS")==0){
	  arg = -2*Pi/(2*N+1);	  
	  
	}
	else if(strcmp( grid,"Gauss")==0){
	  arg = Pi/(N+1);
	}
	else if(strcmp( grid,"Lobatto")==0){
	  arg = Pi/N;	  
	}
	else if(strcmp( grid,"Fourier")==0){
	  arg = 2.*Pi/(2*N+1);
	}
	else{ 
	  printf("Error in get_AngleIncrement: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n",  grid);
	  exit(1);
	}
	
	
	
	return arg;
}
