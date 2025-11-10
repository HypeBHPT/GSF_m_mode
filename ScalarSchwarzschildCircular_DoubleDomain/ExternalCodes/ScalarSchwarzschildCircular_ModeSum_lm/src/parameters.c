#include "Solve_ODE.h"

void set_parameters(parameters *par, double r0, int l, int m, int N){

	//INPUT PHYSICAL PARAMETERS------------------------------------
	sprintf((*par).Equation, "ReggeWheeler"); //Zerilli; BardeenPress;

	(*par).r0_over_M = r0; //Orbital radius in units of M

	(*par).spin = 0; //Field's spin weight
	(*par).q = 1.; //Particle's scalar Charge
	
	
	
	
	
	
	(*par).ell = l; //Angular Mode
	(*par).m = m; //Azimutal Mode

	(*par).rh_over_M = 2.;
	(*par).lambda_over_M = 4.; //Hyperboloidal Length Scale in units of M
	(*par).rho_0 = (*par).rh_over_M/(*par).lambda_over_M; //Hyperboloidal Radial Gauge
	(*par).rho_1 = 0. ; //Hyperboloidal Radial Gauge

	double r0_over_rh = (*par).r0_over_M/(*par).rh_over_M, 
		   rh_over_r0 = 1./r0_over_rh;	
	
	(*par).eta = sqr(r0_over_rh)*sqrt(1.-rh_over_r0)/(1.+ r0_over_rh);; //Size of Excition Region in R = eta rh, with R = radius in particle' s frame
	
	//INPUT NUMERICAL RESOLUTION------------------------------------
	(*par).effective_source_FLAG = 0; // (0) Analytical Expression; (1) Read File

	int N0, N1;
	N0=N;
	N1=N;


	(*par).N[0]    =  N0 ; //Resolution Domain 1
	(*par).N[1]    =  N1 ; //Resolution Domain 2


	sprintf((*par).grid[0], "Lobatto");
	sprintf((*par).grid[1], "Lobatto");

	
	(*par).AnMR_x_boundary[0]= 1.;
	(*par).AnMR_kappa[0] = 0.;

	(*par).AnMR_x_boundary[1]= -1.;
	double Alm = l<100? 0.5 : 1.50;
	(*par).AnMR_kappa[1] = Alm + 0.5*log((*par).r0_over_M);
	//0.;

	// sprintf((*par).SimName, "../../../InputData/r0_over_M_%.5lf/m_%d/ell_%d/", (*par).r0_over_M, (*par).m, (*par).ell);
	sprintf((*par).SimName, "r0_over_M_%.5lf/eta_%lf/m_%d/ell_%d/", (*par).r0_over_M, (*par).eta, (*par).m, (*par).ell);

	//DERIVED PHYSICAL PARAMETERS---------------------------------
	double f0, df0_dr, rh_over_M = 2;
	func_f( (*par).r0_over_M, &f0, &df0_dr);
	(*par).f0 = f0; 
	(*par).E0 = f0/sqrt(1-3./(*par).r0_over_M);
	(*par).L0_over_M = sqrt((*par).r0_over_M)/sqrt(1-3./(*par).r0_over_M);
	(*par).M_Omega0 = sqrt(1./(*par).r0_over_M)/(*par).r0_over_M;

	(*par).r_plus_over_M  = (*par).r0_over_M + sqrt(f0)* (*par).eta*rh_over_M; //Effective Source outter boundary in units of M
	(*par).r_minus_over_M = (*par).r0_over_M - sqrt(f0)* (*par).eta*rh_over_M; //Effective Source inner boundary in units of M

	double sigma0, sigma_minus, sigma_plus, dr_dsigma_0, ddd;
	func_sigma_of_r(*par, (*par).r0_over_M, &sigma0, &dr_dsigma_0);
	func_sigma_of_r(*par, (*par).r_plus_over_M, &sigma_minus, &ddd);
	func_sigma_of_r(*par, (*par).r_minus_over_M, &sigma_plus, &ddd);
	(*par).sigma_minus = sigma_minus;
	(*par).sigma_plus =  sigma_plus;
	
	(*par).sigma[0]=0.; //Location of Scri	

	(*par).sigma[1] = (*par).sigma0 = sigma0; //Location of Particle:

	(*par).sigma[nDom]=1.; //Location of Horizon
	


	//Particle's orbital frequency
	(*par).s = -I* (*par).m * (*par).lambda_over_M * (*par).M_Omega0;


	double c_lm = get_clm(l,m),
	 	   P_lm = plgndr((*par).ell, (*par).m, 0.), Y_lm = c_lm*P_lm,
		   kappa = -4*Pi*(*par).q*Y_lm/((*par).E0*sqr((*par).r0_over_M)),
		   rho_r0, drho_r0_dsigma, d2rho_r0_dsigma, beta_r0, dbeta_r0_dsigma;
	double complex Z_r0, dlnZ_r0_dsigma;

	func_rho(*par, sigma0, &rho_r0, &drho_r0_dsigma, &d2rho_r0_dsigma);
	func_beta(*par, sigma0, &beta_r0, &dbeta_r0_dsigma);
	func_Z( *par, sigma0, &Z_r0, &dlnZ_r0_dsigma);	

	(*par).bar_kappa = (*par).lambda_over_M*kappa*f0*sqr(rho_r0)/(beta_r0*Z_r0);

  	double complex alpha2, dalpha2_dr0, dalpha2_dsig0, jump_phi, jump_dphi;
  	func_alpha2(*par, sigma0, &alpha2, &dalpha2_dsig0, &dalpha2_dr0);

  	jump_phi   = 0.,
    jump_dphi  = (*par).bar_kappa/alpha2;

  
  //Function discontinuity particle [Field 0 - Re(phi)]
  (*par).jump_field[0][0] = creal(jump_phi);

  //Function discontinuity particle [Field 1 - Im(phi)]
  (*par).jump_field[1][0] = cimag(jump_phi);

  //Derivative discontinuity particle [Field 0 - Re(dphi)]
  (*par).jump_derivative[0][0] = creal(jump_dphi);

  //Derivative discontinuity particle [Field 1 - Im(dphi)]
  (*par).jump_derivative[1][0] = cimag(jump_dphi);





	//DERIVED NUMERICAL PARAMETERS------------------------------
	Get_Index_Arrays(par);
	(*par).ntotal = nFields*((*par).mA[nDom]) + nPar;			
    (*par).Ntotal = (*par).ntotal-1;

    char data_dir[500];
    sprintf(data_dir,"data/%s",(*par).SimName);
    create_directory(data_dir);

    (*par).i_omp = omp_get_thread_num();
    char fn[500];
    sprintf(fn, "output_thread_%d.txt",(*par).i_omp);
    (*par).fout = fopen(fn, "a");
	
	return;
}
