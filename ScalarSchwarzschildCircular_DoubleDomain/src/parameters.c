#include "Solve_PDE.h"

void set_parameters(parameters *par, int N, int nbar, int m, double r0){
	//INPUT PHYSICAL PARAMETERS------------------------------------
	

	(*par).r0_over_M = r0; //Orbital radius in units of M
	(*par).rh_over_M = 2; //Horizon radius in units of M
	
	(*par).nbar = nbar; //Puncture order
	(*par).lmax = 60; //Max l modes in the Boundary data for Retarded Field

	double f0, df0_dr, rh_over_r0 = (*par).rh_over_M/(*par).r0_over_M, r0_over_rh = 1/rh_over_r0;
		
	(*par).eta = sqr(r0_over_rh)*sqrt(1.-rh_over_r0)/(1.+ r0_over_rh); //Size of Excition Region in R = eta rh, with R = radius in particle' s frame
		
	func_f( (*par).r0_over_M, &f0, &df0_dr);
	(*par).f0 = f0; 
	(*par).r_plus_over_M  = (*par).r0_over_M + sqrt(f0)* (*par).eta * (*par).rh_over_M; //Effective Source outter boundary in units of M
	(*par).r_minus_over_M = (*par).r0_over_M - sqrt(f0)* (*par).eta * (*par).rh_over_M; //Effective Source inner boundary in units of M	


	(*par).spin = 0; //Field's spin weight
	(*par).q = 1., //Particle's scalar Charge
	
	(*par).m = m; //Azimutal Mode

	(*par).lambda_over_M = 4.; //Hyperboloidal Length Scale in units of M
	(*par).rho_0 = 0.5; //Hyperboloidal Radial Gauge
	(*par).rho_1 = 0. ; //Hyperboloidal Radial Gauge

	double sigma0, dsigma_dr_0, sigma_minus, dsigma_dr_minus, sigma_plus, dsigma_dr_plus;
	func_sigma_of_r(*par, (*par).r0_over_M, &sigma0, &dsigma_dr_0);
	func_sigma_of_r(*par, (*par).r_plus_over_M, &sigma_minus, &dsigma_dr_minus);
	func_sigma_of_r(*par, (*par).r_minus_over_M, &sigma_plus, &dsigma_dr_plus);
	(*par).sigma0 = sigma0;
	(*par).sigma_minus = sigma_minus;
	(*par).sigma_plus =  sigma_plus;

	//RESOLUTION INPUT DATA ------------------------------------
	(*par).N1_PuncSeff = 50;
	(*par).N2_PuncSeff = 50;
	(*par).prec = 63.81836;
	//63.81836;
	//159.54590;//63.81836;//31.90918;////31.90918;//15.95459; //Precision in Mathematica Notebook
	
	
	//INPUT NUMERICAL RESOLUTION------------------------------------
	//-1:Debug Error in function F
	// 0: SVN Decomposition
	// 1: Direct LU Inversion
	// 2: Iterative BiCGStab w/ FinDif Preconditioner
	
	(*par).SOLVER_METHOD = 2;
		
	int N1, N2;
	N1=N/2;
	N2=N;

	(*par).idom_particle = 1;
	(*par).CoordMap_FLAG = 2;
	int iDom;
	//------- Set Up Domain 0--------------------------------
	iDom = 0;
	//DIRECTION 1
	sprintf((*par).grid_1[iDom], "Lobatto"); //Grid 
	(*par).N1[iDom]   =  N1 ; //Resolution 
	
	(*par).AnMR_x_boundary_1[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_1[iDom] = 0.; //AnMR Parameter
	
	//DIRECTION 2
	sprintf((*par).grid_2[iDom], "Lobatto"); //Grid
	(*par).N2[iDom]   =  N2 ; //Resolution 


	(*par).AnMR_x_boundary_2[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_2[iDom] = 0.; //AnMR Parameter
	//------------------------------------------------------

	//------- Set Up Domain 1--------------------------------
	iDom = 1;
	//DIRECTION 1
	sprintf((*par).grid_1[iDom], "Lobatto"); //Grid 
	(*par).N1[iDom]   =  N1 ; //Resolution 

	(*par).AnMR_x_boundary_1[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_1[iDom] = 0.; //AnMR Parameter
	
	//DIRECTION 2
	sprintf((*par).grid_2[iDom], "Radau_RHS"); //Grid
	// sprintf((*par).grid_2[iDom], "Lobatto"); //Grid
	(*par).N2[iDom]   =  N2 ; //Resolution 

	(*par).AnMR_x_boundary_2[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_2[iDom] = 0.; //AnMR Parameter
	//------------------------------------------------------

	sprintf((*par).SimName, "Coord%d/r0_over_M_%3.5f/eta_%2.2f/N_input_%d/Prec_%lf/nbarMax_%d/m_%d/N1_%d_N2_%d/", (*par).CoordMap_FLAG,(*par).r0_over_M, (*par).eta,(*par).N2_PuncSeff, (*par).prec, (*par).nbar, (int)(*par).m, N1, N2 );// Simulation Name
	// sprintf((*par).SimName, "TestField_Reg/nbar_%d", nbar );// Simulation Name
	

	//DERIVED PHYSICAL PARAMETERS---------------------------------
	(*par).E0 = f0/sqrt(1-3./(*par).r0_over_M);
	(*par).L0_over_M = sqrt((*par).r0_over_M)/sqrt(1-3./(*par).r0_over_M);
	(*par).M_Omega0 = sqrt(1./(*par).r0_over_M)/(*par).r0_over_M;

	//Particle's orbital frequency
	(*par).s = -I* (*par).m * (*par).lambda_over_M * (*par).M_Omega0;



	//DERIVED NUMERICAL PARAMETERS------------------------------
	Get_Index_Arrays(par);
	Construct_grid(par);

	(*par).n_2D = nFields*(*par).n1n2[nDom];
	(*par).ntotal = nFields*(*par).n1n2[nDom] + nPar;
    (*par).Ntotal = (*par).ntotal-1;


    char data_dir[500], cp_par[500];
    sprintf(data_dir,"data/%s",(*par).SimName);
    create_directory(data_dir);	
	sprintf(cp_par, "cp src/parameters.c data/%s/",(*par).SimName) ;
	int func_out =  system(cp_par);
	if(func_out){ fprintf(stderr, "Error in set_parameters: Unable to copy paraemeter file\n"); exit(1);}

	(*par).i_omp = omp_get_thread_num();
    char fn[500];
    sprintf(fn, "output_thread_%d.txt",(*par).i_omp);
    // (*par).fout = fopen(fn, "a");

	return;
}
