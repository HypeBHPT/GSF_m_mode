#include "Solve_PDE.h"

void set_parameters(parameters *par, int N, int nbar, int m, double M_Omega0, double a_over_M){
	(*par).TEST_Func_FLAG = 0;	//0: No Test Function; 1: Test Function
	
	//INPUT PHYSICAL PARAMETERS------------------------------------
	(*par).rho_min=0.0;
	(*par).FLAG_Trajec =  (M_Omega0 >= 0) ? +1 : -1;; //Trajectory direction: +1=prograde, -1=retrograde
	(*par).a_over_M = a_over_M; //Black Hole spin in units of M
	(*par).r0_over_M =  pow(1./M_Omega0 - (*par).FLAG_Trajec*(*par).a_over_M,2./3 );//Particle's orbital radius in units of M
	
	(*par).rh_over_M = 1. + sqrt( 1. - sqr((*par).a_over_M) ); //Event Horizon radius in units of M
	(*par).rC_over_M = 1. - sqrt( 1. - sqr((*par).a_over_M)); //Cauchy Horizon radius in units of M
	(*par).kappa = (*par).a_over_M/(*par).rh_over_M; //Black Hole spin in units of rh
	(*par).sigma_hrz = 1.;
	(*par).sigma_Cauchy = 1./sqr((*par).kappa);
	
	(*par).lambda_over_M =(*par).rh_over_M; //Hyperboloidal Length Scale in units of M
	(*par).rho_1 = 0. ; // Radial fixing Gauge
	//	   (*par).rh_over_M/(*par).lambda_over_M * sqr((*par).kappa) // Cauchy fixing Gauge
	(*par).rho_0 = (*par).rh_over_M/(*par).lambda_over_M - (*par).rho_1; //Hyperboloidal Radial Gauge
	if((*par).rho_1 != 0.){
		printf("ATTENTION: CAUCHY FIXING GAUGE NOT FULLY TESTED\n");
		pause();
	}
	
	(*par).spin = 0; //Field's spin weight
	(*par).q = 1.; //Particle's scalar Charge
	(*par).m = m; //Azimutal Mode
	(*par).nbar = nbar; //Puncture order
	double M = 1;
	effsource_init(M, a_over_M);
	

	double f0, df0_dr, sigma0, dsigma_dr_0;
	func_sigma_of_r(*par, (*par).r0_over_M, &sigma0, &dsigma_dr_0);
	(*par).sigma0 = sigma0;
	func_f(*par, (*par).r0_over_M, &f0, &df0_dr);
	(*par).f0 = f0; 
	
	
	// double rh_over_r0 = (*par).rh_over_M/(*par).r0_over_M, r0_over_rh = 1/rh_over_r0;
		
	(*par).eta = (1 - sigma0)/(sigma0*sqrt(f0)*(1+sigma0)); //Size of Excition Region in R = eta rh, with R = radius in particle' s frame
		
	(*par).r_plus_over_M  = (*par).r0_over_M + sqrt(f0)* (*par).eta * (*par).rh_over_M; //Effective Source outter boundary in units of M
	(*par).r_minus_over_M = (*par).r0_over_M - sqrt(f0)* (*par).eta * (*par).rh_over_M; //Effective Source inner boundary in units of M	

	double sigma_minus, dsigma_dr_minus, sigma_plus, dsigma_dr_plus;
	
	func_sigma_of_r(*par, (*par).r_plus_over_M, &sigma_minus, &dsigma_dr_minus);
	func_sigma_of_r(*par, (*par).r_minus_over_M, &sigma_plus, &dsigma_dr_plus);
	
	(*par).sigma_minus = sigma_minus;
	(*par).sigma_plus =  sigma_plus;

	//RESOLUTION INPUT DATA ------------------------------------
	sprintf((*par).grid_1_PuncSeff, "Gauss"); //Grid 	
	(*par).N1_PuncSeff = 50;

	sprintf((*par).grid_2_PuncSeff, "Radau_RHS"); //Grid 
	(*par).N2_PuncSeff = 50;
	(*par).prec = 63.81836;
	//63.81836;
	//159.54590;//63.81836;//31.90918;////31.90918;//15.95459; //Precision in Mathematica Notebook
	
	
	//INPUT NUMERICAL RESOLUTION------------------------------------
	//-1:Debug Error in function F
	// 0: SVN Decomposition
	// 1: Direct LU Inversion
	// 2: Iterative BiCGStab w/ FinDif Preconditioner
	
	(*par).SOLVER_METHOD = 1;
		
	int N1, N2;
	N1=N/2;
	N2=N;

	(*par).Dom_scri = 0;
	(*par).Dom_bulk = 1;
	(*par).Dom_ptcl = 2;
	(*par).Dom_hrzn = 3;

	
	(*par).CoordMap_FLAG = 2;
	
	int iDom;
	//------- Set Up Domain 0--------------------------------
	iDom = (*par).Dom_scri;
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
	iDom = (*par).Dom_bulk;
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

	//------- Set Up Domain 2--------------------------------
	iDom = (*par).Dom_ptcl;
	//DIRECTION 1
	sprintf((*par).grid_1[iDom], "Lobatto"); //Grid 
	(*par).N1[iDom]   =  N1 ; //Resolution 
	
	(*par).AnMR_x_boundary_1[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_1[iDom] = 0.; //AnMR Parameter
	
	//DIRECTION 2
	sprintf((*par).grid_2[iDom], "Radau_RHS"); //Grid
	(*par).N2[iDom]   =  N2 ; //Resolution 
	
	(*par).AnMR_x_boundary_2[iDom]= 1.; //AnMR at Left (-1) or Right (1) Boundary
	(*par).AnMR_kappa_2[iDom] = 0.; //AnMR Parameter
	//------------------------------------------------------
	//------- Set Up Domain 3--------------------------------
	iDom = (*par).Dom_hrzn;
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
		

	if((*par).TEST_Func_FLAG==0)
		sprintf((*par).SimName, "M_Omega0_%3.5f/r0_over_M_%3.5f_a_over_M_%2.2f/N_input_%d/Prec_%lf/nbarMax_%d/m_%d/N1_%d_N2_%d/", M_Omega0, (*par).r0_over_M, (*par).a_over_M, (*par).N2_PuncSeff, (*par).prec, (*par).nbar, (int)(*par).m, N1, N2 );
	else
		sprintf((*par).SimName, "M_Omega0_%3.5f/r0_over_M_%3.5f_a_over_M_%2.2f/TestFunction/m_%d/N1_%d_N2_%d/", M_Omega0,(*par).r0_over_M, (*par).a_over_M,  (int)(*par).m, N1, N2 );
	


	//DERIVED PHYSICAL PARAMETERS---------------------------------
	// Particle's orbital frequency in units of 1/M
	(*par).M_Omega0 = 1.0 / ( pow((*par).r0_over_M, 1.5) + (*par).FLAG_Trajec * (*par).a_over_M );
	if( fabs( (*par).M_Omega0 - M_Omega0) > 1e-10 ){
		printf("Error in set_parameters: M_Omega0 inconsistent with r0_over_M and a_over_M\n");
		printf(" M_Omega0 input = %lf, M_Omega0 computed = %lf\n", M_Omega0, (*par).M_Omega0);
		exit(1);
	}

	// Specific energy E0 (dimensionless)
	(*par).E0 = ( (*par).r0_over_M - 2.0 + (*par).FLAG_Trajec * (*par).a_over_M / sqrt((*par).r0_over_M) )
            / sqrt( (*par).r0_over_M * ( (*par).r0_over_M - 3.0 + 2.0 * (*par).FLAG_Trajec * (*par).a_over_M / sqrt((*par).r0_over_M) ) );

	// Specific angular momentum L0/M (dimensionless)
	(*par).L0_over_M = (*par).FLAG_Trajec * ( sqr((*par).r0_over_M) - 2.0 * (*par).FLAG_Trajec * (*par).a_over_M * sqrt((*par).r0_over_M) + sqr((*par).a_over_M) ) / ( (*par).r0_over_M * sqrt( (*par).r0_over_M - 3.0 + 2.0 * (*par).FLAG_Trajec * (*par).a_over_M / sqrt((*par).r0_over_M) ) );

	struct coordinate xp;
	xp.r = (*par).r0_over_M;
	xp.theta = Pih;
	xp.phi = 0;
	xp.t = 0.;
	double ur = 0; //Particle's 4 velocity u^r
	effsource_set_particle(&xp, (*par).E0, (*par).L0_over_M, ur);
	
	
	
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
    (*par).fout = fopen(fn, "a");

	return;
}
