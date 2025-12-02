#include "Solve_PDE.h"
//-------------------------------------------------------------------------------
void get_ComplexField(parameters par, int indx_Re, int indx_Im, derivs_2D W , complex_derivs_2D *z){

	(*z).d0 = W.d0[indx_Re] +I*W.d0[indx_Im];
	(*z).d1 = W.d1[indx_Re] +I*W.d1[indx_Im];
	(*z).d11 = W.d11[indx_Re] +I*W.d11[indx_Im];
	(*z).d2 = W.d2[indx_Re] +I*W.d2[indx_Im];
	(*z).d22 = W.d22[indx_Re] +I*W.d22[indx_Im];
	(*z).d12 = W.d12[indx_Re] +I*W.d12[indx_Im];

	return;
}
//-------------------------------------------------------------------------------
void FieldEquations(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
  complex_derivs_2D phi;
  double complex A_phi, Source_phi=0.;
	int indx_Re_phi = Index(par, iDom, 0, j1, j2),
		indx_Im_phi = Index(par, iDom, 1, j1, j2);
		  
	get_ComplexField(par, indx_Re_phi, indx_Im_phi, W , &phi);

	A_phi = Diff_Operator(par, phi, iDom, j1, j2, 0);	
	if(par.TEST_Func_FLAG==0)
		Source_phi = HyperboloidalEffectiveSource(par, iDom, j1, j2);
	else
		Source_phi = Test_Effective_Source(par, iDom, j1, j2);
	
		F[0] = creal(A_phi - Source_phi);
		F[1] = cimag(A_phi - Source_phi);


	return;
}
// -------------------------------------------------------------------------------
void LinearFieldEquations(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
 complex_derivs_2D Dphi;
  double complex A_Dphi;
	int indx_Re_phi = Index(par, iDom, 0, j1, j2),
		  indx_Im_phi = Index(par, iDom, 1, j1, j2);
		  
	get_ComplexField(par, indx_Re_phi, indx_Im_phi, DW , &Dphi);
	A_Dphi = Diff_Operator(par, Dphi, iDom, j1, j2, 0);	

	J[0] = creal(A_Dphi);
	J[1] = cimag(A_Dphi);
	

 return;
}
// -------------------------------------------------------------------------------	
// -------------------------------------------------------------------------------	
void BoundaryCondition(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	if(par.rho_min!=0){
		FieldEquations(par, W, Xpar, iDom, j1, j2, F);
		return;
	}
	// return;


	complex_derivs_2D w;
	complex_sigma_derivs phi;

	int indx_Re_phi = Index(par, iDom, 0, j1, j2),
			indx_Im_phi = Index(par, iDom, 1, j1, j2);

	get_ComplexField(par, indx_Re_phi, indx_Im_phi, W , &w);
	get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, w, &phi);

	
	
	switch(iDom){
		case 0: //Dom Scri
			if(iDom!=par.Dom_scri){
				fprintf(stderr, "Error in BoundaryCondition: idom = %d dom_scri=%d\n", iDom, par.Dom_scri);
				exit(-1);
			}

			if(j1 == 0)
				FieldDerivativeJumpEquations_Boundary_chi1(par, W, Xpar, par.Dom_scri, par.Dom_bulk, 0, par.N1[par.Dom_bulk], j2, F);
			else 
				FieldEquations(par, W, Xpar, iDom, j1, j2, F);

			
		break;	

		case 1: //Dom Bulk
			if(iDom!=par.Dom_bulk){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_bulk=%d\n", iDom, par.Dom_bulk);
				exit(-1);
			}
			if(j2 == par.N2[iDom])
				FieldDerivativeJumpEquations_Boundary_chi2(par, W, Xpar, par.Dom_bulk, par.Dom_ptcl, j1, par.N2[iDom], 0, F);
			else if(j1 == 0)
				FieldDerivativeJumpEquations_Boundary_chi1(par, W, Xpar, par.Dom_bulk, par.Dom_hrzn, 0, par.N1[par.Dom_hrzn],j2, F);
			else if(j1 == par.N1[iDom])
		 		FieldJumpEquations_Boundary_chi1(par, W, Xpar, par.Dom_bulk, par.Dom_scri, par.N1[iDom], 0, j2, F);
		 	
		 	else 
				FieldEquations(par, W, Xpar, iDom, j1, j2, F);
				
		break;

		case 2: //Dom Particle
			if(iDom!=par.Dom_ptcl){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_particle=%d\n", iDom, par.Dom_ptcl);
				exit(-1);
			}

			if(j2 == 0)
				FieldJumpEquations_Boundary_chi2(par, W, Xpar, par.Dom_ptcl, par.Dom_bulk, j1, 0, par.N2[par.Dom_bulk],F);
			else	
				FieldEquations(par, W, Xpar, iDom, j1, j2, F);
		break;
			printf("In BoundaryCondition dom %d\n", iDom); 	exit(-1);
		case 3: //Dom Horizon
			if(iDom!=par.Dom_hrzn){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_horizon=%d\n", iDom, par.Dom_hrzn);
				exit(-1);
			}
			
			if(j1 == par.N1[iDom])
				FieldJumpEquations_Boundary_chi1(par, W, Xpar, par.Dom_hrzn, par.Dom_bulk, par.N1[iDom], 0, j2, F);
			else 
				FieldEquations(par, W, Xpar, iDom, j1, j2, F);
		break;

		default:
			fprintf(stderr,"Error in BoundaryCondition: iDom = %d does not exist\n", iDom);
			exit(-1);

	}

	return;
}
// -------------------------------------------------------------------------------	
void LinearBoundaryCondition(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	if(par.rho_min!=0){
		LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
		return;
	}
	
	complex_derivs_2D Dw;
	complex_sigma_derivs Dphi;
		
	int indx_Re_phi = Index(par, iDom, 0, j1, j2),
		  indx_Im_phi = Index(par, iDom, 1, j1, j2);

	get_ComplexField(par, indx_Re_phi, indx_Im_phi, DW , &Dw);
	get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, Dw, &Dphi);

	switch(iDom){
		case 0: //Dom Scri
			if(iDom!=par.Dom_scri){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_scri=%d\n", iDom, par.Dom_scri);
				exit(-1);
			}

			if(j1 == 0)
				LinearFieldDerivativeJumpEquations_Boundary_chi1(par, W, DW, Xpar, DXpar, par.Dom_scri, par.Dom_bulk, 0, par.N1[par.Dom_bulk], j2, J);
			else 
				LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
		break;

		case 1: //Dom Bulk
			if(iDom!=par.Dom_bulk){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_bulk=%d\n", iDom, par.Dom_bulk);
				exit(-1);
			}
			if(j2 == par.N2[iDom])
				LinearFieldDerivativeJumpEquations_Boundary_chi2(par, W, DW, Xpar, DXpar, par.Dom_bulk, par.Dom_ptcl, j1, par.N2[iDom], 0, J);	
			else if(j1 == 0)
				LinearFieldDerivativeJumpEquations_Boundary_chi1(par, W, DW, Xpar, DXpar, par.Dom_bulk, par.Dom_hrzn, 0, par.N1[par.Dom_hrzn], j2, J);
		 	else if(j1 == par.N1[iDom])
				LinearFieldJumpEquations_Boundary_chi1(par, W, DW, Xpar, DXpar, par.Dom_bulk, par.Dom_scri, par.N1[iDom], 0, j2, J);		
		  		
		  else 
				LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);					
		break;

		case 2: //Dom Particle
			if(iDom!=par.Dom_ptcl){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_particle=%d\n", iDom, par.Dom_ptcl);
				exit(-1);
			}

			if(j2 == 0)
				LinearFieldJumpEquations_Boundary_chi2(par, W, DW, Xpar, DXpar, par.Dom_ptcl, par.Dom_bulk, j1, 0, par.N2[par.Dom_bulk], J);
			else 
				LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
		break;

		case 3: //Dom Horizon
			if(iDom!=par.Dom_hrzn){
				fprintf(stderr, "Error in LinearBoundaryCondition: idom = %d dom_horizon=%d\n", iDom, par.Dom_hrzn);
				exit(-1);
			}
			if(j1 == par.N1[iDom])
				LinearFieldJumpEquations_Boundary_chi1(par, W, DW, Xpar, DXpar, par.Dom_hrzn, par.Dom_bulk, par.N1[iDom], 0, j2, J);
			else 
				LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
		break;
		
		default:
			fprintf(stderr,"Error in BoundaryCondition: iDom = %d does not exist\n", iDom);
			exit(-1);

	}

	return;
	
}
//---------------------------------------------------------------------------
void FieldJumpEquations_Boundary_chi1(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *F){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1_a,j2),
			indx_b = Index(par, iDom_b, iF, j1_b,j2);
		F[iF] = W.d0[indx_b]-W.d0[indx_a];			    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void LinearFieldJumpEquations_Boundary_chi1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *J){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1_a,j2),
			  indx_b = Index(par, iDom_b, iF, j1_b,j2);
		J[iF] = DW.d0[indx_b] - DW.d0[indx_a];			    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void FieldDerivativeJumpEquations_Boundary_chi1(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *F){
	int iF;
	complex_derivs_2D U_a, U_b;

	double chi_1_Dom_a = par.grid_chi_1[iDom_a][j1_a],
				 chi_2_Dom_a = par.grid_chi_2[iDom_a][j2],
				 chi_1_Dom_b = par.grid_chi_1[iDom_b][j1_b],
				 chi_2_Dom_b = par.grid_chi_2[iDom_b][j2];

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1_a, j2),
			indx_b = Index(par, iDom_b, iF, j1_b, j2);	

		get_ComplexField(par, indx_a, indx_a, W , &U_a);
		get_ComplexField(par, indx_b, indx_b, W , &U_b);
		MapFirstDerivative_Dom_a_to_Dom_b(par,  iDom_a, chi_1_Dom_a, chi_2_Dom_a, U_a, 
																					  iDom_b, chi_1_Dom_b, chi_2_Dom_b, &U_a);
				
		F[iF] = U_b.d1 - U_a.d1;	    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void LinearFieldDerivativeJumpEquations_Boundary_chi1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1_a, int j1_b, int j2, double *J){
	int iF;
	complex_derivs_2D DU_a, DU_b;

	double chi_1_Dom_a = par.grid_chi_1[iDom_a][j1_a],
				 chi_2_Dom_a = par.grid_chi_2[iDom_a][j2],
				 chi_1_Dom_b = par.grid_chi_1[iDom_b][j1_b],
				 chi_2_Dom_b = par.grid_chi_2[iDom_b][j2];

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1_a, j2),
				indx_b = Index(par, iDom_b, iF, j1_b, j2);	

		get_ComplexField(par, indx_a, indx_a, DW , &DU_a);
		get_ComplexField(par, indx_b, indx_b, DW , &DU_b);
		MapFirstDerivative_Dom_a_to_Dom_b(par,  iDom_a, chi_1_Dom_a, chi_2_Dom_a, DU_a, 
																					  iDom_b, chi_1_Dom_b, chi_2_Dom_b, &DU_a);
				
		J[iF] = DU_b.d1 - DU_a.d1;		    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void FieldJumpEquations_Boundary_chi2(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *F){
	if(iDom_a !=par.Dom_ptcl){
		fprintf(stderr, "Error in function FieldJumpEquations_Boundary_chi2: dom_a = %d does not coincide with particle domain = %d\n", iDom_a, par.Dom_ptcl);
		exit(-1);
	}
	
	int iF;
	double chi_1_prtcl = par.grid_chi_1[par.Dom_ptcl][j1], jump;
	complex_derivs_2D U_a, U_b;
	complex_sigma_derivs phi_Res, phi_Ret;

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1, j2_a),
			indx_b = Index(par, iDom_b, iF, j1, j2_b);
		
		get_ComplexField(par, indx_a, indx_a, W , &U_a);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_a, j1, j2_a, U_a, &phi_Res);

		get_ComplexField(par, indx_b, indx_b, W , &U_b);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_b, j1, j2_b, U_b, &phi_Ret);

		if(par.TEST_Func_FLAG==0){
			double phiP = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc, par.N1_PuncSeff, chi_1_prtcl):
						   		   Clenshaw_Chebyshev(par.Im_cheb_phi_Punc, par.N1_PuncSeff, chi_1_prtcl);	
			jump = -phiP;			
		}
		else{
			func_derivs_2D sigma, y;			
			get_sigma(par,  iDom_a, chi_1_prtcl, par.grid_chi_2[par.Dom_ptcl][j2_a], &sigma);
			get_y(par,  iDom_a, chi_1_prtcl, par.grid_chi_2[par.Dom_ptcl][j2_a], &y);
			func_sigma_derivs_2D f = Test_Func(par, sigma.d0, y.d0);
			jump = 0;;//f.d0;
		}		

		F[iF] = phi_Res.d0 - phi_Ret.d0 - jump;			    	
	}
	
	// return;
}
//---------------------------------------------------------------------------
void LinearFieldJumpEquations_Boundary_chi2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *J){
	int iF;
	complex_derivs_2D DU_a, DU_b;
	complex_sigma_derivs Dphi_Res, Dphi_Ret;

	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1, j2_a),
			indx_b = Index(par, iDom_b, iF, j1, j2_b);
		get_ComplexField(par, indx_a, indx_a, DW , &DU_a);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_a, j1, j2_a, DU_a, &Dphi_Res);

		get_ComplexField(par, indx_b, indx_b, DW , &DU_b);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_b, j1, j2_b, DU_b, &Dphi_Ret);
		
		J[iF] = Dphi_Res.d0 - Dphi_Ret.d0;
	}
	
	// return;
}
//---------------------------------------------------------------------------
void FieldDerivativeJumpEquations_Boundary_chi2(parameters par, derivs_2D W, double *Xpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *F){
	if(iDom_b !=par.Dom_ptcl){
		fprintf(stderr, "Error in function FieldDerivativeJumpEquations_Boundary_chi2: dom_a = %d does not coincide with particle domain = %d\n", iDom_a, par.Dom_ptcl);
		exit(-1);
	}
	
	int iF;

	complex_derivs_2D U_a, U_b;
	double chi_1_prtcl = par.grid_chi_1[par.Dom_ptcl][j1];
	complex_sigma_derivs phi_Res, phi_Ret;

	double chi_1_Dom_a = par.grid_chi_1[iDom_a][j1],
		   chi_2_Dom_a = par.grid_chi_2[iDom_a][j2_a],
		   chi_1_Dom_b = par.grid_chi_1[iDom_b][j1],
		   chi_2_Dom_b = par.grid_chi_2[iDom_b][j2_b];
	
	double phiP_dsigma, phiP_dy, NormalDerv_phi_P, NormalDerv_phi_Res, NormalDerv_phi_Ret,
		   fac_Dom_a, fac_Dom_b, NormalDerv_jump;

	func_derivs_2D sigma_Dom_a, y_Dom_a, sigma_Dom_b, y_Dom_b;
	get_sigma(par, iDom_a, chi_1_Dom_a, chi_2_Dom_a, &sigma_Dom_a);
	get_y(par, iDom_a, chi_1_Dom_a, chi_2_Dom_a, &y_Dom_a);

	get_sigma(par, iDom_b, chi_1_Dom_b, chi_2_Dom_b, &sigma_Dom_b);
	get_y(par, iDom_b, chi_1_Dom_b, chi_2_Dom_b, &y_Dom_b);

	fac_Dom_a = sqrt( sqr(sigma_Dom_a.d1) + sqr(y_Dom_a.d1)  );
	fac_Dom_b = sqrt( sqr(sigma_Dom_b.d1) + sqr(y_Dom_b.d1)  );


	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1, j2_a),
			indx_b = Index(par, iDom_b, iF, j1, j2_b);
		
		get_ComplexField(par, indx_a, indx_a, W , &U_a);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_a, j1, j2_a, U_a, &phi_Ret);
		NormalDerv_phi_Res = ( - y_Dom_a.d1 * phi_Ret.dsigma + sigma_Dom_a.d1 * phi_Ret.dy )/fac_Dom_a;

		get_ComplexField(par, indx_b, indx_b, W , &U_b);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_b, j1, j2_b, U_b, &phi_Res );
		NormalDerv_phi_Ret = ( - y_Dom_b.d1 * phi_Res.dsigma +  sigma_Dom_b.d1*phi_Res.dy)/fac_Dom_b;

		if(par.TEST_Func_FLAG==0){		
			phiP_dsigma = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_sigma, par.N1_PuncSeff, chi_1_prtcl):
								Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_sigma, par.N1_PuncSeff, chi_1_prtcl);

			phiP_dy     = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_y, par.N1_PuncSeff, chi_1_prtcl):
								Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_y, par.N1_PuncSeff, chi_1_prtcl);
			NormalDerv_phi_P = ( -y_Dom_a.d1 * phiP_dsigma + sigma_Dom_a.d1*phiP_dy )/fac_Dom_a;
			NormalDerv_jump = - NormalDerv_phi_P;
		}
		else{
			func_derivs_2D sigma, y;			
			get_sigma(par,  iDom_b, chi_1_prtcl, par.grid_chi_2[par.Dom_ptcl][j2_b], &sigma);
			get_y(par,  iDom_b, chi_1_prtcl, par.grid_chi_2[par.Dom_ptcl][j2_b], &y);
			func_sigma_derivs_2D f = Test_Func(par, sigma.d0, y.d0);

			double NormalDerv_f = ( -y_Dom_a.d1 * f.dsigma + sigma_Dom_a.d1*f.dy )/fac_Dom_a;
			NormalDerv_jump = 0;//NormalDerv_f;

		}

		

		F[iF] = NormalDerv_phi_Res - NormalDerv_phi_Ret + NormalDerv_jump;
	}
	
	return;
}
//---------------------------------------------------------------------------
void LinearFieldDerivativeJumpEquations_Boundary_chi2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom_a, int iDom_b, int j1, int j2_a, int j2_b, double *J){
	int iF;

	complex_derivs_2D DU_a, DU_b;
	complex_sigma_derivs Dphi_Res, Dphi_Ret;

	double chi_1_Dom_a = par.grid_chi_1[iDom_a][j1],
		   chi_2_Dom_a = par.grid_chi_2[iDom_a][j2_a],
		   chi_1_Dom_b = par.grid_chi_1[iDom_b][j1],
		   chi_2_Dom_b = par.grid_chi_2[iDom_b][j2_b];
	
	double NormalDerv_Dphi_Res, NormalDerv_Dphi_Ret,
		   fac_Dom_a, fac_Dom_b;

	func_derivs_2D sigma_Dom_a, y_Dom_a, sigma_Dom_b, y_Dom_b;
	get_sigma(par, iDom_a, chi_1_Dom_a, chi_2_Dom_a, &sigma_Dom_a);
	get_y(par, iDom_a, chi_1_Dom_a, chi_2_Dom_a, &y_Dom_a);

	get_sigma(par, iDom_b, chi_1_Dom_b, chi_2_Dom_b, &sigma_Dom_b);
	get_y(par, iDom_b, chi_1_Dom_b, chi_2_Dom_b, &y_Dom_b);

	fac_Dom_a = sqrt( sqr(sigma_Dom_a.d1) + sqr(y_Dom_a.d1)  );
	fac_Dom_b = sqrt( sqr(sigma_Dom_b.d1) + sqr(y_Dom_b.d1)  );


	for(iF=0; iF<nFields; iF++){
		int indx_a = Index(par, iDom_a, iF, j1, j2_a),
			indx_b = Index(par, iDom_b, iF, j1, j2_b);	

		get_ComplexField(par, indx_a, indx_a, DW , &DU_a);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_a, j1, j2_a, DU_a, &Dphi_Ret);
		NormalDerv_Dphi_Res = ( -y_Dom_a.d1 * Dphi_Ret.dsigma + sigma_Dom_a.d1 * Dphi_Ret.dy )/fac_Dom_a;

		get_ComplexField(par, indx_b, indx_b, DW , &DU_b);
		get_PhysDerv_from_SpecDerv_complex(par, iDom_b, j1, j2_b, DU_b, &Dphi_Res);
		NormalDerv_Dphi_Ret = ( - y_Dom_b.d1 * Dphi_Res.dsigma +  sigma_Dom_b.d1 * Dphi_Res.dy)/fac_Dom_b;

		J[iF] = NormalDerv_Dphi_Res - NormalDerv_Dphi_Ret;
	}
	
	return;
}
//-------------------------------------------------------------------------------
void ParameterEquations(parameters par, derivs_2D w, double *Xpar, double *F){

	return;	
}
//-------------------------------------------------------------------------------
void LinearParameterEquations(parameters par, derivs_2D w, derivs_2D Dw, double *Xpar, double *DXpar, double *J){

 return;					
}