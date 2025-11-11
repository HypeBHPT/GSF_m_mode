#include "Solve_PDE.h"
//-------------------------------------------------------------------------------
void get_ComplexField(parameters par, int indx_Re, int indx_Im, derivs_2D W , complex_derivs_2D *z){

	(*z).d0 = W.d0[indx_Re] +I*W.d0[indx_Im];
	
	(*z).d1 = W.d1[indx_Re] +I*W.d1[indx_Im];
	(*z).d2 = W.d2[indx_Re] +I*W.d2[indx_Im];
	
	(*z).d11 = W.d11[indx_Re] +I*W.d11[indx_Im];
	(*z).d22 = W.d22[indx_Re] +I*W.d22[indx_Im];
	(*z).d12 = W.d12[indx_Re] +I*W.d12[indx_Im];

	(*z).d111 = W.d111[indx_Re] +I*W.d111[indx_Im];
	(*z).d112 = W.d112[indx_Re] +I*W.d112[indx_Im];
	(*z).d122 = W.d122[indx_Re] +I*W.d122[indx_Im];
	(*z).d222 = W.d222[indx_Re] +I*W.d222[indx_Im];

	(*z).d1111 = W.d1111[indx_Re] +I*W.d1111[indx_Im];
	(*z).d1112 = W.d1112[indx_Re] +I*W.d1112[indx_Im];
	(*z).d1122 = W.d1122[indx_Re] +I*W.d1122[indx_Im];
	(*z).d1222 = W.d1222[indx_Re] +I*W.d1222[indx_Im];
	(*z).d2222 = W.d2222[indx_Re] +I*W.d2222[indx_Im];

	(*z).d11222 = W.d11222[indx_Re] +I*W.d11222[indx_Im];

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
	Source_phi =HyperboloidalEffectiveSource(par, iDom, j1, j2);
	//Test_Effective_Source(par, iDom, j1, j2);
	//

	F[0] = creal( A_phi - Source_phi);
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
void BoundaryCondition(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	int indx_Re_phi = Index(par, iDom, 0, j1, j2),
			indx_Im_phi = Index(par, iDom, 1, j1, j2);

	double chi_1 = par.grid_chi_1[iDom][j1],
		   chi_2 = par.grid_chi_2[iDom][j2];

	complex_derivs_2D w;
	complex_sigma_derivs phi;
	func_derivs_2D sigma, y;

	get_ComplexField(par, indx_Re_phi, indx_Im_phi, W , &w);
	get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, w, &phi);
	get_sigma(par, iDom, chi_1, chi_2, &sigma);
	get_y(par, iDom, chi_1, chi_2, &y);

	

	switch(iDom){
		case 0:
			if(j1 == par.N1[iDom]){
				F[0] = creal(phi.d0 - par.phi_ret_minus[j2]);
				F[1] = cimag(phi.d0 - par.phi_ret_minus[j2]);
			}
			else if(j1 == 0){
				F[0] = creal(phi.d0 - par.phi_ret_plus[j2]);
				F[1] = cimag(phi.d0 - par.phi_ret_plus[j2]);
			}
			else if(j2 == par.N2[iDom]){
				FieldJumpEquations_Boundary2(par, W, Xpar, iDom, j1, j2, F);	
			}
			else if(j2 == 0){
				// if(par.m==0){ 
					FieldEquations(par, W, Xpar, iDom, j1, j2, F);
				// }
				// else{
				// 	F[0] = creal(phi.d0);
				// 	F[1] = cimag(phi.d0);
				// }
			}		
			break;

		case 1:
			if(j2 == 0){
				FieldDerivativeJumpEquations_Boundary2(par, W, Xpar, iDom, j1, j2, F);
			}
			else	FieldEquations(par, W, Xpar, iDom, j1, j2, F);
			break;

		default:
			fprintf(stderr,"Error in BoundaryCondition: iDom = %d does not exist\n", iDom);
			exit(-1);

	}

	return;
}
// -------------------------------------------------------------------------------	
void LinearBoundaryCondition(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	complex_derivs_2D Dw;
		complex_sigma_derivs Dphi;
		func_derivs_2D sigma, y;
		int indx_Re_phi = Index(par, iDom, 0, j1, j2),
		    indx_Im_phi = Index(par, iDom, 1, j1, j2);
		double chi_1 = par.grid_chi_1[iDom][j1],
					 chi_2 = par.grid_chi_2[iDom][j2];

		get_ComplexField(par, indx_Re_phi, indx_Im_phi, DW , &Dw);
		get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, Dw, &Dphi);
		get_sigma(par, iDom, chi_1, chi_2, &sigma);
		get_y(par, iDom, chi_1, chi_2, &y);


	switch(iDom){
		case 0:
			if(j1 == par.N1[iDom]){
				J[0] = creal(Dphi.d0);
				J[1] = cimag(Dphi.d0);
			}
			else if(j1 == 0){
				J[0] = creal(Dphi.d0);
				J[1] = cimag(Dphi.d0);
			}
			else if(j2 == par.N2[iDom]){
				LinearFieldJumpEquations_Boundary2(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
			}
			else if(j2 == 0){
				// if(par.m==0){ 
					LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
				// }
				// else{
				// 	J[0] = creal(Dphi.d0);
				// 	J[1] = cimag(Dphi.d0);
				// }
			}						
			
				
			break;

		case 1:
			if(j2 == 0){
				LinearFieldDerivativeJumpEquations_Boundary2(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
			}
			else{
				LinearFieldEquations(par, W, DW, Xpar, DXpar, iDom, j1, j2, J);
			}

			break;

			default:
				fprintf(stderr,"Error in BoundaryCondition: iDom = %d does not exist\n", iDom);
				exit(-1);

	}

	return;
	
}
//---------------------------------------------------------------------------
void FieldJumpEquations_Boundary1(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
			  indx_PreviousDom = Index(par, iDom-1, iF, 0,j2);
		F[iF] = W.d0[indx_Dom]-W.d0[indx_PreviousDom];
		// - par.jump_field[iF][iDom-1];				    	
	}
	
	// return;
}
//---------------------------------------------------------------------------
void LinearFieldJumpEquations_Boundary1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
			  indx_PreviousDom = Index(par, iDom-1, iF, 0,j2);
		J[iF] = DW.d0[indx_Dom] - DW.d0[indx_PreviousDom];
		// - par.jump_field[iF][iDom-1];				    	
	}
	
	// return;
}
//---------------------------------------------------------------------------
void FieldDerivativeJumpEquations_Boundary1(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	int iF;
	complex_derivs_2D U, U_nextDom;
	complex_sigma_derivs f, f_nextDom; 

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
				indx_NextDom = Index(par, iDom+1, iF, par.N1[iDom],j2);	

		get_ComplexField(par, indx_Dom, indx_Dom, W , &U);
		get_ComplexField(par, indx_NextDom, indx_NextDom, W , &U_nextDom);
				
		get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, U, &f);
		get_PhysDerv_from_SpecDerv_complex(par, iDom+1, par.N1[iDom+1], j2, U_nextDom, &f_nextDom);	

		F[iF] = f_nextDom.dsigma - f.dsigma; 
		//- par.jump_derivative[iF][iDom];			    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void LinearFieldDerivativeJumpEquations_Boundary1(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	int iF;
	complex_derivs_2D DU, DU_nextDom;
	complex_sigma_derivs Df, Df_nextDom; 

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
				indx_NextDom = Index(par, iDom+1, iF, par.N1[iDom],j2);	

		get_ComplexField(par, indx_Dom, indx_Dom, DW , &DU);
		get_ComplexField(par, indx_NextDom, indx_NextDom, DW , &DU_nextDom);
				
		get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, DU, &Df);
		get_PhysDerv_from_SpecDerv_complex(par, iDom+1, par.N1[iDom+1], j2, DU_nextDom, &Df_nextDom);	

		J[iF] = Df_nextDom.dsigma - Df.dsigma; 
		//- par.jump_derivative[iF][iDom];			    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void FieldJumpEquations_Boundary2(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	int iF;
	double chi_1_NextDom = par.grid_chi_1[iDom+1][j1], phiP, jump;

	
	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
				indx_NextDom = Index(par, iDom+1, iF, j1, 0);			

		phiP = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc, par.N1_LoadPunc, chi_1_NextDom):
									  Clenshaw_Chebyshev(par.Im_cheb_phi_Punc, par.N1_LoadPunc, chi_1_NextDom);

		jump = -phiP;


		F[iF] = W.d0[indx_NextDom]-W.d0[indx_Dom] - jump;			    	
	}
	
	
	return;
}
//---------------------------------------------------------------------------
void LinearFieldJumpEquations_Boundary2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
				indx_NextDom = Index(par, iDom+1, iF, j1, 0);	
		J[iF] = DW.d0[indx_NextDom] - DW.d0[indx_Dom];			    	
	}
	
	// return;
}
//---------------------------------------------------------------------------
void FieldDerivativeJumpEquations_Boundary2(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
	int iF;
	complex_derivs_2D U, U_PrevDom;
	func_derivs_2D sigma_Dom, y_Dom, sigma_PrevDom, y_PrevDom;
	complex_sigma_derivs phi_Res, phi_Ret;

	double chi_1_Dom = par.grid_chi_1[iDom][j1],
				 chi_2_Dom = par.grid_chi_2[iDom][j2],
				 
				 chi_1_PrevDom = par.grid_chi_1[iDom-1][j1],
				 chi_2_PrevDom = par.grid_chi_2[iDom-1][par.N2[iDom-1]],
				 
				 phiP_dsigma, phiP_dy, NormalDerv_phi_P, NormalDerv_phi_Res, NormalDerv_phi_Ret,
				 fac_Dom, fac_PrevDom, NormalDerv_jump, OneOver_g_sigsig, OneOver_g_yy;

	get_sigma(par, iDom, chi_1_Dom, chi_2_Dom, &sigma_Dom);
	get_y(par, iDom, chi_1_Dom, chi_2_Dom, &y_Dom);

	get_sigma(par, iDom-1, chi_1_PrevDom, chi_2_PrevDom, &sigma_PrevDom);
	get_y(par, iDom-1, chi_1_PrevDom, chi_2_PrevDom, &y_PrevDom);


	OneOver_g_sigsig = 1;//1./(4*(1+sigma_Dom.d0));
	OneOver_g_yy = 1;//4*y_Dom.d0*(1.-y_Dom.d0);

	fac_Dom = sqrt( OneOver_g_yy*sqr(sigma_Dom.d1) + OneOver_g_sigsig*sqr(y_Dom.d1)  );
	fac_PrevDom = sqrt( OneOver_g_yy*sqr(sigma_PrevDom.d1) + OneOver_g_sigsig*sqr(y_PrevDom.d1)  );

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
			  indx_PreviousDom = Index(par, iDom-1, iF, j1,par.N2[iDom-1]);

		phiP_dsigma = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_sigma, par.N1_LoadPunc, chi_1_Dom):
									  			 Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_sigma, par.N1_LoadPunc, chi_1_Dom);

		phiP_dy     = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_y, par.N1_LoadPunc, chi_1_Dom):
									  			 Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_y, par.N1_LoadPunc, chi_1_Dom);


		NormalDerv_phi_P = ( -OneOver_g_sigsig*y_Dom.d1 * phiP_dsigma + OneOver_g_yy*sigma_Dom.d1*phiP_dy )/fac_Dom;
		NormalDerv_jump = - NormalDerv_phi_P;

		get_ComplexField(par, indx_Dom, indx_Dom, W , &U);
		get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, U, &phi_Res);
		NormalDerv_phi_Res = ( -OneOver_g_sigsig * y_Dom.d1 * phi_Res.dsigma + OneOver_g_yy * sigma_Dom.d1 * phi_Res.dy )/fac_Dom;

		get_ComplexField(par, indx_PreviousDom, indx_PreviousDom, W , &U_PrevDom);
		get_PhysDerv_from_SpecDerv_complex(par, iDom-1, j1, par.N2[iDom-1], U_PrevDom, &phi_Ret);
		NormalDerv_phi_Ret = ( -OneOver_g_sigsig * y_PrevDom.d1 * phi_Ret.dsigma + OneOver_g_yy * sigma_PrevDom.d1*phi_Ret.dy)/fac_PrevDom;
				
		F[iF] = NormalDerv_phi_Res - NormalDerv_phi_Ret - NormalDerv_jump;

	}	
	return;


	// int iF;
	// complex_derivs_2D U, U_PrevDom;
	// func_derivs_2D sigma, y;

	// double chi_1_Dom = par.grid_chi_1[iDom][j1],
	// 			 chi_2_Dom = par.grid_chi_2[iDom][j2],
				 
	// 			 chi_1_PrevDom = par.grid_chi_1[iDom-1][j1],
	// 			 chi_2_PrevDom = par.grid_chi_2[iDom-1][par.N2[iDom-1]],
	// 			 phiP_dsigma, phiP_dy, phi_dx2, jump_dx2;

	// get_sigma(par, iDom, chi_1_Dom, chi_2_Dom, &sigma);
	// get_y(par, iDom, chi_1_Dom, chi_2_Dom, &y);

	// for(iF=0; iF<nFields; iF++){
	// 	int indx_Dom = Index(par, iDom, iF, j1,j2),
	// 		  indx_PreviousDom = Index(par, iDom-1, iF, j1,par.N2[iDom-1]);

	// 	// int indx_Dom = Index(par, iDom, iF, j1,j2),
	// 	// 		indx_NextDom = Index(par, iDom+1, iF, j1, 0);	

	// 	phiP_dsigma = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_sigma, par.N1_LoadPunc, chi_1_Dom):
	// 								  			 Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_sigma, par.N1_LoadPunc, chi_1_Dom);

	// 	phiP_dy     = iF%2==0? Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_y, par.N1_LoadPunc, chi_1_Dom):
	// 								  			 Clenshaw_Chebyshev(par.Im_cheb_phi_Punc_y, par.N1_LoadPunc, chi_1_Dom);

	// 	phi_dx2 = (sigma.d2 * phiP_dsigma + y.d2 * phiP_dy);
		
	// 	jump_dx2 = -phi_dx2;

	// 	get_ComplexField(par, indx_Dom, indx_Dom, W , &U);
	// 	get_ComplexField(par, indx_PreviousDom, indx_PreviousDom, W , &U_PrevDom);

		
	// 	MapFirstDerivative_Dom_a_to_Dom_b(par,  iDom-1, chi_1_PrevDom, chi_2_PrevDom, U_PrevDom, 
	// 																				  iDom, 			chi_1_Dom, 		 chi_2_Dom, &U_PrevDom);
				
	// 	F[iF] = U.d2 - U_PrevDom.d2 - jump_dx2;
		    	
	// }
	
	// return;
}
//---------------------------------------------------------------------------
void LinearFieldDerivativeJumpEquations_Boundary2(parameters par, derivs_2D W, derivs_2D DW, double *Xpar, double *DXpar, int iDom, int j1, int j2, double *J){
	int iF;
	complex_derivs_2D DU, DU_PrevDom;
	func_derivs_2D sigma_Dom, y_Dom, sigma_PrevDom, y_PrevDom;
	complex_sigma_derivs Dphi_Res, Dphi_Ret;

	double chi_1_Dom = par.grid_chi_1[iDom][j1],
				 chi_2_Dom = par.grid_chi_2[iDom][j2],
				 
				 chi_1_PrevDom = par.grid_chi_1[iDom-1][j1],
				 chi_2_PrevDom = par.grid_chi_2[iDom-1][par.N2[iDom-1]],
				 DNormalDerv_phi_Res, DNormalDerv_phi_Ret,
				 fac_Dom, fac_PrevDom, OneOver_g_sigsig, OneOver_g_yy;

	get_sigma(par, iDom, chi_1_Dom, chi_2_Dom, &sigma_Dom);
	get_y(par, iDom, chi_1_Dom, chi_2_Dom, &y_Dom);

	get_sigma(par, iDom-1, chi_1_PrevDom, chi_2_PrevDom, &sigma_PrevDom);
	get_y(par, iDom-1, chi_1_PrevDom, chi_2_PrevDom, &y_PrevDom);

	OneOver_g_sigsig = 1.;//1./(4*(1+sigma_Dom.d0));
	OneOver_g_yy = 1.;//4*y_Dom.d0*(1.-y_Dom.d0);

	fac_Dom = sqrt( OneOver_g_yy*sqr(sigma_Dom.d1) + OneOver_g_sigsig*sqr(y_Dom.d1)  );
	fac_PrevDom = sqrt( OneOver_g_yy*sqr(sigma_PrevDom.d1) + OneOver_g_sigsig*sqr(y_PrevDom.d1)  );

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, j1,j2),
			  indx_PreviousDom = Index(par, iDom-1, iF, j1,par.N2[iDom-1]);

		get_ComplexField(par, indx_Dom, indx_Dom, DW , &DU);
		get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, DU, &Dphi_Res);

		DNormalDerv_phi_Res = ( -OneOver_g_sigsig * y_Dom.d1 * Dphi_Res.dsigma + OneOver_g_yy * sigma_Dom.d1 * Dphi_Res.dy )/fac_Dom;

		get_ComplexField(par, indx_PreviousDom, indx_PreviousDom, DW , &DU_PrevDom);
		get_PhysDerv_from_SpecDerv_complex(par, iDom-1, j1, par.N2[iDom-1], DU_PrevDom, &Dphi_Ret);

		DNormalDerv_phi_Ret = ( -OneOver_g_sigsig * y_PrevDom.d1 * Dphi_Ret.dsigma + OneOver_g_yy * sigma_PrevDom.d1 * Dphi_Ret.dy)/fac_PrevDom;
						
		J[iF] = DNormalDerv_phi_Res - DNormalDerv_phi_Ret ;
		    	
	}
	
	return;

	// int iF;
	// complex_derivs_2D DU, DU_PrevDom;
	// // complex_sigma_derivs Df, Df_nextDom; 

	// double chi_1_Dom = par.grid_chi_1[iDom][j1],
	// 			 chi_2_Dom = par.grid_chi_2[iDom][j2],
	// 			 chi_1_PrevDom = par.grid_chi_1[iDom-1][j1],
	// 			 chi_2_PrevDom = par.grid_chi_2[iDom-1][par.N2[iDom-1]];

	// for(iF=0; iF<nFields; iF++){
	// 	int indx_Dom = Index(par, iDom, iF, j1,j2),
	// 		  indx_PreviousDom = Index(par, iDom-1, iF, j1,par.N2[iDom-1]);	

	// 	get_ComplexField(par, indx_Dom, indx_Dom, DW , &DU);
	// 	get_ComplexField(par, indx_PreviousDom, indx_PreviousDom, DW , &DU_PrevDom);

	// 	MapFirstDerivative_Dom_a_to_Dom_b(par,  iDom-1, chi_1_PrevDom, chi_2_PrevDom,  DU_PrevDom, 
	// 																				  iDom, 			chi_1_Dom, 		 chi_2_Dom, &DU_PrevDom);

		
				
	// 	// get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, DU, &Df);
	// 	// get_PhysDerv_from_SpecDerv_complex(par, iDom+1, j1, 0, DU_nextDom, &Df_nextDom);	

	// 	J[iF] = DU.d2 - DU_PrevDom.d2;
	// 	//Df_nextDom.dsigma - Df.dsigma; 
	// 	//- par.jump_derivative[iF][iDom];			    	
	// }
	
	// return;
}
//-------------------------------------------------------------------------------
void ParameterEquations(parameters par, derivs_2D w, double *Xpar, double *F){

	return;	
}
//-------------------------------------------------------------------------------
void LinearParameterEquations(parameters par, derivs_2D w, derivs_2D Dw, double *Xpar, double *DXpar, double *J){

 return;					
}
