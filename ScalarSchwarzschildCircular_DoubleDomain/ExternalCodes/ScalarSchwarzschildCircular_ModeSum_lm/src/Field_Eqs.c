#include "Solve_ODE.h"
void get_ComplexField(parameters par, int indx_Re, int indx_Im, derivs W , complex_derivs *z){

	(*z).d0 = W.d0[indx_Re] +I*W.d0[indx_Im];
	(*z).d1 = W.d1[indx_Re] +I*W.d1[indx_Im];
	(*z).d11 = W.d11[indx_Re] +I*W.d11[indx_Im];

	return;
}
//---------------------------------------------------------------------------
void Field_Eqs(parameters par, int iDom, int i, derivs W, double *F){
	complex_derivs phi;
	double complex A_phi, Source_phi=0.;
	int indx_Re_phi = Index(par, iDom, 0, i),
		indx_Im_phi = Index(par, iDom, 1, i);

	get_ComplexField(par, indx_Re_phi, indx_Im_phi, W , &phi);

	A_phi = Diff_Operator(par, phi, iDom, i, 0);	
	Source_phi = Get_ModeSum_Source(par, iDom, i);
	

	F[0] = creal(A_phi - Source_phi);
	F[1] = cimag(A_phi - Source_phi);
	

	return;
}
//---------------------------------------------------------------------------
void Boundary_Data_Left(parameters par, int iDom, int i, derivs W, double *F){
	Field_Eqs(par, iDom, i, W, F);	
	return;
}
//---------------------------------------------------------------------------
void Boundary_Data_Right(parameters par, int iDom, int i, derivs W, double *F){
	Field_Eqs(par, iDom, i, W, F);
	return;
}
//---------------------------------------------------------------------------
void Transition_Condition_Left(parameters par, int iDom, int i, derivs W, double *F){
	int iF;

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, i),
			indx_PreviousDom = Index(par, iDom-1, iF, 0);
		F[iF] = W.d0[indx_Dom]-W.d0[indx_PreviousDom] - par.jump_field[iF][iDom-2];				    	
	}
	
	return;
}
//---------------------------------------------------------------------------
void Transition_Condition_Right(parameters par, int iDom, int i, derivs W, double *F){
	int iF;
	complex_derivs U, U_nextDom;
	complex_sigma_derivs f, f_nextDom; 

	for(iF=0; iF<nFields; iF++){
		int indx_Dom = Index(par, iDom, iF, i),
			indx_NextDom = Index(par, iDom+1, iF, par.N[iDom]);	
		
		get_ComplexField(par, indx_Dom, indx_Dom, W , &U);
		get_ComplexField(par, indx_NextDom, indx_NextDom, W , &U_nextDom);

		get_SigmaDerv_from_SpecDerv_complex(par, iDom, i, U, &f);
		get_SigmaDerv_from_SpecDerv_complex(par, iDom+1, par.N[iDom], U_nextDom, &f_nextDom);

		F[iF] = f_nextDom.dsigma - f.dsigma - par.jump_derivative[iF][iDom-1];			    	
	}
	
	return;
}