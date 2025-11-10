#include "Solve_ODE.h"


//--------------------------------------------------------------------------------
void F_of_X(parameters par, double *X, double *F)
{
  
	int i, iDom, iF, N, Ntotal = par.Ntotal, ntotal=Ntotal+1;
	double Faux[nFields];	

	derivs W;
	allocate_derivs(&W, ntotal);
	copy_X_to_W(par, X, W); 
	Get_Derivatives(par, W);

	for(iDom=1; iDom<=nDom; iDom++){
		N=par.N[iDom-1];
		for(i=0; i<=N; i++){
			     if(i==N && iDom==1)    Boundary_Data_Left(par, iDom, i, W, Faux);
			else if(i==N && iDom!=1)    Transition_Condition_Left(par, iDom, i, W, Faux);
			else if(i==0 && iDom!=nDom) Transition_Condition_Right(par, iDom, i, W, Faux);
			else if(i==0 && iDom==nDom) Boundary_Data_Right(par, iDom, i, W, Faux);
			else 											  Field_Eqs(par, iDom, i, W, Faux);
			 
			for(iF=0;iF<nFields; iF++){
			 	int idx = Index(par, iDom, iF, i);
			 	F[idx] = Faux[iF];
			}
		}

	}

	free_derivs(&W, ntotal);
	return;
}
// -------------------------------------------------------------------------------
void copy_X_to_W(parameters par, double *X, derivs W){
  int i, Ntotal=par.Ntotal;
  
//   #pragma omp parallel for
  for(i=0; i<=Ntotal; i++){
    W.d0[i]=X[i];
  }
}
// -------------------------------------------------------------------------------
void Jacobian_FD(parameters par, double *X, double **J)

{
	int j, Ntotal = par.Ntotal;
	double  eps = 5.e-06;
	
// 	#pragma omp parallel for
	for(j=0; j<= Ntotal; j++){
	  double *Xp, *Fp, *Xm, *Fm;
	  Xp = dvector(0, Ntotal);
	  Fp = dvector(0, Ntotal);
	  Xm = dvector(0, Ntotal);
	  Fm = dvector(0, Ntotal);
	
	
	  int l;  
	  for(l=0; l<= Ntotal; l++) Xp[l] = Xm[l] = X[l];

	
	  Xp[j] += eps; 
	  Xm[j] -= eps; 
		
	  F_of_X(par, Xp, Fp);
	  F_of_X(par, Xm, Fm);
	  for(l=0; l<= Ntotal; l++) J[l][j] = 0.5*(Fp[l]-Fm[l])/eps;
		
	  Xp[j] = Xm[j] = X[j];
	
	
	  free_dvector(Xp, 0, Ntotal);
	  free_dvector(Xm, 0, Ntotal);
	  free_dvector(Fm, 0, Ntotal);
	  free_dvector(Fp, 0, Ntotal);
	}
}
