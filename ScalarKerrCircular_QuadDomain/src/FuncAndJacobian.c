#include "Solve_PDE.h"

//-------------------------------------------------------------------------------
void copy_X_to_W(parameters par, double *X, derivs_2D W){
  int i, Ntotal=par.Ntotal;
  
//   #pragma omp parallel for
  for(i=0; i<=Ntotal; i++){
    W.d0[i]=X[i];
  }
}
// -------------------------------------------------------------------------------
void get_W_at_Grid(parameters par, int iDom, int j1, int j2, derivs_2D W, derivs_2D W_atGrid){
  int i;

  for(i=0; i<nFields; i++){
    int Indx = Index(par, iDom, i, j1, j2);
    W_atGrid.d0[i]=W.d0[Indx];
    W_atGrid.d1[i]=W.d1[Indx];
    W_atGrid.d11[i]=W.d11[Indx];
    W_atGrid.d2[i]=W.d2[Indx];
    W_atGrid.d22[i]=W.d22[Indx];
    W_atGrid.d12[i]=W.d12[Indx];
  }


}
//--------------------------------------------------------------------------------
void get_ExtraPar_from_X(parameters par, double *X, double *ExtraParameters){
	int i, Ntotal = par.Ntotal;
	ExtraParameters[nPar] = -1.; //Irrelevant extra componente
	for(i=0; i<nPar; i++){
		ExtraParameters[i] = X[Ntotal-i];
	}
	
	return;
}
// -------------------------------------------------------------------------------
void F_of_X(parameters par, double *X, double *F)
{
	int iDom, iFields, j1, j2, Ntotal = par.Ntotal, nTotal = Ntotal+1, N1, N2 ;
	double *Eq=dvector(0, nFields-1), *ExtraParameters=dvector(0,nPar), *EqPar=dvector(0,nPar);
	derivs_2D W, W_atGrid;
	allocate_derivs_2D(&W, nTotal);
	allocate_derivs_2D(&W_atGrid, nFields);
	
	copy_X_to_W(par, X, W);
	Get_Derivatives(par, W);
	get_ExtraPar_from_X(par, X, ExtraParameters);	
	
    
	for(iDom=0; iDom<nDom; iDom++){
		N1 = par.N1[iDom]; 
		N2 = par.N2[iDom];
		for(j2=0; j2<=N2; j2++){
			for(j1=0; j1<=N1; j1++){
				get_W_at_Grid(par, iDom, j1, j2, W, W_atGrid);
				if(j1*j2*( j1-N1 )*( j2-N2 )!=0 )
					FieldEquations(par, W, ExtraParameters, iDom, j1, j2, Eq);
				else
					BoundaryCondition(par, W, ExtraParameters, iDom, j1, j2, Eq);
	      	
	      for(iFields=0; iFields<nFields; iFields++){
	      	int Indx = Index(par, iDom, iFields, j1, j2);
	      	F[Indx]= Eq[iFields];
	      }
	    }
	  }
	}
	
	
	


	ParameterEquations(par, W, ExtraParameters, EqPar);
	for(iFields=0; iFields<nPar; iFields++){
	 	F[Ntotal-iFields] = EqPar[iFields];
	}
	
	  	  
	free_derivs_2D(&W, nTotal);
	free_derivs_2D(&W_atGrid, nFields);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(ExtraParameters, 0, nPar);
	free_dvector(EqPar, 0, nPar);
	
	return;
}
//-------------------------------------------------------------------------------
void J_times_DX(parameters par, double *X,double *DX, double *JDX)
{
	int iDom, iFields, j1, j2, Ntotal = par.Ntotal,nTotal = Ntotal+1, N1, N2;
	
	double *Eq=dvector(0, nFields-1), *EqPar=dvector(0,nPar),
	*ExtraParameters=dvector(0,nPar), *DExtraParameters=dvector(0,nPar);
	derivs_2D W, DW, W_atGrid, DW_atGrid;
	
	allocate_derivs_2D(&W, nTotal);
	allocate_derivs_2D(&DW, nTotal);
	
	copy_X_to_W(par, X, W);
	copy_X_to_W(par, DX, DW);
	get_ExtraPar_from_X(par, X, ExtraParameters);
	get_ExtraPar_from_X(par, DX, DExtraParameters);
	
	Get_Derivatives(par, W);	
	Get_Derivatives(par, DW);
	
	
	allocate_derivs_2D(&W_atGrid, nFields);
	allocate_derivs_2D(&DW_atGrid, nFields);

	for(iDom=0; iDom<nDom; iDom++){
		N1=par.N1[iDom];
		N2=par.N2[iDom];
		for(j2=0; j2<=N2; j2++){
			for(j1=0; j1<=N1; j1++){
				// get_W_at_Grid(par, iDom, j1, j2, W, W_atGrid);
				// get_W_at_Grid(par, iDom, j1, j2, DW, DW_atGrid);

	      
	      if(j1*j2*( j1-N1 )*( j2-N2 )!=0 )
	      	LinearFieldEquations(par, W, DW, ExtraParameters, DExtraParameters, iDom, j1, j2, Eq);
				else
					LinearBoundaryCondition(par, W, DW, ExtraParameters, DExtraParameters, iDom, j1, j2, Eq);
					
	      for(iFields=0; iFields<nFields; iFields++){
	      	int Indx = Index(par, iDom, iFields, j1, j2);
	      	JDX[Indx]=Eq[iFields];
	      }
	    }
	  }
	}
	  
	LinearParameterEquations(par, W, DW, ExtraParameters, DExtraParameters, EqPar);  
	for(iFields=0; iFields<nPar; iFields++){
		JDX[Ntotal-iFields] = EqPar[iFields];
	}

	free_derivs_2D(&W, nTotal);
	free_derivs_2D(&DW, nTotal);
	free_derivs_2D(&W_atGrid, nFields);
	free_derivs_2D(&DW_atGrid, nFields);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(ExtraParameters, 0, nPar);
	free_dvector(DExtraParameters, 0, nPar);
	free_dvector(EqPar, 0, nPar);	

	return;
}
//-------------------------------------------------------------------------------
void J_times_DX_FD(parameters par, double *X,double *DX, double *JDX)
{
	int iDom, iFields, j1, j2, Ntotal = par.Ntotal,nTotal = Ntotal+1, N1, N2;
	
	double *Eq=dvector(0, nFields-1), *EqPar=dvector(0,nPar),
	*ExtraParameters=dvector(0,nPar), *DExtraParameters=dvector(0,nPar);
	derivs_2D W, DW, W_atGrid, DW_atGrid;
	
	allocate_derivs_2D(&W, nTotal);
	allocate_derivs_2D(&DW, nTotal);
	
	copy_X_to_W(par, X, W);
	copy_X_to_W(par, DX, DW);
	get_ExtraPar_from_X(par, X, ExtraParameters);
	get_ExtraPar_from_X(par, DX, DExtraParameters);
	
	Get_Derivatives(par, W);	
	Get_DerivativesFinDif(par, DW);
	
	allocate_derivs_2D(&W_atGrid, nFields);
	allocate_derivs_2D(&DW_atGrid, nFields);

	for(iDom=0; iDom<nDom; iDom++){
		N1=par.N1[iDom];
		N2=par.N2[iDom];
		for(j2=0; j2<=N2; j2++){
			for(j1=0; j1<=N1; j1++){
				// get_W_at_Grid(par, iDom, j1, j2, W, W_atGrid);
				// get_W_at_Grid(par, iDom, j1, j2, DW, DW_atGrid);

	      
	      if(j1*j2*( j1-N1 )*( j2-N2 )!=0 )
	      	LinearFieldEquations(par, W, DW, ExtraParameters, DExtraParameters, iDom, j1, j2, Eq);
				else
					LinearBoundaryCondition(par, W, DW, ExtraParameters, DExtraParameters, iDom, j1, j2, Eq);
					
	      for(iFields=0; iFields<nFields; iFields++){
	      	int Indx = Index(par, iDom, iFields, j1, j2);
	      	JDX[Indx]=Eq[iFields];
	      }
	    }
	  }
	}
	  
	LinearParameterEquations(par, W, DW, ExtraParameters, DExtraParameters, EqPar);  
	for(iFields=0; iFields<nPar; iFields++){
		JDX[Ntotal-iFields] = EqPar[iFields];
	}

	free_derivs_2D(&W, nTotal);
	free_derivs_2D(&DW, nTotal);
	free_derivs_2D(&W_atGrid, nFields);
	free_derivs_2D(&DW_atGrid, nFields);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(ExtraParameters, 0, nPar);
	free_dvector(DExtraParameters, 0, nPar);
	free_dvector(EqPar, 0, nPar);	

	return;
}
// -------------------------------------------------------------------------------
