#include "Solve_PDE.h"
//-------------------------------------------------------------------------
void get_Solution(parameters par, double *X, derivs_2D Sol){
  int i, Ntotal=par.Ntotal;
  
//   #pragma omp parallel for
  for(i=0; i<=Ntotal; i++){
    Sol.d0[i]=X[i];
  }

  Get_Derivatives(par, Sol);
  return;
}
//-------------------------------------------------------------------------------
int solve_equations(parameters par, double *X){

	int newt_iter=-1;


	if(par.SOLVER_METHOD == -1){ //Debug Error in function F
		newt_iter = newton_error(par, X);
	}
	else if(par.SOLVER_METHOD == 0){ //Direct LU Inversion
		newt_iter = newton_SVN(par, X);
	}
	else if(par.SOLVER_METHOD == 1){ //Direct LU Inversion
		newt_iter = newton_direct(par, X);
	}
	else if(par.SOLVER_METHOD == 2){ //Iterative BiCGStab w/ FinDif Preconditioner
		newt_iter = newton(par, X);
	}
	
	else{
		fprintf(stderr, "Error in solve_equations: Method = %d not known\n", par.SOLVER_METHOD );
	}


	return newt_iter;
}