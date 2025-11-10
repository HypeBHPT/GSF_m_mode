#include "Solve_ODE.h"

// -------------------------------------------------------------------
int newton(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1, *indx, iter=0, j;
	double *F, *DX, **J, d, norm, norm0;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	indx  = ivector(0, Ntotal);
	J     = dmatrix(0, Ntotal, 0, Ntotal);
		
	F_of_X(par, X, F);
	
	norm = norm0 = norm2(F, ntotal);
	norm0 = 1.;
	copy_dvector(DX, F, 0, Ntotal); 
	
	if (Newton_verb == 1){
		fprintf(par.fout," Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm/norm0);
		fprintf(par.fout," ------------------------------------------------------------------\n");
// 		pause();
	}
 	while(iter < Newton_itmin || (norm/norm0 > Newton_tol  && iter < Newton_itmax)){
		iter += 1;

// 		Jacobian(par, X, J);
		Jacobian_FD(par, X, J);
	 //       	PrintMatrix(J, 0, Ntotal, 0, Ntotal);
		// exit(1);
		ludcmp(J, Ntotal, indx, &d, 0);
		lubksb(J, Ntotal, indx, DX, 0);


		for (j = 0; j < ntotal; j++) 
			X[j] -= DX[j];

		F_of_X(par, X, F);
		norm = norm2(F, ntotal);
		copy_dvector(DX, F, 0, Ntotal); 
		if (Newton_verb == 1){
// 			fprintf(par.fout," Newton: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e \n",
			if (isinf(norm) || isnan(norm) || norm/norm0 > 1.0e+10) {
				fprintf(par.fout,"\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
// 				exit(1);
			}
			fprintf(par.fout," Newton: iter = %3d \t |F| = %e \n", iter, norm/norm0);
		}	
	}
	if(norm/norm0 > Newton_tol)
		fprintf(par.fout,
		" Newton Raphson Method failed to converge to prescribed tolerance (%3.3e of %3.3e). Now move on to the next sequence element.\n", norm/norm0,Newton_tol
		);
// 
	free_dvector(F,    0, Ntotal);
	free_dvector(DX,   0, Ntotal);
	free_ivector(indx, 0, Ntotal);
	free_dmatrix(J,    0, Ntotal, 0, Ntotal);

	return iter;
}
//-----------------------------------------------------------------------------

// -------------------------------------------------------------------------------
