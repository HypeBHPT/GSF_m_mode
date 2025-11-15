#include "Solve_PDE.h"

// -------------------------------------------------------------------------------
void Jacobian_FD(parameters par, double *X, double **J)
{
	int j, k, Ntotal = par.Ntotal;
	double *Xp, *Fp, *Xm, *Fm, eps = 5.e-6;
	
	Xp = dvector(0, Ntotal);
	Fp = dvector(0, Ntotal);
	Xm = dvector(0, Ntotal);
	Fm = dvector(0, Ntotal);
	
	for(j=0; j<= Ntotal; j++)
		Xp[j] = Xm[j] = X[j];

	for(j=0; j<= Ntotal; j++){
		Xp[j] += eps; 
		Xm[j] -= eps;
		
		F_of_X(par, Xp, Fp);
		F_of_X(par, Xm, Fm);
		
		for(k=0; k<= Ntotal; k++)			
		  J[k][j] = 0.5*(Fp[k]-Fm[k])/eps;
		
		
		Xp[j] = Xm[j] = X[j];
	}
	
	free_dvector(Xp, 0, Ntotal);
	free_dvector(Xm, 0, Ntotal);
	free_dvector(Fm, 0, Ntotal);
	free_dvector(Fp, 0, Ntotal);
}
//-----------------------------------------------------------------------------
void Jacobian(parameters par, double *X, double **J)
{
	int j, l, Ntotal = par.Ntotal, ntotal=Ntotal+1;
	double *DX, *JDX;

	DX  = dvector(0, Ntotal);
	JDX = dvector(0, Ntotal);
	
	fill0_dvector(DX,   0, ntotal);
	fill0_dvector(JDX,  0, ntotal);

	for(j=0; j <= Ntotal; j++){

		DX[j] = 1.;
		
		// J_times_DX_FD(par, X, DX, JDX);
		J_times_DX(par, X, DX, JDX);
		
		
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDX[l];
		
		DX[j] = 0.;
	}
	free_dvector(DX,  0, Ntotal);
	free_dvector(JDX, 0, Ntotal);	
	
}

//-----------------------------------------------------------------------------
void JacobianBand(parameters par, double *X, double **J)
{
	int j, l, Ntotal = par.Ntotal, ntotal=Ntotal+1;
	double *DX, *JDX;

	DX  = dvector(0, Ntotal);
	JDX = dvector(0, Ntotal);
	
	fill0_dvector(DX,   0, ntotal);
	fill0_dvector(JDX,  0, ntotal);

	for(j=0; j <= Ntotal; j++){

		DX[j] = 1.;
		
		J_times_DX_FD(par, X, DX, JDX);
		
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDX[l];
		DX[j] = 0.;
	}
	free_dvector(DX,  0, Ntotal);
	free_dvector(JDX, 0, Ntotal);	
	
}
// -------------------------------------------------------------------
int newton_direct(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1, *indx, iter=0, j, FLAG=-1;
	double *F, *DX, **J, d, norm;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	
	indx  = ivector(0, Ntotal);
	J     = dmatrix(0, Ntotal, 0, Ntotal);

	F_of_X(par, X, F);
	// PrintVector(par, F, 0, Ntotal);
	// exit(-1);

	norm = norm2(F, ntotal);
	
	copy_dvector(DX, F, 0, ntotal); 
	
	if (Newton_verb == 1){
		fprintf(par.fout, " Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm);
		fprintf(par.fout, " ------------------------------------------------------------------\n");
		
	}
	// exit(-1);
	Jacobian(par, X, J);
	// PrintMatrix(par, "Jacobian.txt", J, 0, Ntotal, 0, Ntotal);
	// exit(-1);

	// double **J_FD;
	// J_FD = dmatrix(0, Ntotal, 0, Ntotal);
	// Jacobian_FD(par, X, J_FD);
	// PrintMatrix(par, "Jacobian_FD.txt", J_FD, 0, Ntotal, 0, Ntotal);
	// exit(-1);
	ludcmp(J, Ntotal, indx, &d, 0);		
	
	
	while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
		if(norm < Newton_tol) FLAG++;
		iter += 1;
		// Jacobian_FD(par, X, J);			
		// PrintMatrix(par, "Jacobian_FD.txt", J, 0, Ntotal, 0, Ntotal);

		
		// PrintMatrix(par, "Jacobian.txt", J, 0, Ntotal, 0, Ntotal);
		// exit(-1);
		// JacobianBand(par, X, J);
		// PrintMatrix(par, "Jacobian_Band.txt", J, 0, Ntotal, 0, Ntotal);
		
		
		// int i;
		// Jacobian_FD(par, X, err_J);			
		// for(i=0; i<=Ntotal; i++){
		//   for(j=0; j<=Ntotal; j++){
		//     err_J[i][j]-=J[i][j];
		//   }
		// }
		// PrintMatrix(par, "Jacobian_Diff.txt", err_J, 0, Ntotal, 0, Ntotal);
		// exit(1);
		
		
		lubksb(J, Ntotal, indx, DX, 0);
		

		for (j = 0; j < ntotal; j++) 
			X[j] -= DX[j];

		F_of_X(par, X, F);

		// PrintVector(par, F, 0, Ntotal);
	    // exit(-1);

		norm = norm2(F, ntotal);
		copy_dvector(DX, F, 0, ntotal); 
		if (Newton_verb == 1){
// 			printf(" Newton: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e \n",
			if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				fprintf(par.fout, "\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
// 				exit(1);
			}
			fprintf(par.fout," Newton: iter = %3d \t |F| = %e (%3.2e), %d \n", iter, norm, Newton_tol, FLAG);
		}	
	}
	if(norm > Newton_tol)
		fprintf(par.fout,
		" Newton Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,    0, Ntotal);
	free_dvector(DX,   0, Ntotal);
	free_ivector(indx, 0, Ntotal);
	free_dmatrix(J,    0, Ntotal, 0, Ntotal);
	return iter;
}
// -------------------------------------------------------------------
int newton_SVN(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1, iter=0, i,j, FLAG=-1;
	double *F, **J, norm, **V, *W, **U, *b, *dx;
	
	F     = dvector(0, Ntotal);
	F_of_X(par, X, F);
	norm = norm2(F, ntotal);	
	
	J     = dmatrix(0, Ntotal, 0, Ntotal);
	
	//NOTE: SVN uses Matrix with index i = 1...n
	U     = dmatrix(1, ntotal, 1, ntotal);
	V     = dmatrix(1, ntotal, 1, ntotal);
	W     = dvector(1, ntotal);
	b     = dvector(1, ntotal);
	dx     = dvector(1, ntotal);

	Jacobian(par, X, J);
	
	
	for(i=0; i<=Ntotal; i++){
		b[i+1] = F[i];
		for(j=0; j<=Ntotal; j++){
			U[i+1][j+1]=J[i][j];
		}
	}
	
	svdcmp(U,ntotal,ntotal,W,V);
	

	double wmax=0, wmin=1.e20, w_small, w_threshold=1.e-8;
	for(i=1; i<=ntotal; i++){		
		if(fabs(W[i]) > fabs(wmax)) wmax = W[i];
		if(fabs(W[i]) < fabs(wmin)) wmin = W[i];
	}
	

	w_small = wmin/wmax;
	printf("%3.15e\n", w_small);
	for(i=1; i<=ntotal; i++){
		if(fabs(W[i]/wmax) < w_threshold) W[i]=0.;
	}
	// printf("%3.15e\n", w_threshold);

	for(i=1; i<=ntotal; i++){
		printf("%d, %3.15e\n", i , W[i]);
	}

	// exit(-1);

	
	
	
	
	if (Newton_verb == 1){
		printf(" Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm);
		printf(" ------------------------------------------------------------------\n");
		
	}


	
	
	while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
		if(norm < Newton_tol) FLAG++;
		iter += 1;

		svbksb(U,W,V,ntotal,ntotal, b, dx);
		// lubksb(J, Ntotal, indx, DX, 0);

		for (j = 0; j < ntotal; j++) X[j] -= dx[j+1];

		F_of_X(par, X, F);
		norm = norm2(F, ntotal);
		for(i=0; i<=Ntotal; i++) b[i+1] = F[i];
		
		if (Newton_verb == 1){
// 			printf(" Newton: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e \n",
			if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				printf("\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
// 				exit(1);
			}
			printf(" Newton: iter = %3d \t |F| = %e (%3.2e), %d \n", iter, norm, Newton_tol, FLAG);
		}	
	}
	if(norm > Newton_tol)
		printf(
		" Newton Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,    0, Ntotal);	
	free_dmatrix(J,    0, Ntotal, 0, Ntotal);
	free_dmatrix(U,    1, ntotal, 1, ntotal);
	free_dmatrix(V,    1, ntotal, 1, ntotal);
	free_dvector(W,    1, ntotal);
	free_dvector(b,    1, ntotal);
	free_dvector(dx,   1, ntotal);
	
	return iter;
}
//-----------------------------------------------------------------------------

