#include "Solve_PDE.h"
// -------------------------------------------------------------------
int newton(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1, iter=0, j, output_bicgstab=0, FLAG=-1;
	double *F, *DX, norm, norm_prev;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);

	
	
	F_of_X(par, X, F);
	// PrintVector(par, F, 0, Ntotal);


	norm = norm2(F, ntotal);///sqrt(ntotal);	
	if (Newton_verb == 1){
		fprintf(par.fout, " Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm);fflush(0);
		fprintf(par.fout, " ------------------------------------------------------------------\n");fflush(0);
		// exit(1);	
	}
 	
	norm_prev = 1.e15;
	// while(iter < Newton_itmin || (/*FLAG < 0*/ norm > Newton_tol  && iter < Newton_itmax)){
	while(iter < Newton_itmin || (FLAG < 0 /*norm < Newton_tol*/  && iter < Newton_itmax)){
		// if(norm < Newton_tol) FLAG++;
		double normres=-1.0;
		if(norm < Newton_tol || norm/norm_prev>10*bicgstab_decr) FLAG++;
		norm_prev = norm;
		iter += 1;
		fill0_dvector(DX, 0, ntotal); 
		

		par.bicgstab_tol = bicgstab_decr*norm;

		output_bicgstab = bicgstab(par, X, DX, F, &normres);
// 		printf("in newton after bicgstab\n"); exit(1);
		
		for (j = 0; j < ntotal; j++) X[j] -= DX[j];

		F_of_X(par, X, F);
		norm = norm2(F, ntotal);///sqrt(ntotal);
		
		
		if (Newton_verb == 1){
			
			if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				fprintf(par.fout, "\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n"); fflush(0);
				exit(1);
			}		
			fprintf(par.fout, " Newton: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e (%3.2e), %d \n", iter, output_bicgstab, norm, Newton_tol, FLAG);
		}	
	}
	if(norm > Newton_tol){
		fprintf(par.fout,
		" Newton Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		); fflush(0);
	}

	// PrintVector(par, F, 0, Ntotal);

	free_dvector(F,    0, Ntotal);
	free_dvector(DX,   0, Ntotal);	
	return iter;
}
//-----------------------------------------------------------------------------
int bicgstab(parameters par, double *X, double *DX, double *F, double *normres)
{	// Biconjugate Stabilized Gradient Method, see pages 1, 2, http://www.netlib.org/templates/
	int Ntotal=par.Ntotal,ntotal=Ntotal+1, iter, n;
	double alpha = 0., beta = 0., rho = 0., rho1 = 1., rhotol = 1e-50, omega = 1.0e-50, omegatol = 1e-50,
		*p, *ph, *q, *r, *rt, *s, *sh, *t;

	p  = dvector(0, ntotal-1);		ph = dvector(0, ntotal-1);
	q  = dvector(0, ntotal-1);		r  = dvector(0, ntotal-1);
	s  = dvector(0, ntotal-1);		sh = dvector(0, ntotal-1);
	rt = dvector(0, ntotal-1);		t  = dvector(0, ntotal-1);
	
	
	JFD_Components JFD;	
	get_BandMatrix(par, X, &JFD);	
	// printf("m1=%d, m2=%d\n", JFD.m1, JFD.m2);
	// exit(-1);


	//-------- SOLVE BAND BRUTAL FORCE
	// double **J = dmatrix(0, Ntotal, 0, Ntotal), res, d;
	// int *indx  = ivector(0, Ntotal);
	// JacobianBand(par, X, J);
	// PrintMatrix(par, "Jacobian_Band.txt", J, 0, Ntotal, 0, Ntotal);

	// Jacobian(par, X, J);
	// PrintMatrix(par, "Jacobian.txt", J, 0, Ntotal, 0, Ntotal);
	

	// // ludcmp(J, Ntotal, indx, &d, 0);
	
	// printf("%d %d\n",JFD.m1, JFD.m2);
	// exit(-1);
	//---------



	double resprecond;

	// 	compute initial residual rt = r = p = F - J*DX
	J_times_DX(par, X, DX, r);	
	for (n = 0; n < ntotal; n++) rt[n] = r[n] = p[n] = F[n] - r[n];
	

	*normres = norm2(r, ntotal);	
	if (*normres <= par.bicgstab_tol) return 0;

	for (iter = 0; iter < bicgstab_itmax; iter++){
		rho = scalarproduct(rt, r, ntotal);
		if (fabs(rho) < rhotol)	  break;		
		if (iter > 0){ // compute direction vector p
			beta = (rho/rho1)*(alpha/omega);
			for (n = 0; n < ntotal; n++) 
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		
		// compute direction adjusting vector ph and scalar alpha :		
		PreCond(par, JFD,  X, p, ph, &resprecond); //J_FD*ph=p (solve for ph)		
		// PreCondLU(par, J, indx, X, p, ph, &res);
		
		J_times_DX(par, X, ph, q);	// q = J*ph
		alpha = rho/scalarproduct(rt, q, ntotal);
		
		for (n = 0; n < ntotal; n++) s[n] = r[n] - alpha * q[n];
			
		// early check of tolerance:
		*normres = norm2(s, ntotal);
	
		if (*normres <= par.bicgstab_tol){
			for (n = 0; n < ntotal; n++)
				DX[n] += alpha * ph[n];
			break;
		}
		// compute stabilizer vector sh and scalar omega:
		PreCond(par, JFD, X, s, sh, &resprecond);//J_FD*sh=s (solve for sh)
		// PreCondLU(par, J, indx, X, s, sh, &res);

		J_times_DX(par, X, sh, t);	// t = J*sh
		omega = scalarproduct(t, s, ntotal) / scalarproduct(t, t, ntotal);
		
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++) {
			if(n<ntotal) DX[n] += alpha * ph[n] + omega * sh[n];
			r[n]   = s[n] - omega * t[n];
		}		
		rho1 = rho;

		// are we done? 
		*normres = norm2(r, ntotal);
		if (bicgstab_verb == 1) {fprintf(par.fout, "BiCGStab: iter = %2d  norm = %6.1e\r", iter+1, *normres);fflush(0); }
		if (*normres <= par.bicgstab_tol)  break;
		
		if (fabs(omega) < omegatol)  break;
	}
	

	free_bandMatrix(par,&JFD);
	
	free_dvector(p,  0, ntotal-1); 		free_dvector(ph, 0, ntotal-1); 
	free_dvector(q,  0, ntotal-1); 		free_dvector(r,  0, ntotal-1);
	free_dvector(s,  0, ntotal-1); 		free_dvector(sh, 0, ntotal-1); 
	free_dvector(rt, 0, ntotal-1); 		free_dvector(t,  0, ntotal-1);
		
	/* iteration failed */
	if (iter > bicgstab_itmax) return  -1;
	
	/* breakdown */
	if (fabs(rho) < rhotol) return -10;
	if (fabs(omega) < omegatol) return -11;
	
	/* success! */
	if (bicgstab_verb == 1){fprintf(par.fout, "BiCGStab: iter = %2d  norm = %6.1e\n", iter+1, *normres);fflush(0); }

	// free_dmatrix(J,    0, Ntotal, 0, Ntotal);
	return iter;
}
//--------------------------------------------------------------------------------------------------------------------
void get_BandMatrix(parameters par, double *X, JFD_Components *JFD){

    int n=par.ntotal - nPar, N = n-1, NPar = nPar-1,       
        maxcol = nFields*STENCILSIZE;
    
//  allocating the components of JFD
    (*JFD).J_UL         = dmatrix(0, N, 0, maxcol-1);	fill0_dmatrix((*JFD).J_UL,        0, n,   0, maxcol);
    (*JFD).cols_J_UL    = imatrix(0, N, 0, maxcol-1);	fill0_imatrix((*JFD).cols_J_UL,   0, n,   0, maxcol);
    (*JFD).ncols_J_UL   = ivector(0, N);				fill0_ivector((*JFD).ncols_J_UL,  0, n);

    if(nPar>0){
      (*JFD).J_UR         = dmatrix(0, N, 0, NPar);		fill0_dmatrix((*JFD).J_UR,        0, n,   0, nPar);
      (*JFD).J_LL         = dmatrix(0, NPar, 0, N);		fill0_dmatrix((*JFD).J_LL,        0, nPar,0, n);
      (*JFD).J_LR         = dmatrix(0, NPar, 0, NPar);	fill0_dmatrix((*JFD).J_LR,        0, nPar,      0, nPar);
    }   

    Get_JFD_Matrix(par, X, JFD);

    (*JFD).K_UL         = dmatrix(1, n, 1, (*JFD).m1+(*JFD).m2+1); fill0_dmatrix((*JFD).K_UL,    1, n+1, 1, (*JFD).m1+(*JFD).m2+2);    
    (*JFD).Kl_UL        = dmatrix(1, n, 1, (*JFD).m1);			   fill0_dmatrix((*JFD).Kl_UL,   1, n+1, 1, (*JFD).m1+1);
    (*JFD).iK_UL        = ivector(1, n);						   fill0_ivector((*JFD).iK_UL,   1, n+1);
//     
    if(nPar>0){
      (*JFD).K_UR        = dmatrix(1, n,    1, nPar);			   fill0_dmatrix((*JFD).K_UR,    1, n+1,   1, nPar+1);
      (*JFD).K_LR        = dmatrix(1, nPar,    1, nPar);  		   fill0_dmatrix((*JFD).K_LR,    1, nPar+1,1, nPar+1);    
      (*JFD).indx_LR     = ivector(1, nPar); 				       fill0_ivector((*JFD).indx_LR, 1, nPar+1);
    }

    Get_JFD_Components(par,X, JFD);
   return;
}
//----------------------------------------------------------
void free_bandMatrix(parameters par,JFD_Components *JFD){

  int n=par.ntotal - nPar, N = n-1, NPar = nPar-1,       
      maxcol = nFields*STENCILSIZE;

  free_dmatrix((*JFD).J_UL, 0, N, 0, maxcol-1);
  if(nPar>0){
    free_dmatrix((*JFD).J_UR, 0, N, 0, NPar);
    free_dmatrix((*JFD).J_LL, 0, NPar, 0, N);
    free_dmatrix((*JFD).J_LR, 0, NPar, 0, NPar);
    
    free_dmatrix((*JFD).K_UR, 1, n,    1, nPar);
    free_dmatrix((*JFD).K_LR, 1, nPar,    1, nPar);
    free_ivector((*JFD).indx_LR, 1, nPar);
  }
  
  
  free_imatrix((*JFD).cols_J_UL, 0, N, 0, maxcol-1);
  free_ivector((*JFD).ncols_J_UL , 0, N);	  
  free_dmatrix((*JFD).K_UL , 1, n,    1, (*JFD).m1+(*JFD).m2+1);
  free_dmatrix((*JFD).Kl_UL, 1, n,    1, (*JFD).m1);
  free_ivector((*JFD).iK_UL, 1, n);

  return;
}
// -------------------------------------------------------------------------------
void Get_boundary_i1(parameters par, int jdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1, int *mdom)
{
	int d = FD_ORDER/2 + 1;
	adom[*mdom] = jdom;
	ai0[*mdom]  = par.N1[jdom]-d;
	ai1[*mdom]  = par.N1[jdom];
	aj0[*mdom]  = aj0[0];
	aj1[*mdom]  = aj1[0];
	*mdom += 1;
}
// -------------------------------------------------------------------------------
void Get_boundary_i2(parameters par, int jdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1, int *mdom)
{
	int d = FD_ORDER/2 + 1;
	adom[*mdom] = jdom;
	ai0[*mdom]  = ai0[0]; //par.N1[jdom]-d;
	ai1[*mdom]  = ai1[0]; //par.N1[jdom];
	aj0[*mdom]  = par.N2[jdom]-d;
	aj1[*mdom]  = par.N2[jdom];
	*mdom += 1;
}
// -------------------------------------------------------------------------------
void Get_gridpoints(parameters par, int n, int *mdom, int *adom, int *ai0, int *ai1, int *aj0, int *aj1)
{	// Getting the relevant gridpoints at which for a given unit vector DX the evaluation of the
	// field equations and boundary conditions yields a non-vanishing entry for J*DX
	// At the exit:
	//      mdom: number of domains within which relevant gridpoints are located
	//      for (jdom = 0; jdom < mdom; jdom ++) :
	//          adom[jdom]: contains the number of the domain (jdom), ranging from 0 to 1
	//           ai0[jdom]: contains the minimal x1-index of gridpoints in the domain adom[jdom]
	//           ai1[jdom]: contains the maximal x1-index of gridpoints in the domain adom[jdom]
	//           aj0[jdom]: contains the minimal x2-index of gridpoints in the domain adom[jdom]
	//           aj1[jdom]: contains the maximal x2-index of gridpoints in the domain adom[jdom]

	int j, N1, N2, d = FD_ORDER/2 + 1;
	
	// if(n < n_2D){
		int idom, iF, i;
		// 	n_v   = par.n_v_of_n[n], 
		// 	num_v = par.num_v; // = n_2D - par.nt[nDom-1]*NPOT; X[0]..X[num_v-1]: non-trivial v-values

		// Get_Indices_From_n_v(par, n_v, &idom, &iF, &i, &j);
		Get_Indices_From_Index(par, n, &idom, &iF, &i, &j);
		
		N1 = par.N1[idom];              
		N2 = par.N2[idom];

		*mdom   = 1;
		adom[0] = idom;

		// if (n < num_v){
			ai0[0] = maximum2(0,i-d);
			ai1[0] = minimum2(i+d,N1);
		// }
		// else{
			// ai0[0] = 0;
			// ai1[0] = N1;
		// }
		aj0[0] = maximum2(0,j-d);
		aj1[0] = minimum2(j+d,N2);

		// if (idom == 0 && ai1[0] == N1) Get_boundary_i1(par, 1, adom, ai0, ai1, aj0, aj1, mdom);
		// if (idom == 1 && ai1[0] == N1) Get_boundary_i1(par, 0, adom, ai0, ai1, aj0, aj1, mdom);

		if (idom == 0 && aj1[0] == N2){ 
			// Get_boundary_i2(par, 1, adom, ai0, ai1, aj0, aj1, mdom);
				int jdom=1;
				adom[*mdom] = jdom;
				ai0[*mdom]  = ai0[0]; //par.N1[jdom]-d;
				ai1[*mdom]  = ai1[0]; //par.N1[jdom];
				aj0[*mdom]  = 0;
				aj1[*mdom]  = d;
				*mdom += 1;
		}
		if (idom == 1 && aj0[0] == 0){ 
			// Get_boundary_i2(par, 0, adom, ai0, ai1, aj0, aj1, mdom);
				int jdom=0;
				adom[*mdom] = 0;
				ai0[*mdom]  = ai0[0]; //par.N1[jdom]-d;
				ai1[*mdom]  = ai1[0]; //par.N1[jdom];
				aj0[*mdom]  = par.N2[jdom]-d;
				aj1[*mdom]  = par.N2[jdom];
				*mdom += 1;
		}
	// }
	// else{
	// 	exit(-1);
	// 	N2    = par.N2[nDom-1];
	// 	// n2    = N2+1;
		
	// 	j     = n - n_2D;
	// 	*mdom = nDom;
	// 	for(jdom=0; jdom<nDom; jdom++){
	// 		adom[jdom] = jdom;
	// 		ai0[jdom]  = 0;
	// 		ai1[jdom]  = par.N1[jdom];
	// 		if(j>0 && j<N2){
	// 			int j = n - n_2D;
	// 			aj0[0] = maximum2(0,j-d);
	// 			aj1[0] = minimum2(j+d,N2);
	// 		}
	// 		else{ // taking ALL gridpoints for extra parameters [ and 1D functions(not implemented) ]
	// 			aj0[jdom]  = 0;
	// 			aj1[jdom]  = par.N2[jdom];
	// 		}
	// 	}
	// }
}
//----------------------------------------------------------------------------------------------------------
void Get_JFD_Matrix(parameters par, double *X, JFD_Components *JFD){ // Calculating JFD
	int iDom, N1, N2,
	    Ntotal = par.Ntotal, ntotal=par.ntotal, 
	    n=ntotal - nPar,
	    iDom_read, iF, i, j, row, column, mcol, ii, jj, J, i0, i1, j0, j1, jdom,
	    mdom, adom[nDom], ai0[nDom], aj0[nDom], ai1[nDom], aj1[nDom];
	
	
	double *JDX, *Eq=dvector(0, nFields-1), *EqPar=dvector(0,nPar),
	*ExtraParameters=dvector(0,nPar), *DExtraParameters=dvector(0,nPar);

	derivs_2D W, DW, W_atGrid, DW_atGrid;	
	allocate_derivs_2D(&W, ntotal);
	
	copy_X_to_W(par, X, W);
	Get_Derivatives(par, W);	
	get_ExtraPar_from_X(par, X, ExtraParameters);
	
	allocate_derivs_2D(&DW, ntotal);
	fill0_derivs_2D(DW, ntotal);	
	fill0_dvector(DExtraParameters,0,nPar);
	
	allocate_derivs_2D(&W_atGrid, nFields);
	allocate_derivs_2D(&DW_atGrid, nFields);


	JDX = dvector(0, Ntotal);
	fill0_dvector(JDX,  0, ntotal);

	(*JFD).m1 = (*JFD).m2 = 0;

	for(column = 0; column < n; column ++){	
		DW.d0[column]=1.;
		

		Get_Indices_From_Index(par, column, &iDom_read, &iF, &i, &j);
		// N1=par.N1[iDom_read]; 
		// N2=par.N2[iDom_read]; 
		// iDom = iDom_read;	

		// int d = FD_ORDER/2 + 1; 
	 //      i0 = maximum2(0,i-d);
	 //      i1 = minimum2(i+d,N1);
	 //      j0 = maximum2(0,j-d);
	      // j1 = minimum2(j+d,N2);		
		
		Get_gridpoints(par, column, &mdom, adom, ai0, ai1, aj0, aj1);
		// printf("%d (%d) --- mdom = %d, iDom = %d, iF = %d, i1 = %d, i2 = %d\n", column,n, mdom, iDom_read, iF, i , j);

		
		for(jdom=0; jdom < mdom; jdom++){
			iDom = adom[jdom];
			i0   = ai0[jdom];	i1 = ai1[jdom];
			j0   = aj0[jdom]; j1 = aj1[jdom]; 
			N1=par.N1[iDom]; 
			N2=par.N2[iDom]; 

			for(ii=i0; ii<=i1; ii++) 
				for(jj=j0; jj<=j1; jj++)
					Get_DerivativesFinDif_grid(par, ii, jj, iF, iDom, DW);
 		
	    for(jj=j0; jj<=j1; jj++){
	    	for(ii=i0; ii<=i1; ii++){
	    		if(ii*jj*( ii-N1 )*( jj-N2 )!=0 )
	    			LinearFieldEquations(par, W, DW, ExtraParameters, DExtraParameters, iDom, ii, jj, Eq);    		
	    		else
	    			LinearBoundaryCondition(par, W, DW, ExtraParameters, DExtraParameters, iDom, ii, jj, Eq);

	    		for(J=0;J<nFields; J++){
	    			row = Index(par, iDom, J, ii, jj);
	    			JDX[row]=Eq[J];
	  //   		}
	  //   	}
	  //   }
	  // }
	  // for(jdom=0; jdom < mdom; jdom++){
			// iDom = adom[jdom];
			// i0   = ai0[jdom];	i1 = ai1[jdom];
			// j0   = aj0[jdom]; j1 = aj1[jdom];
			// for(J=0;J<nFields; J++){
			// 	for(jj=j0; jj<=j1; jj++){
	  //   		for(ii=i0; ii<=i1; ii++){
	  //   			row = Index(par, iDom, J, ii, jj);

	    			if (fabs(JDX[row]) > TINY){
	    				mcol = (*JFD).ncols_J_UL[row];
	    				(*JFD).cols_J_UL[row][mcol] =  column;
	    				(*JFD).J_UL[row][mcol]      =  JDX[row];
	    				(*JFD).ncols_J_UL[row]     +=  1;

	    				if(column > row) // determining the numbers m1 and m2
	    					(*JFD).m2=maximum2((*JFD).m2, column-row);
	    				else
	    					(*JFD).m1=maximum2((*JFD).m1, row-column);
	    			}
	    			// JDX[row] =  0.;
	    		}
	    	}
    	}
    }

    	for(J=0; J<nPar; J++){
    		LinearParameterEquations(par, W, DW, ExtraParameters, DExtraParameters, EqPar);
    		(*JFD).J_LL[J][column]=EqPar[J];
    	}
    fill0_derivs_2D(DW, ntotal);
  }
	
	for(column = 0; column < nPar; column ++){
		DExtraParameters[column]=1.;
		
		for(iDom=0; iDom<nDom; iDom++){
			N1=par.N1[iDom];
			N2=par.N2[iDom];
			for(jj=0; jj<=N2; jj++){
				for(ii=0; ii<=N1; ii++){
					if(ii*jj*( ii-N1 )*( jj-N2 )!=0 )
						LinearFieldEquations(par, W, DW, ExtraParameters, DExtraParameters, iDom, ii, jj, Eq);
					else
						LinearBoundaryCondition(par, W, DW, ExtraParameters, DExtraParameters, iDom, ii, jj, Eq);

					for(J=0; J<nFields; J++){
						int Indx = Index(par, iDom, J, ii, jj);
						(*JFD).J_UR[Indx][column] = Eq[J];
					}
				}
			}
		}

		for(J=0; J<nPar; J++){
			LinearParameterEquations(par, W, DW, ExtraParameters, DExtraParameters, EqPar);
			(*JFD).J_LR[J][column] = EqPar[J];
		}
		DExtraParameters[column]=0.;
	}
	

	free_derivs_2D(&W, ntotal);
	free_derivs_2D(&DW, ntotal);
	free_derivs_2D(&W_atGrid, nFields);
	free_derivs_2D(&DW_atGrid, nFields);
	free_dvector(JDX, 0, ntotal-1);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(ExtraParameters, 0, nPar);
	free_dvector(DExtraParameters, 0, nPar);
	free_dvector(EqPar, 0, nPar);

	return;
}
// //-----------------------------------------------------------------------------
void Get_JFD_Components(parameters par, double *X, JFD_Components *JFD){
// Calculates the remaining components of JFD
 // Note: For the band-matrix operations we use the 'Numerical Recipes in C'
 // which provide the corresponding routines for matrices and vectors with indices
 // starting from 1.
	int n=par.ntotal - nPar, N = n-1,
	    m1=(*JFD).m1, m2=(*JFD).m2,
		j, k, row_K, col_K, col_J;
	double d=-1.;

	for(j=0; j<n; j++){
		row_K = j+1; // see Note
		for(k=0; k < (*JFD).ncols_J_UL[j]; k++){
			col_J = (*JFD).cols_J_UL[j][k];
			col_K = col_J+1; // see Note
			(*JFD).K_UL[row_K][m1+1+col_K-row_K] = (*JFD).J_UL[j][k];	// writing K
		}
	}

	// PrintMatrix(par, "K_UL_Band.txt", (*JFD).K_UL, 1, n, 1, (*JFD).m1+(*JFD).m2+1);
	
	// Calculating the LU-decomposition of K; introducing Kl and iK,
	// see 'Numerical Recipes in C' pages 51-54
	bandec((*JFD).K_UL, n, m1, m2, (*JFD).Kl_UL, (*JFD).iK_UL, &d);
	
	if(nPar>0){
	  
	  double *col_J_UR=dvector(0,N);
	  
	  int iF, jF;
	  for(iF=0; iF<nPar; iF++){
	    for(j=0; j<n; j++) col_J_UR[j] = (*JFD).J_UR[j][iF];
	    
	    banbks((*JFD).K_UL, n, (*JFD).m1, (*JFD).m2, (*JFD).Kl_UL, (*JFD).iK_UL, col_J_UR);
	    
	    for(j=0; j<n; j++) (*JFD).K_UR[j+1][iF+1] = col_J_UR[j];
	  }
	  
	  for(iF=0; iF<nPar; iF++){
	    for(jF=0; jF<nPar; jF++){
	      (*JFD).K_LR[iF+1][jF+1] = (*JFD).J_LR[iF][jF];
	      
	      for(j=0; j<n; j++)
		(*JFD).K_LR[iF+1][jF+1] -= (*JFD).J_LL[iF][j] * (*JFD).K_UR[j+1][jF+1];	      
	    }
	  }
	  
	  double d_LR;
	  
	  ludcmp((*JFD).K_LR, nPar, (*JFD).indx_LR, &d_LR, 1);		  
	  
	  free_dvector(col_J_UR,0,N);	  
	}
}
// -------------------------------------------------------------------
void PreCond(parameters par, JFD_Components JFD, double *X, double *b, double *Y, double *res){	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, n=par.ntotal - nPar, N = n-1, jF, j;
	double *c_U, *c_L;
	
	c_U=dvector(0, N);
	copy_dvector(c_U, b, 0, n);
	
	banbks(JFD.K_UL, n, JFD.m1, JFD.m2, JFD.Kl_UL, JFD.iK_UL, c_U);
	
	for(j=0; j<n; j++)
	  Y[j] = c_U[j];
	
	
	
	if(nPar>0){
	  
	  c_L=dvector(1, nPar);
	  
	  for(jF=0; jF<nPar; jF++){
	    c_L[jF+1]=b[Ntotal-jF];
	    
	    for(j=0; j<n; j++)
	      c_L[jF+1]-= JFD.J_LL[jF][j]*c_U[j];
	  }
	
	  lubksb(JFD.K_LR, nPar, JFD.indx_LR, c_L, 1);
	
	  for(jF=0; jF<nPar; jF++)	  
	    Y[Ntotal-jF] = c_L[jF+1];
	
	
	  
	  for(j=0; j<n; j++)
	    for(jF=0; jF<nPar; jF++) 
	      Y[j] -= JFD.K_UR[j+1][jF+1]*c_L[jF+1];	  
	
	  free_dvector(c_L, 1, nPar);
	 }
	
	free_dvector(c_U, 0, N);
	
	
	

}
// -------------------------------------------------------------------
void PreCondLU(parameters par, double **J, int *indx, double *X, double *b, double *Y, double *res)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1;
	
	copy_dvector(Y, b, 0, ntotal); 
	lubksb(J, Ntotal, indx, Y, 0);	
}
// -------------------------------------------------------------------
int newton_error(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int Ntotal = par.Ntotal, ntotal=Ntotal+1, iter=-1;
	double *F, norm;
	
	
	

	
	F     = dvector(0, Ntotal);
		
	F_of_X(par, X, F);
// 	PrintVector(par, F, 0, Ntotal);

	norm = norm2(F, ntotal);
	
	

		printf(" Newton Raphson Method: Initial Residual: \t |F| = %e \n", norm);fflush(0);
		printf(" ------------------------------------------------------------------\n");fflush(0);

	free_dvector(F,    0, Ntotal);
	return iter;
}
//-----------------------------------------------------------------------------
