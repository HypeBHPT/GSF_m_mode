#include "Solve_ODE.h"
#define  j1_max 100
//-------------------------------------------------------------------------------
void output_cheb(parameters par, double *X){
  
  FILE *fp;
  char fn[500];

  
  int iField,iDom,
      i, N1, n1;

	for(iDom=1; iDom<=nDom; iDom++){
	N1=par.N[iDom-1]; n1=N1+1;
  
    	for(iField=0; iField<nFields; iField++){
    		double *p=dvector(0, n1-1), *cp=dvector(0, n1-1);
    		sprintf(fn, "data/%s/cheb_Dom_%d_Field_%d.dat", par.SimName, iDom, iField);
    		fp=fopen(fn,"w");
    		for(i=0; i<=N1; i++){
    			int indx = Index(par, iDom,iField,i);
	  			p[i]=X[indx];
			}
			Chebyshev_Coefficients(p, cp, N1, par.grid[iDom-1]);
	
			fprintf(fp, "#1: i\t2:c_i\n");
			for(i=0; i<=N1; i++)
				fprintf(fp, "%d\t%3.15e\n", i, cp[i]);
			
			fprintf(fp, "\n");
			free_dvector(p,0,N1); free_dvector(cp, 0, N1);
			fclose(fp);
		}
	}

	return;
}
//-------------------------------------------------------------------------------
void output_Solution(parameters par, double *X){
  
  int iDom, iF, j1, N1, n1;
  FILE *fp;
  char fn[500];	
  
  double *Sol, *c;
  
  for(iDom=1; iDom<=nDom; iDom++){
  N1=par.N[iDom-1]; n1=N1+1;
  Sol=dvector(0, N1); c=dvector(0, N1);
  
  
  	for(iF=0; iF<nFields; iF++){
    
    	sprintf(fn, "data/%s/Solution_Dom_%d_fields_%d.dat",par.SimName, iDom, iF);
    	fp=fopen(fn ,"w");
    
    	for(j1=0; j1<n1; j1 ++) {
        	int indx = Index(par, iDom, iF, j1);
			Sol[j1]=X[indx];      
    	}
    	Chebyshev_Coefficients(Sol, c, N1, par.grid[iDom-1]);
    
    	fprintf(fp, "#1:sigma\t 2: sigma_bar\t 3:Conformal field\n");
    	for(j1=0; j1<=j1_max; j1++){
    		double chi = -1. + 2.*j1/j1_max, sigma_bar, x, sigma, dsigma=par.sigma[iDom]-par.sigma[iDom-1];
    		
    		x = get_x_from_chi(par.AnMR_x_boundary[iDom-1], par.AnMR_kappa[iDom-1], chi);
    		sigma_bar = par.sigma[iDom-1] + dsigma*(1.+ chi)*0.5;
    		sigma = par.sigma[iDom-1] + dsigma*(1.+ x)*0.5;

    		double out = Clenshaw_Chebyshev(c, N1, chi);
    		fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma,  sigma_bar, out);
    	}
    	fprintf(fp, "\n");
    	fclose(fp);
    }
    
    double par_sol;
  	for(iF=0; iF<nPar; iF++){
	
		sprintf(fn, "data/%s/Solution_Parameter_%d.dat", par.SimName, iF);
		fp=fopen(fn ,"w");
	
	    int indx = par.Ntotal-iF;
	    par_sol= X[indx];
	
	   	fprintf(fp, "%3.15e\n", par_sol);
	   	fclose(fp);
	}
	free_dvector(Sol,0, N1);
	free_dvector(c,0, N1);
   }

	
	return;
}//----------------------------------------------------------------------------------------------------------------------------------------
void output_Ret_Field_Boundary_data(parameters par, double *X){
	derivs W;
	allocate_derivs(&W, par.ntotal);
	copy_X_to_W(par, X, W); 
	Get_Derivatives(par, W);	


	double sigma_minus_output = par.sigma_minus, sigma_plus_output = par.sigma_plus;
	int idom_minus_output = 1, idom_plus_output = 2;

	char fn[500], dir_name[500];
	
	sprintf(dir_name, "../../InputData/r0_over_M_%.5lf/eta_%3.5lf/m_%d", par.r0_over_M, par.eta, par.m);
	create_directory(dir_name);

	sprintf(fn, "%s/phi_ret_boundary_ell_%d.dat", dir_name,   par.ell);

	
	FILE *fp=fopen(fn,"w");
	fprintf(fp, "eta\t%3.15e\tsigma_minus\t%3.15e\tsigma_plus\t%3.15e\n", par.eta, sigma_minus_output, sigma_plus_output);

	int idom, i, N;
	double *phi_real, *c_real, *phi_imag, *c_imag,
				 *dphidsig_real, *cd_real, *dphidsig_imag, *cd_imag,
				 *d2phidsig2_real, *cd2_real, *d2phidsig2_imag, *cd2_imag,
				 real_phi, imag_phi,
				 // real_dphidsig, imag_dphidsig,
				 // real_d2phid2sig, imag_d2phid2sig,
				 x, chi;

	complex_derivs f;
	complex_sigma_derivs phi;

	idom = idom_minus_output;
	x = get_x_from_sigma(par, idom, sigma_minus_output);
	chi = get_chi_from_x(par.AnMR_x_boundary[idom-1], par.AnMR_kappa[idom-1], x);

	N=par.N[idom-1];
	
	phi_real = dvector(0,N);c_real = dvector(0,N);
	phi_imag = dvector(0,N);c_imag = dvector(0,N);

	dphidsig_real = dvector(0,N);cd_real = dvector(0,N);
	dphidsig_imag = dvector(0,N);cd_imag = dvector(0,N);

	d2phidsig2_real = dvector(0,N);cd2_real = dvector(0,N);
	d2phidsig2_imag = dvector(0,N);cd2_imag = dvector(0,N);

	for(i=0; i<=N; i++){
		int indx_re = Index(par, idom, 0, i),
				indx_im = Index(par, idom, 1, i);

				f.d0  = W.d0[indx_re] + I*W.d0[indx_im];
				f.d1  = W.d1[indx_re] + I*W.d1[indx_im];
				f.d11 = W.d11[indx_re] + I*W.d11[indx_im];

				get_SigmaDerv_from_SpecDerv_complex(par, idom, i, f, &phi);

				phi_real[i] = creal(phi.d0);
				phi_imag[i] = cimag(phi.d0);

				dphidsig_real[i]  = creal(phi.dsigma);
				dphidsig_imag[i]  = cimag(phi.dsigma);

				d2phidsig2_real[i] = creal(phi.d2sigma);
				d2phidsig2_imag[i] = cimag(phi.d2sigma);
	}
	Chebyshev_Coefficients(phi_real, c_real, N, par.grid[idom-1]); Chebyshev_Coefficients(phi_imag, c_imag, N, par.grid[idom-1]);
	Chebyshev_Coefficients(dphidsig_real, cd_real, N, par.grid[idom-1]); Chebyshev_Coefficients(dphidsig_imag, cd_imag, N, par.grid[idom-1]);
	Chebyshev_Coefficients(d2phidsig2_real, cd2_real, N, par.grid[idom-1]); Chebyshev_Coefficients(d2phidsig2_imag, cd2_imag, N, par.grid[idom-1]);

	real_phi = Clenshaw_Chebyshev(c_real, N, chi);
	imag_phi = Clenshaw_Chebyshev(c_imag, N, chi);

	// real_dphidsig = Clenshaw_Chebyshev(cd_real, N, chi);
	// imag_dphidsig = Clenshaw_Chebyshev(cd_imag, N, chi);

	// real_d2phid2sig = Clenshaw_Chebyshev(cd2_real, N, chi);
	// imag_d2phid2sig = Clenshaw_Chebyshev(cd2_imag, N, chi);

	fprintf(fp, "real phi_minus\t%3.15e\t imag phi_minus\t%3.15e\n", real_phi, imag_phi);
	// fprintf(fp, "real dphidsig_minus\t%3.15e\t imag dphidsig_minus\t%3.15e\n", real_dphidsig, imag_dphidsig);
	// fprintf(fp, "real d2phidsig2_minus\t%3.15e\t imag d2phidsig2_minus\t%3.15e\n", real_d2phid2sig, imag_d2phid2sig);

	free_dvector(phi_real, 0, N); free_dvector(c_real, 0, N);
	free_dvector(phi_imag, 0, N); free_dvector(c_imag, 0, N);

	free_dvector(dphidsig_real, 0, N); free_dvector(cd_real, 0, N);
	free_dvector(dphidsig_imag, 0, N); free_dvector(cd_imag, 0, N);
	
	free_dvector(d2phidsig2_real, 0, N); free_dvector(cd2_real, 0, N);
	free_dvector(d2phidsig2_imag, 0, N); free_dvector(cd2_imag, 0, N);

//-------------------
	idom = idom_plus_output;
	x = get_x_from_sigma(par, idom, sigma_plus_output);
	chi = get_chi_from_x(par.AnMR_x_boundary[idom-1], par.AnMR_kappa[idom-1], x);

	N=par.N[idom-1];
	phi_real = dvector(0,N);c_real = dvector(0,N);
	phi_imag = dvector(0,N);c_imag = dvector(0,N);

	dphidsig_real = dvector(0,N);cd_real = dvector(0,N);
	dphidsig_imag = dvector(0,N);cd_imag = dvector(0,N);

	d2phidsig2_real = dvector(0,N);cd2_real = dvector(0,N);
	d2phidsig2_imag = dvector(0,N);cd2_imag = dvector(0,N);

	for(i=0; i<=N; i++){
		int indx_re = Index(par, idom, 0, i),
				indx_im = Index(par, idom, 1, i);

				f.d0  = W.d0[indx_re] + I*W.d0[indx_im];
				f.d1  = W.d1[indx_re] + I*W.d1[indx_im];
				f.d11 = W.d11[indx_re] + I*W.d11[indx_im];

				get_SigmaDerv_from_SpecDerv_complex(par, idom, i, f, &phi);

				phi_real[i] = creal(phi.d0);
				phi_imag[i] = cimag(phi.d0);

				dphidsig_real[i]  = creal(phi.dsigma);
				dphidsig_imag[i]  = cimag(phi.dsigma);

				d2phidsig2_real[i] = creal(phi.d2sigma);
				d2phidsig2_imag[i] = cimag(phi.d2sigma);
	}
	Chebyshev_Coefficients(phi_real, c_real, N, par.grid[idom-1]); Chebyshev_Coefficients(phi_imag, c_imag, N, par.grid[idom-1]);
	Chebyshev_Coefficients(dphidsig_real, cd_real, N, par.grid[idom-1]); Chebyshev_Coefficients(dphidsig_imag, cd_imag, N, par.grid[idom-1]);
	Chebyshev_Coefficients(d2phidsig2_real, cd2_real, N, par.grid[idom-1]); Chebyshev_Coefficients(d2phidsig2_imag, cd2_imag, N, par.grid[idom-1]);

	real_phi = Clenshaw_Chebyshev(c_real, N, chi);
	imag_phi = Clenshaw_Chebyshev(c_imag, N, chi);

	// real_dphidsig = Clenshaw_Chebyshev(cd_real, N, chi);
	// imag_dphidsig = Clenshaw_Chebyshev(cd_imag, N, chi);

	// real_d2phid2sig = Clenshaw_Chebyshev(cd2_real, N, chi);
	// imag_d2phid2sig = Clenshaw_Chebyshev(cd2_imag, N, chi);

	fprintf(fp, "real phi_plus\t%3.15e\t imag phi_plus\t%3.15e\n", real_phi, imag_phi);
	// fprintf(fp, "real dphidsig_plus\t%3.15e\t imag dphidsig_plus\t%3.15e\n", real_dphidsig, imag_dphidsig);
	// fprintf(fp, "real d2phidsig2_plus\t%3.15e\t imag d2phidsig2_plus\t%3.15e\n", real_d2phid2sig, imag_d2phid2sig);

	free_dvector(phi_real, 0, N); free_dvector(c_real, 0, N);
	free_dvector(phi_imag, 0, N); free_dvector(c_imag, 0, N);

	free_dvector(dphidsig_real, 0, N); free_dvector(cd_real, 0, N);
	free_dvector(dphidsig_imag, 0, N); free_dvector(cd_imag, 0, N);
	
	free_dvector(d2phidsig2_real, 0, N); free_dvector(cd2_real, 0, N);
	free_dvector(d2phidsig2_imag, 0, N); free_dvector(cd2_imag, 0, N);
	

	fclose(fp);

	free_derivs(&W, par.ntotal);
	return;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void output_phi_ell(parameters par, int ell_min, int ell_max, double *Sum_DomLeft, double *Sum_DomRight){
  
  	  char fn[500];
	  sprintf(fn, "data/%s/../../Regularization.dat", par.SimName);
	  
	  FILE *fp;
	  fp = fopen(fn, "w");
	  fprintf(fp, "#1: ell\t #2:Phi_ell_Dom1@r0\t #3#2:Phi_ell_Dom2@r0\t\n");	    

	  
	  int l;
	  for(l=ell_min; l<=ell_max; l++){
	  	    
	    fprintf(fp, "%d\t%3.15e\t%3.15e\n", l, Sum_DomLeft[l], Sum_DomRight[l]);
	   
	  }
	  fclose(fp);
  
  
  return;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void output_phi(parameters par, int ell_min, int ell_max, double *Sum_Dom){
  
  	  char fn[500];
	  sprintf(fn, "data/%s/../../SelfField.dat", par.SimName);
	  
	  FILE *fp;
	  fp = fopen(fn, "w");
	  fprintf(fp, "#1: ell\t #2:Phi_DomLeft@r0\t #3#2:Phi_ell_DomRight@r0\t\n");	    

	  double PartialSum=0., SumRef=0.;
	  
	  int l;
	  
	  for(l=ell_min; l<=ell_max; l++)
	    SumRef+=Sum_Dom[l];
  
	  
	  
	  for(l=ell_min; l<=ell_max; l++){
	    PartialSum+=Sum_Dom[l];
	    
	  	    
	    fprintf(fp, "%d\t%3.15e\t%3.15e\n", l, PartialSum, fabs(PartialSum-SumRef));
	   
	  }
	  fclose(fp);
  
  
  return;
}

////----------------------------------------------------------------------------------------------------------------------------------------
//void output_Convergence(parameters par, double *X){
//  static int count = 0;
//  int iF, j1, N=par.N;
//  double Sol_error;
//  FILE *fp,*fr;
//  char fn[200], fnr[200];
//  
//  if(nDom>1){printf("Error in output_Convergence: output not implemented for nDom>1\n"); exit(1);}
//  
//    
//    
//    
//    
//    double *Sol=dvector(0, N), *c=dvector(0, N);
//    for(iF=0; iF<nFields; iF++){
//	
//	sprintf(fnr, "data/%s/Solution_fields_%d.dat",par.SimName, iF);
//	sprintf(fn, "data/%s/Convergence_fields_%d.dat",par.SimName, iF);
//	
//	
//	fr=fopen(fnr ,"r");
//	if(count==0){
//	  fp=fopen(fn, "w");
//	  fprintf(fp, "#1:N\t2:Error (Nref = %d)\n", N);	  
//	}
//	else{
//	  fp=fopen(fn, "a");
//	}
//	  
//	  for(j1=0; j1<=N; j1 ++){
//	      int Ind = Index(par, iF,  j1);
//	      Sol[j1]=X[Ind];
//	    
//	  }
//	  Chebyshev_Coefficients(Sol, c, N, par.grid);
//	  
//	  
//	  Sol_error=0.;
//	  
//	  for(j1=0; j1<=j1_max; j1++){
//	      double s = -1. + 2.*j1/j1_max, dz=z1-z0, zeta=z0 + dz*(1.+s)*0.5, z = z_of_zeta(zeta);
//
//	      double Sol = Clenshaw_Chebyshev(c, N, s), Sol_ref, z_read, zeta_read;
//	      fscanf(fr, "%lf\t%lf\t%lf\n", &zeta_read, &z_read, &Sol_ref);
//	      
//	      if(fabs(zeta_read-zeta)>1.0e-16){
//		printf("\n Error in output_Convergence: read grid z=%lf, evaluating at z = %lf (dz=%3.15e)\n\n", z_read, z, fabs(z_read-z));
//		exit(1);
//	      }
//	      Sol_error=dmaximum2( fabsl(Sol - Sol_ref) , Sol_error) ;
//	      
//	    }
//	    fprintf(fp, "%d\t%3.15e\n", N, Sol_error);
//	fclose(fp);
//	fclose(fr);	    
//	}
//
//
//	free_dvector(Sol,0, N);
//	free_dvector(c,0, N);
//
//  
//  double par_sol, par_error;
//  for(iF=0; iF<nPar; iF++){
//	par_error=0.;
//	sprintf(fnr, "data/%s/Solution_Parameter_%d.dat",par.SimName, iF);
//	fr=fopen(fnr ,"r");
//	
//	sprintf(fn, "data/%s/Convergence_Parameter_%d.dat",par.SimName, iF);
//	if(count==0){
//	  fp=fopen(fn ,"w");
//	  fprintf(fp, "#1:N\t2:Error parameter (Nref = %d)\n", N);
//	}
//	else{
//	  fp=fopen(fn ,"a");
//	
//	    int Ind = par.Ntotal-iF;
//	    double par_ref;
//	    par_sol= X[Ind];
//	    fscanf(fr, "%lf\n", &par_ref);
//	    par_error = dmaximum2( fabsl(1. - par_sol/par_ref) , par_error) ;
//// 	    par_error = dmaximum2( fabsl(par_ref - par_sol) , par_error) ;
//
//	    
//	    fprintf(fp, "%d\t%3.15e\n", N, par_error);	    
//	}
//	fclose(fp);
//	fclose(fr);
//	
//	
//	   	
//  }
//  count ++;
//  return;
//}
//
//
//
