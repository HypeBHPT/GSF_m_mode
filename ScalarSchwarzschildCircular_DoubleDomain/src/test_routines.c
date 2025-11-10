#include "Solve_PDE.h"

func_sigma_derivs_2D Test_Func(parameters par, double sigma, double y){
  func_sigma_derivs_2D f;
  
  double sig, sig3, sig2, sig4, y2, y3, y4, den, den_2, den_3;
  
  sig   = sigma - par.sigma0;
  sig2  = sqr(sig);
  sig3  = sig2*sig;
  sig4  = sqr(sig2);
  y2    = sqr(y);
  y3    = y2*y;
  y4    = sqr(y2);
  den   = sig2 + y2;
  den_2 = sqr(den);
  den_3 = den_2*den;
  


  f.d0 = sig2*y2/den;

  f.dsigma = 2*y4*sig/den_2;
  f.d2sigma = 2*y4*(y2 - 3.*sig2)/den_3;

  f.dy = 2*y*sig4/den_2;
  f.d2y = 2*sig4*(sig2 - 3.*y2)/den_3;

  f.d2sigmay = 8*y3*sig3/den_3;


  return f;
}
// -------------------------------------------------------------------------------
double complex Test_Effective_Source(parameters par, int iDom, int i1, int i2){
  
  double complex HypSource, a2, a1, a0, g2, g1, da0, da1, da2_dsig, da2;
  double chi_1, chi_2;
  func_derivs_2D sigma, y;
  func_sigma_derivs_2D f;

  chi_1 = par.grid_chi_1[iDom][i1];
  chi_2 = par.grid_chi_2[iDom][i2];

  get_sigma(par, iDom, chi_1, chi_2, &sigma);
  get_y(par, iDom, chi_1, chi_2, &y);
  f = Test_Func(par, sigma.d0, y.d0);

  func_alpha0(par, sigma.d0, &a0, &da0);
  func_alpha1(par, sigma.d0, &a1, &da1);
  func_alpha2(par, sigma.d0, &a2, &da2_dsig, &da2);
  func_gamma2(par, y.d0, &g2);
  func_gamma1(par, y.d0, &g1);

  HypSource = a2*f.d2sigma + a1*f.dsigma + a0*f.d0 + g2*f.d2y+ g1*f.dy;

  return HypSource;
}
//--------------------------------
void output_test_function(parameters par){
	char fn[200], fn_chi[200], fn_sig[200], fn_y[200];
	FILE *fp, *fp_chi, *fp_sig, *fp_y;
	func_sigma_derivs_2D f;

	int iDom;
	for(iDom=0; iDom<nDom; iDom++){
    sprintf(fn, "data/%s/TestFunc_dom%d.dat", par.SimName, iDom);
    sprintf(fn_chi, "data/%s/TestFunc_chi_dom%d.dat", par.SimName, iDom);
    sprintf(fn_sig, "data/%s/Coord_sigma_dom%d.dat", par.SimName, iDom);
    sprintf(fn_y, "data/%s/Coord_y_dom%d.dat", par.SimName, iDom);


    fp = fopen(fn, "w");
    fp_chi = fopen(fn_chi, "w");
    fp_sig = fopen(fn_sig, "w");
    fp_y = fopen(fn_y, "w");
    

    int j1, j2, N1=par.N1[iDom], N2=par.N2[iDom];
	    for(j1=0; j1<=N1; j1++){
	      for(j2=0; j2<=N2; j2++){
	        double chi_1, chi_2;
	        func_derivs_2D sigma, y;

	        chi_1 = par.grid_chi_1[iDom][j1];
	        chi_2 = par.grid_chi_2[iDom][j2];

	        get_sigma( par, iDom, chi_1, chi_2, &sigma);
	        get_y( par, iDom, chi_1, chi_2, &y);

	        f = Test_Func(par, sigma.d0, y.d0);

	        fprintf(fp,"%3.15e %3.15e %3.15e\n", sigma.d0, y.d0, f.d0 );
	        fprintf(fp_chi,"%3.15e %3.15e %3.15e\n", chi_1, chi_2, f.d0 );

	        fprintf(fp_sig,"%3.15e %3.15e %3.15e\n", chi_1, chi_2, sigma.d0 );
	        fprintf(fp_y,"%3.15e %3.15e %3.15e\n", chi_1, chi_2,  y.d0);
	      }
	      fprintf(fp, "\n" );
	      fprintf(fp_chi, "\n" );
	      fprintf(fp_sig, "\n" );
	      fprintf(fp_y, "\n" );

	    }
  }
  fclose(fp);
  fclose(fp_chi);
  fclose(fp_sig);
  fclose(fp_y);

  return;
}
//--------------------------------
void output_test_Effective_Source(parameters par){
	char fn_real[200], fn_real_chi[200], fn_imag[200], fn_imag_chi[200];
	FILE *fp_real, *fp_real_chi, *fp_imag, *fp_imag_chi;
	double complex f;

	int iDom;
	for(iDom=0; iDom<nDom; iDom++){
    sprintf(fn_real, "data/%s/TestSeff_Real_dom%d.dat", par.SimName, iDom);
    sprintf(fn_real_chi, "data/%s/TestSeff_Real_chi_dom%d.dat", par.SimName, iDom);
    sprintf(fn_imag, "data/%s/TestSeff_Imag_dom%d.dat", par.SimName, iDom);
    sprintf(fn_imag_chi, "data/%s/TestSeff_Imag_chi_dom%d.dat", par.SimName, iDom);

    fp_real = fopen(fn_real, "w");
    fp_real_chi = fopen(fn_real_chi, "w");
    fp_imag = fopen(fn_imag, "w");
    fp_imag_chi = fopen(fn_imag_chi, "w");

    int j1, j2, N1=par.N1[iDom], N2=par.N2[iDom];
	    for(j1=0; j1<=N1; j1++){
	      for(j2=0; j2<=N2; j2++){
	        double chi_1, chi_2;
	        func_derivs_2D sigma, y;

	        chi_1 = par.grid_chi_1[iDom][j1];
	        chi_2 = par.grid_chi_2[iDom][j2];

	        get_sigma( par, iDom, chi_1, chi_2, &sigma);
	        get_y( par, iDom, chi_1, chi_2, &y);

	        f = Test_Effective_Source(par, iDom, j1, j2);

	        fprintf(fp_real,"%3.15e %3.15e %3.15e\n", sigma.d0, y.d0, creal(f) );
	        fprintf(fp_real_chi,"%3.15e %3.15e %3.15e\n", chi_1, chi_2, creal(f) );
	        fprintf(fp_imag,"%3.15e %3.15e %3.15e\n", sigma.d0, y.d0, cimag(f) );
	        fprintf(fp_imag_chi,"%3.15e %3.15e %3.15e\n", chi_1, chi_2, cimag(f) );
	      }
	      fprintf(fp_real, "\n" );
	      fprintf(fp_real_chi, "\n" );
	      fprintf(fp_imag, "\n" );
	      fprintf(fp_imag_chi, "\n" );

	    }
  }
  fclose(fp_real);
  fclose(fp_real_chi);
  fclose(fp_imag);
  fclose(fp_imag_chi);

  return;
}
//-------------------------------------------------------------------------------
void output_test_function_error(parameters par, double *X){
  
  int iDom, iF, j1, j2, j1_max=50, j2_max=50, N1, N2;
  FILE *fp;
  char fn[200];
  double **Sol, **c;
  func_sigma_derivs_2D f;
	
  for(iDom=0; iDom<nDom; iDom++){
    N1=par.N1[iDom];
    N2=par.N2[iDom];
    Sol=dmatrix(0, N1, 0, N2); c=dmatrix(0,N1, 0, N2);

    for(iF=0; iF<nFields; iF++){
    
      for(j1=0; j1<=N1; j1 ++) {
        for(j2=0; j2<=N2; j2 ++){
          int indx = Index(par, iDom, iF, j1, j2);
          Sol[j1][j2]=X[indx];
        }
      }
      Chebyshev_Coefficients_2D(Sol, c, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);


    
      sprintf(fn, "data/%s/ErrorSolution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
      fp=fopen(fn ,"w");
      fprintf(fp, "#1: sigma\t 2:y\t 3:Solution \n");

      for(j1=0; j1<=j1_max; j1++){
        double chi_1 = -1. + 2.*j1/j1_max ;
        func_derivs_2D sigma;

        for(j2=0; j2<=j2_max; j2++){
          double chi_2 = -1. + 2.*j2/j2_max;
          func_derivs_2D y;

          get_sigma(par, iDom, chi_1, chi_2, &sigma);
          get_y(par, iDom, chi_1, chi_2, &y);
          f = Test_Func(par, sigma.d0, y.d0);

          
          double out = Clenshaw_Chebyshev_2D(c, N1, N2, chi_1, chi_2), error = iF == 0? fabs(f.d0 - out):fabs(out);
          fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, error );
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
    free_dmatrix(Sol,0, N1, 0, N2);
    free_dmatrix(c, 0,N1, 0, N2);
  }

  for(iF=0; iF<nPar; iF++){
    sprintf(fn, "data/%s/Solution_parameter_%d.dat",par.SimName, iF);
    fp = fopen(fn, "w");
    double out = X[par.Ntotal-iF];
    fprintf(fp, "%3.15e\t\n", out );
    fclose(fp);
  }

  return;
}







//------------------------
void Func(parameters par, int iDom, double chi_1, double chi_2, 
	double *f, double *f_sig, double *f_sigsig, double *f_y, double *f_yy, double *f_sigy){
  
  double sig, sig2, y2, d_sig2, d2_sig2, d_y2, d2_y2;
  func_derivs_2D sigma, y;
  get_sigma( par, iDom, chi_1, chi_2, &sigma);
  get_y( par, iDom, chi_1, chi_2, &y);

  sig = sigma.d0 - par.sigma0;
  sig2=sqr(sig);
  d_sig2 = 2*sig;
  d2_sig2 = 2.;

  y2 = sqr(y.d0);
  d_y2 = 2*y.d0;
  d2_y2 = 2.;


 *f = sig2/(sig2 + y2);
 *f_sig = d_sig2/(sig2 + y2) - sig2*d_sig2/sqr(sig2 + y2);
 *f_sigsig = d2_sig2/(sig2 + y2) - 2*sqr(d_sig2)/sqr(sig2 + y2)  - sig2*d2_sig2/sqr(sig2 + y2) +2*sig2*sqr(d_sig2)/pow(sig2 + y2,3);

 *f_y = -sig2*d_y2/sqr(sig2 + y2);
 *f_yy = -sig2*d2_y2/sqr(sig2 + y2) + 2*sig2*sqr(d_y2)/pow(sig2 + y2,3);
 *f_sigy = -d_sig2*d_y2/sqr(sig2 + y2) + 2*sig2*d_sig2*d_y2/pow(sig2 + y2,3);

  return;
}
//------------------------
void Check_Derivative(parameters par, char *fn, int j1, int j2, derivs_2D W, int FLAG){
	
	
	int indx_0 = Index(par, 0, 0, j1, j2);
	int indx_1 = Index(par, 0, 1, j1, j2); 
	
	complex_derivs_2D w;
	get_ComplexField(par, indx_0, indx_1, W , &w);

	double chi_1 = par.grid_chi_1[0][j1];
	double chi_2 = par.grid_chi_2[0][j2];
	
	func_derivs_2D sigma, y;
	get_sigma(par, 0, chi_1, chi_2, &sigma);
	get_y(par, 0, chi_1, chi_2, &y);
	
	double f, f_sig, f_sigsig, f_y, f_yy, f_sigy;
	Func(par, 0, chi_1, chi_2, &f, &f_sig, &f_sigsig, &f_y, &f_yy, &f_sigy);
	
	complex_sigma_derivs phi;
	get_PhysDerv_from_SpecDerv_complex(par, 0, j1, j2, w, &phi);
	

	double err_fsig = fabs(f_sig - creal(phi.dsigma))+1.e-20;
	double err_fsigsig = fabs(f_sigsig - creal(phi.d2sigma))+1.e-20;
	double err_fy = fabs(f_y - creal(phi.dy))+1.e-20;
	double err_fyy = fabs(f_yy - creal(phi.d2y))+1.e-20;
	double err_fsigy = fabs(f_sigy - creal(phi.d2sigmay))+1.e-20;

	FILE *fp;
	fp = FLAG==0? fopen(fn,"w") : fopen(fn,"a");
	fprintf(fp, "%d\t %3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\n", 
		par.N1[0],err_fsig,err_fsigsig, err_fy, err_fyy, err_fsigy);
	fclose(fp);
}
