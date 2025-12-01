#include "Solve_PDE.h"

double get_clm(int l, int m){
  int x=0;
  double c_lm = 1.;

  while(x<= 2*m-1){

    c_lm/=sqrt(l + m - x);
    x++;
  }
  c_lm *=sqrt( (2.*l +1.)/(4*Pi));

  return c_lm;
}
//----------------------------------------------------------------------------------------
double complex Get_phi_From_hyp_phi(parameters par, double sigma, double y, int FLAG_NS,  double complex hyp_phi){
  
  double complex phi, Z, dlnZdsigma, dlnZ;  
  func_Z(par, sigma, y, FLAG_NS,  &Z, &dlnZdsigma, &dlnZ);
  
  phi = Z *  hyp_phi;  
  
  return phi;
}
//----------------------------------------------------------------------------------------
complex_derivs_2D get_complex_dervs_conjugate(complex_derivs_2D z){
  complex_derivs_2D conj_z;

  conj_z.d0  = conj(z.d0);
  conj_z.d1  = conj(z.d1);
  conj_z.d11 = conj(z.d11);
  conj_z.d12 = conj(z.d12);
  conj_z.d2  = conj(z.d2);
  conj_z.d22 = conj(z.d22);

  return conj_z;
}
//----------------------------------------------------------------------------------------
void Get_SelfForce_Ft_m(parameters par, double *X, double complex *Ft_m_bulk, double complex *Ft_m_equator){
  //NOTE: EXPRESSION FOR A SINGLE m>0
  //NOTE: OUTPUT WILL TAKE INTO ACCOUNT THE CONTRIBUTIONS FROM m AND -m
  int idom_part = par.Dom_ptcl, N1, N2,  j1, j2;
  double **Re_Sol, **Re_c, **Im_Sol, **Im_c, q =par.q, lambda=par.lambda_over_M;
  double complex phi_hyp_bulk, phi_hyp_equator, phi_bulk,  s=par.s;


  N1=par.N1[idom_part];
  N2=par.N2[idom_part];
  
  Re_Sol=dmatrix(0, N1, 0, N2); 
  Re_c=dmatrix(0,N1, 0, N2);


  Im_Sol=dmatrix(0, N1, 0, N2); 
  Im_c=dmatrix(0,N1, 0, N2);

  
  for(j1=0; j1<=N1; j1 ++) {
    for(j2=0; j2<=N2; j2 ++){
      int indx_Real = Index(par, idom_part, 0, j1, j2),
          indx_Imag = Index(par, idom_part, 1, j1, j2);
      
      Re_Sol[j1][j2]=X[indx_Real];
      
      Im_Sol[j1][j2]=X[indx_Imag];
    }
  }

  Chebyshev_Coefficients_2D(Re_Sol, Re_c, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);
  Chebyshev_Coefficients_2D(Im_Sol, Im_c, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);
  
  phi_hyp_bulk = Clenshaw_Chebyshev_2D(Re_c, N1, N2, 0., -1.) + I*Clenshaw_Chebyshev_2D(Im_c, N1, N2, 0., -1.);
  phi_hyp_equator = Clenshaw_Chebyshev_2D(Re_c, N1, N2, 1., -1.) + I*Clenshaw_Chebyshev_2D(Im_c, N1, N2, 1., -1.);

  phi_bulk = Get_phi_From_hyp_phi(par, par.sigma0, 0., 1, phi_hyp_bulk);
  phi_hyp_equator = Get_phi_From_hyp_phi(par, par.sigma0, 0., 1, phi_hyp_equator);

  *Ft_m_bulk = q/lambda * s * phi_bulk ;
  *Ft_m_equator = q/lambda * s * phi_hyp_equator;
  
  free_dmatrix(Re_Sol,0, N1, 0, N2); free_dmatrix(Re_c, 0,N1, 0, N2);
  free_dmatrix(Im_Sol,0, N1, 0, N2); free_dmatrix(Im_c, 0,N1, 0, N2);

  return;

}
//----------------------------------------------------------------------------------------
void Get_SelfForce_Fr_m(parameters par, double *X, double complex *Fr_m_bulk, double complex *Fr_m_equator){
  //NOTE: EXPRESSION FOR A SINGLE m>0
  //NOTE: OUTPUT WILL TAKE INTO ACCOUNT THE CONTRIBUTIONS FROM m AND -m
  int idom_part = par.Dom_ptcl, N1, N2,  j1, j2;
  double **Re_Sol, **Re_c, **Im_Sol, **Im_c, 
         **Re_Sol_sigma, **Re_c_sigma, **Im_Sol_sigma, **Im_c_sigma, H, dH_dsigma,
         q =par.q, sigma0 =par.sigma0, sigma0_2 = sqr(sigma0), rh = par.rh_over_M;
  double complex s=par.s, Z, dlnZ_dsig, dlnZ_dy, phi_hyp_bulk, phi_hyp_equator, phi_hyp_sigma_bulk, phi_hyp_sigma_equator;  
  func_Z(par, sigma0, 0., 1, &Z, &dlnZ_dsig, &dlnZ_dy);
  func_Height(par, sigma0, &H, &dH_dsigma);
  
  derivs_2D W;
	allocate_derivs_2D(&W, par.ntotal);
  copy_X_to_W(par, X, W);	
	Get_Derivatives(par, W);

  complex_derivs_2D Sol_x1x2;
  complex_sigma_derivs Sol_sigma_y;

  N1=par.N1[idom_part];
  N2=par.N2[idom_part];
  
  Re_Sol=dmatrix(0, N1, 0, N2);  Re_c=dmatrix(0,N1, 0, N2);
  Re_Sol_sigma=dmatrix(0, N1, 0, N2);  Re_c_sigma=dmatrix(0,N1, 0, N2);

  Im_Sol=dmatrix(0, N1, 0, N2);  Im_c=dmatrix(0,N1, 0, N2);
  Im_Sol_sigma=dmatrix(0, N1, 0, N2);  Im_c_sigma=dmatrix(0,N1, 0, N2);
  
  for(j1=0; j1<=N1; j1 ++) {
    for(j2=0; j2<=N2; j2 ++){
      int indx_Real = Index(par, idom_part, 0, j1, j2),
          indx_Imag = Index(par, idom_part, 1, j1, j2);
      
      Sol_x1x2.d0 =  W.d0[indx_Real] + I* W.d0[indx_Imag];
      Sol_x1x2.d1 =  W.d1[indx_Real] + I* W.d1[indx_Imag];
      Sol_x1x2.d11 = W.d11[indx_Real] + I* W.d11[indx_Imag];
      Sol_x1x2.d2 = W.d2[indx_Real] + I* W.d2[indx_Imag];
      Sol_x1x2.d22 = W.d22[indx_Real] + I* W.d22[indx_Imag];
      Sol_x1x2.d12 = W.d12[indx_Real] + I* W.d12[indx_Imag];

      get_PhysDerv_from_SpecDerv_complex(par, idom_part, j1, j2, Sol_x1x2, &Sol_sigma_y);
      
      Re_Sol[j1][j2]=creal(Sol_sigma_y.d0); 
      Im_Sol[j1][j2]=cimag(Sol_sigma_y.d0);

      Re_Sol_sigma[j1][j2]=creal(Sol_sigma_y.dsigma); 
      Im_Sol_sigma[j1][j2]=cimag(Sol_sigma_y.dsigma);
    }
  }

  Chebyshev_Coefficients_2D(Re_Sol, Re_c, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);
  Chebyshev_Coefficients_2D(Im_Sol, Im_c, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);

  Chebyshev_Coefficients_2D(Re_Sol_sigma, Re_c_sigma, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);
  Chebyshev_Coefficients_2D(Im_Sol_sigma, Im_c_sigma, N1, par.grid_1[idom_part], N2, par.grid_2[idom_part]);

  phi_hyp_bulk = Clenshaw_Chebyshev_2D(Re_c, N1, N2, 0., -1.) + I*Clenshaw_Chebyshev_2D(Im_c, N1, N2, 0., -1.);
  phi_hyp_equator = Clenshaw_Chebyshev_2D(Re_c, N1, N2, 1., -1.) + I*Clenshaw_Chebyshev_2D(Im_c, N1, N2, 1., -1.);

  phi_hyp_sigma_bulk = Clenshaw_Chebyshev_2D(Re_c_sigma, N1, N2, 0., -1.) + I*Clenshaw_Chebyshev_2D(Im_c_sigma, N1, N2, 0., -1.);
  phi_hyp_sigma_equator = Clenshaw_Chebyshev_2D(Re_c_sigma, N1, N2, 1., -1.) + I*Clenshaw_Chebyshev_2D(Im_c_sigma, N1, N2, 1., -1.);

  *Fr_m_bulk    = -q * sigma0_2/rh * Z * ( s * dH_dsigma * phi_hyp_bulk    + phi_hyp_sigma_bulk    +  phi_hyp_bulk/sigma0  );
  *Fr_m_equator = -q * sigma0_2/rh * Z * ( s * dH_dsigma * phi_hyp_equator + phi_hyp_sigma_equator +  phi_hyp_equator/sigma0  );
 
  free_dmatrix(Re_Sol,0, N1, 0, N2); free_dmatrix(Re_c, 0,N1, 0, N2);
  free_dmatrix(Re_Sol_sigma,0, N1, 0, N2); free_dmatrix(Re_c_sigma, 0,N1, 0, N2);
  free_dmatrix(Im_Sol,0, N1, 0, N2); free_dmatrix(Im_c, 0,N1, 0, N2);
  free_dmatrix(Im_Sol_sigma,0, N1, 0, N2); free_dmatrix(Im_c_sigma, 0,N1, 0, N2);
  free_derivs_2D(&W, par.ntotal);

  return;

}
//----------------------------------------------------------------------------------------
void get_Puncture_EffectiveSource(parameters *par){
 int idom_part = (*par).Dom_ptcl, i1, i2, N1 = (*par).N1_PuncSeff, N2 = (*par).N2_PuncSeff, m = (*par).m;
 
 
 double chi_1, chi_2, r, dr_dsigma, theta,
 Punc_BL_Bound[2], dPunc_BL_Bound_dx[8], d2Punc_BL_Bound_dx2[20], Seff_BL[2],
 *Re_Hyp_Punc_Bound = dvector(0, N1), *Im_Hyp_Punc_Bound = dvector(0, N1),
 *Re_dHyp_Punc_dsigma_Bound = dvector(0, N1), *Im_dHyp_Punc_dsigma_Bound = dvector(0, N1),
 *Re_dHyp_Punc_dy_Bound = dvector(0, N1), *Im_dHyp_Punc_dy_Bound = dvector(0, N1),
 **Re_Hyp_Seff = dmatrix(0, N1, 0, N2), **Im_Hyp_Seff = dmatrix(0, N1, 0, N2);
 func_derivs_2D sigma, y;

 double complex Hyp_phi_Punc_Bound_complex, dHyp_phi_Punc_dsig_Bound_complex, dHyp_phi_Punc_dy_Bound_complex,
                Seff_BL_complex,  Hyp_Seff_complex;

     
 double **Re_Hyp_Punc = dmatrix(0, N1, 0, N2), **Im_Hyp_Punc = dmatrix(0, N1, 0, N2);
 
 FILE *fp = fopen("Hyp_Seff_m.dat", "w");
 
 for(i1=0; i1<=N1; i1++){
   chi_1 = get_grid( (*par).grid_1_PuncSeff, i1, N1 );
   for(i2=0; i2<=N2; i2++){
    chi_2 = get_grid( (*par).grid_2_PuncSeff, i2, N2 );

    get_sigma(*par, idom_part, chi_1, chi_2, &sigma);    
    func_r_of_sigma(*par, sigma.d0, &r, &dr_dsigma); 

    get_y(*par, idom_part, chi_1, chi_2, &y);
    theta = acos(sqrt(y.d0)) ;
    struct coordinate xBL_coord = {r, theta, 0., 0.};

    effsource_calc_m(m, &xBL_coord, Punc_BL_Bound, dPunc_BL_Bound_dx, d2Punc_BL_Bound_dx2, Seff_BL); 

    chi_2 = get_grid( (*par).grid_2_PuncSeff, i2, N2 );
     get_sigma(*par, idom_part, chi_1, chi_2, &sigma);    
    func_r_of_sigma(*par, sigma.d0, &r, &dr_dsigma); 

    get_y(*par, idom_part, chi_1, chi_2, &y);
    theta = acos(sqrt(y.d0)) ;   
    
    if( chi_2 == 1. ){ //Boundary points
      double complex Punc_BL_Bound_complex = Punc_BL_Bound[0] + I* Punc_BL_Bound[1],
        dPunc_BL_dr_Bound_complex = dPunc_BL_Bound_dx[2] + I* dPunc_BL_Bound_dx[3],
        dPunc_BL_dtheta_Bound_complex = dPunc_BL_Bound_dx[4] + I* dPunc_BL_Bound_dx[5], 
        dPunc_BL_dsigma_Bound_complex, dPunc_BL_dy_Bound_complex;   

      Hyp_phi_Punc_Bound_complex = Get_HypFunc_From_BLFunc(*par, sigma.d0, y.d0, 1, Punc_BL_Bound_complex);  
      get_HypDerv_from_BLDerv(*par, sigma.d0, y.d0, 1,  dPunc_BL_dr_Bound_complex, dPunc_BL_dtheta_Bound_complex, &dPunc_BL_dsigma_Bound_complex, &dPunc_BL_dy_Bound_complex);
      get_Derivatives_HypFunc_From_BLFunc(*par, sigma.d0, y.d0, 1, Punc_BL_Bound_complex,dPunc_BL_dsigma_Bound_complex, dPunc_BL_dy_Bound_complex, &dHyp_phi_Punc_dsig_Bound_complex, &dHyp_phi_Punc_dy_Bound_complex);
      
      Re_Hyp_Punc_Bound[i1] = creal(Hyp_phi_Punc_Bound_complex); Im_Hyp_Punc_Bound[i1] = cimag(Hyp_phi_Punc_Bound_complex);
      Re_dHyp_Punc_dsigma_Bound[i1] = creal(dHyp_phi_Punc_dsig_Bound_complex); Im_dHyp_Punc_dsigma_Bound[i1] = cimag(dHyp_phi_Punc_dsig_Bound_complex);
      Re_dHyp_Punc_dy_Bound[i1] = creal(dHyp_phi_Punc_dy_Bound_complex); Im_dHyp_Punc_dy_Bound[i1] = cimag(dHyp_phi_Punc_dy_Bound_complex);
    }
    Seff_BL_complex = Seff_BL[0] + I* Seff_BL[1];
    Hyp_Seff_complex = Get_HypFunc_From_BLFunc(*par, sigma.d0, y.d0, 1, Seff_BL_complex); 
    Hyp_Seff_complex = Rescale_Source(*par, sigma.d0, y.d0, Hyp_Seff_complex);

    Re_Hyp_Seff[i1][i2] = creal(Hyp_Seff_complex); Im_Hyp_Seff[i1][i2] = cimag(Hyp_Seff_complex);

    double complex Punc_BL_Bound_complex = Punc_BL_Bound[0] + I* Punc_BL_Bound[1],
                    Hyp_phi_Punc_complex = Get_HypFunc_From_BLFunc(*par, sigma.d0, y.d0, 1, Punc_BL_Bound_complex); 
     
    Re_Hyp_Punc[i1][i2] = creal(Hyp_phi_Punc_complex); Im_Hyp_Punc[i1][i2] = cimag(Hyp_phi_Punc_complex);

    if(fabs(Seff_BL[0])>1.e1){
    printf("dr=%3.15e, theta/pi=%lf, Seff=%lf, Hyp_Seff=%lf \n", fabs(r-par->r0_over_M), theta/Pi, Seff_BL[0],creal(Hyp_Seff_complex) );
    }

    fprintf(fp, "%3.15e %3.15e %3.15e %3.15e %3.15e %3.15e \n",chi_1, chi_2, sigma.d0, y.d0, creal(Hyp_Seff_complex), cimag(Hyp_Seff_complex) );
        

  }
  fprintf(fp, "\n");
 }
 fclose(fp);

 (*par).Re_cheb_phi_Punc = dvector(0, N1); Chebyshev_Coefficients(Re_Hyp_Punc_Bound, (*par).Re_cheb_phi_Punc, N1, (*par).grid_1_PuncSeff); 
 (*par).Im_cheb_phi_Punc = dvector(0, N1); Chebyshev_Coefficients(Im_Hyp_Punc_Bound, (*par).Im_cheb_phi_Punc, N1, (*par).grid_1_PuncSeff);
  
 (*par).Re_cheb_phi_Punc_sigma = dvector(0, N1); Chebyshev_Coefficients(Re_dHyp_Punc_dsigma_Bound, (*par).Re_cheb_phi_Punc_sigma, N1, (*par).grid_1_PuncSeff);
 (*par).Im_cheb_phi_Punc_sigma = dvector(0, N1); Chebyshev_Coefficients(Im_dHyp_Punc_dsigma_Bound, (*par).Im_cheb_phi_Punc_sigma, N1, (*par).grid_1_PuncSeff);
 
 (*par).Re_cheb_phi_Punc_y = dvector(0, N1); Chebyshev_Coefficients(Re_dHyp_Punc_dy_Bound, (*par).Re_cheb_phi_Punc_y, N1, (*par).grid_1_PuncSeff);
 (*par).Im_cheb_phi_Punc_y = dvector(0, N1); Chebyshev_Coefficients(Im_dHyp_Punc_dy_Bound, (*par).Im_cheb_phi_Punc_y, N1, (*par).grid_1_PuncSeff);
 
 (*par).Re_cheb_Seff = dmatrix(0, N1, 0, N2); Chebyshev_Coefficients_2D(Re_Hyp_Seff, (*par).Re_cheb_Seff, N1, (*par).grid_1_PuncSeff, N2, (*par).grid_2_PuncSeff); output_Seff_SpecCoef(*par, Re_Hyp_Seff, "Re_Hyp_Seff");
 (*par).Im_cheb_Seff = dmatrix(0, N1, 0, N2); Chebyshev_Coefficients_2D(Im_Hyp_Seff, (*par).Im_cheb_Seff, N1, (*par).grid_1_PuncSeff, N2, (*par).grid_2_PuncSeff); output_Seff_SpecCoef(*par, Im_Hyp_Seff, "Im_Hyp_Seff");
 
 
 

 free_dvector(Re_Hyp_Punc_Bound, 0, N1); free_dvector(Im_Hyp_Punc_Bound, 0, N1);
 free_dvector(Re_dHyp_Punc_dsigma_Bound, 0, N1); free_dvector(Im_dHyp_Punc_dsigma_Bound, 0, N1); 
 free_dvector(Re_dHyp_Punc_dy_Bound, 0, N1); free_dvector(Im_dHyp_Punc_dy_Bound, 0, N1);
 free_dmatrix(Re_Hyp_Seff, 0, N1, 0, N2); free_dmatrix(Im_Hyp_Seff, 0, N1, 0, N2);
 
 free_dmatrix(Re_Hyp_Punc, 0, N1, 0, N2); free_dmatrix(Im_Hyp_Punc, 0, N1, 0, N2);
 


 return;
}
//----------------------------------------------------------------------------------------
double complex Get_HypFunc_From_BLFunc(parameters par, double sigma, double y, int FLAG_NS, double complex BL_Func){
  

  double complex HypFunc, Z, dlnZ_dsigma, dlnZ_dy;
  func_Z(par, sigma, y, FLAG_NS,  &Z, &dlnZ_dsigma, &dlnZ_dy);
 
  HypFunc = BL_Func/Z;
  
  return HypFunc;
 }
 //----------------------------------------------------------------------------------------
void get_HypDerv_from_BLDerv(parameters par, double sigma, double y, int FLAG_NS,  double complex d_phi_dr, double complex d_phi_dtheta, double complex *d_phi_dsigma, double complex *d_phi_dy){
  double r, dr_dsigma, theta = FLAG_NS*acos(sqrt(y)), cos_theta = cos(theta), sin_theta = sin(theta);
  
  func_r_of_sigma(par, sigma, &r, &dr_dsigma);

  *d_phi_dsigma = dr_dsigma * d_phi_dr;  

  double dy_dtheta = -2* sin_theta * cos_theta;
  *d_phi_dy = 1/(dy_dtheta) * d_phi_dtheta;
  return;
}
 //----------------------------------------------------------------------------------------
 void get_Derivatives_HypFunc_From_BLFunc(parameters par, double sigma, double y, int FLAG_NS,
                     double complex BL_Func, double complex dBL_Func_dsigma, double complex dBL_Func_dy,
                                          double complex *dHypFunc_dsigma, double complex *dHypFunc_dy){   
 
   double complex Z, dlnZ_dsigma, dlnZ_dy;
   func_Z(par, sigma, y, FLAG_NS,  &Z, &dlnZ_dsigma, &dlnZ_dy);  

  
   *dHypFunc_dsigma = dBL_Func_dsigma/Z - BL_Func * dlnZ_dsigma/Z;
   *dHypFunc_dy     = dBL_Func_dy/Z - BL_Func * dlnZ_dy/Z;
   
   return;
  }
  //----------------------------------------------------------------------------------------
 double complex Rescale_Source(parameters par, double sigma, double y, double complex SourceIn){
   
 
   double complex SourceOut;
   double r, dr_dsigma, Sigma;
   func_r_of_sigma(par, sigma, &r, &dr_dsigma);
   func_Sigma(par, r, y, &Sigma);
  
   SourceOut = Sigma*SourceIn;
   
   return SourceOut;
  }
