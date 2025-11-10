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
double complex Get_phi_From_v(parameters par, double sigma, double y, double complex hyp_phi){
  
  double complex phi, Z, dlnZ;  
  func_Z(par, sigma, y,  &Z, &dlnZ);
  
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
double get_ToyScndOrdSource(parameters par, complex_sigma_derivs phi_Phys, complex_sigma_derivs conjugate_phi_Phys){
  double m = 1.*par.m, m2=sqr(m);

  double complex Toy2ndSource;

  Toy2ndSource = m2 * phi_Phys.d0*conjugate_phi_Phys.d0;// +
                // m  * (phi_Phys.d0*conjugate_phi_Phys.dsigma + phi_Phys.d0*conjugate_phi_Phys.dy) +
                // phi_Phys.dsigma*conjugate_phi_Phys.dsigma + phi_Phys.dsigma*conjugate_phi_Phys.dy + phi_Phys.dy*conjugate_phi_Phys.dy +
                // phi_Phys.d0*conjugate_phi_Phys.d2sigma + phi_Phys.d0*conjugate_phi_Phys.d2sigmay  + phi_Phys.d0*conjugate_phi_Phys.d2y;

  return 2*creal(Toy2ndSource); 
}
//----------------------------------------------------------------------------------------
void get_ToyScndOrdSource_L2Norm(parameters par, int iDom, derivs_2D Sol, double *L2_Toy2ndSource_RR, double *L2_Toy2ndSource_SR){
  int i1, i2,  N1 = par.N1[iDom], N2 = par.N2[iDom];
  double  Jac,  **S_2nd_RR = dmatrix(0,N1, 0,N2), **cS_2nd_RR = dmatrix(0,N1, 0,N2),
                **S_2nd_SR = dmatrix(0,N1, 0,N2), **cS_2nd_SR = dmatrix(0,N1, 0,N2);
  
  complex_derivs_2D phi_reg, conjugate_phi_reg;  
  complex_sigma_derivs phi_reg_Phys, conjugate_phi_reg_Phys, phi_punc_Phys, conjugate_phi_punc_Phys;
  double Toy2ndSource_RR, Toy2ndSource_SR;
  

  for(i1=0; i1<=N1; i1++){
    for(i2=0; i2<=N2; i2++){
      int indx_real = Index(par, iDom, 0, i1, i2),
          indx_imag = Index(par, iDom, 1, i1, i2);

      Jac = get_Jacobian(par, iDom, i1, i2);

      get_ComplexField(par, indx_real, indx_imag, Sol , &phi_reg); 
      conjugate_phi_reg = get_complex_dervs_conjugate(phi_reg); 

      get_PhysDerv_from_SpecDerv_complex(par, iDom, i1, i2, phi_reg, &phi_reg_Phys); 
      get_PhysDerv_from_SpecDerv_complex(par, iDom, i1, i2, conjugate_phi_reg, &conjugate_phi_reg_Phys);
      Toy2ndSource_RR = get_ToyScndOrdSource(par, phi_reg_Phys, conjugate_phi_reg_Phys);
      S_2nd_RR[i1][i2] = sqr(Toy2ndSource_RR) * Jac;
      
      double chi_1 = par.grid_chi_1[iDom][i1], chi_2 = par.grid_chi_2[iDom][i2]; 
      phi_punc_Phys.d0 = Clenshaw_Chebyshev_2D(par.Re_cheb_PuncField, par.N1_LoadSeff, par.N2_LoadSeff, chi_1, chi_2) + I * Clenshaw_Chebyshev_2D(par.Im_cheb_PuncField, par.N1_LoadSeff, par.N2_LoadSeff, chi_1, chi_2);
      conjugate_phi_punc_Phys.d0 = Clenshaw_Chebyshev_2D(par.Re_cheb_PuncField, par.N1_LoadSeff, par.N2_LoadSeff, chi_1, chi_2) - I * Clenshaw_Chebyshev_2D(par.Im_cheb_PuncField, par.N1_LoadSeff, par.N2_LoadSeff, chi_1, chi_2);
      Toy2ndSource_SR = get_ToyScndOrdSource(par, phi_reg_Phys, conjugate_phi_punc_Phys);
      S_2nd_SR[i1][i2] = sqr(Toy2ndSource_SR) * Jac;

    }
  }

  Chebyshev_Coefficients_2D(S_2nd_RR, cS_2nd_RR, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
  Chebyshev_Coefficients_2D(S_2nd_SR, cS_2nd_SR, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
  

  *L2_Toy2ndSource_RR = sqrt( Chebyshev_Definite_2DIntegral_Coefficients(cS_2nd_RR, N1, N2) ) ;
  *L2_Toy2ndSource_SR = iDom == par.idom_particle?  sqrt( Chebyshev_Definite_2DIntegral_Coefficients(cS_2nd_SR, N1, N2) ) : sqrt(-1.) ;



  free_dmatrix(S_2nd_RR, 0 , N1, 0, N2);   free_dmatrix(cS_2nd_RR, 0 , N1,0, N2);
  free_dmatrix(S_2nd_SR, 0 , N1, 0, N2);   free_dmatrix(cS_2nd_SR, 0 , N1,0, N2);

  return;
}
//----------------------------------------------------------------------------------------
void Get_SelfForce_Ft_m(parameters par, double *X, double complex *Ft_m_bulk, double complex *Ft_m_equator){
  //NOTE: EXPRESSION FOR A SINGLE m>0
  //NOTE: OUTPUT WILL TAKE INTO ACCOUNT THE CONTRIBUTIONS FROM m AND -m
  int idom_part = par.idom_particle, N1, N2,  j1, j2;
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

  phi_bulk = Get_phi_From_v(par, par.sigma0, 0., phi_hyp_bulk);
  phi_hyp_equator = Get_phi_From_v(par, par.sigma0, 0., phi_hyp_equator);

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
  int idom_part = par.idom_particle, N1, N2,  j1, j2;
  double **Re_Sol, **Re_c, **Im_Sol, **Im_c, 
         **Re_Sol_sigma, **Re_c_sigma, **Im_Sol_sigma, **Im_c_sigma, H, dH_dsigma,
         q =par.q, sigma0 =par.sigma0, sigma0_2 = sqr(sigma0), rh = par.rh_over_M;
  double complex s=par.s, Z, dlnZ, phi_hyp_bulk, phi_hyp_equator, phi_hyp_sigma_bulk, phi_hyp_sigma_equator;  
  func_Z(par, sigma0, 0.,  &Z, &dlnZ);
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
