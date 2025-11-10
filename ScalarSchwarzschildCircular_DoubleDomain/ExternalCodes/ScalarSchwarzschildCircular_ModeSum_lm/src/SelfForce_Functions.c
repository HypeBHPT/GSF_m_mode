#include "Solve_ODE.h"

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
double complex Get_phi_From_v(parameters par, double sigma, double complex V){

  double complex phi, Z, dlnZ;
  func_Z(par, sigma, &Z, &dlnZ);

  phi = Z*V;


  return phi;
}

//----------------------------------------------------------------------------------------
double complex Get_psi_From_w(parameters par, double sigma, double complex V, double complex W){

  double complex psi, Z, dlnZ_dsigma, dlnZ_dr0;

  double x, dx_dsigma, h, dh_dsigma;
  func_tortoise_x(par, sigma, &x, &dx_dsigma);
  func_height(par, sigma, &h, &dh_dsigma);
  func_Z(par, sigma, &Z, &dlnZ_dsigma);

  dlnZ_dr0 = par.ds_dr0*(x+h);

  psi = Z*(W+dlnZ_dr0*V);
  return psi;
}
//----------------------------------------------------------------------------------------
double complex Get_psi_From_w_at_particle(parameters par, double complex V, double complex W){

  double complex psi, Z, dlnZ_dsigma, dlzZ_dr0;

  double x, dx_dsigma, h, dh_dsigma, sig_0, dsigma_dr0; // = - 0.5*sqr(par.sigma0);
  func_sigma_of_r(par, par.r0_over_M, &sig_0, &dsigma_dr0);
  func_tortoise_x(par, par.sigma0, &x, &dx_dsigma);
  func_height(par, par.sigma0, &h, &dh_dsigma);
  func_Z(par, par.sigma0, &Z, &dlnZ_dsigma);

  dlzZ_dr0 = par.ds_dr0*(x+h) + dsigma_dr0*dlnZ_dsigma;

  psi = Z*(W+dlzZ_dr0*V);


  return psi;
}
//----------------------------------------------------------------------------------------
void get_Ft_SSF_lm(parameters par, double *X, double complex *Ft_SSF_lm_dom1, double complex *Ft_SSF_lm_dom2){
  double sigma0=par.sigma0, lambda=par.lambda_over_M, l=1.*par.ell, m=1.*par.m, q = par.q,
         
         c_lm = sqrt( (2.* l+1.)*factrl(l-m)/(4*Pi*factrl(l+m))),
         P_lm = plgndr(l, m, 0.), 
         Y_lm = c_lm*P_lm;

  double complex s = par.s, v_dom1, phi_dom1, v_dom2, phi_dom2, dphi_dt_dom1, dphi_dt_dom2;

  int I_real_particle_dom1 = Index(par, 1, 0, 0),
      I_imag_particle_dom1 = Index(par, 1, 1, 0),
      I_real_particle_dom2 = Index(par, 2, 0, par.N[1]),
      I_imag_particle_dom2 = Index(par, 2, 1, par.N[1]);

      v_dom1 = X[I_real_particle_dom1] + I*X[I_imag_particle_dom1];
      
      phi_dom1 = Get_phi_From_v(par, sigma0, v_dom1);
      dphi_dt_dom1 = s*phi_dom1;

      v_dom2 = X[I_real_particle_dom2] + I*X[I_imag_particle_dom2];
      phi_dom2 = Get_phi_From_v(par, sigma0, v_dom2);
      dphi_dt_dom2 = s*phi_dom2;


  *Ft_SSF_lm_dom1 = q*dphi_dt_dom1 * Y_lm/lambda;
  *Ft_SSF_lm_dom2 = q*dphi_dt_dom2 * Y_lm/lambda;


  return;
}
//----------------------------------------------------------------------------------------
void get_Flux_lm(parameters par, double *X, double *Flux_lm_hrzn, double *Flux_lm_scri, double *Flux_lm){
  int I_real_scri = Index(par, 1, 0, par.N[0]),
      I_imag_scri = Index(par, 1, 1, par.N[0]),
      I_real_hrz = Index(par, 2, 0, 0),
      I_imag_hrz = Index(par, 2, 1, 0);
  double lambda = par.lambda_over_M;
  double complex s=par.s, v_scri, v_hrz, dv_dt_scri, dv_dt_hrz;

  v_scri = X[I_real_scri] + I*X[I_imag_scri];
  dv_dt_scri = v_scri*s;

  v_hrz  = X[I_real_hrz] + I*X[I_imag_hrz];
  dv_dt_hrz = v_hrz*s;

  *Flux_lm_hrzn = sqr(cabs(dv_dt_hrz)) /(4*Pi*sqr(2*lambda));///(4*Pi*Pi*sqr(2*lambda))/(sqr(0.5)) * 4;
  *Flux_lm_scri = sqr(cabs(dv_dt_scri)) /(4*Pi*sqr(2*lambda));///(4*Pi*Pi*sqr(2*lambda))/(sqr(0.5)) * 4;
  

  *Flux_lm = *Flux_lm_hrzn + *Flux_lm_scri;


  return;
}
//----------------------------------------------------------------------------------------
double get_Ft_FLux_lm(parameters par, double *X){
  double Flux_lm_hrzn, Flux_lm_scri, Flux_lm, ut , Ft_lm_flux;

  get_Flux_lm(par, X, &Flux_lm_hrzn, &Flux_lm_scri, &Flux_lm);
  
  ut=par.E0/par.f0; 
  
  Ft_lm_flux = Flux_lm * ut;

  return Ft_lm_flux;
}
//----------------------------------------------------------------------------------------
double get_Ft_PN(parameters par){
  double r0 = par.r0_over_M, V = sqrt(1./r0),
         V2 = sqr(V), V3=V2*V, V4=V3*V, V5 = V4*V,
         Ft_PN;

  Ft_PN = V4*(1./3  - V2/6. + 2*Pi*V3/3. - 77.*V4/24. + 9*Pi*V5/5.)/sqr(r0);
  return Ft_PN;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void Get_Bfield(parameters par, double *Bfield, double complex *Hyp_Bfield){
  double q=par.q,
       r0=par.r0_over_M, r0_2=sqr(r0),
       L0 = r0/Sqrt(r0-3.), L0_2=sqr(L0),
       f0, df0, ak, K;

  double complex Z0, dlnZ0;

  func_f(r0, &f0, &df0);
  ak=sqrt(1./(r0*f0));
  K=ellf(Pih, ak);


  *Bfield = 2*q*K/(Pi*Sqrt(L0_2 + r0_2));

  func_Z(par, par.sigma0, &Z0, &dlnZ0);
  *Hyp_Bfield = (*Bfield)*r0/Z0;
}

//----------------------------------------------------------------------------------------------------------------------------------------
void Sum_over_m(parameters par, double *X, double *sum_ell_Dom1, double *sum_ell_Dom2){

  int Idx_Particle_Re_DomLeft = Index(par, 1, 0, 0),Idx_Particle_Im_DomLeft = Index(par, 1, 1, 0),
      Idx_Particle_Re_DomRight = Index(par, 2, 0, par.N[0]),Idx_Particle_Im_DomRight = Index(par, 2, 1, par.N[0]);




  double complex V_DomLeft, V_DomRight, phi_DomLeft, phi_DomRight;

      V_DomLeft  = X[Idx_Particle_Re_DomLeft]+I*X[Idx_Particle_Im_DomLeft];
      V_DomRight = X[Idx_Particle_Re_DomRight]+I*X[Idx_Particle_Im_DomRight];

      phi_DomLeft= Get_phi_From_v(par, par.sigma0, V_DomLeft);
      phi_DomRight= Get_phi_From_v(par, par.sigma0, V_DomRight);

      double c_lm = sqrt( (2.*par.ell+1.)*factrl(par.ell-par.m)/(4*Pi*factrl(par.ell+par.m))),
       P_lm = plgndr(par.ell, par.m, 0.);

      if (par.m == 0 ){
  *sum_ell_Dom1 += creal(phi_DomLeft)*c_lm*P_lm;
  *sum_ell_Dom2 += creal(phi_DomRight)*c_lm*P_lm;
      }
      else{
  *sum_ell_Dom1 += 2.*creal(phi_DomLeft)*c_lm*P_lm;
  *sum_ell_Dom2 += 2.*creal(phi_DomRight)*c_lm*P_lm;
      }



}
