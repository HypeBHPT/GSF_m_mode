#include "Solve_PDE.h"

void Spectral_derivative(char *grid, double *f, double *df, double *d2f, int N){
  
  
  if(strcmp( grid,"Radau_RHS" ) ==0){
    double c[N+1], d1c[N+1], d2c[N+1];
    Chebyshev_Coefficients_Radau_RHS(f, c, N);
    Chebyshev_Coefficients_Derivative(c, d1c, N);
    Chebyshev_Collocations_Radau_RHS(df, d1c, N);
    Chebyshev_Coefficients_Derivative(d1c, d2c, N);
    Chebyshev_Collocations_Radau_RHS(d2f, d2c, N);
  }
  else if(strcmp( grid,"Radau_LHS" ) ==0){
    double c[N+1], d1c[N+1], d2c[N+1];
    Chebyshev_Coefficients_Radau_LHS(f, c, N);
    Chebyshev_Coefficients_Derivative(c, d1c, N);
    Chebyshev_Collocations_Radau_LHS(df, d1c, N);
    Chebyshev_Coefficients_Derivative(d1c, d2c, N);
    Chebyshev_Collocations_Radau_LHS(d2f, d2c, N);
  }
  else if(strcmp( grid,"Lobatto" ) ==0){
    double c[N+1], d1c[N+1], d2c[N+1];
    Chebyshev_Coefficients_Lobatto(f, c, N);
    Chebyshev_Coefficients_Derivative(c, d1c, N);
    Chebyshev_Collocations_Lobatto(df, d1c, N);
    Chebyshev_Coefficients_Derivative(d1c, d2c, N);
    Chebyshev_Collocations_Lobatto(d2f, d2c, N);
  }
  else if(strcmp( grid,"Gauss" ) ==0){
    double c[N+1], d1c[N+1], d2c[N+1];
    Chebyshev_Coefficients_Gauss(f, c, N);
    Chebyshev_Coefficients_Derivative(c, d1c, N);
    Chebyshev_Collocations_Gauss(df, d1c, N);
    Chebyshev_Coefficients_Derivative(d1c, d2c, N);
    Chebyshev_Collocations_Gauss(d2f, d2c, N);
  }
  else if(strcmp( grid,"Fourier" ) ==0){
    double alpha[N+1], d1alpha[N+1], d2alpha[N+1],
	   beta[N+1], d1beta[N+1], d2beta[N+1];
    Fourier_Coefficients(f, alpha, beta, N);
    Fourier_Coefficients_Derivative(alpha, beta, d1alpha, d1beta, N);
    Fourier_Collocations(df, d1alpha, d1beta, N);

    Fourier_Coefficients_Derivative(d1alpha, d1beta, d2alpha, d2beta, N);
    Fourier_Collocations(d2f, d2alpha, d2beta, N);   
  }
  else{
    printf("Error in Spectral_derivative: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n", grid);
    exit(1);
  }
  
  
}

void Get_Derivatives(parameters par, derivs_2D W)
{
	double 	*Y, *d1Y, *d2Y, 
	        *Z, *d1Z, *d2Z,
			*A, *d1A, *d2A,
			*B, *d1B, *d2B,
			*C, *d1C, *d2C,
			*D, *d1D, *d2D,
			*E, *d1E, *d2E;
	int iDom, iF, j1, j2, N1, N2, n1, n2, N;

	for(iDom=0; iDom<nDom; iDom++){
		N1=par.N1[iDom]; n1=par.n1[iDom];
		N2=par.N2[iDom]; n2=par.n2[iDom];
		N = maximum2(N1, 2*N2);
		Y = dvector(0, N);   d1Y = dvector(0, N);   d2Y = dvector(0, N);	
		Z = dvector(0, N);   d1Z = dvector(0, N);   d2Z = dvector(0, N);
		A = dvector(0, N);   d1A = dvector(0, N);   d2A = dvector(0, N);
		B = dvector(0, N);   d1B = dvector(0, N);   d2B = dvector(0, N);
		C = dvector(0, N);   d1C = dvector(0, N);   d2C = dvector(0, N);
		D = dvector(0, N);   d1D = dvector(0, N);   d2D = dvector(0, N);
		E = dvector(0, N);   d1E = dvector(0, N);   d2E = dvector(0, N);
		

		for(iF=0; iF < nFields; iF++){
			
			// Computing 1st, 2nd and 3rd x1-Derivatives
			for(j2=0; j2<n2; j2 ++){

				for(j1=0; j1<n1; j1 ++) {
					int k = Index(par, iDom, iF, j1, j2);
					Y[j1] = W.d0[k];
				}
				Spectral_derivative(par.grid_1[iDom], Y, d1Y, d2Y, N1);
				Spectral_derivative(par.grid_1[iDom], d2Y, d1A, d2A, N1);		
		
				for(j1=0; j1<n1; j1 ++){
					int k = Index(par, iDom, iF, j1, j2);
					W.d1[k]    = d1Y[j1];
					W.d11[k]   = d2Y[j1];
					W.d111[k]  = d1A[j1];
					W.d1111[k] = d2A[j1];
				}
			}
			// Computing 1st and 2nd x2-Derivatives as well as mixed (x1,x2)-Derivatives
			for(j1=0; j1<n1; j1 ++){
				for(j2=0; j2<n2; j2 ++){
					int k = Index(par, iDom, iF, j1, j2);
					Y[j2] =   W.d0[k];
					Z[j2] =   W.d1[k];
					A[j2] =   W.d11[k];
					C[j2] =   W.d111[k];
				}
				Spectral_derivative(par.grid_2[iDom], Y, d1Y, d2Y, N2);
				Spectral_derivative(par.grid_2[iDom], Z, d1Z, d2Z, N2);
				Spectral_derivative(par.grid_2[iDom], A, d1A, d2A, N2);
				Spectral_derivative(par.grid_2[iDom], d2Y, d1B, d2B, N2);
				Spectral_derivative(par.grid_2[iDom], C, d1C, d2C, N2);
				Spectral_derivative(par.grid_2[iDom], d2Z, d1D, d2D, N2);
				Spectral_derivative(par.grid_2[iDom], d2A, d1E, d2E, N2);

				for(j2=0; j2<n2; j2 ++){
					int k = Index(par, iDom, iF, j1, j2);
					W.d2[k]  = d1Y[j2];
					W.d22[k] = d2Y[j2];					
					W.d12[k]  = d1Z[j2];

					W.d112[k] = d1A[j2];
					W.d122[k] = d2Z[j2];
					W.d222[k] = d1B[j2];

					W.d1112[k]  = d1C[j2];
					W.d1122[k] = d2A[j2];
					W.d1222[k] = d1D[j2];
					W.d2222[k] = d2B[j2];

					W.d11222[k] = d1E[j2];
				}
			}


		}
		free_dvector(Y, 0, N);  free_dvector(d1Y, 0, N); free_dvector(d2Y, 0, N);
		free_dvector(Z, 0, N);  free_dvector(d1Z, 0, N); free_dvector(d2Z, 0, N);
		free_dvector(A, 0, N);  free_dvector(d1A, 0, N); free_dvector(d2A, 0, N);
		free_dvector(B, 0, N);  free_dvector(d1B, 0, N); free_dvector(d2B, 0, N);
		free_dvector(C, 0, N);  free_dvector(d1C, 0, N); free_dvector(d2C, 0, N);
		free_dvector(D, 0, N);  free_dvector(d1D, 0, N); free_dvector(d2D, 0, N);
		free_dvector(E, 0, N);  free_dvector(d1E, 0, N); free_dvector(d2E, 0, N);
	}
	return;  
}
// -------------------------------------------------------------------------------
void Get_DerivativesFinDif(parameters par, derivs_2D v)
{
	int n1, n2;
	    

	int i, j,i_field, idom;	 
	for(idom=0;idom<nDom; idom++){
		n1=par.n1[idom]; n2=par.n2[idom];
		for(i_field=0; i_field<nFields; i_field++)
			for(j=0; j<n2; j++)
				for(i=0; i < n1; i++)
					Get_DerivativesFinDif_grid( par,  i,  j,  i_field, idom, v);

	}
}
// -------------------------------------------------------------------------------
void Get_DerivativesFinDif_grid(parameters par, int i, int j, int i_field, int idom, derivs_2D v)
{
	int N1=par.N1[idom], N2=par.N2[idom]; 

	  	// Calculation of Derivatives w.r.t. x1-Direction (here called direction a)  (Finite Dif)
		int indx  = Index(par, idom, i_field, i, j),
		indx_a_P1 = Index(par, idom, i_field, i+1, j),
		indx_a_M1 = Index(par, idom, i_field, i-1, j),
		indx_a_P2 = Index(par, idom, i_field, i+2, j),
		indx_a_M2 = Index(par, idom, i_field, i-2, j),
		indx_a_P3 = Index(par, idom, i_field, i+3, j),
		indx_a_M3 = Index(par, idom, i_field, i-3, j),
		indx_a_P4 = Index(par, idom, i_field, i+4, j),
		indx_a_M4 = Index(par, idom, i_field, i-4, j),

		
		indx_b_P1 = Index(par, idom, i_field, i, j+1),
		indx_b_M1 = Index(par, idom, i_field, i, j-1),
		indx_b_P2 = Index(par, idom, i_field, i, j+2),
		indx_b_M2 = Index(par, idom, i_field, i, j-2),
		indx_b_P3 = Index(par, idom, i_field, i, j+3),
		indx_b_M3 = Index(par, idom, i_field, i, j-3),
		indx_b_P4 = Index(par, idom, i_field, i, j+4),
		indx_b_M4 = Index(par, idom, i_field, i, j-4),
		
		indx_a_P3_b_M3 = Index(par, idom, i_field, i+3, j-3),
		indx_a_P3_b_M2 = Index(par, idom, i_field, i+3, j-2),
		indx_a_P3_b_M1 = Index(par, idom, i_field, i+3, j-1),
		indx_a_P3_b_P1 = Index(par, idom, i_field, i+3, j+1),
		indx_a_P3_b_P2 = Index(par, idom, i_field, i+3, j+2),
		indx_a_P3_b_P3 = Index(par, idom, i_field, i+3, j+3),
		
		indx_a_P2_b_M3 = Index(par, idom, i_field, i+2, j-3),
		indx_a_P2_b_M2 = Index(par, idom, i_field, i+2, j-2),
		indx_a_P2_b_M1 = Index(par, idom, i_field, i+2, j-1),
		indx_a_P2_b_P1 = Index(par, idom, i_field, i+2, j+1),
		indx_a_P2_b_P2 = Index(par, idom, i_field, i+2, j+2),
		indx_a_P2_b_P3 = Index(par, idom, i_field, i+2, j+3),
		
		indx_a_P1_b_M3 = Index(par, idom, i_field, i+1, j-3),
		indx_a_P1_b_M2 = Index(par, idom, i_field, i+1, j-2),
		indx_a_P1_b_M1 = Index(par, idom, i_field, i+1, j-1),
		indx_a_P1_b_P1 = Index(par, idom, i_field, i+1, j+1),
		indx_a_P1_b_P2 = Index(par, idom, i_field, i+1, j+2),
		indx_a_P1_b_P3 = Index(par, idom, i_field, i+1, j+3),

		indx_a_M3_b_M3 = Index(par, idom, i_field, i-3, j-3),
		indx_a_M3_b_M2 = Index(par, idom, i_field, i-3, j-2),
		indx_a_M3_b_M1 = Index(par, idom, i_field, i-3, j-1),
		indx_a_M3_b_P1 = Index(par, idom, i_field, i-3, j+1),
		indx_a_M3_b_P2 = Index(par, idom, i_field, i-3, j+2),
		indx_a_M3_b_P3 = Index(par, idom, i_field, i-3, j+3),
		
		indx_a_M2_b_M3 = Index(par, idom, i_field, i-2, j-3),
		indx_a_M2_b_M2 = Index(par, idom, i_field, i-2, j-2),
		indx_a_M2_b_M1 = Index(par, idom, i_field, i-2, j-1),
		indx_a_M2_b_P1 = Index(par, idom, i_field, i-2, j+1),
		indx_a_M2_b_P2 = Index(par, idom, i_field, i-2, j+2),
		indx_a_M2_b_P3 = Index(par, idom, i_field, i-2, j+3),
		
		indx_a_M1_b_M3 = Index(par, idom, i_field, i-1, j-3),
		indx_a_M1_b_M2 = Index(par, idom, i_field, i-1, j-2),
		indx_a_M1_b_M1 = Index(par, idom, i_field, i-1, j-1),
		indx_a_M1_b_P1 = Index(par, idom, i_field, i-1, j+1),
		indx_a_M1_b_P2 = Index(par, idom, i_field, i-1, j+2),
		indx_a_M1_b_P3 = Index(par, idom, i_field, i-1, j+3);			

		double da, d2a, d4a,ha, ha_2, ha_4, ca, ca2, sa, sa2, //ca3, sa3,sa4,
		       db, d2b, d4b,hb, hb_2, hb_4, cb, cb2, sb, sb2, //cb3, sb3,sb4,
		       dab, d2ab, d2ba, d2a2b;
		
		// Calculation of Derivatives w.r.t. x1-Direction (here called direction a)  (Finite Dif)
		ha =get_AngleIncrement(par.grid_1[idom], N1); 
		ha_2=sqr(ha); 
		ha_4 = sqr(ha_2);
		
		ca=get_grid(par.grid_1[idom], i, N1); 
		ca2=sqr(ca); 
		// ca3=ca*ca2; 
		sa=sqrt(1.-ca2); 
		sa2=sqr(sa);
		// sa3=sa*sa2;
		// sa4=sqr(sa2);
		
		if( FD_ORDER ==2 ){
		    da =  (0.5*v.d0[indx_a_P1]- 0.5*v.d0[indx_a_M1])/(ha),
		    d2a = ( v.d0[indx_a_M1] - 2*v.d0[indx] + v.d0[indx_a_P1] )/(ha_2);
		    d4a = (v.d0[indx_a_M2] -4*v.d0[indx_a_M1] + 6*v.d0[indx] - 4*v.d0[indx_a_P1] + v.d0[indx_a_P2])/(ha_4);
		}		
		else if (FD_ORDER == 4){
		da =( (1./12)*v.d0[indx_a_M2] -(2./3)*v.d0[indx_a_M1] + (2./3)*v.d0[indx_a_P1] - (1./12)*v.d0[indx_a_P2])/(ha);
		d2a=(-(1./12)*v.d0[indx_a_M2] +(4./3)*v.d0[indx_a_M1] -(5./2)*v.d0[indx] + (4./3)*v.d0[indx_a_P1] - (1./12)*v.d0[indx_a_P2])/(ha_2);
		d4a=(-(1./6)*v.d0[indx_a_M3] + 2*v.d0[indx_a_M2] - (13./2)*v.d0[indx_a_M1] +(28./3)*v.d0[indx] - (13./2)*v.d0[indx_a_P1] + 2*v.d0[indx_a_P2] - (1./6)*v.d0[indx_a_P3])/(ha_4);
		}
		else if (FD_ORDER == 6){
		da =( -(1./60)*v.d0[indx_a_M3] + (3./20)*v.d0[indx_a_M2] -(3./4)*v.d0[indx_a_M1] + (3./4)*v.d0[indx_a_P1] - (3./20)*v.d0[indx_a_P2]+ (1./60)*v.d0[indx_a_P3])/(ha);
		d2a=((1./90)*v.d0[indx_a_M3]-(3./20)*v.d0[indx_a_M2] + (3./2)*v.d0[indx_a_M1] -(49./18)*v.d0[indx] + (3./2)*v.d0[indx_a_P1] - (3./20)*v.d0[indx_a_P2] + (1./90)*v.d0[indx_a_P3])/(ha_2);
		d4a=((7./240)*v.d0[indx_a_M4]-(2./5)*v.d0[indx_a_M3] + (169./60)*v.d0[indx_a_M2] - (122./15)*v.d0[indx_a_M1] + (91./8)*v.d0[indx] - (122./15)*v.d0[indx_a_P1] + (169./60)*v.d0[indx_a_P2] - (2./5)*v.d0[indx_a_P3]+(7./240)*v.d0[indx_a_P4])/(ha_4);
		}
		
		if(strcmp( par.grid_1[idom],"Fourier" ) ==0){
		    v.d1[indx]=da;
		    v.d11[indx]=d2a;
		}
		else{ 
		    if( fabs(fabs(ca)-1.) < TINY){
		    v.d1[indx] = -d2a/ca;
		    v.d11[indx] = (d4a + d2a)/3.;
		    }
		    else
		    {
		    v.d1[indx] = -da/sa;
		    v.d11[indx] = (d2a + ca*v.d1[indx])/sa2;
		    }
		}
		
		// Calculation of Derivatives w.r.t. x2-Direction (here called direction b) (Finite Dif)
		hb=get_AngleIncrement(par.grid_2[idom], N2); 
		hb_2=sqr(hb); 
		hb_4=sqr(hb_2);
		
		cb=get_grid(par.grid_2[idom], j, N2); 
		cb2=sqr(cb); 
		// cb3=cb*cb2; 
		sb=sqrt(1.-cb2); 
		sb2=sqr(sb); 
		// sb3=sb*sb2; 
		// sb4=sqr(sb2);
		

		if( FD_ORDER ==2 ){
		    db =  (0.5*v.d0[indx_b_P1]- 0.5*v.d0[indx_b_M1])/(hb),
		    d2b = ( v.d0[indx_b_M1] - 2*v.d0[indx] + v.d0[indx_b_P1] )/(hb_2);
		    d4b = (v.d0[indx_b_M2] -4*v.d0[indx_b_M1] + 6*v.d0[indx] - 4*v.d0[indx_b_P1] + v.d0[indx_b_P2])/(hb_4);
		}		
		else if (FD_ORDER == 4){
		db =( (1./12)*v.d0[indx_b_M2] -(2./3)*v.d0[indx_b_M1] + (2./3)*v.d0[indx_b_P1] - (1./12)*v.d0[indx_b_P2])/(hb);
		d2b=(-(1./12)*v.d0[indx_b_M2] +(4./3)*v.d0[indx_b_M1] -(5./2)*v.d0[indx] + (4./3)*v.d0[indx_b_P1] - (1./12)*v.d0[indx_b_P2])/(hb_2);
		d4b=(-(1./6)*v.d0[indx_b_M3] + 2*v.d0[indx_b_M2] - (13./2)*v.d0[indx_b_M1] +(28./3)*v.d0[indx] - (13./2)*v.d0[indx_b_P1] + 2*v.d0[indx_b_P2] - (1./6)*v.d0[indx_b_P3])/(hb_4);
		}
		else if (FD_ORDER == 6){
		db =( -(1./60)*v.d0[indx_b_M3] + (3./20)*v.d0[indx_b_M2] -(3./4)*v.d0[indx_b_M1] + (3./4)*v.d0[indx_b_P1] - (3./20)*v.d0[indx_b_P2]+ (1./60)*v.d0[indx_b_P3])/(hb);
		d2b=((1./90)*v.d0[indx_b_M3]-(3./20)*v.d0[indx_b_M2] + (3./2)*v.d0[indx_b_M1] -(49./18)*v.d0[indx] + (3./2)*v.d0[indx_b_P1] - (3./20)*v.d0[indx_b_P2] + (1./90)*v.d0[indx_b_P3])/(hb_2);
		d4b=((7./240)*v.d0[indx_b_M4]-(2./5)*v.d0[indx_b_M3] + (169./60)*v.d0[indx_b_M2] - (122./15)*v.d0[indx_b_M1] + (91./8)*v.d0[indx] - (122./15)*v.d0[indx_b_P1] + (169./60)*v.d0[indx_b_P2] - (2./5)*v.d0[indx_b_P3]+(7./240)*v.d0[indx_b_P4])/(hb_4);
		}
		
		if(strcmp( par.grid_2[idom],"Fourier" ) ==0){
		    v.d2[indx]=db;
		    v.d22[indx]=d2b;
		}
		else{ 
		    if( fabs(fabs(cb)-1.) < TINY){
		    v.d2[indx] = -d2b/cb;
		    v.d22[indx] = (d4b + d2b)/3.;
		    }
		    else
		    {
		    v.d2[indx] = -db/sb;
		    v.d22[indx] = (d2b + cb*v.d2[indx])/sb2;
		    }
		}
		
	  // Calculation of mixed Derivatives w.r.t. x1 and x2-Directions (here called directions a and b) (Finite Dif)
		if( FD_ORDER ==2 ){
		  dab = (
		          0.5*(0.5*v.d0[indx_a_P1_b_P1]- 0.5*v.d0[indx_a_P1_b_M1])/(hb) 
		        - 0.5*(0.5*v.d0[indx_a_M1_b_P1]- 0.5*v.d0[indx_a_M1_b_M1])/(hb) 
			)/(ha);
		  
		  d2ab = ( 
			  (0.5*v.d0[indx_a_M1_b_P1]- 0.5*v.d0[indx_a_M1_b_M1])/(hb)
			- 2*(0.5*v.d0[indx_b_P1]- 0.5*v.d0[indx_b_M1])/(hb) 
			+ (0.5*v.d0[indx_a_P1_b_P1]- 0.5*v.d0[indx_a_P1_b_M1])/(hb)
			)/(ha_2);
			
		  d2ba =  (
			    0.5*( v.d0[indx_a_P1_b_M1] - 2*v.d0[indx_a_P1] + v.d0[indx_a_P1_b_P1] )/(hb_2)
			  - 0.5*( v.d0[indx_a_M1_b_M1] - 2*v.d0[indx_a_M1] + v.d0[indx_a_M1_b_P1] )/(hb_2)
			  )/(ha);
		  
		  d2a2b = ( 
			     ( v.d0[indx_a_M1_b_M1] - 2*v.d0[indx_a_M1] + v.d0[indx_a_M1_b_P1] )/(hb_2) 
			    -2*( v.d0[indx_b_M1] - 2*v.d0[indx] + v.d0[indx_b_P1] )/(hb_2)
			    +( v.d0[indx_a_P1_b_M1] - 2*v.d0[indx_a_P1] + v.d0[indx_a_P1_b_P1] )/(hb_2)
			  )/(ha_2);
		

		}		
		else if (FD_ORDER == 4){
		dab =( 
		      (1./12)*( (1./12)*v.d0[indx_a_M2_b_M2] -(2./3)*v.d0[indx_a_M2_b_M1] + (2./3)*v.d0[indx_a_M2_b_P1] - (1./12)*v.d0[indx_a_M2_b_P2])/(hb)
		     -(2./3)*( (1./12)*v.d0[indx_a_M1_b_M2] -(2./3)*v.d0[indx_a_M1_b_M1] + (2./3)*v.d0[indx_a_M1_b_P1] - (1./12)*v.d0[indx_a_M1_b_P2])/(hb)
		    + (2./3)*( (1./12)*v.d0[indx_a_P1_b_M2] -(2./3)*v.d0[indx_a_P1_b_M1] + (2./3)*v.d0[indx_a_P1_b_P1] - (1./12)*v.d0[indx_a_P1_b_P2])/(hb)
		    - (1./12)*( (1./12)*v.d0[indx_a_P2_b_M2] -(2./3)*v.d0[indx_a_P2_b_M1] + (2./3)*v.d0[indx_a_P2_b_P1] - (1./12)*v.d0[indx_a_P2_b_P2])/(hb)
		    )/(ha);
		    
		d2ab=(
		      -(1./12)*( (1./12)*v.d0[indx_a_M2_b_M2] -(2./3)*v.d0[indx_a_M2_b_M1] + (2./3)*v.d0[indx_a_M2_b_P1] - (1./12)*v.d0[indx_a_M2_b_P2])/(hb)
		      +(4./3)*( (1./12)*v.d0[indx_a_M1_b_M2] -(2./3)*v.d0[indx_a_M1_b_M1] + (2./3)*v.d0[indx_a_M1_b_P1] - (1./12)*v.d0[indx_a_M1_b_P2])/(hb)
		      -(5./2)*( (1./12)*v.d0[indx_b_M2] -(2./3)*v.d0[indx_b_M1] + (2./3)*v.d0[indx_b_P1] - (1./12)*v.d0[indx_b_P2])/(hb) 
		      +(4./3)*( (1./12)*v.d0[indx_a_P1_b_M2] -(2./3)*v.d0[indx_a_P1_b_M1] + (2./3)*v.d0[indx_a_P1_b_P1] - (1./12)*v.d0[indx_a_P1_b_P2])/(hb)
		      - (1./12)*( (1./12)*v.d0[indx_a_P2_b_M2] -(2./3)*v.d0[indx_a_P2_b_M1] + (2./3)*v.d0[indx_a_P2_b_P1] - (1./12)*v.d0[indx_a_P2_b_P2])/(hb)
		    )/(ha_2);
		d2ba =( 
			(1./12)*(-(1./12)*v.d0[indx_a_M2_b_M2] +(4./3)*v.d0[indx_a_M2_b_M1] -(5./2)*v.d0[indx_a_M2] + (4./3)*v.d0[indx_a_M2_b_P1] - (1./12)*v.d0[indx_a_M2_b_P2])/(hb_2)
			-(2./3)*(-(1./12)*v.d0[indx_a_M1_b_M2] +(4./3)*v.d0[indx_a_M1_b_M1] -(5./2)*v.d0[indx_a_M1] + (4./3)*v.d0[indx_a_M1_b_P1] - (1./12)*v.d0[indx_a_M1_b_P2])/(hb_2)
			+(2./3)*(-(1./12)*v.d0[indx_a_P1_b_M2] +(4./3)*v.d0[indx_a_P1_b_M1] -(5./2)*v.d0[indx_a_P1] + (4./3)*v.d0[indx_a_P1_b_P1] - (1./12)*v.d0[indx_a_P1_b_P2])/(hb_2)
			-(1./12)*(-(1./12)*v.d0[indx_a_P2_b_M2] +(4./3)*v.d0[indx_a_P2_b_M1] -(5./2)*v.d0[indx_a_P2] + (4./3)*v.d0[indx_a_P2_b_P1] - (1./12)*v.d0[indx_a_P2_b_P2])/(hb_2)
		      )/(ha);
		      
		d2a2b = ( 
			(-(1./12)*v.d0[indx_a_M1_b_M2] +(4./3)*v.d0[indx_a_M1_b_M1] -(5./2)*v.d0[indx_a_M1] + (4./3)*v.d0[indx_a_M1_b_P1] - (1./12)*v.d0[indx_a_M1_b_P2])/(hb_2) 
		     - 2*(-(1./12)*v.d0[indx_b_M2] +(4./3)*v.d0[indx_b_M1] -(5./2)*v.d0[indx] + (4./3)*v.d0[indx_b_P1] - (1./12)*v.d0[indx_b_P2])/(hb_2) 
		     + (-(1./12)*v.d0[indx_a_P1_b_M2] +(4./3)*v.d0[indx_a_P1_b_M1] -(5./2)*v.d0[indx_a_P1] + (4./3)*v.d0[indx_a_P1_b_P1] - (1./12)*v.d0[indx_a_P1_b_P2])/(hb_2) 
		      )/(ha_2);

		}
		else if (FD_ORDER == 6){
		dab =( 
			-(1./60)*( -(1./60)*v.d0[indx_a_M3_b_M3] + (3./20)*v.d0[indx_a_M3_b_M2] -(3./4)*v.d0[indx_a_M3_b_M1] + (3./4)*v.d0[indx_a_M3_b_P1] - (3./20)*v.d0[indx_a_M3_b_P2]+ (1./60)*v.d0[indx_a_M3_b_P3])/(hb) 
		        +(3./20)*( -(1./60)*v.d0[indx_a_M2_b_M3] + (3./20)*v.d0[indx_a_M2_b_M2] -(3./4)*v.d0[indx_a_M2_b_M1] + (3./4)*v.d0[indx_a_M2_b_P1] - (3./20)*v.d0[indx_a_M2_b_P2]+ (1./60)*v.d0[indx_a_M2_b_P3])/(hb)
		        -(3./4)*( -(1./60)*v.d0[indx_a_M1_b_M3] + (3./20)*v.d0[indx_a_M1_b_M2] -(3./4)*v.d0[indx_a_M1_b_M1] + (3./4)*v.d0[indx_a_M1_b_P1] - (3./20)*v.d0[indx_a_M1_b_P2]+ (1./60)*v.d0[indx_a_M1_b_P3])/(hb) 
		        +(3./4)*( -(1./60)*v.d0[indx_a_P1_b_M3] + (3./20)*v.d0[indx_a_P1_b_M2] -(3./4)*v.d0[indx_a_P1_b_M1] + (3./4)*v.d0[indx_a_P1_b_P1] - (3./20)*v.d0[indx_a_P1_b_P2]+ (1./60)*v.d0[indx_a_P1_b_P3])/(hb) 
		        -(3./20)*( -(1./60)*v.d0[indx_a_P2_b_M3] + (3./20)*v.d0[indx_a_P2_b_M2] -(3./4)*v.d0[indx_a_P2_b_M1] + (3./4)*v.d0[indx_a_P2_b_P1] - (3./20)*v.d0[indx_a_P2_b_P2]+ (1./60)*v.d0[indx_a_P2_b_P3])/(hb)
		        +(1./60)*( -(1./60)*v.d0[indx_a_P3_b_M3] + (3./20)*v.d0[indx_a_P3_b_M2] -(3./4)*v.d0[indx_a_P3_b_M1] + (3./4)*v.d0[indx_a_P3_b_P1] - (3./20)*v.d0[indx_a_P3_b_P2]+ (1./60)*v.d0[indx_a_P3_b_P3])/(hb)
		    )/(ha);
		
		
		
		
		d2ab=(
		       (1./90)*( -(1./60)*v.d0[indx_a_M3_b_M3] + (3./20)*v.d0[indx_a_M3_b_M2] -(3./4)*v.d0[indx_a_M3_b_M1] + (3./4)*v.d0[indx_a_M3_b_P1] - (3./20)*v.d0[indx_a_M3_b_P2]+ (1./60)*v.d0[indx_a_M3_b_P3])/(hb)
		      -(3./20)*( -(1./60)*v.d0[indx_a_M2_b_M3] + (3./20)*v.d0[indx_a_M2_b_M2] -(3./4)*v.d0[indx_a_M2_b_M1] + (3./4)*v.d0[indx_a_M2_b_P1] - (3./20)*v.d0[indx_a_M2_b_P2]+ (1./60)*v.d0[indx_a_M2_b_P3])/(hb) 
		      +(3./2)*( -(1./60)*v.d0[indx_a_M1_b_M3] + (3./20)*v.d0[indx_a_M1_b_M2] -(3./4)*v.d0[indx_a_M1_b_M1] + (3./4)*v.d0[indx_a_M1_b_P1] - (3./20)*v.d0[indx_a_M1_b_P2]+ (1./60)*v.d0[indx_a_M1_b_P3])/(hb) 
		      -(49./18)*( -(1./60)*v.d0[indx_b_M3] + (3./20)*v.d0[indx_b_M2] -(3./4)*v.d0[indx_b_M1] + (3./4)*v.d0[indx_b_P1] - (3./20)*v.d0[indx_b_P2]+ (1./60)*v.d0[indx_b_P3])/(hb)
		      +(3./2)*( -(1./60)*v.d0[indx_a_P1_b_M3] + (3./20)*v.d0[indx_a_P1_b_M2] -(3./4)*v.d0[indx_a_P1_b_M1] + (3./4)*v.d0[indx_a_P1_b_P1] - (3./20)*v.d0[indx_a_P1_b_P2]+ (1./60)*v.d0[indx_a_P1_b_P3])/(hb) 
		      -(3./20)*( -(1./60)*v.d0[indx_a_P2_b_M3] + (3./20)*v.d0[indx_a_P2_b_M2] -(3./4)*v.d0[indx_a_P2_b_M1] + (3./4)*v.d0[indx_a_P2_b_P1] - (3./20)*v.d0[indx_a_P2_b_P2]+ (1./60)*v.d0[indx_a_P2_b_P3])/(hb) 
		      +(1./90)*( -(1./60)*v.d0[indx_a_P3_b_M3] + (3./20)*v.d0[indx_a_P3_b_M2] -(3./4)*v.d0[indx_a_P3_b_M1] + (3./4)*v.d0[indx_a_P3_b_P1] - (3./20)*v.d0[indx_a_P3_b_P2]+ (1./60)*v.d0[indx_a_P3_b_P3])/(hb)
		    )/(ha_2);
		    
	       d2ba =( 
			-(1./60)*((1./90)*v.d0[indx_a_M3_b_M3]-(3./20)*v.d0[indx_a_M3_b_M2] + (3./2)*v.d0[indx_a_M3_b_M1] -(49./18)*v.d0[indx_a_M3] + (3./2)*v.d0[indx_a_M3_b_P1] - (3./20)*v.d0[indx_a_M3_b_P2] + (1./90)*v.d0[indx_a_M3_b_P3])/(hb_2) 
			+(3./20)*((1./90)*v.d0[indx_a_M2_b_M3]-(3./20)*v.d0[indx_a_M2_b_M2] + (3./2)*v.d0[indx_a_M2_b_M1] -(49./18)*v.d0[indx_a_M2] + (3./2)*v.d0[indx_a_M2_b_P1] - (3./20)*v.d0[indx_a_M2_b_P2] + (1./90)*v.d0[indx_a_M2_b_P3])/(hb_2) 
			-(3./4)*((1./90)*v.d0[indx_a_M1_b_M3]-(3./20)*v.d0[indx_a_M1_b_M2] + (3./2)*v.d0[indx_a_M1_b_M1] -(49./18)*v.d0[indx_a_M1] + (3./2)*v.d0[indx_a_M1_b_P1] - (3./20)*v.d0[indx_a_M1_b_P2] + (1./90)*v.d0[indx_a_M1_b_P3])/(hb_2)
			+(3./4)*((1./90)*v.d0[indx_a_P1_b_M3]-(3./20)*v.d0[indx_a_P1_b_M2] + (3./2)*v.d0[indx_a_P1_b_M1] -(49./18)*v.d0[indx_a_P1] + (3./2)*v.d0[indx_a_P1_b_P1] - (3./20)*v.d0[indx_a_P1_b_P2] + (1./90)*v.d0[indx_a_P1_b_P3])/(hb_2) 
			-(3./20)*((1./90)*v.d0[indx_a_P2_b_M3]-(3./20)*v.d0[indx_a_P2_b_M2] + (3./2)*v.d0[indx_a_P2_b_M1] -(49./18)*v.d0[indx_a_P2] + (3./2)*v.d0[indx_a_P2_b_P1] - (3./20)*v.d0[indx_a_P2_b_P2] + (1./90)*v.d0[indx_a_P2_b_P3])/(hb_2)
			+(1./60)*((1./90)*v.d0[indx_a_P3_b_M3]-(3./20)*v.d0[indx_a_P3_b_M2] + (3./2)*v.d0[indx_a_P3_b_M1] -(49./18)*v.d0[indx_a_P3] + (3./2)*v.d0[indx_a_P3_b_P1] - (3./20)*v.d0[indx_a_P3_b_P2] + (1./90)*v.d0[indx_a_P3_b_P3])/(hb_2)
		     )/(ha);
		
		d2a2b=(
			(1./90)*((1./90)*v.d0[indx_a_M3_b_M3]-(3./20)*v.d0[indx_a_M3_b_M2] + (3./2)*v.d0[indx_a_M3_b_M1] -(49./18)*v.d0[indx_a_M3] + (3./2)*v.d0[indx_a_M3_b_P1] - (3./20)*v.d0[indx_a_M3_b_P2] + (1./90)*v.d0[indx_a_M3_b_P3])/(hb_2)
		       -(3./20)*((1./90)*v.d0[indx_a_M2_b_M3]-(3./20)*v.d0[indx_a_M2_b_M2] + (3./2)*v.d0[indx_a_M2_b_M1] -(49./18)*v.d0[indx_a_M2] + (3./2)*v.d0[indx_a_M2_b_P1] - (3./20)*v.d0[indx_a_M2_b_P2] + (1./90)*v.d0[indx_a_M2_b_P3])/(hb_2) 
		       +(3./2)*((1./90)*v.d0[indx_a_M1_b_M3]-(3./20)*v.d0[indx_a_M1_b_M2] + (3./2)*v.d0[indx_a_M1_b_M1] -(49./18)*v.d0[indx_a_M1] + (3./2)*v.d0[indx_a_M1_b_P1] - (3./20)*v.d0[indx_a_M1_b_P2] + (1./90)*v.d0[indx_a_M1_b_P3])/(hb_2) 
		       -(49./18)*((1./90)*v.d0[indx_b_M3]-(3./20)*v.d0[indx_b_M2] + (3./2)*v.d0[indx_b_M1] -(49./18)*v.d0[indx] + (3./2)*v.d0[indx_b_P1] - (3./20)*v.d0[indx_b_P2] + (1./90)*v.d0[indx_b_P3])/(hb_2) 
		       +(3./2)*((1./90)*v.d0[indx_a_P1_b_M3]-(3./20)*v.d0[indx_a_P1_b_M2] + (3./2)*v.d0[indx_a_P1_b_M1] -(49./18)*v.d0[indx_a_P1] + (3./2)*v.d0[indx_a_P1_b_P1] - (3./20)*v.d0[indx_a_P1_b_P2] + (1./90)*v.d0[indx_a_P1_b_P3])/(hb_2) 
		       -(3./20)*((1./90)*v.d0[indx_a_P2_b_M3]-(3./20)*v.d0[indx_a_P2_b_M2] + (3./2)*v.d0[indx_a_P2_b_M1] -(49./18)*v.d0[indx_a_P2] + (3./2)*v.d0[indx_a_P2_b_P1] - (3./20)*v.d0[indx_a_P2_b_P2] + (1./90)*v.d0[indx_a_P2_b_P3])/(hb_2) 
		       +(1./90)*((1./90)*v.d0[indx_a_P3_b_M3]-(3./20)*v.d0[indx_a_P3_b_M2] + (3./2)*v.d0[indx_a_P3_b_M1] -(49./18)*v.d0[indx_a_P3] + (3./2)*v.d0[indx_a_P3_b_P1] - (3./20)*v.d0[indx_a_P3_b_P2] + (1./90)*v.d0[indx_a_P3_b_P3])/(hb_2)
		    )/(ha_2);

		}
		
		if(strcmp( par.grid_1[idom],"Fourier" ) ==0){		  
		    if(strcmp( par.grid_2[idom],"Fourier" ) ==0){
		      v.d12[indx]=dab;
		    }
		    else{
		      if( fabs(fabs(cb)-1.) < TINY){
						v.d12[indx] = -d2ba/cb;
		      }
		      else{
						v.d12[indx] = -dab/sb;
		      }
		    }
		}
		else if(strcmp( par.grid_2[idom],"Fourier" ) ==0){
		  
		    if(strcmp( par.grid_1[idom],"Fourier" ) ==0){
		      v.d12[indx]=dab;
		    }
		    else{
		      if( fabs(fabs(ca)-1.) < TINY){
			v.d12[indx] = -d2ab/ca;
		      }
		      else{
			v.d12[indx] = -dab/sa;
		      }
		    }
		}
		else{
		    if( fabs(fabs(ca)-1.) > TINY && fabs(fabs(cb)-1.) > TINY ){
		      v.d12[indx] = dab/(sa*sb);
		    }
		    else if( fabs(fabs(ca)-1.) < TINY && fabs(fabs(cb)-1.) > TINY ){
		      v.d12[indx] = d2ab/(ca*sb);
		    }
		    else if( fabs(fabs(ca)-1.) > TINY && fabs(fabs(cb)-1.) < TINY){
		      v.d12[indx] = d2ba/(cb*sa);
		    }
		    else
		    {
		      v.d12[indx] = d2a2b/(ca*cb);		    
		    }
		}
		
		return;
}
// // -------------------------------------------------------------------------------
