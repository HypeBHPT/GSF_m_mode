#include "Solve_ODE.h"

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
//-------------------------------------------------------------------------------------------------------
void Get_Derivatives(parameters par, derivs W)
{
	double 	*Y, *d1Y, *d2Y;
	int iDom, iField, j1, N, n;
	  
	
	for(iDom=1; iDom<=nDom; iDom++){
		N=par.N[iDom-1]; n = N+1;
		Y = dvector(0, N);   d1Y = dvector(0, N);   d2Y = dvector(0, N);
	

		for(iField=0; iField < nFields; iField++){
		
			// Computing 1st and 2nd x1-Derivatives
				for(j1=0; j1<n; j1 ++) {
					int k = Index(par, iDom, iField, j1);
					Y[j1] = W.d0[k];
				}

				Spectral_derivative(par.grid[iDom-1], Y, d1Y, d2Y, N);
			
	
				for(j1=0; j1<n; j1 ++){
						int k = Index(par, iDom, iField, j1);
						W.d1[k]  = d1Y[j1];
						W.d11[k] = d2Y[j1];
				}
		}
	
		free_dvector(Y, 0, N);  free_dvector(d1Y, 0, N); free_dvector(d2Y, 0, N); 
	}  
}