#include "Solve_ODE.h"


//---------------------------------------------------------------
void Initial_Guess(parameters par, double *X){
  
  int iDom, i, N;

  for(iDom=1; iDom<=nDom; iDom++){
  	N=par.N[iDom-1];
  	for(i=0; i<=N; i++){
    	  double sigma, dx_dsigma;
  		  get_sigma(par, iDom, i, &sigma, &dx_dsigma);
   
    	int indx_Re_phi = Index(par, iDom, 0,  i), 
	    	  indx_Im_phi = Index(par, iDom, 1,  i);
	    

   		
				X[indx_Re_phi]=-1.;//(1+sigma);
				X[indx_Im_phi]=1.;//sqr(sigma);   


  	}
  }

}
