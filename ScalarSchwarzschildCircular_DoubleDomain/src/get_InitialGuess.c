#include "Solve_PDE.h"

void get_InitialGuess(parameters par, double *X){

//-------ZERO VECTOR-----------------------------------
	fill0_dvector(X, 0, par.ntotal);

//-------READ FROM FILE-----------------------------------
 // input_Solution(par, X);


  // int iDom;
  // for(iDom=0; iDom<nDom; iDom++){
  //   int j1, j2, N1=par.N1[iDom], N2=par.N2[iDom];
  //   for(j1=0; j1<=N1; j1++){
  //     for(j2=0; j2<=N2; j2++){
  //       int indx_Re_phi = Index(par, iDom, 0,  j1,j2), indx_Im_phi = Index(par, iDom, 1,  j1,j2);
  //       double chi_1 = par.grid_chi_1[iDom][j1], 
  //              chi_2 = par.grid_chi_2[iDom][j2];
  //       func_derivs_2D sigma;

  //       get_sigma(par, iDom, chi_1, chi_2, &sigma);      

  //       X[indx_Re_phi]=0.;
  //       X[indx_Im_phi]=0.;
  //     }
  //   }
  // }

  return;
}
//-------------------------------------------------------------------------------
void input_Solution(parameters par, double *X){
  
  
  //------------REFERENCE SIMULATION DATA-------------------------------------------------
  char RefSim[200], fn[200];
  int m = 0;
  double r0_over_M = 10., eta = 2.; 
  
  sprintf(RefSim, "StaticMode_r0%lf_m%d_eta%lf", r0_over_M, m, eta);
  //-------------------------------------------------------------------------------------
  
  
  int iF, j1, j2, j1read, j2read, N1read, N2read, iDom, N1, N2;
  double **c;
  FILE *fr;
  char grid1_read[50], grid2_read[50];
  
  for(iDom=0; iDom<nDom; iDom++){
    N1 = par.N1[iDom];
    N2 = par.N2[iDom];
    for(iF=0; iF<nFields; iF++){
    
      sprintf(fn, "InputData/%s/ChebSolution_dom_%d_fields_%d.dat",RefSim, iDom, iF); 
      
        
      fr=fopen(fn ,"r");
      if(fr==NULL){
        printf("Error in input_Solution:\nFile doesn't exist\n%s\n", fn);
        exit(1);
      }
      
      
      fscanf(fr, "#N1 = %d\tgrid1 = %s\n", &N1read, grid1_read);
      fscanf(fr, "#N2 = %d\tgrid2 = %s\n", &N2read, grid2_read);
      fscanf(fr, "#1:j1\t2:j2\t3:c[j1][j2]\n");
      
      c=dmatrix(0,N1read, 0, N2read);
          
      
      for(j1=0; j1<=N1read; j1 ++){ 
        for(j2=0; j2<=N2read; j2 ++){
            fscanf(fr, "%d\t%d\t%lf\n", &j1read, &j2read, &(c[j1][j2]) );

  	
  	
  	         if(j1!=j1read ||j2!=j2read ){
              fprintf(stderr, "Error in input_Solution: j1 = %d j1read = %d j2 = %d j2read = %d\n", j1, j1read, j2, j2read);
              exit(1);
            }
        }
      }         
      
      for(j1=0; j1<=N1; j1 ++){
        double chi_1 = par.grid_chi_1[iDom][j1];
  	     
        for(j2=0; j2<=N2; j2 ++){
          double chi_2 = par.grid_chi_2[iDom][j2];

          int Indx = Index(par, iDom, iF, j1, j2);
          X[Indx]= Clenshaw_Chebyshev_2D(c, N1read, N2read, chi_1, chi_2);	
        }
      }
       
      fclose(fr);
      free_dmatrix(c, 0,N1read, 0, N2read);      
    }
  }
  
//   X[par.Ntotal-1] = Omega_read;
//   X[par.Ntotal] = 1.-gamma;
// //   printf("%f, %f\n", X[par.Ntotal-1], X[par.Ntotal]);
// //   exit(1);

}
//--------------------------------------------------------------------
void get_InitialGuess_test(parameters par, double *X){
    
    //-------test function VECTOR-----------------------------------
    int j1, j2;
    int iDom;
    double chi_1, chi_2;
    func_derivs_2D sigma, y;
    func_sigma_derivs_2D f,g;
   
    fill0_dvector(X, 0, par.ntotal);
    
    for(iDom=0; iDom<nDom; iDom++){
        int N1=par.N1[iDom], N2=par.N2[iDom];
        for(j1=0; j1<=N1; j1++){
            for(j2=0; j2<=N2; j2++){
                chi_1 = par.grid_chi_1[iDom][j1];
                chi_2 = par.grid_chi_2[iDom][j2];
                
                get_sigma( par, iDom, chi_1, chi_2, &sigma);
                get_y( par, iDom, chi_1, chi_2, &y);

                f = Test_Func(par, sigma.d0, y.d0);
                g = Test_Func(par, sigma.d0, y.d0);
                
                //printf("j1=%d, j2=%d, iDom=%d, %f, %f,%f, %f\n", j1,j2,iDom,creal(f.d0),cimag(f.d0),creal(g.d0),cimag(g.d0));
                
                int indx_Re_Even = Index(par, iDom, 0, j1, j2),
                    indx_Im_Even = Index(par, iDom, 1, j1, j2);
                    // indx_Re_Odd  = Index(par, iDom, 2, j1, j2),
                    // indx_Im_Odd  = Index(par, iDom, 3, j1, j2);
                X[indx_Re_Even]=creal(f.d0);
                X[indx_Im_Even]=cimag(f.d0);
                // X[indx_Re_Odd] =creal(g.d0);
                // X[indx_Im_Odd] =cimag(g.d0);
            }
        }
    }
    
}
//------------
