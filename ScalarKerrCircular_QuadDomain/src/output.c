#include "Solve_PDE.h"
//-------------------------------------------------------------------------------
void output_Solution(parameters par, double *X){
  
  int iDom, iF, j1, j2, j1_max=50, j2_max=50, N1, N2;
  FILE *fp, *fp_SpcCoord;
  char fn[200], fn_SpcCoord[200];
  double **Sol, **c;
	
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

      
      sprintf(fn, "data/%s/Solution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
      sprintf(fn_SpcCoord, "data/%s/Solution_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
      fp=fopen(fn ,"w");
      fp_SpcCoord=fopen(fn_SpcCoord ,"w");

      fprintf(fp, "#1: sigma\t 2:y\t 3:Solution \n");
      fprintf(fp_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution \n");

      for(j1=0; j1<=j1_max; j1++){
        double chi_1 = -1. + 2.*j1/j1_max ;
        func_derivs_2D sigma;

        for(j2=0; j2<=j2_max; j2++){
          double chi_2 = -1. + 2.*j2/j2_max;
          func_derivs_2D y;

          get_sigma(par, iDom, chi_1, chi_2, &sigma);
          get_y(par, iDom, chi_1, chi_2, &y);
          
          double out = Clenshaw_Chebyshev_2D(c, N1, N2, chi_1, chi_2);
          fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out );
          fprintf(fp_SpcCoord, "%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out );
        }
        fprintf(fp, "\n");
        fprintf(fp_SpcCoord, "\n");
      }
      fclose(fp);
      fclose(fp_SpcCoord);
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
//-------------------------------------------------------------------------------
void output_SpecCoef(parameters par, double *X){
  
  FILE *fp;
  char fn[200];  
  int iDom, iField, j1, j2, N1, N2;
  double *p, *cp, **c;

  for(iDom=0; iDom<nDom; iDom++){
    N1=par.N1[iDom];
    N2=par.N2[iDom];
    c=dmatrix(0,N1, 0, N2);
    for(iField=0; iField<nFields; iField++){
      p=dvector(0, N1); cp=dvector(0, N1);

      for(j2=0; j2<=N2; j2++){        
        for(j1=0; j1<=N1; j1++){
          int indx = Index(par, iDom, iField,j1,j2);
          p[j1]=X[indx];
        }
        Chebyshev_Coefficients(p, cp, N1, par.grid_1[iDom]);
        for(j1=0; j1<=N1; j1++) c[j1][j2] = cp[j1];
      }
      free_dvector(p,0,N1); free_dvector(cp, 0, N1);

      sprintf(fn, "data/%s/ChebCoef_x1_dom_%d_field_%d.dat", par.SimName, iDom, iField);
      fp=fopen(fn,"w");
      fprintf(fp, "#1:i1\t");
      for(j2=0; j2<=N2; j2++){
        fprintf(fp, "#%d: j2=%d (chi_2=%3.3e)\t", j2+2, j2, par.grid_chi_2[iDom][j2] );
      }

      for(j1=0; j1<=N1; j1++){
        fprintf(fp, "\n%d\t", j1);
        for(j2=0; j2<=N2; j2++) fprintf(fp, "%3.15e\t", c[j1][j2] );
      }
      fclose(fp);

      p=dvector(0,N2); cp=dvector(0,N2);
      for(j1=0; j1<=N1; j1++){
        for(j2=0; j2<=N2; j2++){
          int indx = Index(par, iDom, iField,j1,j2);
          p[j2]=X[indx];
        }
        Chebyshev_Coefficients(p, cp, N2, par.grid_2[iDom]);
        for(j2=0; j2<=N2; j2++) c[j1][j2] = cp[j2];
      }
      free_dvector(p,0,N2); free_dvector(cp, 0, N2);

      sprintf(fn, "data/%s/ChebCoef_x2_dom_%d_field_%d.dat", par.SimName, iDom, iField);
      fp=fopen(fn,"w");
      fprintf(fp, "#1:i2\t");
      for(j1=0; j1<=N1; j1++){       
        fprintf(fp, "#%d: j1=%d (chi_1=%3.3e)\t", j1+2, j1,  par.grid_chi_1[iDom][j1]);
      }

      for(j2=0; j2<=N2; j2++){
        fprintf(fp, "\n%d\t", j2);
        for(j1=0; j1<=N1; j1++) fprintf(fp, "%3.15e\t", c[j1][j2] );
      }
      fclose(fp);
    }
    free_dmatrix(c, 0,N1, 0, N2);
  }

  return;
}
//-------------------------------------------------------------------------------
void output_SolutionChebyshevPostProc(parameters par, double *X){
  
  int iDom, iF, j1, j2,  N1, N2;
  FILE *fp;
  char fn[200];
  double **Sol, **c;
	
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

      //OUTPUT 2D - CHEBYSCHEV COEFFICIENTS FOR POST PROCESSING DATA
      // NOTE: Include AnMR when Implemented
      sprintf(fn, "data/%s/ChebSolution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
      fp=fopen(fn ,"w");
      fprintf(fp, "#N1 = %d\tgrid1 = %s\n", par.N1[iDom], par.grid_1[iDom]);
      fprintf(fp, "#N2 = %d\tgrid2 = %s\n", par.N2[iDom], par.grid_2[iDom]);
      fprintf(fp, "#1:j1\t2:j2\t3:c[j1][j2]\n");

      for(j1=0; j1<=N1; j1 ++)
        for(j2=0; j2<=N2; j2 ++)	fprintf(fp, "%d\t%d\t%3.15e\n", j1, j2, c[j1][j2]);

      fclose(fp);
      //------------------------------------------------------------------------------
    
    
    }
    free_dmatrix(Sol,0, N1, 0, N2);
    free_dmatrix(c, 0,N1, 0, N2);
  }


  return;
}
//-------------------------------------------------------------------------------
// void output_Solution_derivatives(parameters par, double *X){
  
//   int iDom, iF, j1, j2, j1_max=50, j2_max=50, N1, N2;
//   FILE *fp, *fp_SpcCoord;
//   char fn[200], fn_SpcCoord[200];
//   double **Sol, **c;

  
	
//   for(iDom=0; iDom<nDom; iDom++){
//     N1=par.N1[iDom];
//     N2=par.N2[iDom];
//     Sol=dmatrix(0, N1, 0, N2); c=dmatrix(0,N1, 0, N2);

//     for(iF=0; iF<nFields; iF++){
    
//       for(j1=0; j1<=N1; j1 ++) {
//         for(j2=0; j2<=N2; j2 ++){
//           int indx = Index(par, iDom, iF, j1, j2);
//           Sol[j1][j2]=Sol_data.d22[indx];
//         }
//       }
//       Chebyshev_Coefficients_2D(Sol, c, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

//       //OUTPUT 2D - CHEBYSCHEV COEFFICIENTS FOR POST PROCESSING DATA
//       // NOTE: Include AnMR when Implemented
//       sprintf(fn, "data/%s/ChebSolution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       fp=fopen(fn ,"w");
//       fprintf(fp, "#N1 = %d\tgrid1 = %s\n", par.N1[iDom], par.grid_1[iDom]);
//       fprintf(fp, "#N2 = %d\tgrid2 = %s\n", par.N2[iDom], par.grid_2[iDom]);
//       fprintf(fp, "#1:j1\t2:j2\t3:c[j1][j2]\n");

//       for(j1=0; j1<=N1; j1 ++)
//         for(j2=0; j2<=N2; j2 ++)	fprintf(fp, "%d\t%d\t%3.15e\n", j1, j2, c[j1][j2]);

//       fclose(fp);
//       //------------------------------------------------------------------------------
    
//       sprintf(fn, "data/%s/Solution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_SpcCoord, "data/%s/Solution_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
//       fp=fopen(fn ,"w");
//       fp_SpcCoord=fopen(fn_SpcCoord ,"w");

//       fprintf(fp, "#1: sigma\t 2:y\t 3:Solution \n");
//       fprintf(fp_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution \n");

//       for(j1=0; j1<=j1_max; j1++){
//         double chi_1 = -1. + 2.*j1/j1_max ;
//         func_derivs_2D sigma;

//         for(j2=0; j2<=j2_max; j2++){
//           double chi_2 = -1. + 2.*j2/j2_max;
//           func_derivs_2D y;

//           get_sigma(par, iDom, chi_1, chi_2, &sigma);
//           get_y(par, iDom, chi_1, chi_2, &y);
          
//           double out = Clenshaw_Chebyshev_2D(c, N1, N2, chi_1, chi_2);
//           fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out );
//           fprintf(fp_SpcCoord, "%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out );
//         }
//         fprintf(fp, "\n");
//         fprintf(fp_SpcCoord, "\n");
//       }
//       fclose(fp);
//       fclose(fp_SpcCoord);
//     }
//     free_dmatrix(Sol,0, N1, 0, N2);
//     free_dmatrix(c, 0,N1, 0, N2);
//   }

//   for(iF=0; iF<nPar; iF++){
//     sprintf(fn, "data/%s/Solution_parameter_%d.dat",par.SimName, iF);
//     fp = fopen(fn, "w");
//     double out = X[par.Ntotal-iF];
//     fprintf(fp, "%3.15e\t\n", out );
//     fclose(fp);
//   }

//   return;
// }
//-------------------------------------------------------------------------------
// void output_FieldEquator(parameters par){
  
//   int iDom, iF, j1, j2, j2_max=50, N1, N2;
//   FILE *fp, *fp_physical;
//   char fn[200], fn_physical[200];
//   double *X;
//   X=dvector(0, par.Ntotal);
//   load_X_from_ChebSolution(par, X);

//   derivs_2D W;
//   allocate_derivs_2D(&W, par.ntotal);
// 	copy_X_to_W(par, X, W);
// 	Get_Derivatives(par,  W);

//   complex_derivs_2D sol_local;
//   complex_sigma_derivs sol_local_Phys;

  
//   double **Sol,  **Sol_d1, **Sol_d11, **Sol_d12, **Sol_d2, **Sol_d22,
//                  **Sol_dsigma, **Sol_d2sigma, **Sol_d2sigmay, **Sol_dy, **Sol_d2y,
//          **cSol, **cSol_d1, **cSol_d11, **cSol_d12, **cSol_d2, **cSol_d22,
//                  **cSol_dsigma, **cSol_d2sigma, **cSol_d2sigmay, **cSol_dy, **cSol_d2y;
	
//   // for(iDom=0; iDom<nDom; iDom++){
//     iDom = par.idom_particle;
//     N1=par.N1[iDom];
//     N2=par.N2[iDom];
//     Sol    =dmatrix(0, N1, 0, N2); cSol=dmatrix(0,N1, 0, N2);
//     Sol_d1 =dmatrix(0, N1, 0, N2); cSol_d1=dmatrix(0,N1, 0, N2);
//     Sol_d11 =dmatrix(0, N1, 0, N2); cSol_d11=dmatrix(0,N1, 0, N2);
//     Sol_d12 =dmatrix(0, N1, 0, N2); cSol_d12=dmatrix(0,N1, 0, N2);
//     Sol_d2  =dmatrix(0, N1, 0, N2); cSol_d2=dmatrix(0,N1, 0, N2);
//     Sol_d22 =dmatrix(0, N1, 0, N2); cSol_d22=dmatrix(0,N1, 0, N2);

//     Sol_dsigma=dmatrix(0, N1, 0, N2); cSol_dsigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigma=dmatrix(0, N1, 0, N2); cSol_d2sigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigmay=dmatrix(0, N1, 0, N2); cSol_d2sigmay=dmatrix(0,N1, 0, N2);
//     Sol_dy=dmatrix(0, N1, 0, N2); cSol_dy=dmatrix(0,N1, 0, N2);
//     Sol_d2y=dmatrix(0, N1, 0, N2); cSol_d2y=dmatrix(0,N1, 0, N2);

//     for(iF=0; iF<nFields; iF++){
    
      
//       for(j1=0; j1<=N1; j1 ++){
//         for(j2=0; j2<=N2; j2 ++){
//           int indx = Index(par, iDom, iF, j1, j2);
//           Sol[j1][j2]     = sol_local.d0  = W.d0[indx];
//           Sol_d1[j1][j2]  = sol_local.d1  = W.d1[indx];
//           Sol_d11[j1][j2] = sol_local.d11  = W.d11[indx];
//           Sol_d12[j1][j2] = sol_local.d12  = W.d12[indx];
//           Sol_d2[j1][j2]  = sol_local.d2  = W.d2[indx];
//           Sol_d22[j1][j2] = sol_local.d22  = W.d22[indx];

//           get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, sol_local, &sol_local_Phys);

//           Sol_dsigma[j1][j2]   = creal(sol_local_Phys.dsigma);
//           Sol_d2sigma[j1][j2]  = creal(sol_local_Phys.d2sigma);
//           Sol_d2sigmay[j1][j2] = creal(sol_local_Phys.d2sigmay);
//           Sol_dy[j1][j2]       = creal(sol_local_Phys.dy);
//           Sol_d2y[j1][j2]      = creal(sol_local_Phys.d2y);

//         }
//       }
//       Chebyshev_Coefficients_2D(Sol, cSol, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d1, cSol_d1, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d11, cSol_d11, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d12, cSol_d12, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2, cSol_d2, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d22, cSol_d22, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

//       Chebyshev_Coefficients_2D(Sol_dsigma, cSol_dsigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigma, cSol_d2sigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigmay, cSol_d2sigmay, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_dy, cSol_dy, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2y, cSol_d2y, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

          
//       sprintf(fn, "data/%s/EquatorSpecGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_physical, "data/%s/EquatorPhysicalGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
      
//       fp=fopen(fn ,"w");
//       fp_physical=fopen(fn_physical ,"w");
      
      
//       fprintf(fp, "#1: chi_1\t 2:chi_2\t 3:Solution\t 4:Solution,1\t 5:Solution,11\t 6:Solution,12\t 7:Solution,2\t 8:Solution,22 \n");
//       fprintf(fp_physical, "#1: chi_1\t 2:chi_2\t 3:Solution\t 4:Solution,sigma\t 5:Solution,sigmasigma\t 6:Solution,sigmay\t 7:Solution,y\t 8:Solution,yy \n");
            
//       // for(j1=0; j1<=j1_max; j1++){
//         double chi_1 = -1.; //-1. + 2.*j1/j1_max ;
//         func_derivs_2D sigma;

//         for(j2=0; j2<=j2_max; j2++){
//           // double chi_2 = -1.;
//           double chi_2 = -1. + 2.*j2/j2_max;
//           func_derivs_2D y;

//           get_sigma(par, iDom, chi_1, chi_2, &sigma);
//           get_y(par, iDom, chi_1, chi_2, &y);
          
//           double out    = Clenshaw_Chebyshev_2D(cSol, N1, N2, chi_1, chi_2),
//                  out_d1 = Clenshaw_Chebyshev_2D(cSol_d1, N1, N2, chi_1, chi_2),
//                  out_d11 = Clenshaw_Chebyshev_2D(cSol_d11, N1, N2, chi_1, chi_2),
//                  out_d12 = Clenshaw_Chebyshev_2D(cSol_d12, N1, N2, chi_1, chi_2),
//                  out_d2 = Clenshaw_Chebyshev_2D(cSol_d2, N1, N2, chi_1, chi_2),
//                  out_d22 = Clenshaw_Chebyshev_2D(cSol_d22, N1, N2, chi_1, chi_2),
                 
//                  out_dsigma = Clenshaw_Chebyshev_2D(cSol_dsigma, N1, N2, chi_1, chi_2),
//                  out_d2sigma = Clenshaw_Chebyshev_2D(cSol_d2sigma, N1, N2, chi_1, chi_2),
//                  out_d2sigmay = Clenshaw_Chebyshev_2D(cSol_d2sigmay, N1, N2, chi_1, chi_2),
//                  out_dy = Clenshaw_Chebyshev_2D(cSol_dy, N1, N2, chi_1, chi_2),
//                  out_d2y = Clenshaw_Chebyshev_2D(cSol_d2y, N1, N2, chi_1, chi_2);
                 
                
          
//           fprintf(fp, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_d1, out_d11, out_d12, out_d2, out_d22 );
//           fprintf(fp_physical, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_dsigma, out_d2sigma, out_d2sigmay, out_dy, out_d2y);
          
//         }

        
//       // }
//       fclose(fp);
//       fclose(fp_physical);
      
      
//     }
//     free_dmatrix(Sol,0, N1, 0, N2); free_dmatrix(cSol, 0,N1, 0, N2);
//     free_dmatrix(Sol_d1,0, N1, 0, N2); free_dmatrix(cSol_d1, 0,N1, 0, N2);
//     free_dmatrix(Sol_d11,0, N1, 0, N2); free_dmatrix(cSol_d11, 0,N1, 0, N2);
//     free_dmatrix(Sol_d12,0, N1, 0, N2); free_dmatrix(cSol_d12, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2,0, N1, 0, N2); free_dmatrix(cSol_d2, 0,N1, 0, N2);
//     free_dmatrix(Sol_d22,0, N1, 0, N2); free_dmatrix(cSol_d22, 0,N1, 0, N2);

//     free_dmatrix(Sol_dsigma,0, N1, 0, N2); free_dmatrix(cSol_dsigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigma,0, N1, 0, N2); free_dmatrix(cSol_d2sigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigmay,0, N1, 0, N2); free_dmatrix(cSol_d2sigmay, 0,N1, 0, N2);
//     free_dmatrix(Sol_dy,0, N1, 0, N2); free_dmatrix(cSol_dy, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2y,0, N1, 0, N2); free_dmatrix(cSol_d2y, 0,N1, 0, N2);
    
//   // }


//   free_derivs_2D(&W, par.ntotal);
//   free_dvector(X, 0, par.Ntotal);
//   return;
// }
//-------------------------------------------------------------------------------
// void output_FieldParticle(parameters par){
  
//   int iDom, iF, j1, j2, j1_max=50, N1, N2;
//   FILE *fp, *fp_physical, *fp_reg2derv;
//   char fn[200], fn_physical[200], fn_reg2derv[200];
//   double *X;
//   X=dvector(0, par.Ntotal);
//   load_X_from_ChebSolution(par, X);

//   derivs_2D W;
//   allocate_derivs_2D(&W, par.ntotal);
// 	copy_X_to_W(par, X, W);
// 	Get_Derivatives(par,  W);

//   complex_derivs_2D sol_local;
//   complex_sigma_derivs sol_local_Phys;

//   double complex Eq_Reg_dsigma, Eq_Reg_d2sigma, Eq_Reg_dy, Eq_Reg_d2y, Eq_Reg_d2sigmay;
//   double sigma0=par.sigma0, sigma0_2=sqr(sigma0),
//          **Sol,  **Sol_d1, **Sol_d11, **Sol_d12, **Sol_d2, **Sol_d22,
//                  **Sol_dsigma, **Sol_d2sigma, **Sol_d2sigmay, **Sol_dy, **Sol_d2y,
//                  **Reg_dsigma, **Reg_d2sigma, **Reg_dy, **Reg_d2y, **Reg_d2sigmay,
//          **cSol, **cSol_d1, **cSol_d11, **cSol_d12, **cSol_d2, **cSol_d22,
//                  **cSol_dsigma, **cSol_d2sigma, **cSol_d2sigmay, **cSol_dy, **cSol_d2y,
//         **cReg_dsigma, **cReg_d2sigma, **cReg_dy, **cReg_d2y, **cReg_d2sigmay;
	
//   // for(iDom=0; iDom<nDom; iDom++){
//     iDom = par.idom_particle;
//     N1=par.N1[iDom];
//     N2=par.N2[iDom];
//     Sol    =dmatrix(0, N1, 0, N2); cSol=dmatrix(0,N1, 0, N2);
//     Sol_d1 =dmatrix(0, N1, 0, N2); cSol_d1=dmatrix(0,N1, 0, N2);
//     Sol_d11 =dmatrix(0, N1, 0, N2); cSol_d11=dmatrix(0,N1, 0, N2);
//     Sol_d12 =dmatrix(0, N1, 0, N2); cSol_d12=dmatrix(0,N1, 0, N2);
//     Sol_d2  =dmatrix(0, N1, 0, N2); cSol_d2=dmatrix(0,N1, 0, N2);
//     Sol_d22 =dmatrix(0, N1, 0, N2); cSol_d22=dmatrix(0,N1, 0, N2);

//     Sol_dsigma=dmatrix(0, N1, 0, N2); cSol_dsigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigma=dmatrix(0, N1, 0, N2); cSol_d2sigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigmay=dmatrix(0, N1, 0, N2); cSol_d2sigmay=dmatrix(0,N1, 0, N2);
//     Sol_dy=dmatrix(0, N1, 0, N2); cSol_dy=dmatrix(0,N1, 0, N2);
//     Sol_d2y=dmatrix(0, N1, 0, N2); cSol_d2y=dmatrix(0,N1, 0, N2);

//     Reg_dsigma =dmatrix(0, N1, 0, N2);   cReg_dsigma = dmatrix(0, N1, 0, N2);
//     Reg_d2sigma =dmatrix(0, N1, 0, N2);  cReg_d2sigma = dmatrix(0, N1, 0, N2);
//     Reg_dy =dmatrix(0, N1, 0, N2);       cReg_dy = dmatrix(0, N1, 0, N2);
//     Reg_d2y =dmatrix(0, N1, 0, N2);      cReg_d2y = dmatrix(0, N1, 0, N2);
//     Reg_d2sigmay =dmatrix(0, N1, 0, N2); cReg_d2sigmay = dmatrix(0, N1, 0, N2);

//     for(iF=0; iF<nFields; iF++){
    
      
//       for(j1=0; j1<=N1; j1 ++){
//         for(j2=0; j2<=N2; j2 ++){
//           int indx = Index(par, iDom, iF, j1, j2);
//           get_ComplexField(par, indx, indx, W , &sol_local);
//           get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, sol_local, &sol_local_Phys);

//           Sol[j1][j2]     = creal(sol_local.d0);  
//           Sol_d1[j1][j2]  = creal(sol_local.d1); 
//           Sol_d11[j1][j2] = creal(sol_local.d11); 
//           Sol_d12[j1][j2] = creal(sol_local.d12);
//           Sol_d2[j1][j2]  = creal(sol_local.d2);  
//           Sol_d22[j1][j2] = creal(sol_local.d22);           

//           Sol_dsigma[j1][j2]   = creal(sol_local_Phys.dsigma);
//           Sol_d2sigma[j1][j2]  = creal(sol_local_Phys.d2sigma);
//           Sol_d2sigmay[j1][j2] = creal(sol_local_Phys.d2sigmay);
//           Sol_dy[j1][j2]       = creal(sol_local_Phys.dy);
//           Sol_d2y[j1][j2]      = creal(sol_local_Phys.d2y);

//           double x1=par.grid_chi_1[iDom][j1];
//           Eq_Reg_dsigma = sol_local.d1;
//           Eq_Reg_d2sigma = sol_local.d11;
//           Eq_Reg_dy = sol_local.d2 - x1*sol_local.d12;
//           Eq_Reg_d2y = sol_local.d122- x1*sol_local.d1122;          
//           Eq_Reg_d2sigmay = sol_local.d222 - (x1*sol_local.d1222 - 1./3. * x1*x1*sol_local.d11222);

//           Reg_dsigma[j1][j2] = creal(Eq_Reg_dsigma);
//           Reg_d2sigma[j1][j2] = creal(Eq_Reg_d2sigma);
//           Reg_dy[j1][j2] = creal(Eq_Reg_dy);
//           Reg_d2y[j1][j2] = creal(Eq_Reg_d2y);
//           Reg_d2sigmay[j1][j2] = creal(Eq_Reg_d2sigmay);
//         }
//       }
//       Chebyshev_Coefficients_2D(Sol, cSol, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d1, cSol_d1, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d11, cSol_d11, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d12, cSol_d12, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2, cSol_d2, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d22, cSol_d22, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

//       Chebyshev_Coefficients_2D(Sol_dsigma, cSol_dsigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigma, cSol_d2sigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigmay, cSol_d2sigmay, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_dy, cSol_dy, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2y, cSol_d2y, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

//       Chebyshev_Coefficients_2D(Reg_dsigma, cReg_dsigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Reg_d2sigma, cReg_d2sigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Reg_dy, cReg_dy, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Reg_d2y, cReg_d2y, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Reg_d2sigmay, cReg_d2sigmay, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

  
    
//       sprintf(fn, "data/%s/ParticleSpecGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_physical, "data/%s/ParticlePhysicalGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_reg2derv, "data/%s/ParticleReg2Derv_dom_%d_fields_%d.dat",par.SimName, iDom, iF);

//       fp=fopen(fn ,"w");
//       fp_physical=fopen(fn_physical ,"w");
//       fp_reg2derv=fopen(fn_reg2derv ,"w");
      
//       fprintf(fp, "#1: chi_1\t 2:chi_2\t 3:Solution\t 4:Solution,1\t 5:Solution,11\t 6:Solution,12\t 7:Solution,2\t 8:Solution,22 \n");
//       fprintf(fp_physical, "#1: chi_1\t 2:chi_2\t 3:Solution\t 4:Solution,sigma\t 5:Solution,sigmasigma\t 6:Solution,sigmay\t 7:Solution,y\t 8:Solution,yy \n");
//       fprintf(fp_reg2derv, "#1: chi_1\t 2:chi_2\t 3:Reg d_sigma\t 4:Reg d2_sigma\t 5:Reg d_y\t 6:Reg d2_y\t 7:Reg d_2ysigma\n");
      
//       for(j1=0; j1<=j1_max; j1++){
//         double chi_1 = -1. + 2.*j1/j1_max ;
//         func_derivs_2D sigma;

//         // for(j2=0; j2<=j2_max; j2++){
//           double chi_2 = -1.;
//           func_derivs_2D y;

//           get_sigma(par, iDom, chi_1, chi_2, &sigma);
//           get_y(par, iDom, chi_1, chi_2, &y);
          
//           double out    = Clenshaw_Chebyshev_2D(cSol, N1, N2, chi_1, chi_2),
//                  out_d1 = Clenshaw_Chebyshev_2D(cSol_d1, N1, N2, chi_1, chi_2),
//                  out_d11 = Clenshaw_Chebyshev_2D(cSol_d11, N1, N2, chi_1, chi_2),
//                  out_d12 = Clenshaw_Chebyshev_2D(cSol_d12, N1, N2, chi_1, chi_2),
//                  out_d2 = Clenshaw_Chebyshev_2D(cSol_d2, N1, N2, chi_1, chi_2),
//                  out_d22 = Clenshaw_Chebyshev_2D(cSol_d22, N1, N2, chi_1, chi_2),
                 
//                  out_dsigma = Clenshaw_Chebyshev_2D(cSol_dsigma, N1, N2, chi_1, chi_2),
//                  out_d2sigma = Clenshaw_Chebyshev_2D(cSol_d2sigma, N1, N2, chi_1, chi_2),
//                  out_d2sigmay = Clenshaw_Chebyshev_2D(cSol_d2sigmay, N1, N2, chi_1, chi_2),
//                  out_dy = Clenshaw_Chebyshev_2D(cSol_dy, N1, N2, chi_1, chi_2),
//                  out_d2y = Clenshaw_Chebyshev_2D(cSol_d2y, N1, N2, chi_1, chi_2),
                 
//                  out_Reg_dsigma=Clenshaw_Chebyshev_2D(cReg_dsigma, N1, N2, chi_1, chi_2),
//                  out_Reg_d2sigma=Clenshaw_Chebyshev_2D(cReg_d2sigma, N1, N2, chi_1, chi_2),
//                  out_Reg_dy=Clenshaw_Chebyshev_2D(cReg_dy, N1, N2, chi_1, chi_2),
//                  out_Reg_d2y=Clenshaw_Chebyshev_2D(cReg_d2y, N1, N2, chi_1, chi_2),
//                  out_Reg_d2sigmay=Clenshaw_Chebyshev_2D(cReg_d2sigmay, N1, N2, chi_1, chi_2);
          
//           fprintf(fp, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_d1, out_d11, out_d12, out_d2, out_d22 );
//           fprintf(fp_physical, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_dsigma, out_d2sigma, out_d2sigmay, out_dy, out_d2y);
//           fprintf(fp_reg2derv, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out_Reg_dsigma, out_Reg_d2sigma, out_Reg_dy, out_Reg_d2y,  out_Reg_d2sigmay);
//         // }

        
//       }
//       fclose(fp);
//       fclose(fp_physical);
//       fclose(fp_reg2derv);
      
//     }
//     free_dmatrix(Sol,0, N1, 0, N2); free_dmatrix(cSol, 0,N1, 0, N2);
//     free_dmatrix(Sol_d1,0, N1, 0, N2); free_dmatrix(cSol_d1, 0,N1, 0, N2);
//     free_dmatrix(Sol_d11,0, N1, 0, N2); free_dmatrix(cSol_d11, 0,N1, 0, N2);
//     free_dmatrix(Sol_d12,0, N1, 0, N2); free_dmatrix(cSol_d12, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2,0, N1, 0, N2); free_dmatrix(cSol_d2, 0,N1, 0, N2);
//     free_dmatrix(Sol_d22,0, N1, 0, N2); free_dmatrix(cSol_d22, 0,N1, 0, N2);

//     free_dmatrix(Sol_dsigma,0, N1, 0, N2); free_dmatrix(cSol_dsigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigma,0, N1, 0, N2); free_dmatrix(cSol_d2sigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigmay,0, N1, 0, N2); free_dmatrix(cSol_d2sigmay, 0,N1, 0, N2);
//     free_dmatrix(Sol_dy,0, N1, 0, N2); free_dmatrix(cSol_dy, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2y,0, N1, 0, N2); free_dmatrix(cSol_d2y, 0,N1, 0, N2);
    
//     free_dmatrix(Reg_dsigma,0, N1, 0, N2); free_dmatrix(cReg_dsigma, 0,N1, 0, N2);
//     free_dmatrix(Reg_d2sigma,0, N1, 0, N2); free_dmatrix(cReg_d2sigma, 0,N1, 0, N2);
//     free_dmatrix(Reg_dy,0, N1, 0, N2); free_dmatrix(cReg_dy, 0,N1, 0, N2);
//     free_dmatrix(Reg_d2y,0, N1, 0, N2); free_dmatrix(cReg_d2y, 0,N1, 0, N2);
//     free_dmatrix(Reg_d2sigmay,0, N1, 0, N2); free_dmatrix(cReg_d2sigmay, 0,N1, 0, N2);
//   // }


//   free_derivs_2D(&W, par.ntotal);
//   free_dvector(X, 0, par.Ntotal);
//   return;
// }
//-------------------------------------------------------------------------------
// void output_LoadedData_SpecGrid(parameters par){
  
//   int iDom, iF, j1, j2, j1_max=50, j2_max=50, N1, N2;
//   FILE *fp, *fp_SpcCoord, *fp_physical, *fp_physical_SpcCoord;
//   char fn[200], fn_SpcCoord[200], fn_physical[200], fn_physical_SpcCoord[200];
//   double *X;
//   X=dvector(0, par.Ntotal);
//   load_X_from_ChebSolution(par, X);

//   derivs_2D W;
//   allocate_derivs_2D(&W, par.ntotal);
// 	copy_X_to_W(par, X, W);
// 	Get_Derivatives(par,  W);

//   complex_derivs_2D sol_local;
//   complex_sigma_derivs sol_local_Phys;


//   double **Sol,  **Sol_d1, **Sol_d11, **Sol_d12, **Sol_d2, **Sol_d22,
//                  **Sol_dsigma, **Sol_d2sigma, **Sol_d2sigmay, **Sol_dy, **Sol_d2y,
//          **cSol, **cSol_d1, **cSol_d11, **cSol_d12, **cSol_d2, **cSol_d22,
//                  **cSol_dsigma, **cSol_d2sigma, **cSol_d2sigmay, **cSol_dy, **cSol_d2y;
	
//   for(iDom=0; iDom<nDom; iDom++){
//     N1=par.N1[iDom];
//     N2=par.N2[iDom];
//     Sol=dmatrix(0, N1, 0, N2); cSol=dmatrix(0,N1, 0, N2);
//     Sol_d1=dmatrix(0, N1, 0, N2); cSol_d1=dmatrix(0,N1, 0, N2);
//     Sol_d11=dmatrix(0, N1, 0, N2); cSol_d11=dmatrix(0,N1, 0, N2);
//     Sol_d12=dmatrix(0, N1, 0, N2); cSol_d12=dmatrix(0,N1, 0, N2);
//     Sol_d2=dmatrix(0, N1, 0, N2); cSol_d2=dmatrix(0,N1, 0, N2);
//     Sol_d22=dmatrix(0, N1, 0, N2); cSol_d22=dmatrix(0,N1, 0, N2);

//     Sol_dsigma=dmatrix(0, N1, 0, N2); cSol_dsigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigma=dmatrix(0, N1, 0, N2); cSol_d2sigma=dmatrix(0,N1, 0, N2);
//     Sol_d2sigmay=dmatrix(0, N1, 0, N2); cSol_d2sigmay=dmatrix(0,N1, 0, N2);
//     Sol_dy=dmatrix(0, N1, 0, N2); cSol_dy=dmatrix(0,N1, 0, N2);
//     Sol_d2y=dmatrix(0, N1, 0, N2); cSol_d2y=dmatrix(0,N1, 0, N2);

//     for(iF=0; iF<nFields; iF++){
    
//       for(j1=0; j1<=N1; j1 ++) {
//         for(j2=0; j2<=N2; j2 ++){
//           int indx = Index(par, iDom, iF, j1, j2);
//           Sol[j1][j2]     = sol_local.d0  = W.d0[indx];
//           Sol_d1[j1][j2]  = sol_local.d1  = W.d1[indx];
//           Sol_d11[j1][j2] = sol_local.d11  = W.d11[indx];
//           Sol_d12[j1][j2] = sol_local.d12  = W.d12[indx];
//           Sol_d2[j1][j2]  = sol_local.d2  = W.d2[indx];
//           Sol_d22[j1][j2] = sol_local.d22  = W.d22[indx];

//           get_PhysDerv_from_SpecDerv_complex(par, iDom, j1, j2, sol_local, &sol_local_Phys);

//           Sol_dsigma[j1][j2]   = creal(sol_local_Phys.dsigma);
//           Sol_d2sigma[j1][j2]  = creal(sol_local_Phys.d2sigma);
//           Sol_d2sigmay[j1][j2] = creal(sol_local_Phys.d2sigmay);
//           Sol_dy[j1][j2]       = creal(sol_local_Phys.dy);
//           Sol_d2y[j1][j2]      = creal(sol_local_Phys.d2y);
//         }
//       }
//       Chebyshev_Coefficients_2D(Sol, cSol, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d1, cSol_d1, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d11, cSol_d11, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d12, cSol_d12, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2, cSol_d2, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d22, cSol_d22, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

//       Chebyshev_Coefficients_2D(Sol_dsigma, cSol_dsigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigma, cSol_d2sigma, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2sigmay, cSol_d2sigmay, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_dy, cSol_dy, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);
//       Chebyshev_Coefficients_2D(Sol_d2y, cSol_d2y, N1, par.grid_1[iDom], N2, par.grid_2[iDom]);

  
    
//       sprintf(fn, "data/%s/LoadedSolutionSpecGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_physical, "data/%s/LoadedSolutionPhysicalGrid_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
//       sprintf(fn_SpcCoord, "data/%s/LoadedSolutionSpecGrid_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF); 
//       sprintf(fn_physical_SpcCoord, "data/%s/LoadedSolutionPhysicalGrid_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
//       fp=fopen(fn ,"w");
//       fp_SpcCoord=fopen(fn_SpcCoord ,"w");
//       fp_physical=fopen(fn_physical ,"w");
//       fp_physical_SpcCoord=fopen(fn_physical_SpcCoord ,"w");

//       fprintf(fp, "#1: chi1\t 2:chi2\t 3:Solution\t 4:Solution,1\t 5:Solution,11\t 6:Solution,12\t 7:Solution,2\t 8:Solution,22 \n");
//       fprintf(fp_physical, "#1: sigma\t 2:y\t 3:Solution\t 4:Solution,sigma\t 5:Solution,sigmasigma\t 6:Solution,sigmay\t 7:Solution,y\t 8:Solution,yy \n");
//       fprintf(fp_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution\t 4:Solution,1\t 5:Solution,11\t 6:Solution,12\t 7:Solution,2\t 8:Solution,22 \n");
//       fprintf(fp_physical_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution\t 4:Solution,sigma\t 5:Solution,sigmasigma\t 6:Solution,sigmay\t 7:Solution,y\t 8:Solution,yy \n");

//       int d_j = 0;
//       for(j1=0; j1<=j1_max; j1++){
//         double chi_1 = -1. + 2.*j1/j1_max ;
//         func_derivs_2D sigma;

//         for(j2=0+d_j; j2<=j2_max-d_j; j2++){
//           double chi_2 = -1. + 2.*j2/j2_max;
//           func_derivs_2D y;

//           get_sigma(par, iDom, chi_1, chi_2, &sigma);
//           get_y(par, iDom, chi_1, chi_2, &y);
          
//           double out    = Clenshaw_Chebyshev_2D(cSol, N1, N2, chi_1, chi_2),
//                  out_d1 = Clenshaw_Chebyshev_2D(cSol_d1, N1, N2, chi_1, chi_2),
//                  out_d11 = Clenshaw_Chebyshev_2D(cSol_d11, N1, N2, chi_1, chi_2),
//                  out_d12 = Clenshaw_Chebyshev_2D(cSol_d12, N1, N2, chi_1, chi_2),
//                  out_d2 = Clenshaw_Chebyshev_2D(cSol_d2, N1, N2, chi_1, chi_2),
//                  out_d22 = Clenshaw_Chebyshev_2D(cSol_d22, N1, N2, chi_1, chi_2),
                 
//                  out_dsigma = Clenshaw_Chebyshev_2D(cSol_dsigma, N1, N2, chi_1, chi_2),
//                  out_d2sigma = Clenshaw_Chebyshev_2D(cSol_d2sigma, N1, N2, chi_1, chi_2),
//                  out_d2sigmay = Clenshaw_Chebyshev_2D(cSol_d2sigmay, N1, N2, chi_1, chi_2),
//                  out_dy = Clenshaw_Chebyshev_2D(cSol_dy, N1, N2, chi_1, chi_2),
//                  out_d2y = Clenshaw_Chebyshev_2D(cSol_d2y, N1, N2, chi_1, chi_2);
          
//           fprintf(fp, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out, out_d1, out_d11, out_d12, out_d2, out_d22 );
//           fprintf(fp_SpcCoord, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_d1, out_d11, out_d12, out_d2, out_d22 );

//           fprintf(fp_physical, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out, out_dsigma, out_d2sigma, out_d2sigmay, out_dy, out_d2y);
//           fprintf(fp_physical_SpcCoord, "%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out, out_dsigma, out_d2sigma, out_d2sigmay, out_dy, out_d2y);
//         }
//         fprintf(fp, "\n");
//         fprintf(fp_SpcCoord, "\n");
//         fprintf(fp_physical, "\n");
//         fprintf(fp_physical_SpcCoord, "\n");
//       }
//       fclose(fp);
//       fclose(fp_SpcCoord);
//       fclose(fp_physical);
//       fclose(fp_physical_SpcCoord);
//     }
//     free_dmatrix(Sol,0, N1, 0, N2); free_dmatrix(cSol, 0,N1, 0, N2);
//     free_dmatrix(Sol_d1,0, N1, 0, N2); free_dmatrix(cSol_d1, 0,N1, 0, N2);
//     free_dmatrix(Sol_d11,0, N1, 0, N2); free_dmatrix(cSol_d11, 0,N1, 0, N2);
//     free_dmatrix(Sol_d12,0, N1, 0, N2); free_dmatrix(cSol_d12, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2,0, N1, 0, N2); free_dmatrix(cSol_d2, 0,N1, 0, N2);
//     free_dmatrix(Sol_d22,0, N1, 0, N2); free_dmatrix(cSol_d22, 0,N1, 0, N2);

//     free_dmatrix(Sol_dsigma,0, N1, 0, N2); free_dmatrix(cSol_dsigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigma,0, N1, 0, N2); free_dmatrix(cSol_d2sigma, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2sigmay,0, N1, 0, N2); free_dmatrix(cSol_d2sigmay, 0,N1, 0, N2);
//     free_dmatrix(Sol_dy,0, N1, 0, N2); free_dmatrix(cSol_dy, 0,N1, 0, N2);
//     free_dmatrix(Sol_d2y,0, N1, 0, N2); free_dmatrix(cSol_d2y, 0,N1, 0, N2);
//   }

//   for(iF=0; iF<nPar; iF++){
//     sprintf(fn, "data/%s/LoadedSolution_parameter_%d.dat",par.SimName, iF);
//     fp = fopen(fn, "w");
//     double out = X[par.Ntotal-iF];
//     fprintf(fp, "%3.15e\t\n", out );
//     fclose(fp);
//   }

//   free_derivs_2D(&W, par.ntotal);
//   free_dvector(X, 0, par.Ntotal);
//   return;
// }

//-------------------------------------------------------------------------------
void output_F_of_X(parameters par, double *X){
  
  int iDom, iF, j1, j2, N1, N2;
  FILE *fp, *fp_SpcCoord;
  char fn[200], fn_SpcCoord[200];
  double *F=dvector(0,par.Ntotal), **Sol, **c;
  F_of_X(par, X, F);
  
  for(iDom=0; iDom<nDom; iDom++){
    N1=par.N1[iDom];
    N2=par.N2[iDom];
    Sol=dmatrix(0, N1, 0, N2); c=dmatrix(0,N1, 0, N2);

    for(iF=0; iF<nFields; iF++){

      sprintf(fn, "data/%s/F_of_X_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
      sprintf(fn_SpcCoord, "data/%s/F_of_X_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
      fp=fopen(fn ,"w");
      fp_SpcCoord=fopen(fn_SpcCoord ,"w");

      fprintf(fp, "#1: sigma\t 2:y\t 3:Solution \n");
      fprintf(fp_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution \n");
    
      for(j1=0; j1<=N1; j1++){
        double chi_1 = par.grid_chi_1[iDom][j1];//-1. + 2.*j1/j1_max ;
        func_derivs_2D sigma;

        for(j2=0; j2<=N2; j2++){
          int indx = Index(par, iDom, iF, j1, j2);
          double chi_2 = par.grid_chi_2[iDom][j2];//-1. + 2.*j2/j2_max;
          func_derivs_2D y;

          get_sigma(par, iDom, chi_1, chi_2, &sigma);
          get_y(par, iDom, chi_1, chi_2, &y);
          
          double out = F[indx];
          fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out );
          fprintf(fp_SpcCoord, "%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out );
        }
        fprintf(fp, "\n");
        fprintf(fp_SpcCoord, "\n");
      }
      fclose(fp);
      fclose(fp_SpcCoord);
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

  free_dvector(F, 0 , par.Ntotal);
  return;
}
//-------------------------------------------------------------------------------
void output_EffectiveSource(parameters par){
  int N1 = par.N1_PuncSeff,
      N2=par.N2_PuncSeff,
      idom_part = par.Dom_ptcl, i1, i2, 
      N1_out=100, //par.N1[idom_part], 
      N2_out=100; //par.N2[idom_part];
  FILE *fp;
  char fn[500];
  func_derivs_2D sigma, y;
  double x1, x2, Seff_real, Seff_imag;

  

  sprintf(fn, "data/%s/Seff_N1_%d_N2_%d.dat",par.SimName, N1, N2);
  fp = fopen(fn, "w");
  fprintf(fp, "#1:x1\t 2:x2\t 3:sigma\t 4:y\t 5:Re Seff\t 6:Im Seff\n");

  for(i1=0; i1<=N1_out; i1++){
    // x1 = get_grid( par.grid_1_PuncSeff, i1, N1 );//par.grid_chi_1[idom_part][i1];
    x1 = -1. + 2.*i1/N1_out;
    for(i2=0; i2<=N2_out; i2++){
      // x2=get_grid( par.grid_2_PuncSeff, i2, N2 );//par.grid_chi_2[idom_part][i2];
      x2 = -1. + 2.*i2/N2_out;
      
      get_sigma(par, idom_part, x1, x2, &sigma);
      get_y(par, idom_part, x1, x2, &y); 

      Seff_real=Clenshaw_Chebyshev_2D(par.Re_cheb_Seff, N1, N2, x1, x2);
      Seff_imag=Clenshaw_Chebyshev_2D(par.Im_cheb_Seff, N1, N2, x1, x2);

      fprintf(fp, "%3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\n", x1, x2, sigma.d0, y.d0, Seff_real,Seff_imag);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);


  return;

}
//----------------------------------------------------------------------------
void output_Seff_SpecCoef(parameters par, double **Seff, char *Seff_sector){
  //NOTE: OUTPUTS THE CHEBYSHEV COEFFICIENTS FOR THE EVEN AND ODD SECTORS
  
  FILE *fp;
  char fn[900];  
  int j1, j2, N1, N2;
  
  N1=par.N1_PuncSeff;
  N2=par.N2_PuncSeff;
  
  double *p, *cp, **c;
  c=dmatrix(0,N1, 0, N2);

  //Coefficients along x1 direction
  p=dvector(0, N1); cp=dvector(0, N1);
  
  for(j2=0; j2<=N2; j2++){        
        
    for(j1=0; j1<=N1; j1++){
      p[j1]=Seff[j1][j2];
    }
    Chebyshev_Coefficients(p, cp, N1, par.grid_1_PuncSeff);
    for(j1=0; j1<=N1; j1++) c[j1][j2] = cp[j1];
  }
  free_dvector(p,0,N1); free_dvector(cp, 0, N1);

  sprintf(fn, "data/%s/%s_ChebCoef_x1.dat", par.SimName, Seff_sector);
  fp=fopen(fn,"w");
  fprintf(fp, "#1:i1\t");
  for(j2=0; j2<=N2; j2++){
    double chi_2 = get_grid( par.grid_2_PuncSeff, j2, N2 );   
    fprintf(fp, "#%d: j2=%d (chi_2=%3.3e)\t", j2+2, j2, chi_2 );
  }

  for(j1=0; j1<=N1; j1++){
    fprintf(fp, "\n%d\t", j1);
    for(j2=0; j2<=N2; j2++) fprintf(fp, "%3.15e\t", c[j1][j2] );
  }
  fclose(fp);

  //Coefficients along x2 direction
  p=dvector(0,N2); cp=dvector(0,N2);
  for(j1=0; j1<=N1; j1++){

    for(j2=0; j2<=N2; j2++){
      p[j2]=Seff[j1][j2];
    }
    Chebyshev_Coefficients(p, cp, N2, par.grid_2_PuncSeff);
    for(j2=0; j2<=N2; j2++) c[j1][j2] = cp[j2];
  }
  free_dvector(p,0,N2); free_dvector(cp, 0, N2);

  sprintf(fn, "data/%s/%s_ChebCoef_x2.dat", par.SimName, Seff_sector);
  fp=fopen(fn,"w");
  fprintf(fp, "#1:i2\t");
  for(j1=0; j1<=N1; j1++){    
    double chi_1 = get_grid( par.grid_1_PuncSeff, j1, N1 );   
    fprintf(fp, "#%d: j1=%d (chi_1=%3.3e)\t", j1+2, j1,  chi_1);
  }
  for(j2=0; j2<=N2; j2++){
    fprintf(fp, "\n%d\t", j2);
    for(j1=0; j1<=N1; j1++) fprintf(fp, "%3.15e\t", c[j1][j2] );
  }
  fclose(fp);

  free_dmatrix(c, 0,N1, 0, N2);
  return;

}

//--------------------------------
void output_Puncture_at_Boundary(parameters par){
	FILE *fp;
	char fn[500];
  
  sprintf(fn, "data/%s/Cheb_PunctureField_at_Boundary.dat",par.SimName);
  fp = fopen(fn,"w");
  fprintf(fp, "#1:i1\t #2:cheb(Re@phi_punc)\t #3:cheb(Re@phi_punc,sigma)\t #4:cheb(Re@phi_punc,y)\n");
  
  int i1, N1=par.N1_PuncSeff;
  for(i1=0; i1<=N1; i1++){
    fprintf(fp, "%d\t %3.15e\t %3.15e\t %3.15e\n", i1, par.Re_cheb_phi_Punc[i1],par.Re_cheb_phi_Punc_sigma[i1],par.Re_cheb_phi_Punc_y[i1]);
  }
  fclose(fp);
  
  
  
  
  sprintf(fn, "data/%s/PunctureField_at_Boundary.dat",par.SimName);
  fp = fopen(fn,"w");
  fprintf(fp, "#1:x1\t #2:phi_punc\t #3:phi_punc,sigma\t #4:phi_punc,y\n");
  int i1_max=100;
  for(i1=0; i1<=i1_max; i1++){
    double chi_1 = -1. + 2.*i1/i1_max, phi, phi_sigma, phi_y;
    phi = Clenshaw_Chebyshev(par.Re_cheb_phi_Punc, N1, chi_1);
    phi_sigma = Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_sigma, N1, chi_1);
    phi_y = Clenshaw_Chebyshev(par.Re_cheb_phi_Punc_y, N1, chi_1);

    fprintf(fp, "%lf\t %3.15e\t %3.15e\t %3.15e\n", chi_1, phi, phi_sigma, phi_y);
  }
  fclose(fp);

	return;
}
//-------------------------------------------------------------------------------
void output_PunctureField(parameters par){
  int N1 = par.N1_PuncSeff, N1_out=100,
      N2=par.N2_PuncSeff, N2_out=100,
      idom_part = par.Dom_ptcl, i1, i2;
  FILE *fp;
  char fn[500];
  func_derivs_2D sigma, y;
  double x1, x2, PuncField_real, PuncField_imag;

  sprintf(fn, "data/%s/PunctureField_N1_%d_N2_%d.dat",par.SimName, N1, N2);
  fp = fopen(fn, "w");
  fprintf(fp, "#1:x1\t 2:x2\t 3:sigma\t 4:y\t 5:Re Seff\t 6:Im Seff\n");

  for(i1=0; i1<=N1_out; i1++){
     x1 = -1. + 2.*i1/N1_out;
    for(i2=0; i2<=N2_out; i2++){
      x2 = -1. + 2.*i2/N2_out;
      
      get_sigma(par, idom_part, x1, x2, &sigma);
      get_y(par, idom_part, x1, x2, &y); 

      PuncField_real=Clenshaw_Chebyshev_2D(par.Re_cheb_PuncField, N1, N2, x1, x2);
      PuncField_imag=Clenshaw_Chebyshev_2D(par.Im_cheb_PuncField, N1, N2, x1, x2);

      fprintf(fp, "%3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\t %3.15e\n", x1, x2, sigma.d0, y.d0, PuncField_real,PuncField_imag);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);


  return;

}


//-------------------------------------------------------------------------------
void output_SelfForce(parameters par, double *X){
  //NOTE: OUTPUT TAKES INTO ACCOUNT THE CONTRIBUTIONS FROM m AND -m
  FILE *fp;
  char fn[500];
  sprintf(fn, "data/%s/Ft_m.dat",par.SimName);
  fp = fopen(fn, "w");
  double output_bulk, output_equator;

  double complex Ft_m_bulk, Ft_m_equator;
  Get_SelfForce_Ft_m(par, X, &Ft_m_bulk, &Ft_m_equator);
  output_bulk = 2*creal(Ft_m_bulk);
  output_equator = 2*creal(Ft_m_equator);

  fprintf(fp, "#1: Re(Ft_m_bulk)\t 2: Re(Ft_m_equator)\n");
  fprintf(fp, "%3.15e\t %3.15e\n", output_bulk, output_equator);
  fclose(fp);

  sprintf(fn, "data/%s/Fr_m.dat",par.SimName);
  fp = fopen(fn, "w");

  double complex Fr_m_bulk, Fr_m_equator;
  Get_SelfForce_Fr_m(par, X, &Fr_m_bulk, &Fr_m_equator);
  output_bulk    = par.m ==0?  creal(Fr_m_bulk) : 2*creal(Fr_m_bulk);
  output_equator = par.m ==0?  creal(Fr_m_equator) : 2*creal(Fr_m_equator);

  fprintf(fp, "#1: Re(Fr_m_bulk)\t 2: Re(Fr_m_equator)\n");
  fprintf(fp, "%3.15e\t %3.15e\n", output_bulk, output_equator);
  fclose(fp);

  return;
}
//-------------------------------------------------------------------------------
void output_DifferentialOperator(parameters par, double *X){
  double *AX=dvector(0, par.Ntotal);
  A_of_X(par, X, AX);
  output_SpecCoef_AX(par, AX);
  output_Solution_AX(par, AX);

  free_dvector(AX, 0 , par.Ntotal);
  return;
}
// -------------------------------------------------------------------------------
void A_of_X(parameters par, double *X, double *F)
{
  int iDom, iFields, j1, j2, Ntotal = par.Ntotal, nTotal = Ntotal+1, N1, N2 ;
  
  double *Eq=dvector(0, nFields-1), *ExtraParameters=dvector(0,nPar), *EqPar=dvector(0,nPar);
  derivs_2D W, W_atGrid;
  allocate_derivs_2D(&W, nTotal);
  allocate_derivs_2D(&W_atGrid, nFields);
  
  copy_X_to_W(par, X, W);
  get_ExtraPar_from_X(par, X, ExtraParameters);

  Get_Derivatives(par, W);

  for(iDom=0; iDom<nDom; iDom++){
    N1 = par.N1[iDom]; 
    N2 = par.N2[iDom];
    for(j2=0; j2<=N2; j2++){
      for(j1=0; j1<=N1; j1++){
          FieldEquations_AX(par, W, ExtraParameters, iDom, j1, j2, Eq);
        
        for(iFields=0; iFields<nFields; iFields++){
          int Indx = Index(par, iDom, iFields, j1, j2);
          F[Indx]=Eq[iFields];
        }
      }
    }
  }
  


  ParameterEquations(par, W, ExtraParameters, EqPar);
  for(iFields=0; iFields<nPar; iFields++){
    F[Ntotal-iFields] = EqPar[iFields];
  }
  
        
  free_derivs_2D(&W, nTotal);
  free_derivs_2D(&W_atGrid, nFields);
  free_dvector(Eq, 0, nFields-1);
  free_dvector(ExtraParameters, 0, nPar);
  free_dvector(EqPar, 0, nPar);
  
  return;
}
//-------------------------------------------------------------------------------
void FieldEquations_AX(parameters par, derivs_2D W, double *Xpar, int iDom, int j1, int j2, double *F){
  complex_derivs_2D phi;
  double complex A_phi;
  int indx_Re_phi = Index(par, iDom, 0, j1, j2),
      indx_Im_phi = Index(par, iDom, 1, j1, j2);
      
  get_ComplexField(par, indx_Re_phi, indx_Im_phi, W , &phi);

  A_phi = Diff_Operator(par, phi, iDom, j1, j2, 0); 
  
  //Test_Effective_Source(par, iDom, j1, j2);
  //Get_ModeSum_Source(par, iDom, j1, j2);

  F[0] = creal(A_phi);
  F[1] = cimag(A_phi); 

  return;
}
//-------------------------------------------------------------------------------
void output_SpecCoef_AX(parameters par, double *X){
  
  FILE *fp;
  char fn[200];  
  int iDom, iField, j1, j2, N1, N2;
  double *p, *cp, **c;

  for(iDom=0; iDom<nDom; iDom++){
    N1=par.N1[iDom];
    N2=par.N2[iDom];
    c=dmatrix(0,N1, 0, N2);
    for(iField=0; iField<nFields; iField++){
      p=dvector(0, N1); cp=dvector(0, N1);

      for(j2=0; j2<=N2; j2++){        
        for(j1=0; j1<=N1; j1++){
          int indx = Index(par, iDom, iField,j1,j2);
          p[j1]=X[indx];
        }
        Chebyshev_Coefficients(p, cp, N1, par.grid_1[iDom]);
        for(j1=0; j1<=N1; j1++) c[j1][j2] = cp[j1];
      }
      free_dvector(p,0,N1); free_dvector(cp, 0, N1);

      sprintf(fn, "data/%s/DiffOp_ChebCoef_x1_dom_%d_field_%d.dat", par.SimName, iDom, iField);
      fp=fopen(fn,"w");
      fprintf(fp, "#1:i1\t");
      for(j2=0; j2<=N2; j2++){
        fprintf(fp, "#%d: j2=%d (chi_2=%3.3e)\t", j2+2, j2, par.grid_chi_2[iDom][j2] );
      }

      for(j1=0; j1<=N1; j1++){
        fprintf(fp, "\n%d\t", j1);
        for(j2=0; j2<=N2; j2++) fprintf(fp, "%3.15e\t", c[j1][j2] );
      }
      fclose(fp);

      p=dvector(0,N2); cp=dvector(0,N2);
      for(j1=0; j1<=N1; j1++){
        for(j2=0; j2<=N2; j2++){
          int indx = Index(par, iDom, iField,j1,j2);
          p[j2]=X[indx];
        }
        Chebyshev_Coefficients(p, cp, N2, par.grid_2[iDom]);
        for(j2=0; j2<=N2; j2++) c[j1][j2] = cp[j2];
      }
      free_dvector(p,0,N2); free_dvector(cp, 0, N2);

      sprintf(fn, "data/%s/DiffOp_ChebCoef_x2_dom_%d_field_%d.dat", par.SimName, iDom, iField);
      fp=fopen(fn,"w");
      fprintf(fp, "#1:i2\t");
      for(j1=0; j1<=N1; j1++){       
        fprintf(fp, "#%d: j1=%d (chi_1=%3.3e)\t", j1+2, j1,  par.grid_chi_1[iDom][j1]);
      }

      for(j2=0; j2<=N2; j2++){
        fprintf(fp, "\n%d\t", j2);
        for(j1=0; j1<=N1; j1++) fprintf(fp, "%3.15e\t", c[j1][j2] );
      }
      fclose(fp);
    }
    free_dmatrix(c, 0,N1, 0, N2);
  }

  return;
}
//-------------------------------------------------------------------------------
void output_Solution_AX(parameters par, double *X){
  
  int iDom, iF, j1, j2, j1_max=50, j2_max=50, N1, N2;
  FILE *fp, *fp_SpcCoord;
  char fn[200], fn_SpcCoord[200];
  double **Sol, **c;
  
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

  
      sprintf(fn, "data/%s/DiffOp_Solution_dom_%d_fields_%d.dat",par.SimName, iDom, iF);
      sprintf(fn_SpcCoord, "data/%s/DiffOp_Solution_chi_dom_%d_fields_%d.dat",par.SimName, iDom, iF);    
      fp=fopen(fn ,"w");
      fp_SpcCoord=fopen(fn_SpcCoord ,"w");

      fprintf(fp, "#1: sigma\t 2:y\t 3:Solution \n");
      fprintf(fp_SpcCoord, "#1: chi1\t 2:chi2\t 3:Solution \n");


      j1_max=N1;
      j2_max=N2;
      for(j1=0; j1<=j1_max; j1++){
        // double chi_1 = -1. + 2.*j1/j1_max ;
        double chi_1 = par.grid_chi_1[iDom][j1];;
        func_derivs_2D sigma;

        for(j2=0; j2<=j2_max; j2++){
          // double chi_2 = -1. + 2.*j2/j2_max;
          double chi_2 = par.grid_chi_2[iDom][j2];
          func_derivs_2D y;

          get_sigma(par, iDom, chi_1, chi_2, &sigma);
          get_y(par, iDom, chi_1, chi_2, &y);
          
          double out = Clenshaw_Chebyshev_2D(c, N1, N2, chi_1, chi_2);
          fprintf(fp, "%3.15e\t%3.15e\t%3.15e\n", sigma.d0, y.d0, out );
          fprintf(fp_SpcCoord, "%3.15e\t%3.15e\t%3.15e\n", chi_1, chi_2, out );
        }
        fprintf(fp, "\n");
        fprintf(fp_SpcCoord, "\n");
      }
      fclose(fp);
      fclose(fp_SpcCoord);
    }
    free_dmatrix(Sol,0, N1, 0, N2);
    free_dmatrix(c, 0,N1, 0, N2);
  }

  for(iF=0; iF<nPar; iF++){
    sprintf(fn, "data/%s/DiffOp_Solution_parameter_%d.dat",par.SimName, iF);
    fp = fopen(fn, "w");
    double out = X[par.Ntotal-iF];
    fprintf(fp, "%3.15e\t\n", out );
    fclose(fp);
  }

  return;
}
//-------------------------------------------------------------------------------
// void output_convergence(parameters par, double **c, double **c_Ref){
//   int iDom, iFields;




//   return;
// }


