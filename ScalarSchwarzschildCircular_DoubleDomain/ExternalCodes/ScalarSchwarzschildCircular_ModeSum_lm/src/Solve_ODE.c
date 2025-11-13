#include "Solve_ODE.h"

//---------------------------------------------------------------
int main() 
{
	printf("\n");
	printf("*** Hyperboloidal Scalar Self-force: Radial orbit Schwarzschild - 14 Jun 21  *** \n");
	printf("--------------------------------------- \n");
	printf("\n");


	int n_omp = omp_get_max_threads();
	int m, 
		m_max = 2, 
		N = 150, 
		l_max = 60;
	
	double r0_over_M = 10.;	
	

	double start_time, final_time;

	start_time = omp_get_wtime();

	
	for(m=2; m<=m_max; m++){
	
	
		int l, l_min=m;
		#pragma omp parallel for
		for(l=l_min; l<=l_max; l+=2){			
			double *X;
			parameters par;	

			//Read Parameters-----------------------
			set_parameters(&par, r0_over_M, l, m, N);
			// par.fout = stdout;
			//--------------------------------------   
			printf("Solving eta = %lf, ell = %d, m = %d on thread = %d (%d)\n", par.eta, par.ell, par.m, par.i_omp, n_omp ); 
			fprintf(par.fout, "Solving eta = %lf, ell = %d, m = %d on thread = %d (%d)\n", par.eta, par.ell, par.m, par.i_omp, n_omp ); 
			
			//Solve individual angular mode-----------------------	    
			X = dvector(0, par.Ntotal);	
			Initial_Guess(par, X);

			int Newton_iter = newton(par, X);
			fprintf(par.fout, "Newton exit on thread = %d after %d iterations\n\n", par.i_omp, Newton_iter); 
			//--------------------------------------  

			//Output Solution------------------------------
			// output_cheb(par, X);
			// output_Solution(par, X);
			output_Ret_Field_Boundary_data(par, X);
			//----------------------------------------------

			free_dvector(X, 0, par.Ntotal);
			fclose(par.fout);
		}
	}
	final_time = omp_get_wtime();
	printf("number of threads = %d, Total time = %3.15e s\n", n_omp, final_time-start_time);
	 	          
    return 0;
}

