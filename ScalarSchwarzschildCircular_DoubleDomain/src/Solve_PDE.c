#include "Solve_PDE.h"

int main() 
{
	Check_OS();
	
	printf("\n");fflush(0);
	printf("***2D Hyperboloidal Frequency Domain Self Force (10 Dez 21)  *** \n");fflush(0);
	printf("---------------------------------------------- \n");fflush(0);
	printf("\n");fflush(0);	
	
	
	double  start_time, final_time;	
	double r0 = 10.; //Particle's orbital radius in units of M

	start_time = omp_get_wtime();
	int n_omp = omp_get_max_threads();

	int nbar = 2; //Order of Puncture Scheme
	int N = 40; //Order of numerical resolution (assuming same N in all domains and all directions)

	
	int	m, m_min =2 , m_max =2, delta_m=1;			
	
	#pragma omp parallel for 
	for(m=m_min; m<=m_max; m+=delta_m){		
				
	
				double  *X;
				// //Read Parameters----------------------------------
				parameters par;
				set_parameters(&par, N, nbar, m, r0);
				par.fout = stdout;
				
				//------------------------------------------------- 
				printf("Solving r0_over_M = %lf, m = %d, nbar = %d on thread = %d (%d)\n", par.r0_over_M, par.m, par.nbar, par.i_omp, n_omp ); 
				fprintf(par.fout, "m = %d\t nbar_max =%d\t Data for N=%d\n", m, nbar, N);

				//Load External Data-----------------------------
				load_Retarded_at_Boundary(&par);
				// output_RetardedField_Boundary(par);
				// output_Cheb_RetardedField_Boundary(par);
				// output_Legendre_Retarded_at_Boundary(par);
				
				
				load_Puncture_at_Boundary(&par);
				load_EffectiveSource(&par);
				
				// load_PunctureField(&par);								
				// output_PunctureField(par);
				// output_Puncture_at_Boundary(par);	
				// exit(-1);			
				//----------------------------------------------
				

				//Solve individual angular mode-----------------------   
				X  = dvector(0, par.Ntotal);
				get_InitialGuess(par, X); //initialise solution field
				//------------
				
				double  start_time_solver, final_time_solver;
				start_time_solver = clock(); //Start measuring time	
				int Newton_iter = solve_equations(par, X);
				final_time_solver = clock(); //Stop measuring time					
				fprintf(par.fout, "Time= %3.3e s (Newton Iteration = %d)\n\n", (final_time_solver-start_time_solver)/CLOCKS_PER_SEC, Newton_iter);
				//-----------------------------------------------------------------------------------------------	

				//OUTPUT SOLUTION-------------------------------------------------------------------------------	
				output_Solution(par, X);
				output_SpecCoef(par, X);	
				output_SelfForce(par, X);				
				
				// output_SolutionChebyshevPostProc(par, X);
				// output_Toy2ndSource(par, X, N, N_min, N_max, delta_N, m, m_min, m_max);
				// output_FieldParticle(par, X);
				// output_SolutionDerivatives(par, X);
				// output_DifferentialOperator(par, X);
				// output_F_of_X(par, X);
				
				
				

				
				
				//-----------------------------------------------------------------------------------------------
				
				free_dvector(X, 0, par.Ntotal);	
				free_grid(&par);
				free_external_data(&par);
				// fclose(par.fout);
				
				
		
		
	}
	final_time = omp_get_wtime();
	printf("number of threads = %d, Total time = %3.5e s\n", n_omp, final_time-start_time);


	return 1;
}
