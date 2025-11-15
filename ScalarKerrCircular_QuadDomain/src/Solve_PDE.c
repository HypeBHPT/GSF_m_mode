#include "Solve_PDE.h"

int main() 
{
	Check_OS();
	
	printf("\n");fflush(0);
	printf("***2D Hyperboloidal Frequency Domain Self Force (10 Dez 21)  *** \n");fflush(0);
	printf("---------------------------------------------- \n");fflush(0);
	printf("\n");fflush(0);	
	
	
	double start_time, final_time;	
	double r0_over_M_Schwarzschild = 10, 
		   M_Omega0 = pow(r0_over_M_Schwarzschild, -3./2), 
		   a_over_M = 0.;

	start_time = omp_get_wtime();
	int n_omp = omp_get_max_threads();

	int nbar = 2; //Order of Puncture Scheme
	int N = 16; //Order of numerical resolution (assuming same N in all domains and all directions)
	int	m=2; //Azimutal Mode
	
	
	// int m_min =2 , m_max =2, delta_m=1;	
	// double a_min, a_max, delta_a;
	// int i_a, N_a=10;

	// a_min = 0.; 
	// a_max=1.; 
	// delta_a= (a_max - a_min)/N_a;
	
	// int N_min=10, N_max=32, delta_N=2;
	
	// for(N=N_max; N>=N_min; N-=delta_N){
		
	// 	#pragma omp parallel for 
	// 	for(i_a=0; i_a<=N_a; i_a++){
		parameters par;
	// 	a_over_M = a_min + i_a*delta_a;
			

			// #pragma omp parallel for 
			// for(m=m_min; m<=m_max; m+=delta_m){
						double  *X;
						// //Read Parameters----------------------------------
						
						set_parameters(&par, N, nbar, m, M_Omega0, a_over_M);

						//------------------------------------------------
						double dr=0.1, dtheta=0.1, Punc[2], Seff[2], dPunc[8],d2Punc[20];
						int Nr = 101, Ntheta=101, ir, itheta;

						FILE *fp=fopen("Test_EffSource.dat", "w");
						fprintf(fp, "# r/M \t theta/pi \t Re(Seff) \t Im(Seff) \n");
						

						for(ir=0;ir<=Nr; ir++){
							double r = par.r0_over_M+(1 -2.*ir/Nr) * dr;  
							for(itheta=0; itheta<=Ntheta; itheta++){
								double theta = Pih+(1 -2.*itheta/Ntheta) * dtheta;
								struct coordinate xBL_coord = {r, theta, 0., 0.};
								effsource_calc_m(m, &xBL_coord, Punc, dPunc, d2Punc, Seff);
								fprintf(fp, "%3.15e %3.15e %3.15e %3.15e \n", r-par.r0_over_M, theta-Pih, Seff[0], Seff[1] );  

							}
							fprintf(fp, "\n");

						}
						fclose(fp);
						exit(-1);
						//------------------------------------------------


						
						// par.fout = stdout;
						
						//------------------------------------------------- 
						printf("Solving a_over_M = %lf, M_Omega0 = %lf, r0_over_M = %lf, m = %d, nbar = %d, N=%d, on thread = %d (%d)\n", par.a_over_M, par.M_Omega0, par.r0_over_M, par.m, par.nbar, N, par.i_omp, n_omp ); 
						fprintf(par.fout, "a_over_M = %lf, m = %d\t nbar_max =%d\t Data for N=%d\n", par.a_over_M, m, nbar, N);

						
					
						if(par.TEST_Func_FLAG==0){
							get_Puncture_EffectiveSource(&par);
							
							// load_EffectiveSource(&par);
							// load_Puncture_at_Boundary(&par);
							// load_PunctureField(&par);								
							// output_PunctureField(par);
							output_Puncture_at_Boundary(par);	
							output_EffectiveSource(par);
							
							exit(-1);		
						}							
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
						// -----------------------------------------------------------------------------------------------	

						//OUTPUT SOLUTION-------------------------------------------------------------------------------	
						output_Solution(par, X);
						output_SpecCoef(par, X);	
						output_SelfForce(par, X);
						if(par.TEST_Func_FLAG==1){
							output_test_function(par);
							output_test_function_error(par, X);
						}	
						
						// output_SolutionChebyshevPostProc(par, X);
						// output_Toy2ndSource(par, X, N, N_min, N_max, delta_N, m, m_min, m_max);
						// output_FieldParticle(par, X);
						// output_SolutionDerivatives(par, X);
						// output_DifferentialOperator(par, X);
						// output_F_of_X(par, X);
						
						
						

						
						
						//-----------------------------------------------------------------------------------------------
						
						free_dvector(X, 0, par.Ntotal);	
						free_grid(&par);
						if(par.TEST_Func_FLAG==0) free_external_data(&par);
			// }
			fclose(par.fout);
		// }		
	// }
	final_time = omp_get_wtime();
	printf("number of threads = %d, Total time = %3.5e s\n", n_omp, final_time-start_time);


	return 1;
}
