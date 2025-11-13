#include "Solve_PDE.h"

void PrintVector(parameters par, double *J, int na, int nb){
 int i;
 
 for(i=na; i<=nb; i++){     
     double out =  fabs(J[i ])<1.e-11? 0. : J[i ]; 
     printf("vec[%d]=%3.3e\n",i, out);
     // if(fabs( J[i])>5.0e-1 /*isnan(J[i])!=0*/ ){
            // int j1_grid, j2_grid, iF, iDom;
     //    Get_Indices_From_Index(par, i, &iDom, &iF, &j1_grid, &j2_grid);
     //     printf("\t\tvec[%d]=%3.15e  (iF = %d, j1 = %d, j2 = %d)\n",i, J[i], iF, j1_grid, j2_grid);
       
     //   pause();
     // }
 }
}
//------------------------------------------------------------------------------
void PrintMatrix(parameters par, char *fn, double **J, int la, int lb, int ca, int cb){
 //static int count = 0;
 int i, j; 
 
 FILE *fp = fopen(fn, "w");

 
 for(i=la; i<=lb; i++){
    // int j1_grid, j2_grid, iF, iDom;
        // Get_Indices_From_Index(par, i, &iDom, &iF, &j1_grid, &j2_grid);
        // printf("I=%d, idom=%d, iF=%d, j1=%d, j2=%d\n", i, iDom, iF, j1_grid, j2_grid);
   for(j=ca; j<=cb; j++){
     double out = fabs(J[i][j])<1.0e-12? 0.:J[i][j];
//      printf(fp,"%3.5e ", out);
     fprintf(fp,"%3.2lf ", out);
       
     // if(fabs( J[i][j])>5.0e-8 /*isnan(J[i])!=0*/ ){

        
     //    fprintf(fp, "\t\tJ[%d][%d]=%3.15e : Equation (iDom = %d, iF = %d, j1 = %d, j2 = %d)\n", i, j, out, iDom, iF, j1_grid, j2_grid);
       
     //   // pause();
     // }    
     
   }
//    printf("\n");
   fprintf(fp, "\n");
 }
 fclose(fp);
  // count++;
}
//------------------------------------------------------------------------------
void pause() {

    printf("\nPress any key to continue\n(or 'q' to quit)\n");

    int pause = getchar(); // Read a character from input

    if (pause == EOF) {
        printf("Error in 'Pause'\n");
        return;
    }

    if (pause == 'q' || pause == 'Q') {
        printf("Exiting...\n");
        exit(0); // Exit the program
    }

    // int c;
    // while ((c = getchar()) != '\n' && c != EOF); // Consume characters until newline or EOF

    printf("\n");
    return;
}

// //------------------------------------------------------------------------------
// void pause(){
// char pause;
// printf("\nPress any key to continue\n(or 'q' to quit)\n");
//   __fpurge(stdin);  //FOR LINUX
// //   fpurge(stdin);  //FOR MAC		
// int ret;
// ret = scanf("%c", &pause);
// if (ret == EOF)
// {
//   printf("...Error in 'Pause'\n");
//   exit(1);
// }
// if(pause=='q')
// 	exit(1);
//   __fpurge(stdin); //FOR LINUX
// //  fpurge(stdin); //FOR MAC
// printf("\n");
// }
//------------------------------------------------------------------------------
void create_directory(char *dir_name){
  DIR *dir = opendir(dir_name);
  
  if (dir) {/* Directory exists. */    
    closedir(dir);
  }
  else if (ENOENT == errno){/* Directory does not exist. */
    char mkdir[500];
    sprintf(mkdir, "mkdir -p %s", dir_name) ;
    int func_out =  system(mkdir);
    if(func_out){ fprintf(stderr, "Error in create_directory: Unable to create directory %s\n", dir_name); exit(1);}
    
  }
  else {
    fprintf(stderr, "Error in create_directory:\n");
    exit(-1);
  }
  return;
}
//------------------------------------------------------------------------------
void Check_OS(){
  struct utsname unameData;
  uname(&unameData);

  // Check if the operating system is Mac (Darwin)
  if (strcmp(unameData.sysname, "Darwin") == 0) {
      printf("Running on Mac\n\n");
  }
  // Check if the operating system is Linux
  else if (strcmp(unameData.sysname, "Linux") == 0) {
      printf("Running on Linux\n\n");
  }
  // If it's neither Mac nor Linux
  else {
      printf("Unknown OS\n\n");
  }

  return;
}
