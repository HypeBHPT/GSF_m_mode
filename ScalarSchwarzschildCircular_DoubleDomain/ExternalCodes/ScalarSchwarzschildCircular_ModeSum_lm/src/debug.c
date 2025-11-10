#include "Solve_ODE.h"

void PrintVector(double *J, int na, int nb){
 int i;
 
 for(i=na; i<=nb; i++){
   printf("vec[%d]=%3.15e\n",i, J[i]);
//    pause();
 }
}
//------------------------------------------------------------------------------
void PrintMatrix(double **J, int la, int lb, int ca, int cb){
 int i, j;
 
 for(i=la; i<=lb; i++){
   for(j=ca; j<=cb; j++){
     printf("%3.2e ", J[i][j]);
   }
   printf("\n");
 }
}
//------------------------------------------------------------------------------
void pause(){
char pause;
printf("\nPress any key to continue\n(or 'q' to quit)\n");
// __fpurge(stdin);
fpurge(stdin);
int ret;
ret = scanf("%c", &pause);
if (ret == EOF)
{
  printf("...Error in 'Pause'\n");
  exit(1);
}
if(pause=='q')
	exit(1);
// __fpurge(stdin);
 fpurge(stdin);
printf("\n");
}//------------------------------------------------------------------------------
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
