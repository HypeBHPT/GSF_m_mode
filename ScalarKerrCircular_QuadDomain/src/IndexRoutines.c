#include "Solve_PDE.h"

// -------------------------------------------------------------------------------
int Index_Boundary(char grid[], int i, int n)
{
	int N = n -1, i1 = i;
	if(strcmp( grid,"Fourier" ) ==0){
	  if (i1 < 0)  i1 = i1 + n;
	  if (i1 > N) i1 = i1 - n;
	}
	else if(strcmp( grid,"Lobatto")==0){
	  if (i1 < 0)  i1 = -i1;
	  if (i1 > N) i1 = 2*N-i1;
	}
	else if(strcmp( grid,"Gauss")==0){
	  if (i1 < 0)  i1 = -(i1+1);
	  if (i1 > N) i1 = 2*N+1-i1;
	}
	else if(strcmp( grid,"Radau_RHS" ) ==0){
	  if (i1 < 0)  i1 = -i1;
	  if (i1 > N) i1 = 2*N+1-i1;
	}
	else if(strcmp( grid,"Radau_LHS" ) ==0){
	  if (i1 < 0)  i1 = -i1;
	  if (i1 > N) i1 = 2*N+1-i1;
	}
	else{
	  printf("Error in Index_Boundary: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n", grid);
	  exit(1);
	}
	
	return i1;
}
// -------------------------------------------------------------------------------
void Get_Index_Arrays(parameters *par)
{
	int iDom;

	(*par).n1n2[0] = 0;
	for(iDom=0; iDom<nDom; iDom++){
		(*par).n1[iDom] = total_grid_points((*par).grid_1[iDom], (*par).N1[iDom]);
		(*par).n2[iDom] = total_grid_points((*par).grid_2[iDom], (*par).N2[iDom]);
		(*par).n1n2[iDom+1] = (*par).n1n2[iDom] + (*par).n1[iDom]*(*par).n2[iDom];
	}

	return;
}
//---------------------------------------------------------------
int Index(parameters par, int iDom, int iField, int j1, int j2){
  
  int n1=par.n1[iDom], n2=par.n2[iDom], n1n2=par.n1n2[iDom],
  		i1= Index_Boundary(par.grid_1[iDom], j1, n1),
  		i2= Index_Boundary(par.grid_2[iDom], j2, n2);

  		int II;
  		// if (n1<=n2) 
  			II = i1 +i2*n1 + n1n2;
  		// else  II = i2 + i1*n2 + n1n2;

  return iField + nFields*II;
}
//-------------------------------------------------------------------------------------------
void Get_Indices_From_Index(parameters par, int indx, int *idom, int *If, int *i, int *j)
{
  
	int na, nb, n;
  *idom=nDom;
  while(indx < nFields*par.n1n2[*idom])	*idom -=1;  

  if(*idom<nDom){
  	int n1=par.n1[*idom];//,n2=par.n2[*idom];
  	na = indx - nFields*par.n1n2[*idom];
  	
  	// if(n1> n2){  		
  	// 	n=n2;
  	// 	*i = na/(nFields*n);
  	// 	nb = na - nFields*n*(*i);
  	// 	*j = nb/nFields;
  	// 	*If = nb - nFields*(*j);
  	// }
  	// else{
  		n=n1;
  		*j    = na/(nFields*n);
  		nb    = na - nFields*n*(*j);
  		*i    = nb/(nFields);
  		*If = nb - nFields*(*i);
  	// }


  }

 

  return;      	
}