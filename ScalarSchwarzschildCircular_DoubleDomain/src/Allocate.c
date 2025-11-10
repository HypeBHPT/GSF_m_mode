#include "Solve_PDE.h"

// -------------------------------------------------------------------------------
void allocate_derivs_2D(derivs_2D *v, int n)
{	
	(*v).d0 = dvector(0, n-1);
	(*v).d1 = dvector(0, n-1);
	(*v).d2 = dvector(0, n-1); 
	
	(*v).d11= dvector(0, n-1);	
	(*v).d22= dvector(0, n-1); 
	(*v).d12= dvector(0, n-1);

	(*v).d111= dvector(0, n-1);	
	(*v).d112= dvector(0, n-1);	
	(*v).d122= dvector(0, n-1);
	(*v).d222= dvector(0, n-1); 

	(*v).d1111= dvector(0, n-1);	
	(*v).d1112= dvector(0, n-1);
	(*v).d1122= dvector(0, n-1);
	(*v).d1222= dvector(0, n-1);
	(*v).d2222= dvector(0, n-1); 

	(*v).d11222= dvector(0, n-1);
	
}
// -------------------------------------------------------------------------------
void fill0_derivs_2D(derivs_2D v, int n)
{	
	fill0_dvector(v.d0,  0, n);
	fill0_dvector(v.d1,  0, n);
	fill0_dvector(v.d2,  0, n);
	
	fill0_dvector(v.d11, 0, n);	 
	fill0_dvector(v.d22, 0, n); 
	fill0_dvector(v.d12, 0, n); 

	fill0_dvector(v.d111, 0, n);	 
	fill0_dvector(v.d112, 0, n); 
	fill0_dvector(v.d122, 0, n);
	fill0_dvector(v.d222, 0, n); 

	fill0_dvector(v.d1111, 0, n);	 
	fill0_dvector(v.d1112, 0, n); 
	fill0_dvector(v.d1122, 0, n);
	fill0_dvector(v.d1222, 0, n);
	fill0_dvector(v.d2222, 0, n); 

	fill0_dvector(v.d11222, 0, n);
}
// -------------------------------------------------------------------------------
void free_derivs_2D(derivs_2D *v, int n)
{	
	free_dvector((*v).d0,  0, n-1);
	
	free_dvector((*v).d1,  0, n-1);
	free_dvector((*v).d2,  0, n-1); 

	free_dvector((*v).d11, 0, n-1);
	free_dvector((*v).d12, 0, n-1); 
	free_dvector((*v).d22, 0, n-1); 

	free_dvector((*v).d111, 0, n-1);
	free_dvector((*v).d112, 0, n-1); 
	free_dvector((*v).d122, 0, n-1);
	free_dvector((*v).d222, 0, n-1);

	free_dvector((*v).d1111, 0, n-1);
	free_dvector((*v).d1112, 0, n-1); 
	free_dvector((*v).d1122, 0, n-1); 
	free_dvector((*v).d1222, 0, n-1);
	free_dvector((*v).d2222, 0, n-1);

	free_dvector((*v).d11222, 0, n-1);
}
