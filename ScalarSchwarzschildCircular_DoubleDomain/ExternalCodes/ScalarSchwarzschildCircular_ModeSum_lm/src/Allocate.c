#include "Solve_ODE.h"

// -------------------------------------------------------------------------------
void allocate_derivs(derivs *v, int n)
{	
	(*v).d0 = dvector(0, n-1);
	(*v).d1 = dvector(0, n-1);
	(*v).d11= dvector(0, n-1);
	
}
// -------------------------------------------------------------------------------
void fill0_derivs(derivs v, int n)
{	
	fill0_dvector(v.d0,  0, n);
	fill0_dvector(v.d1,  0, n);
	fill0_dvector(v.d11, 0, n);
	
}
// -------------------------------------------------------------------------------
void free_derivs(derivs *v, int n)
{	
	free_dvector((*v).d0,  0, n-1);
	free_dvector((*v).d1,  0, n-1);
	free_dvector((*v).d11, 0, n-1);
	
}