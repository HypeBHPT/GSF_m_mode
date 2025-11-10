#include "Solve_ODE.h"

// -------------------------------------------------------------------------------
int Index_i(parameters par, int iDom, int i)
{
// 	int NA = par.N[iDom], i1 = i;
// 	
// 	if (i1 < 0)  i1 = -i1;
// 	if (i1 > NA) i1 = 2*NA-i1;
// 	
	return i; //i1;
}
// -------------------------------------------------------------------------------
int Index(parameters par, int iDom, int iField, int i)
{	// mA[iDom] = Sum[ nA[jDom] , {jDom, 0, iDom-1} ]; mA[0] = 0

	int mA = par.mA[iDom-1], ii = Index_i(par, iDom, i);

	return  iField + nFields*(mA + ii);
}
// -------------------------------------------------------------------------------
void Get_Index_Arrays(parameters *par)
{
	int iDom;

	(*par).mA[0] = 0;
	for(iDom=1; iDom<=nDom; iDom++)
		(*par).mA[iDom] = (*par).mA[iDom-1]+((*par).N[iDom-1]+1);


}