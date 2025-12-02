#include "Solve_PDE.h"

void read_header(parameters par, FILE *fr, int *N1, int *N2){

	double r0_read, eta_read, small = 1.e-15;
	int m_read, nbar_read, N1_read, N2_read;
	fscanf(fr, "r0\t%lf\n", &r0_read);
	fscanf(fr, "m\t%d\n", &m_read);
	fscanf(fr, "eta\t%lf\n", &eta_read);
	fscanf(fr, "nbarmax\t%d\n\n", &nbar_read);

	fscanf(fr, "N1\t%d\n", &N1_read);
	fscanf(fr, "N2\t%d\n", &N2_read);
	// 

	if( (r0_read-par.r0_over_M)!=0 || (m_read-par.m)!=0 || fabs(1.- eta_read/par.eta) > small || (nbar_read-par.nbar) !=0.){
		fprintf(stderr,"Error in read_header: input parameters do not match simulation\n");
		fprintf(stderr,"r0 = %3.20e (%3.20e)\n", r0_read, par.r0_over_M );
		fprintf(stderr,"m = %d (%d)\n", m_read, par.m );
		fprintf(stderr,"eta = %3.15e (%3.15e)\n", eta_read, par.eta );
		fprintf(stderr,"nbar = %d (%d)\n", nbar_read, par.nbar );
		exit(-1);
	}



	*N1 = N1_read;
	*N2 = N2_read;
	return;
}
//--------------------------------
void load_Puncture_at_Boundary(parameters *par){
	FILE *fr;

	char fn[500];
	sprintf(fn,"InputData/Coord%d/r0_over_M_%.5lf/eta_%3.5lf/m_%d/nbarMax_%d/N_%d/Prec%3.5f/phi_punc_boundary_nbar%d.dat", (*par).CoordMap_FLAG,(*par).r0_over_M, (*par).eta, (*par).m, (*par).nbar,  (*par).N1_PuncSeff, (*par).prec,(*par).nbar);

	if( (fr = fopen(fn,"r")) == NULL){
		fprintf(stderr, "Error in load_Puncture_at_Boundary: file %s not found\n", fn);
		// return;
		exit(-1);
	}

	int i1, N1, N2;
	read_header(*par, fr, &N1, &N2);

	(*par).N1_LoadPunc = N1;

	(*par).Re_cheb_phi_Punc = dvector(0, N1);
	(*par).Im_cheb_phi_Punc = dvector(0, N1);

	(*par).Re_cheb_phi_Punc_sigma = dvector(0, N1);
	(*par).Im_cheb_phi_Punc_sigma = dvector(0, N1);

	(*par).Re_cheb_phi_Punc_y = dvector(0, N1);
	(*par).Im_cheb_phi_Punc_y = dvector(0, N1);

	fscanf(fr,"#1:i1	 #2:Real(Cheb_PhiP)	 #3:Imag(Cheb_PhiP)	 #4:Real(Cheb_PhiP,sigma)	 #5:Imag(Cheb_PhiP,sigma)	 #6:Real(Cheb_PhiP,y)	 #7:Imag(Cheb_PhiP,y)\n");

	for(i1=0; i1<=N1; i1++){
		int i1_read;
		double Re_c_phi, Re_c_phi_sigma, Re_c_phi_y,
			   Im_c_phi, Im_c_phi_sigma, Im_c_phi_y;

		fscanf(fr,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &i1_read, &Re_c_phi, &Im_c_phi, &Re_c_phi_sigma, &Im_c_phi_sigma, &Re_c_phi_y, &Im_c_phi_y);
		(*par).Re_cheb_phi_Punc[i1] = Re_c_phi; 
		(*par).Im_cheb_phi_Punc[i1] = Im_c_phi;

		(*par).Re_cheb_phi_Punc_sigma[i1] = Re_c_phi_sigma;
		(*par).Im_cheb_phi_Punc_sigma[i1] = Im_c_phi_sigma;

		(*par).Re_cheb_phi_Punc_y[i1] = Re_c_phi_y; 
		(*par).Im_cheb_phi_Punc_y[i1] = Im_c_phi_y;
	}
	fclose(fr);

	return;
}
//--------------------------------
void load_EffectiveSource(parameters *par){
	FILE *fr;

	char fn[500];
	sprintf(fn,"InputData/Coord%d/r0_over_M_%.5lf/eta_%3.5lf/m_%d/nbarMax_%d/N_%d/Prec%3.5f/S_eff_nbar%d.dat", (*par).CoordMap_FLAG, (*par).r0_over_M, (*par).eta, (*par).m, (*par).nbar, (*par).N1_PuncSeff, (*par).prec,(*par).nbar);

	if( (fr = fopen(fn,"r")) == NULL){
		fprintf(stderr, "Error in load_EffectiveSource: file %s not found\n", fn);
		// return;
		exit(-1);
	}

	int i1, i2, N1, N2;
	read_header(*par, fr, &N1, &N2);

	(*par).N1_LoadSeff = N1;
	(*par).N2_LoadSeff = N2;

	(*par).Re_cheb_Seff = dmatrix(0, N1, 0 , N2);
	(*par).Im_cheb_Seff = dmatrix(0, N1, 0 , N2);

	fscanf(fr,"#1:i1\t #2:i2\t #3:Real(Cheb S_eff)\t #4:Imag(Cheb Seff)\n");

	for(i1=0; i1<=N1; i1++){
		for(i2=0; i2<=N2; i2++){
			int i1_read, i2_read;
			double Re_c_Seff, Im_c_Seff;

			fscanf(fr,"%d\t%d\t%lf\t%lf\n", &i1_read, &i2_read, &Re_c_Seff, &Im_c_Seff);
			(*par).Re_cheb_Seff[i1][i2] = Re_c_Seff;
			(*par).Im_cheb_Seff[i1][i2] = Im_c_Seff;
		}
		fscanf(fr,"\n");		
	}
	

	fclose(fr);
	return;
}
//--------------------------------
void load_PunctureField(parameters *par){
	FILE *fr;

	char fn[500];
	sprintf(fn,"InputData/Coord%d/r0_over_M_%.5lf/eta_%3.5lf/m_%d/nbarMax_%d/N_%d/Prec%3.5f/phi_punc_nbar%d.dat", (*par).CoordMap_FLAG, (*par).r0_over_M, (*par).eta, (*par).m, (*par).nbar, (*par).N1_PuncSeff, (*par).prec,(*par).nbar);

	if( (fr = fopen(fn,"r")) == NULL){
		fprintf(stderr, "Error in load_EffectiveSource: file %s not found\n", fn);
		// return;
		exit(-1);
	}

	int i1, i2, N1, N2;
	read_header(*par, fr, &N1, &N2);

	(*par).N1_LoadSeff = N1;
	(*par).N2_LoadSeff = N2;

	(*par).Re_cheb_PuncField = dmatrix(0, N1, 0 , N2);
	(*par).Im_cheb_PuncField = dmatrix(0, N1, 0 , N2);

	fscanf(fr,"#1:i1	 #2:i2	 #3:Real(Cheb phi_punc)	 #4:Imag(Cheb phi_punc)\n");

	for(i1=0; i1<=N1; i1++){
		for(i2=0; i2<=N2; i2++){
			int i1_read, i2_read;
			double Re_c_Seff, Im_c_Seff;

			fscanf(fr,"%d\t%d\t%lf\t%lf\n", &i1_read, &i2_read, &Re_c_Seff, &Im_c_Seff);
			(*par).Re_cheb_PuncField[i1][i2] = Re_c_Seff;
			(*par).Im_cheb_PuncField[i1][i2] = Im_c_Seff;
		}
		fscanf(fr,"\n");
	}
	

	fclose(fr);
	return;
}
//--------------------------------
void free_external_data(parameters *par){
	
	free_dvector( (*par).Re_cheb_phi_Punc, 0, (*par).N1_PuncSeff);
	free_dvector( (*par).Im_cheb_phi_Punc, 0, (*par).N1_PuncSeff);

	free_dvector( (*par).Re_cheb_phi_Punc_sigma, 0, (*par).N1_PuncSeff);
	free_dvector( (*par).Im_cheb_phi_Punc_sigma, 0, (*par).N1_PuncSeff);

	free_dvector( (*par).Re_cheb_phi_Punc_y, 0, (*par).N1_PuncSeff);
	free_dvector( (*par).Im_cheb_phi_Punc_y, 0, (*par).N1_PuncSeff);

	free_dmatrix( (*par).Re_cheb_Seff, 0, (*par).N1_PuncSeff, 0, (*par).N2_PuncSeff);
	free_dmatrix( (*par).Im_cheb_Seff, 0, (*par).N1_PuncSeff, 0, (*par).N2_PuncSeff);
	return;
}

//--------------------------------
void free_puncture_domain(parameters *par){
		
	free_dmatrix( (*par).Re_cheb_PuncField, 0, (*par).N1_PuncSeff, 0, (*par).N2_PuncSeff);
	free_dmatrix( (*par).Im_cheb_PuncField, 0, (*par).N1_PuncSeff, 0, (*par).N2_PuncSeff);
	return;
}