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
void load_Retarded_at_Boundary(parameters *par){
	FILE *fr = NULL;
	char fn[500];
	int l, m = (*par).m, l_min = m, lmax=(*par).lmax, iDom = 0, i2, N2_dom0 = (*par).N2[iDom];
	double eta_read, sigma_minus_read, sigma_plus_read,
		   real_phi_ret_lm, imag_phi_ret_lm, small = 1.e-15;
	double complex phi_ret_plus, phi_ret_minus;

	(*par).phi_ret_plus  = cvector(0, N2_dom0); fill0_cvector((*par).phi_ret_plus,  0, N2_dom0+1);
	(*par).phi_ret_minus = cvector(0, N2_dom0); fill0_cvector((*par).phi_ret_minus, 0, N2_dom0+1);

	for(l=l_min; l<=lmax; l+=2){
		sprintf(fn,"InputData/r0_over_M_%.5lf/eta_%3.5lf/m_%d/phi_ret_boundary_ell_%d.dat",(*par).r0_over_M, (*par).eta, (*par).m, l);
		

		if( (fr = fopen(fn,"r")) == NULL){
				fprintf(stderr, "Error in load_Retarded_at_Boundary: file %s not found\n", fn);
				exit(-1);
		}

		fscanf(fr, "eta\t%lf\tsigma_minus\t%lf\tsigma_plus\t%lf\n", &eta_read, &sigma_minus_read, &sigma_plus_read);
		

		if( fabs(1.- eta_read/(*par).eta) > small || fabs(sigma_minus_read-(*par).sigma_minus) > small || fabs(sigma_plus_read-(*par).sigma_plus) > small  ){
			fprintf(stderr,"Error in load_Retarded_at_Boundary: input parameters do not match simulation\n");
			fprintf(stderr,"eta = %3.20e (%3.20e)\n", eta_read, (*par).eta );
			fprintf(stderr,"sigma_minus = %3.15e (%3.15e)\n", sigma_minus_read, (*par).sigma_minus );
			fprintf(stderr,"sigma_plus = %3.15e (%3.15e)\n", sigma_plus_read, (*par).sigma_plus );
			
			exit(-1);
		}
		fscanf(fr, "real phi_minus\t%lf\t imag phi_minus\t%lf\n", &real_phi_ret_lm, &imag_phi_ret_lm);
		
		phi_ret_minus = real_phi_ret_lm + I*imag_phi_ret_lm;


		fscanf(fr, "real phi_plus\t%lf\t imag phi_plus\t%lf\n", &real_phi_ret_lm, &imag_phi_ret_lm);
		phi_ret_plus = real_phi_ret_lm + I*imag_phi_ret_lm;
		
		for(i2=0; i2<=N2_dom0; i2++){
			func_derivs_2D y_minus, y_plus;
			
			double Pl_axis_norm_minus, Ylm_norm_minus,
				   Pl_axis_norm_plus,  Ylm_norm_plus,
				   
				   chi_2 = (*par).grid_chi_2[iDom][i2],
				   
				   CosTheta_minus, CosTheta_plus, 
				   c_lm = get_clm(l,m);

			get_y(*par, iDom, -1., chi_2, &y_minus);
			
			CosTheta_minus = sqrt(y_minus.d0);

			get_y(*par, iDom,  1., chi_2, &y_plus);
			CosTheta_plus = sqrt(y_plus.d0);	

			// Pl_norm = Pl/sin(theta)^m;

			Pl_axis_norm_minus = plgndr(l, m, CosTheta_minus )*pow((1-y_minus.d0), m/2.);//plgndr_axis_normalised(l, m, CosTheta_minus );
			Pl_axis_norm_plus  = plgndr(l, m, CosTheta_plus )*pow((1-y_plus.d0), m/2.);;//plgndr_axis_normalised(l, m, CosTheta_plus );
			
			Ylm_norm_minus = c_lm*Pl_axis_norm_minus;
			Ylm_norm_plus = c_lm*Pl_axis_norm_plus;

			(*par).phi_ret_minus[i2] += phi_ret_minus * Ylm_norm_minus;
			(*par).phi_ret_plus[i2]  += phi_ret_plus  * Ylm_norm_plus;
		}
		fclose(fr);
	}
	
	return;
}
//--------------------------------
void free_external_data(parameters *par){
	
	free_dvector( (*par).Re_cheb_phi_Punc, 0, (*par).N1_LoadPunc);
	free_dvector( (*par).Im_cheb_phi_Punc, 0, (*par).N1_LoadPunc);

	free_dvector( (*par).Re_cheb_phi_Punc_sigma, 0, (*par).N1_LoadPunc);
	free_dvector( (*par).Im_cheb_phi_Punc_sigma, 0, (*par).N1_LoadPunc);

	free_dvector( (*par).Re_cheb_phi_Punc_y, 0, (*par).N1_LoadPunc);
	free_dvector( (*par).Im_cheb_phi_Punc_y, 0, (*par).N1_LoadPunc);

	free_dmatrix( (*par).Re_cheb_Seff, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);
	free_dmatrix( (*par).Im_cheb_Seff, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);

	// free_dmatrix( (*par).Re_cheb_PuncField, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);
	// free_dmatrix( (*par).Im_cheb_PuncField, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);

	free_cvector( (*par).phi_ret_plus  , 0,  (*par).N2[0] );
	free_cvector( (*par).phi_ret_minus , 0,  (*par).N2[0] );
	return;
}

//--------------------------------
void free_puncture_domain(parameters *par){
		
	free_dmatrix( (*par).Re_cheb_PuncField, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);
	free_dmatrix( (*par).Im_cheb_PuncField, 0, (*par).N1_LoadSeff, 0, (*par).N2_LoadSeff);
	return;
}