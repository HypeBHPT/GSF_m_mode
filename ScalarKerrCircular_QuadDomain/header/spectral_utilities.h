#define ChebyCoeff_cut 1.e-50

double Clenshaw_cos(double *alpha, int N, double phi);
double Clenshaw_sin(double *beta, int N, double phi);
double Clenshaw_Fourier(double *alpha, double *beta, int N, double phi);
double Clenshaw_Chebyshev(double *c, int N, double x);

void Fourier_Coefficients(double *psi, double *alpha, double *beta, int N);
void Chebyshev_Coefficients_Radau_RHS(double *psi, double *c, int N);
void Chebyshev_Coefficients_Radau_LHS(double *psi, double *c, int N);
void Chebyshev_Coefficients_Gauss(double *psi, double *c, int N);
void Chebyshev_Coefficients_Lobatto(double *psi, double *c, int N);
void Chebyshev_Coefficients(double *psi, double *c, int N, char *grid);


void Fourier_Collocations(double *psi, double *alpha, double *beta, int N);
void Chebyshev_Collocations_Radau_RHS(double *psi, double *c, int N);
void Chebyshev_Collocations_Radau_LHS(double *psi, double *c, int N);
void Chebyshev_Collocations_Gauss(double *psi, double *c, int N);
void Chebyshev_Collocations_Lobatto(double *psi, double *c, int N);
void Chebyshev_Collocations(double *psi, double *c, int N, char *grid);

void Fourier_Coefficients_Derivative(double *alpha, double *beta, double *d_alpha, double *d_beta, int N);
void Chebyshev_Coefficients_Derivative(double *c, double *dc, int N);
void Chebyshev_Coefficients_Integral(double *c, double *C, int N);
double Chebyshev_Definite_Integral_Coefficients(double *c, int N);

void Chebyshev_Integration_Vector_Gauss(int N, double *Int);
double Chebyshev_Definite_Integral_Collocations(double *psi, double *Int, int N);

void Chebyshev_Differentiation_Matrices_Radau_RHS(int N, double **D1, double **D2);
void Chebyshev_Differentiation_Matrices_Radau_LHS(int N, double **D1, double **D2);
void Chebyshev_Differentiation_Matrices_Gauss(int N, double **D1, double **D2);
void Chebyshev_Differentiation_Matrices_Lobatto(int N, double **D1, double **D2);

void Chebyshev_Collocations_Derivatives(double *psi, double *d1psi, double *d2psi, double **D1, double **D2, int N);

void Chebyshev_Coefficients_2D(double **X, double **c2D, int N1, char *grid1, int N2, char *grid2);
void ChebyshevFourier_Coefficients_2D(double **X, double **alpha2D, double **beta2D, int N1, char *grid1, int N2);
double Clenshaw_Chebyshev_2D(double **c2D, int N1, int N2, double x1, double x2);
double Clenshaw_ChebyshevFourier(double **alpha2D, double **beta2D, int N1, int N2, double x1, double x2);
double Chebyshev_Definite_2DIntegral_Coefficients(double **c, int N1, int N2);

