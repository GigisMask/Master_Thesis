#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const double SQRT_PI = sqrt(M_PI);

void corr_pot(double *V, const double V0, double sigma, int Nx, int Ny, int Vix, int Viy, int Vnx, int Vny, double dl);

void random_values(double *wn, const double sigma, int Nx, int Ny, int Vix, int Viy, int Vnx, int Vny);