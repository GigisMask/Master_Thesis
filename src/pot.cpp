#include "pot.h"

void corr_pot(double *V, const double V0, double sigma, int Nx, int Ny, double dl)
{
    double wn_sig = 1;
    double *wn = (double *)malloc(sizeof(double) * Nx * Ny);

    random_values(wn, wn_sig, Nx, Ny);

    // Define the correlated potential in matricial form
    double g; // correlation function
    // k and l are the indices for the position of the final potential
    for (int k = 0; k < Nx; k++)
    {
        for (int l = 0; l < Ny; l++)
        {
            // i and j are the indices for used to do the correlation between the white noise values
            g = 0;
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    g += exp(-(((k - i) * dl) * ((k - i) * dl) + ((l - j) * dl) * ((l - j) * dl)) / (2 * sigma * sigma)) * wn[i + j * Ny];
            g *= V0 / (SQRT_PI * sigma);

            V[k + l * Ny] = g;
        }
    }
    free(wn);
}

void random_values(double *wn, const double sigma, int Nx, int Ny)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            wn[i + j * Ny] = gsl_ran_gaussian(r, sigma); // we suppose the random variate to be centered
    // gsl_matrix_set(wn, i, j, gsl_ran_gaussian(r, sigma));
    gsl_rng_free(r);
}