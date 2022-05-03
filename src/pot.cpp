#include "pot.h"

void corr_pot(double *V, const double V0, double ksi, int Nx, int Ny, double dl)
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
                    g += exp(-((i - (Nx + (Nx % 2)) / 2 - k) * (i - (Nx + (Nx % 2)) / 2 - k) + (j - (Ny + (Ny % 2)) / 2 - l) * (j - (Ny + (Ny % 2)) / 2 - l)) * dl * dl / (2 * ksi * ksi)) * wn[i + j * Nx];
            g *= V0 / (SQRT_PI * ksi);

            V[k + l * Nx] = g;
        }
    }
    free(wn);
}

void random_values(double *wn, const double sigma, int Nx, int Ny, int V_Nx, int V_Ny)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    int iV_x = (int)(Nx - V_Nx) / 2.;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
        {
            if (i > (int)(Nx - V_Nx) / 2. && i < (int)(Nx + V_Nx) / 2.)
                wn[i + j * Nx] = gsl_ran_gaussian(r, sigma); // we suppose the random variate to be centered
            else
                wn[i + j * Nx] = 0.;
        }
    gsl_rng_free(r);
}