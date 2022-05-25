#include "pot.h"

void corr_pot(double *V, const double V0, double ksi, int Nx, int Ny, int Vix, int Viy, int Vnx, int Vny, double dl)
{
    double wn_sig = 1;
    double *wn = (double *)malloc(sizeof(double) * Nx * Ny);

    random_values(wn, wn_sig, Nx, Ny, Vix, Viy, Vnx, Vny);

    // Define the correlated potential in matricial form
    double g; // correlation function
    // k and l are the indices for the position of the final potential
    int i, j, k, l, m, n = 0;
    for (k = 0; k < Nx; k++)
    {
        for (l = 0; l < Ny; l++)
        {
            // i and j are the indices for used to do the correlation between the white noise values
            g = 0;
            for (i = k - (Nx - Nx % 2) / 2; i < k + (Nx + Nx % 2) / 2; i++)
            {
                if (i < 0)
                    m = i + Nx;
                else if (i >= Nx)
                    m = i - Nx;
                else
                    m = i;

                for (j = l - (Ny - Ny % 2) / 2; j < l + (Ny + Ny % 2) / 2; j++)
                {
                    if (j < 0)
                        n = j + Ny;
                    else if (j >= Ny)
                        n = j - Ny;
                    else
                        n = j;

                    g += exp(-((k - i) * (k - i) + (l - j) * (l - j)) * dl * dl / (2 * ksi * ksi)) * wn[m * Ny + n];
                }
            }
            g *= V0 / (SQRT_PI * ksi);

            V[k * Ny + l] = g;
        }
    }
    free(wn);
}

void random_values(double *wn, const double wn_sigma, int Nx, int Ny, int Vix, int Viy, int Vnx, int Vny)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    bool xCond, yCond;

    for (int i = 0; i < Nx; i++)
    {
        xCond = i >= Vix && i < Vix + Vnx;
        for (int j = 0; j < Ny; j++)
        {
            yCond = j >= Viy && j < Viy + Vny;
            if (xCond && yCond)
                wn[i * Ny + j] = gsl_ran_gaussian(r, wn_sigma); // we suppose the random variate to be centered
            else
                wn[i * Ny + j] = 0.;
        }
    }
    gsl_rng_free(r);
}