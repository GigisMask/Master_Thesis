#include "pot.h"

void corr_pot(double *V, const double V0, double ksi, int Nx, int Ny, int Vix, int Viy, int Vnx, int Vny, double dl)
{
    double wn_sig = 1;
    double *wn = (double *)malloc(sizeof(double) * Nx * Ny);

    double corr_radius = 5*ksi;

    random_values(wn, wn_sig, Nx, Ny, Vix, Viy, Vnx, Vny);

    // Define the correlated potential in matricial form
    double g; // correlation function
    // k and l are the indices for the position of the final potential
    int i, j, k, l, n, m = 0;
    for (k = 0; k < Nx; k++)
    {
        for (l = 0; l < Ny; l++)
        {
            // i and j are the indices for used to do the correlation between the white noise values
            g = 0;
            for (i = k - corr_radius; i < k + corr_radius; i++)
            {
                int dist_x = dist(k, i, Nx);
                if (i < 0)
                    n = i + Nx;
                else if (i > Nx)
                    n = i - Nx;
                else
                    n = i;

                for (j = l - corr_radius; j < l + corr_radius; j++)
                {
                    int dist_y = dist(l, j, Ny);
                    double dist2 = dist_x * dist_x + dist_y * dist_y;

                    if (j < 0)
                        m = j + Ny;
                    else if (j > Ny)
                        m = j - Ny;
                    else
                        m = j;

                    if (dist2 > corr_radius * corr_radius)
                        g += 0;
                    else
                        g += exp(-dist2 * dl * dl / (2. * ksi * ksi)) * wn[n * Ny + m];
                }
            }
            g *= V0 / (SQRT_PI * ksi);

            V[k * Ny + l] = g;
        }
    }
    free(wn);
}
int dist(int i, int n, int N)
{
    int mod = (N + N % 2) / 2;
    double d;
    int d_in = abs(n - i);

    if (d_in == 0)
        d = 0;
    else if (d_in > mod)
        d = (N - N % 2) / 2 - d_in % mod;
    else if (d_in < mod)
        d = d_in % mod;
    else
        d = d_in - N % 2;

    return d;
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