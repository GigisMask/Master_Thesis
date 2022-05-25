#include "const.h"
#include <fftw3.h>
#include <math.h>

void plane_wave_momentum(fftw_complex *arr, double *k, int Nx, int Ny, double dl)
{
    int n0[2] = {(int)(Nx * dl * k[0] / (2 * M_PI) + Nx / 2.), (int)(Ny * dl * k[1] / (2 * M_PI) + Ny / 2.)};
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            arr[i * Ny + j][REAL] = (n0[0] == i && n0[1] == j);
            arr[i * Ny + j][IMAG] = 0.;
        }
    }
}

void gaussian_wavepacket(fftw_complex *arr, double *k0, double *r0, double sigma, int Nx, int Ny, double dl)
{
    int i, j; // indices for the i-th and j-th elements of the matrix
    double kx, ky;
    double a, theta;
    double Coeff = sigma / sqrt(2. * M_PI);
    i = 0;
    for (kx = -M_PI / dl; kx < M_PI / dl; kx += 2. * M_PI / (Nx * dl))
    {
        j = 0;
        for (ky = -M_PI / dl; ky < M_PI / dl; ky += 2. * M_PI / (Ny * dl))
        {
            // Amplitude and angle of the imaginary number
            a = Coeff * exp(-1. / 4. * sigma * sigma * ((kx - k0[0]) * (kx - k0[0]) + (ky - k0[1]) * (ky - k0[1])));
            theta = ((kx - k0[0]) * r0[0] + (ky - k0[1]) * r0[1]);
            // Enter the values in the matrix
            arr[i * Ny + j][REAL] = a * cos(theta);
            arr[i * Ny + j][IMAG] = a * sin(theta);
            j++;
        }
        i++;
    }
}