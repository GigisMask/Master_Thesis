#include "const.h"
#include <fftw3.h>
#include <math.h>

void plane_wave(fftw_complex *arr, int Nx, int Ny, double A, double *k, double dl, double *r0, int n_wavefronts) //    k = 2Pi/lambda
{
    int ix_wave = (int)r0[0] / dl;
    int iy_wave = (int)r0[1] / dl;
    int Nx_wave = (int)k[0] * n_wavefronts / (2 * M_PI);
    int Ny_wave = (int)k[1] * n_wavefronts / (2 * M_PI);

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            arr[i + j * Nx][REAL] = A * cos((k[0] * (i - ix_wave) + k[1] * (j - iy_wave)) * dl);
            arr[i + j * Nx][IMAG] = A * sin((k[0] * (i - ix_wave) + k[1] * (j - iy_wave)) * dl);
        }
    }
}

void plane_wave_momentum(fftw_complex *arr, double *k, int Nx, int Ny, double dl)
{
    int n0[2] = {(int)(Nx * dl * k[0] / (2 * M_PI) + Nx / 2.), (int)(Ny * dl * k[1] / (2 * M_PI) + Ny / 2.)};
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            arr[i + j * Nx][REAL] = (n0[0] == i && n0[1] == j);
            arr[i + j * Nx][IMAG] = 0.;
        }
    }
}

void gaussian_wavepacket(fftw_complex *arr, double *k0, double *r0, double sigma, int Nx, int Ny, double dl)
{
    int i, j; // indices for the i-th and j-th elements of the matrix
    double kx, ky;
    double a, theta;
    double Coeff = sqrt(2. / M_PI);
    i = 0;
    for (kx = -M_PI / dl; kx < M_PI / dl; kx += 2 * M_PI / (Nx * dl))
    {
        j = 0;
        for (ky = -M_PI / dl; ky < M_PI / dl; ky += 2 * M_PI / (Ny * dl))
        {
            // Amplitude and angle of the imaginary number
            a = Coeff * exp(-1. / 4. * sigma * sigma * ((kx - k0[0]) * (kx - k0[0]) + (ky - k0[1]) * (ky - k0[1])));
            theta = ((kx - k0[0]) * r0[0] + (ky - k0[1]) * r0[1]);
            // Enter the values in the matrix
            arr[i + j * Nx][REAL] = a * cos(theta);
            arr[i + j * Nx][IMAG] = a * sin(theta);
            j++;
        }
        i++;
    }
}