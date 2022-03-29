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
            arr[i + j * Ny][REAL] = A * cos((k[0] * (i - ix_wave) + k[1] * (j - iy_wave)) * dl);
            arr[i + j * Ny][IMAG] = A * sin((k[0] * (i - ix_wave) + k[1] * (j - iy_wave)) * dl);
        }
    }
}

void plane_wave_momentum(fftw_complex *arr, int *k, int Nx, int Ny)
{
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            arr[i + j * Ny][REAL] = (k[0]==i && k[1] == j);
            arr[i + j * Ny][IMAG] = 0.;
        }
        
    }
    
    arr[k[0] + k[1] * Ny][REAL] = 1.;
}

void gaussian_wavepacket(fftw_complex *arr, int Nx, int Ny, double dl, double *r0, double a, double *k0)
{
    int ix_wave = (int)r0[0] / dl;
    int iy_wave = (int)r0[1] / dl;
    double coeff;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            coeff = pow(2 * a * a * M_1_PI, 3. / 4.) * 1 / (a * a * a);
            arr[i + j * Ny][REAL] =
                coeff *
                exp(-1 / (a * a) * ((i - ix_wave) * (i - ix_wave) + (j - iy_wave) * (j - iy_wave)) * dl * dl) *
                cos(k0[0] * (i - ix_wave) * dl + k0[1] * (j - iy_wave) * dl);
            arr[i + j * Ny][IMAG] =
                coeff *
                exp(-1 / (a * a) * ((i - ix_wave) * (i - ix_wave) + (j - iy_wave) * (j - iy_wave)) * dl * dl) *
                sin(k0[0] * (i - ix_wave) * dl + k0[1] * (j - iy_wave) * dl);
        }
    }
}