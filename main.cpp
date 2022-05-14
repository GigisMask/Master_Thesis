#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <fftw3.h>

#include "src/init_state.h"
#include "src/pot.h"
#include "src/const.h"
#include "src/ops.h"

#include <boost/filesystem.hpp>
//----------GSL libraries---------//

// Random number generator
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main()
{
    int Nx = 100;
    int Ny = 100;
    int Vix = 5;
    int Viy = 5;
    int Vnx = 1;
    int Vny = 1;
    double dl = 1;
    double *V = (double *)malloc(sizeof(double) * Nx * Ny);
    fftw_complex *Psi = fftw_alloc_complex(Nx * Ny); // Wavefunction matrix in position space
    fftw_complex *Psi_p = fftw_alloc_complex(Nx * Ny);
    fftw_plan pl_momentum_position = fftw_plan_dft_2d(Nx, Ny, Psi_p, Psi, FFTW_BACKWARD, FFTW_MEASURE);

    double r0[2] = {0., 10};
    double k0[2] = {0., M_PI/(4* Ny * dl)};
    double sigma = 10;
    // random_values(V, 1, Nx, Ny,Vix, Viy, Vnx, Vny);
    // corr_pot(V, 1, 3, Nx, Ny, Vix, Viy, Vnx, Vny, dl);
    gaussian_wavepacket(Psi_p, k0, r0, sigma, Nx, Ny, dl);
    int i, j;
    int sign;
    /*
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            if ((i + j) % 2 == 0)
                sign = 1;
            else
                sign = -1;

            Psi_p[i + j * Nx][REAL] *= sign;
            Psi_p[i + j * Nx][IMAG] *= sign;
        }
    }
    fftw_execute(pl_momentum_position);*/
    for (int i = 0; i < Nx; i++) // Normalize
        for (int j = 0; j < Ny; j++)
        {
            Psi[i + j * Nx][REAL] = 1. / (Nx * Ny) * Psi[i + j * Nx][REAL];
            Psi[i + j * Nx][IMAG] = 1. / (Nx * Ny) * Psi[i + j * Nx][IMAG];
        }
    std::ofstream pot_out("data/pot_out.dat", std::ios::trunc);
    if (!pot_out.is_open())
    {
        free(V);
        fftw_free(Psi);
        fftw_free(Psi_p);
        fftw_destroy_plan(pl_momentum_position);
        return -2;
    }
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            //pot_out << V[i * Ny + j] << "\t";
            pot_out << Psi_p[i*Nx + j][REAL]*Psi_p[i*Nx + j][REAL] + Psi_p[i*Nx + j][IMAG]*Psi_p[i*Nx + j][IMAG] << "\t";
            //pot_out << Psi[i*Nx + j][REAL]*Psi[i*Nx + j][REAL] + Psi[i*Nx + j][IMAG]*Psi[i*Nx + j][IMAG] << "\t";
        }
        pot_out << std::endl;
    }
    free(V);
    fftw_free(Psi);
    fftw_free(Psi_p);
    fftw_destroy_plan(pl_momentum_position);
    return 0;
}