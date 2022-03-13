#include <iostream>
#include <fstream>
#include <math.h>
#include <fftw3.h>

#include "src/init_state.h"
#include "src/pot.h"
#include "src/const.h"
#include "src/ops.h"

//----------GSL libraries---------//

// Random number generator
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//--------------Main--------------//
int main(int, char **)
{
    double Lx = 100.;
    double Ly = 100.;
    double dl = 1.;

    double tf = 1e-2, dt = 1e6;

    int Nx = (int)Lx / dl;
    int Ny = (int)Ly / dl;

    fftw_complex *Psi = fftw_alloc_complex(Nx * Ny);   // Wavefunction matrix in space
    fftw_complex *Psi_p = fftw_alloc_complex(Nx * Ny); // Wavefunction matrix in the momentum space
    fftw_plan p_space_momentum = fftw_plan_dft_2d(Nx, Ny, Psi, Psi_p, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan p_momentum_space = fftw_plan_dft_2d(Nx, Ny, Psi_p, Psi, FFTW_BACKWARD, FFTW_MEASURE);

    double sigma = 5; // Indirectly defines the correlation length between sites
    double V0 = 1;

    double *V = (double *)malloc(Nx * Ny * sizeof(double));

    double r0[2] = {50, 50};

    // plane_wave(Psi, Nx, Ny, 1, k, dl, r0, 100);
    double k0[2] = {0.1, 0.1};
    gaussian_wavepacket(Psi, Nx, Ny, dl, r0, 10, k0);

    corr_pot(V, 0, sigma, Nx, Ny, dl);

    setvar(Nx, Ny, dl, dt);
    

    U_ev(Psi, V, Psi_p, p_space_momentum, p_momentum_space);
    normalize(Psi);

    /*
        for (double t = 0; t < tf; t += dt)
        {
            U_ev(Psi, V, Psi_p, p_space_momentum, p_momentum_space);
            if ((int)(tf - t) % 5 == 0)
                normalize(Psi);
        }
    */

    std::ofstream pot_file("data/pot.dat", std::ios::trunc);
    if (pot_file.is_open() == false)
    {
        std::cout << "cannot create the file" << std::endl;
        return -5;
    }

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
            pot_file << V[i + j * Ny] << "\t";
        pot_file << "\n";
    }
    pot_file.close();

    std::ofstream pw_file("data/2d_planewave.dat", std::ios::trunc);
    if (pw_file.is_open())
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                pw_file << sqrt(Psi[i + j * Ny][REAL] * Psi[i + j * Ny][REAL] + Psi[i + j * Ny][IMAG] * Psi[i + j * Ny][IMAG]) << "\t";
                //pw_file << Psi[i + j * Ny][REAL] << "\t" ;
            }
            pw_file << std::endl;
        }
    }
    else
        std::cout << "Error creating Psi file" << std::endl;

    free(V);
    fftw_free(Psi);
    fftw_free(Psi_p);
    fftw_destroy_plan(p_space_momentum);
    fftw_destroy_plan(p_momentum_space);
    return 0;
}