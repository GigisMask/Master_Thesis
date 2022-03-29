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
    double Lx = 100.; // voir pq lorsque L augmente il faut un tf plus grand
    double Ly = 100.;
    double dl = 1.;

    double tf = 50, dt = 1;

    int Nx = (int)Lx / dl;
    int Ny = (int)Ly / dl;

    fftw_complex *Psi = fftw_alloc_complex(Nx * Ny);   // Wavefunction matrix in space
    fftw_complex *Psi_p = fftw_alloc_complex(Nx * Ny); // Wavefunction matrix in the momentum space
    fftw_plan pl_position_momentum = fftw_plan_dft_2d(Nx, Ny, Psi, Psi_p, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan pl_momentum_position = fftw_plan_dft_2d(Nx, Ny, Psi_p, Psi, FFTW_BACKWARD, FFTW_MEASURE);

    double sigma = 2; // Indirectly defines the correlation length between sites
    double V0 = 0.1;

    double *V = (double *)malloc(Nx * Ny * sizeof(double));

    double r0[2] = {50, 50};

    // plane_wave(Psi, Nx, Ny, 1, k, dl, r0, 100);
    //Corriger et fixer k phys et non pas k num
     /// ==> En impulsion Point qui vaut 1 et le reste 0 en impulsion pour Eikr

    double k[2] = {1,0};
    int n0[2] = {(int) (Lx * k[0] / (2*M_PI) + Nx/2.), (int) (Ly * k[1] / (2*M_PI) + Ny/2.)};
    plane_wave_momentum(Psi_p, n0, Nx, Ny);

    corr_pot(V, V0, sigma, Nx, Ny, dl);

    setvar(Nx, Ny, dl, dt);

    // U_ev(Psi, V, Psi_p, p_position_momentum, p_momentum_position);

    /* --------File Output----------*/
    std::string filename;
    std::ofstream file;

    /*---------Video output parameters------*/
    int fps = 24;

    /*
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
    */

    /*

     std::ofstream ph_esp_file("data/phase_space_ev.dat", std::ios::trunc);
     if (ph_esp_file.is_open())
     {
         ph_esp_file << Nx << "\t" << Ny << "\t" << dt << "\t" << tf << std::endl;

         for (double t = 0; t < tf; t += dt)
         {
             for (int i = 0; i < Nx; i++)
                 for (int j = 0; j < Ny; j++)
                     ph_esp_file << sqrt(Psi_p[i + j * Ny][REAL] * Psi_p[i + j * Ny][REAL] + Psi_p[i + j * Ny][IMAG] * Psi_p[i + j * Ny][IMAG]) << "\t";
             ph_esp_file << std::endl;

             Uev_p(Psi_p, Psi, V, pl_position_momentum, pl_momentum_position);
         }

         Uev_p(Psi_p, Psi, V, pl_position_momentum, pl_momentum_position);
         ph_esp_file << sqrt(Psi_p[i + j * Ny][REAL] * Psi_p[i + j * Ny][REAL] + Psi_p[i + j * Ny][IMAG] * Psi_p[i + j * Ny][IMAG]) << "\t";
         ph_esp_file << std::endl;
     }
     else
         std::cout << "Error creating Psi file" << std::endl;
     */

    //
    std::ofstream ph_esp_file("data/2d_planewave.dat", std::ios::trunc);
    for (double t = 0; t < tf; t+=dt)
    {
        Uev_p(Psi_p, Psi, V, pl_position_momentum, pl_momentum_position);
    }
    
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
            ph_esp_file << sqrt(Psi_p[i + j * Ny][REAL] * Psi_p[i + j * Ny][REAL] + Psi_p[i + j * Ny][IMAG] * Psi_p[i + j * Ny][IMAG]) << "\t";
        ph_esp_file << std::endl;
    }

    free(V);
    fftw_free(Psi);
    fftw_free(Psi_p);
    fftw_destroy_plan(pl_position_momentum);
    fftw_destroy_plan(pl_momentum_position);
    return 0;
}