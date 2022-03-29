#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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
int main(int argc, char **argv)
{
    if (argc != 9)
    {
        printf("Please introduce the arguments: x_size y_size dl_size tot_time dt_time corr_len V0 seed\n");
        return 0;
    }
    /*--------------Physical space parameters-----------------*/
    double Lx = atof(argv[1]); // voir pq lorsque L augmente il faut un tf plus grand
    double Ly = atof(argv[2]);
    double dl = atof(argv[3]);
    int Nx = (int)Lx / dl;
    int Ny = (int)Ly / dl;
    if (Lx <= dl || Ly <= dl)
    {
        printf("x_size and y_size must be greater than dl_size");
        return -1;
    }
    /*-------------------Time parameters----------------------*/
    double tf = atof(argv[4]);
    double dt = atof(argv[5]);
    if (tf <= dt)
    {
        printf("tot_time must be greater than dt_time");
        return -1;
    }
    /*-----------------Potential parameters-------------------*/
    double corr_len = atof(argv[6]); // 2; // Indirectly defines the correlation length between sites
    double V0 = atof(argv[7]);       // 0.1;
    //int seed = atof(argv[8]);
    double *V = (double *)malloc(Nx * Ny * sizeof(double));
    gsl_rng_default_seed = atof(argv[8]);
    std::cout << Lx << "\t" << Ly << "\t" << dl << "\t" << tf << "\t" << dt << "\t" << corr_len << "\t" << V0 << std::endl;
    
    corr_pot(V, V0, corr_len, Nx, Ny, dl);
    /*------------------Prepare for simulation---------------*/
    setvar(Nx, Ny, dl, dt);
    fftw_complex *Psi = fftw_alloc_complex(Nx * Ny);   // Wavefunction matrix in position space
    fftw_complex *Psi_p = fftw_alloc_complex(Nx * Ny); // Wavefunction matrix in the momentum space
    fftw_plan pl_position_momentum = fftw_plan_dft_2d(Nx, Ny, Psi, Psi_p, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan pl_momentum_position = fftw_plan_dft_2d(Nx, Ny, Psi_p, Psi, FFTW_BACKWARD, FFTW_MEASURE);
    /*----------------Initial wavefunction------------------*/
    double k[2] = {1, 0};
    int n0[2] = {(int)(Lx * k[0] / (2 * M_PI) + Nx / 2.), (int)(Ly * k[1] / (2 * M_PI) + Ny / 2.)};
    plane_wave_momentum(Psi_p, n0, Nx, Ny);
    // plane_wave(Psi, Nx, Ny, 1, k, dl, r0, 100);
    // Corriger et fixer k phys et non pas k num

    /* --------------File output parameters-----------------*/
    std::string filename;
    std::ofstream file;
    int fps = 24;

    /*---------------------File output---------------------*/
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
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                ph_esp_file << sqrt(Psi_p[i + j * Ny][REAL] * Psi_p[i + j * Ny][REAL] + Psi_p[i + j * Ny][IMAG] * Psi_p[i + j * Ny][IMAG]) << "\t";
        ph_esp_file << std::endl;
    }
    else
        std::cout << "Error creating output file" << std::endl;

    free(V);
    fftw_free(Psi);
    fftw_free(Psi_p);
    fftw_destroy_plan(pl_position_momentum);
    fftw_destroy_plan(pl_momentum_position);
    return 0;
}