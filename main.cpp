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

//--------------Main--------------//
int main(int argc, char **argv)
{
    double Lx; // voir pq lorsque L augmente il faut un tf plus grand
    double Ly;
    double dl;
    double tf;
    double dt;
    double corr_len; // 2; // Indirectly defines the correlation length between sites
    double V0;
    int tot_frames;
    if (argc < 9)
    {
        printf("Please introduce the arguments: x_size y_size dl_size tot_time dt_time corr_len V0 seed tot_frames(optional int) name(optional int)\n");
        Lx = 100; // voir pq lorsque L augmente il faut un tf plus grand
        Ly = 100;
        dl = 1;
        tf = 50;
        dt = 1;
        corr_len = 2; // Indirectly defines the correlation length between sites
        V0 = 0.1;
        tot_frames = tf / dt;
    }
    else
    {
        /*--------------Physical space parameters-----------------*/
        Lx = atof(argv[1]); // voir pq lorsque L augmente il faut un tf plus grand
        Ly = atof(argv[2]);
        dl = atof(argv[3]);
        /*-------------------Time parameters----------------------*/
        tf = atof(argv[4]);
        dt = atof(argv[5]);
        /*-----------------Potential parameters-------------------*/
        corr_len = atof(argv[6]); // Indirectly defines the correlation length between sites
        V0 = atof(argv[7]);
        /*-----------------------Seed-----------------------------*/
        gsl_rng_default_seed = atoi(argv[8]);
        /*----------------------Total frames----------------------*/
        tot_frames = atoi(argv[9]);
    }

    int Nx = (int)Lx / dl; // Number of sites on the x-axis
    int Ny = (int)Ly / dl; // Number of sites on the y-axis
    if (Lx <= dl || Ly <= dl)
    {
        printf("x_size and y_size must be greater than dl_size");
        return -1;
    }

    if (tf <= dt)
    {
        printf("tot_time must be greater than dt_time");
        return -1;
    }

    double *V = (double *)malloc(Nx * Ny * sizeof(double));

    // std::cout << Lx << "\t" << Ly << "\t" << dl << "\t" << tf << "\t" << dt << "\t" << corr_len << "\t" << V0 << std::endl;

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
    std::ostringstream dirName;
    dirName << Lx << "x" << Ly << "_V0_" << V0 << "_" << corr_len;
    std::string fileName;
    if (atof(argv[10]) == 0)
        fileName = "debug";
    else
        fileName = argv[10];

    if (boost::filesystem::exists("data/" + dirName.str()) == false)
    {
        bool mkdir = boost::filesystem::create_directories("data/" + dirName.str());
        if (mkdir == false)
            return -2;
    }

    std::ofstream ph_esp_file("data/" + dirName.str() + "/" + fileName + ".dat", std::ios::trunc);

    /*---------------------File output---------------------*/
    if (ph_esp_file.is_open())
    {
        for (double t = 0; t < tf; t += dt)
        {
            if ((int)(t / dt) % (int)tf/tot_frames == 0)
            {
                for (int i = 0; i < Nx; i++)
                    for (int j = 0; j < Ny; j++)
                        ph_esp_file << sqrt(Psi_p[i + j * Ny][REAL] * Psi_p[i + j * Ny][REAL] + Psi_p[i + j * Ny][IMAG] * Psi_p[i + j * Ny][IMAG]) << "\t";
                ph_esp_file << std::endl;
            }

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