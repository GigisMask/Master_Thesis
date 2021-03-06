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
#include <boost/lexical_cast.hpp>

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
        dl = 0.1;
        tf = 2;
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

    corr_pot(V, V0, corr_len, Nx, Ny, 0, 0, Nx, Ny, dl);
    /*------------------Prepare for simulation---------------*/
    setvar(Nx, Ny, dl, dt);

    fftw_complex *Psi = fftw_alloc_complex(Nx * Ny);   // Wavefunction matrix in position space
    fftw_complex *Psi_p = fftw_alloc_complex(Nx * Ny); // Wavefunction matrix in the momentum space
    fftw_plan pl_position_momentum = fftw_plan_dft_2d(Nx, Ny, Psi, Psi_p, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan pl_momentum_position = fftw_plan_dft_2d(Nx, Ny, Psi_p, Psi, FFTW_BACKWARD, FFTW_MEASURE);

    /*----------------Initial wavefunction------------------*/
    double k0[2] = {2./corr_len, 0};
    double r0[2] = {0, 0};
    double sigma = 2;
    // plane_wave(Psi, Nx, Ny, 1, k0, dl, r0, 100);
    //plane_wave_momentum(Psi_p, k0, Nx, Ny, dl);
    
    gaussian_wavepacket(Psi_p, k0, r0, sigma, Nx, Ny, dl);
    //  Corriger et fixer k phys et non pas k num

    /* --------------File output parameters-----------------*/
    std::ostringstream dirName;
    dirName << int(Lx/dl) << "x" << int(Ly/dl) << "_" << dl << "_V0_" << V0 << "_" << corr_len;
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
    char a[32], b[32]; // used to output values in the form of a+b*j with j sqrt(-1) used for python postprocessing
    if (ph_esp_file.is_open())
    {
        for (double t = 0; t <= tf; t += dt)
        {
            if (((int)(t / dt)) % ((int)tf / tot_frames) == 0)
            {
                for (int i = 0; i < Nx; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        sprintf(a, "%.8f", Psi_p[i * Ny + j][REAL]); //10
                        sprintf(b, "%+.8f", Psi_p[i * Ny + j][IMAG]);
                        ph_esp_file << a << b << "j"
                                    << "\t";
                    }
                }
                ph_esp_file << std::endl;
                // printf("%f\n", t);
            }
            if (t < tf)
                Uev_p(Psi_p, Psi, V, pl_position_momentum, pl_momentum_position);
        }
        ph_esp_file.close();
    }
    else
        std::cout << "Error creating output file" << std::endl;
    
    ops_clear();
    free(V);
    fftw_destroy_plan(pl_position_momentum);
    fftw_destroy_plan(pl_momentum_position);
    fftw_free(Psi);
    fftw_free(Psi_p);
    return 0;
}