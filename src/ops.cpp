#include "ops.h"

int Nx, Ny;
int sign;
double *kx, *ky, kx_i2, kij2;
double dt, dl;

void setvar(int i_Nx, int i_Ny, double i_dl, double i_dt)
{
    Nx = i_Nx;
    Ny = i_Ny;
    dt = i_dt;
    dl = i_dl;
    kx = (double *)malloc(sizeof(double) * Nx);
    ky = (double *)malloc(sizeof(double) * Ny);

    double pi_dl = M_PI / dl;
    for (int n = 0; n < Nx; n++)
        kx[n] = (-1 + 2. * n / Nx) * pi_dl;

    for (int m = 0; m < Ny; m++)
        ky[m] = (-1 + 2. * m / Ny) * pi_dl;
}

void Uk(fftw_complex *arr, fftw_complex *proc, fftw_plan pl_f, fftw_plan pl_b) // Applies the Uk operator to the arr matrix and gives the output out. This function modifies the input matrix and gives a non normalized version of it
{
    // Transformer en impulsion
    int i, j;
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            if ((i + j) % 2 == 0)
                sign = 1;
            else
                sign = -1;

            arr[i + j * Nx][REAL] *= sign;
            arr[i + j * Nx][IMAG] *= sign;
        }
    }
    fftw_execute(pl_f);

    double a, b;

    for (i = 0; i < Nx; i++)
    {
        kx_i2 = kx[i] * kx[i];
        for (j = 0; j < Ny; j++)
        {
            a = proc[i + j * Nx][REAL];
            b = proc[i + j * Nx][IMAG];
            kij2 = kx_i2 + ky[j] * ky[j];
            proc[i + j * Nx][REAL] = a * cos(-dt / 4. * kij2) - b * sin(-dt / 4. * kij2);
            proc[i + j * Nx][IMAG] = a * sin(-dt / 4. * kij2) + b * cos(-dt / 4. * kij2);
        }
    }
    fftw_execute(pl_b);
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            if ((i + j) % 2 == 0)
                sign = 1;
            else
                sign = -1;

            arr[i + j * Nx][REAL] *= sign;
            arr[i + j * Nx][IMAG] *= sign;
        }
    }
}

void Uv(fftw_complex *arr, double *V)
{
    int i, j;
    double a, b;
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            a = arr[i + j * Nx][REAL];
            b = arr[i + j * Nx][IMAG];
            arr[i + j * Nx][REAL] = a * cos(-V[i + j * Nx] * dt) - b * sin(-V[i + j * Nx] * dt);
            arr[i + j * Nx][IMAG] = a * sin(-V[i + j * Nx] * dt) + b * cos(-V[i + j * Nx] * dt);
        }
    }
}

void U_ev(fftw_complex *arr, fftw_complex *proc, double *V, fftw_plan pl_f, fftw_plan pl_b)
{
    Uk(arr, proc, pl_f, pl_b);
    Uv(arr, V);
    Uk(arr, proc, pl_f, pl_b);

    // normalize with Nx*Ny
}

/*========Momentum space operators =====*/

void Uk_p(fftw_complex *momentum_space) // Applies the Uk operator to the arr matrix and gives the output out. This function modifies the input matrix and gives a non normalized version of it
{
    int i, j;
    double a, b;

    for (i = 0; i < Nx; i++)
    {
        kx_i2 = kx[i] * kx[i];
        for (j = 0; j < Ny; j++)
        {
            a = momentum_space[i + j * Nx][REAL];
            b = momentum_space[i + j * Nx][IMAG];
            kij2 = kx_i2 + ky[j] * ky[j];
            momentum_space[i + j * Nx][REAL] = a * cos(-dt / 4. * kij2) - b * sin(-dt / 4. * kij2);
            momentum_space[i + j * Nx][IMAG] = a * sin(-dt / 4. * kij2) + b * cos(-dt / 4. * kij2);
        }
    }
}

void Uv_p(fftw_complex *real_space, fftw_complex *momentum_space, double *V, fftw_plan pl_position_momentum, fftw_plan pl_momentum_position)
{
    int i, j;
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            if ((i + j) % 2 == 0)
                sign = 1;
            else
                sign = -1;

            momentum_space[i + j * Nx][REAL] *= sign;
            momentum_space[i + j * Nx][IMAG] *= sign;
        }
    }
    fftw_execute(pl_momentum_position); // transformer en position

    Uv(real_space, V); // Apply Uv operator

    fftw_execute(pl_position_momentum); // transformer en impulsion

    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            if ((i + j) % 2 == 0)
                sign = 1;
            else
                sign = -1;

            momentum_space[i + j * Nx][REAL] *= sign;
            momentum_space[i + j * Nx][IMAG] *= sign;
        }
    }
}

void Uev_p(fftw_complex *momentum_space, fftw_complex *position_space, double *V, fftw_plan pl_position_momentum, fftw_plan pl_momentum_position)
{
    Uk_p(momentum_space);
    Uv_p(position_space, momentum_space, V, pl_position_momentum, pl_momentum_position);
    Uk_p(momentum_space);

    for (int i = 0; i < Nx; i++) // Normalize
        for (int j = 0; j < Ny; j++)
        {
            momentum_space[i + j * Nx][REAL] = 1. / (Nx * Ny) * momentum_space[i + j * Nx][REAL];
            momentum_space[i + j * Nx][IMAG] = 1. / (Nx * Ny) * momentum_space[i + j * Nx][IMAG];
        }
}
