#include <fftw3.h>
#include <cmath>
#include "const.h"

void setvar(int i_Nx, int i_Ny, double i_dl, double i_dt);
void Uk(fftw_complex *arr, fftw_complex *proc, fftw_plan pl_f, fftw_plan pl_b);
void Uv(fftw_complex *arr, double *V);
void normalize(fftw_complex *arr);
void U_ev(fftw_complex *arr, double *V, fftw_complex *proc, fftw_plan pl_f, fftw_plan pl_b);