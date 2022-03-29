
void plane_wave(fftw_complex *arr, int Nx, int Ny, double A, double *k, double dl, double *r0 ,int n_wavefronts);
void gaussian_wavepacket(fftw_complex *arr, int Nx, int Ny, double dl, double *r0, double sigma, double* k0);
void plane_wave_momentum(fftw_complex *arr, int *k, int Nx, int Ny);
