
void plane_wave(fftw_complex *arr, int Nx, int Ny, double A, double *k, double dl, double *r0 ,int n_wavefronts);
void plane_wave_momentum(fftw_complex *arr, double *k, int Nx, int Ny, double dl);
void gaussian_wavepacket(fftw_complex *arr, double *k0, double *r0, double sigma, int Nx, int Ny, double dl);
