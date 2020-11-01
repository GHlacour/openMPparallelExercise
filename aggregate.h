// Diagonalize with LAPACK (destructive version)
void diagonalizeLPD(float* H, float* v, int N);
// Make spectral histogram
void histogram(float min,float max,int bins,int N,float *E,float *c,float *array);
void hamiltonian(float *H,int N,float J,float w0,float sigma);
