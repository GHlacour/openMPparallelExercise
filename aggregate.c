#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "lapack.h"
#include "aggregate.h"
#include "randomlib.h"

/* Diagonalize with LAPACK (destructive version) */
void diagonalizeLPD(float* H, float* v, int N) {
     int INFO, lwork;
     float *work, *Hcopy;
     int i, j;
     /* Find lwork */
     lwork = -1;
     work = (float *)calloc(1, sizeof(float));
     ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
     lwork = work[0];
     /*  printf("LAPACK work dimension %d\n",lwork); */
     free(work);
     work = (float *)calloc(lwork, sizeof(float));
     Hcopy = (float *)calloc(N * N, sizeof(float));
     /* Copy Hamiltonian */
     for (i = 0; i < N; i++) {
         for (j = 0; j < N; j++) {
             Hcopy[i * N + j] = H[i * N + j];
         }
     }
     /* Call Lapack Routine */
     ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
     /*  printf("LAPACK opt. %f %f\n",work[0],work[0]/N); */
     if (INFO != 0) {
         printf("Something went wrong trying to diagonalize a matrix...\n");
         exit(0);
     }
     /* Move eigenvectors */
     for (i = 0; i < N; i++) {
         for (j = 0; j < N; j++) {
             H[i*N+j] = Hcopy[j*N+i]; /* Converting from FORTRAN format */
         }
     }
     /* Free memory */
     free(Hcopy), free(work);
     return;
}

/* Calculate contributions to spectral histogram */
void histogram(float min,float max,int bins,int N,float *E,float *c,float *array){
    int i,j,bin;
    float binwidth;
    float intensity;
    binwidth=(max-min)/bins;
    /* Loop over states */
    for (i=0;i<N;i++){
        /* Find bin number */
        if (E[i]>min && E[i]<max){
            bin=round((E[i]-min)/binwidth);
            /* Add intensity to the histogram */
            intensity=0;
            for (j=0;j<N;j++){
               intensity+=c[j*N + i];
            }
            array[bin]+=intensity*intensity; 
        }        
    }
    return;
}

/* Build Aggregate Hamiltonain */
void hamiltonian(float *H,int N,float J,float w0,float sigma){
  int i,j;
  float idist,idist3;
  /* Generate Random Diagonal part of Hamiltonian */
  for(i=0;i<N;i++){
     H[i*N+i]=w0+RandomGaussian(0,sigma);
  }
  /* Generate off-diagonal part of Hamiltonian */
  for (i=0;i<N;i++){
     for (j=0;j<N;j++){
       if (i!=j){
          idist=1.0/fabs(i-j);
          idist3=idist*idist*idist;
          H[i*N+j]=J*idist3;
       }
     }
  }
}

/* Main Program */
int main(int argc,char *argv[]){
  int length,bins;
  float sigma;
  float w0,w,J;
  float *H,*v;
  float *spectrum;
  int i,r;
  int realizations;
  FILE *output;
  float wmin,wmax;

  length=atoi(argv[1]); /* Number of dye molecules / length of aggregate */
  sigma=atof(argv[2]); /* Magintude of molecular disorder in cm-1 */
  w0=atof(argv[3]); /* Monomer average absorption energy in cm-1 */
  J=atof(argv[4]); /* Nearest Neighbor Coupling in cm-1 */
  realizations=atoi(argv[5]); /* Number of disorder realizations */
  bins=atoi(argv[6]); /* Number of bins in histogram */  

  /* Initialize arrays */
  v=(float *)calloc(length, sizeof(float));
  spectrum=(float *)calloc(bins, sizeof(float));
  H=(float *)calloc(length * length, sizeof(float));
  /* Initialize spectral parameters (guess reasonable spectral range) */
  wmin=w0-3*J-2*sigma;
  wmax=w0+3*J+2*sigma;
  /* Initialize Random Number Generator */
  RandomInitialise(2511,1974);
  
  /* Loop over relaizations */
  for(r=0;r<realizations;r++){ 
    /* Generate Random Diagonal part of Hamiltonian */
    hamiltonian(H,length,J,w0,sigma);

    /* Diagonalize Hamiltonian */
    diagonalizeLPD(H,v, length);
    
    /* Add contribution to the spectrum */
    histogram(wmin,wmax,bins,length,v,H,spectrum);
  } 
  
  /* Save spectrum to file */
  output=fopen("Spectrum.txt","w");
  for(i=0;i<bins;i++){
     w=wmin+(0.5+i)*(wmax-wmin)/bins;
     fprintf(output,"%f %f\n",w,spectrum[i]);
  } 
  fclose(output);
  free(v);
  free(spectrum);
  free(H);
  return 0;
}
