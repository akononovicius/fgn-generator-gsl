/*  This program simulates fractional Gaussian noise or fractional          */
/*  Brownian motion using the approximate circulant algorithm.              */
/*  The C-packages Ranlib and Meschach are used, both available             */
/*  via Netlib (http://www.netlib.org).                                     */

/*  Reference:                                                              */
/*  A.B. Dieker and M. Mandjes (2002),                                      */
/*  On spectral simulation of fractional Brownian motion,                   */
/*  submitted for publication.                                              */

/*  Copyright Ton Dieker                                                    */
/*  Centre of Mathematics and Computer Science (CWI) Amsterdam              */
/*  April 2002                                                              */

/*  ton@cwi.nl                                                              */

/*  Modified to be used with GSL (instead of Netlib libraries) by           */
/*  Aleksejus Kononovicius                                                  */

#include "apprcirc.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

void apprcirc(long n, double Hurst, double L, long seed, double *output) {
    /* function that generates a fractional Brownian motion or fractional  */
    /* Gaussian noise sample using the approximate circulant method.       */
    /* Input:  n      determines the sample size N by N=2^(*n)             */
    /*         Hurst  the Hurst parameter of the trace                     */
    /*         L      the sample is generated on [0,L]                     */
    /*         seed   seed for the random generator                        */
    /* Output: *output the resulting sample is stored in this array        */
  
    long i, N, halfN, generator;
    double scaling, H;
    double aux;
    
    halfN=pow(2,n);
    H=Hurst;
    N=2*halfN;
    
    /* allocate memory */
    double *pow_spec=(double*)malloc((halfN+1)*sizeof(double));
    double *data=(double *)malloc(2*N*sizeof(double));

    /* set random generator and seeds */
    gsl_rng_env_setup();
    gsl_rng * rng=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng,seed);
    
    /* approximate spectral density */
    FGN_spectrum(pow_spec,halfN,H);
    
    REAL(data,0)=sqrt(2*(pow(N,2*H)-pow(N-1,2*H)))*gsl_ran_gaussian(rng,1);
    IMAG(data,0)=0.;
    REAL(data,halfN)=sqrt(2*pow_spec[halfN])*gsl_ran_gaussian(rng,1);
    IMAG(data,halfN)=0.;
    for(i=1;i<halfN;i++) {
        aux=sqrt(pow_spec[i]);
        REAL(data,i)=aux*gsl_ran_gaussian(rng,1);
        IMAG(data,i)=aux*gsl_ran_gaussian(rng,1);
    }
    for(i=halfN+1;i<N;i++) {
        REAL(data,i)=REAL(data,N-i);
        IMAG(data,i)=-IMAG(data,N-i);
    }
    
    /* real part of Fourier transform of data gives sample path */
    gsl_fft_complex_radix2_backward(data,1,N);

    /* rescale to obtain a sample of size 2^(*n) on [0,L] */
    scaling=pow(L/halfN,H)/sqrt(2*N);
    for(i=0;i<halfN;i++) {
        output[i]=scaling*(REAL(data,i));
    }
    
    /* free memory */
    free(pow_spec);
    free(data);
    gsl_rng_free(rng);
}
