/*  Copyright Ton Dieker                                                    */
/*  Centre of Mathematics and Computer Science (CWI) Amsterdam              */
/*  April 2002                                                              */

/*  ton@cwi.nl                                                              */

/*  Modified to be used with GSL (instead of Netlib libraries) by           */
/*  Aleksejus Kononovicius                                                  */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include "spectrum.h"

extern void apprcirc(long n, double H, double L, long seed, double *output);
