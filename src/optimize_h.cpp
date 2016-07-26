///////////////////////////////////////////////////////////////////////////////
// This file is adapted from lda-alpha.c in the LDA-C package. It's updated to
// incorporate R lgamma, digamma, and trigamma functions. It also includes some
// minor changes to the notation.
///////////////////////////////////////////////////////////////////////////////
//
// The original copyright notice is given below:
//
// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)
//
// This file is part of LDA-C.
//
// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

# include "optimize_h.h"


/*
* Objective function and its derivatives
*
*/
double alhood(double a, double ss, int num_samples, int num_dim) {
  return(num_samples * (lgamma(num_dim * a) - num_dim * lgamma(a))
           + (a - 1) * ss);
}

double d_alhood(double a, double ss, int num_samples, int num_dim) {
  return(num_samples * (num_dim * Rf_digamma(num_dim * a)
                          - num_dim * Rf_digamma(a)) + ss);
}

double d2_alhood(double a, int num_samples, int num_dim) {
  return(num_samples * (num_dim * num_dim * Rf_trigamma(num_dim * a)
                          - num_dim * Rf_trigamma(a)));
}


/*
* Newton's method
*
* References:
*   Beyond Newton's Method. Minka (2000). Technical Report. Microsoft.
*
* This function is adapted from the original function
*
* double opt_alpha(double ss, int D, int K)
*
* given in the LDA-C package. All rights are with the author David M. Blei.
*
*/
double opt_hp(double init_a, double ss, int num_samples, int num_dim) {
  double a, log_a;
  double f, df, d2f;
  int iter = 0;

  printf("\nHyperparamter Optimization (ss: %5.5f samples: %d dim: %d init-hp: %.3f)\n",
         ss, num_samples, num_dim, init_a);
  printf("\n(c) Copyright 2004, David M. Blei\n\n");
  log_a = log(init_a);
  do {
    iter++;
    a = exp(log_a);
    if (isnan(a)) {
      init_a = init_a * 10;
      printf("warning : parameter is nan; new init = %5.5f\n", init_a);
      a = init_a;
      log_a = log(a);
    }
    f = alhood(a, ss, num_samples, num_dim);
    df = d_alhood(a, ss, num_samples, num_dim);
    d2f = d2_alhood(a, num_samples, num_dim);
    log_a = log_a - df/(d2f * a + df);
    printf("iter %5d: hp: %5.5f  f: %5.5f  df: %5.5f  d2f: %5.5f\n",
           iter, exp(log_a), f, df, d2f);
  }
  while ((fabs(df) > NEWTON_THRESH) && (iter < MAX_ALPHA_ITER));

  return (exp(log_a));
}
