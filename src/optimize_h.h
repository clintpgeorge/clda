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


# include "utils.h"


#define NEWTON_THRESH 1e-5
#define MAX_ALPHA_ITER 50

/*
 * Objective function and its derivatives
 *
 */
extern
  double alhood(double a, double ss, int num_samples, int num_dim);

extern
  double d_alhood(double a, double ss, int num_samples, int num_dim);

extern
  double d2_alhood(double a, int num_samples, int num_dim);


/*
* Newton's method
*
*/
extern
  double opt_hp(double init_a, double ss, int num_samples, int num_dim);
