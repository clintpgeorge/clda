// /////////////////////////////////////////////////////////////////////////////
//
// Utility Functions:
//
// This file is part of clda.
//
// Copyright (c) 2016  Clint P. George
//
// clda is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// clda is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
//
// /////////////////////////////////////////////////////////////////////////////

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <math.h>
# include <assert.h>

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

template<typename T>
extern
  T mod(T a, int n);

extern
  unsigned int sample_uniform_int (unsigned int K);


extern
  unsigned int sample_multinomial (arma::vec theta);

extern
  arma::vec sample_dirichlet (
      unsigned int num_elements,
      arma::vec alpha
    );

extern
  arma::rowvec sample_dirichlet_row_vec (
      unsigned int num_elements,
      arma::rowvec alpha
    );

extern
  arma::uvec randperm (unsigned int n);

extern
  arma::vec log_gamma_vec (arma::vec x_vec);

extern
  arma::rowvec log_gamma_rowvec(arma::rowvec x_vec);

extern
  arma::vec digamma_vec (arma::vec x_vec);

extern
  arma::rowvec digamma_rowvec(arma::rowvec x_vec);

extern
  arma::vec trigamma_vec (arma::vec x_vec);

extern
  arma::vec tetragamma_vec (arma::vec x_vec);

extern
  arma::vec gamma_col_vec (arma::vec x_vec);

extern
  double sample_antoniak(unsigned int N, double alpha);
