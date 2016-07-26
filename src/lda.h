// /////////////////////////////////////////////////////////////////////////////
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

# include "utils.h"

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

extern
double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    double alpha_h,
    double eta_h
  );

extern
double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    vector < vector < unsigned int > > doc_word_class,
    double alpha_h,
    double eta_h);
