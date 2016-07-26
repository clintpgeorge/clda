// /////////////////////////////////////////////////////////////////////////////
//
// Constriant Newton's method for \tau (Collection-level topic mixtures)
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



# include "optimize_tau.h"

///////////////////////////////////////////////////////////////////////////////
// Newton's method for \tau (Collection-level topic mixtures)
//
// We first split the Dirichlet variable \btau_j into a base measure \bomega_j
// and scale parameter a_j. We then optimize the lower-bound on \btau_j with
// respect to \bomega_j and a_j alternately.
//
// Modified on: May 28, 2016
// Created by:  Clint P. George
///////////////////////////////////////////////////////////////////////////////
mat optimize_tau(
    mat vi_tau,
    mat vi_tau_rho_ss,
    vec collection_size,
    unsigned int num_collections,
    unsigned int num_topics,
    double alpha_h,
    double gamma_h,
    unsigned int max_opt_iter,
    double step_size,
    unsigned int verbose
){

  unsigned int j;
  unsigned int opt_iter;
  double a_j;
  double grad_lb_a_j;
  double hess_lb_a_j;
  double vi_a_delta_j;
  double vi_lb_tau_j;
  double conv_ratio;
  double vi_lb_tau_j_old;
  double conv_thresh = 1e-10;
  vec grad_lb_omega_j;
  vec hess_lb_omega_j;
  vec vi_omega_delta_j;
  vec tmp;
  vec omega_j;

  if (verbose > 1){
    cout << endl << "tau - constraint Newton iterations..." << endl << endl;
    cout << "Reference: Kim, Voelker, and Saul (2013)" << endl << endl;
  }

  for (j = 0; j < num_collections; j++){ // for each collection

    // initialize the concentration parameter and the base measure

    a_j = sum(vi_tau.col(j));
    omega_j = vi_tau.col(j) / a_j;
    vi_lb_tau_j = 0;

    opt_iter = 0;
    conv_ratio = 1.;

    while ((fabs(conv_ratio) > conv_thresh) && (opt_iter < max_opt_iter)) { // for each Newton step

      opt_iter++;
      vi_lb_tau_j_old = vi_lb_tau_j;

      if (verbose > 1){
        cout << "collection #" << (j + 1);
        cout << " iter #" << opt_iter;
      }

      /////////////////////////////////////////////////////////////////////////
      // Constraint Newton update for \omega_j
      /////////////////////////////////////////////////////////////////////////
      tmp = (
        alpha_h
        + collection_size(j)
        - a_j * omega_j
        - ( gamma_h * collection_size(j) * omega_j )
      );

      grad_lb_omega_j = a_j * trigamma_vec( a_j * omega_j ) % tmp;
      grad_lb_omega_j -= gamma_h * collection_size(j) * digamma_vec(a_j * omega_j);
      grad_lb_omega_j += gamma_h * vi_tau_rho_ss.col(j);
      grad_lb_omega_j -= collection_size(j) * (
        gamma_h * digamma_vec(gamma_h * omega_j)
        + (1. / omega_j)
        - gamma_h
        - gamma_h * log(a_j * omega_j)
      );

      hess_lb_omega_j = a_j * a_j * tetragamma_vec(a_j * omega_j) % tmp;
      hess_lb_omega_j -= a_j * trigamma_vec(a_j * omega_j) * (a_j + 2. * gamma_h * collection_size(j));
      hess_lb_omega_j -= collection_size(j) * (
        gamma_h * gamma_h * trigamma_vec(gamma_h * omega_j)
        - ( 1. / (omega_j % omega_j) )
        - (gamma_h / omega_j)
      );


      hess_lb_omega_j = (1. / hess_lb_omega_j); // inv(diag(H))
      assert(hess_lb_omega_j.is_finite());
      vi_omega_delta_j = (
        (sum(grad_lb_omega_j % hess_lb_omega_j) / sum(hess_lb_omega_j)) * hess_lb_omega_j
        - grad_lb_omega_j % hess_lb_omega_j
      );
      if (verbose > 1){
        cout << " tau-delta-sum: " << sum(vi_omega_delta_j);
      }
      // assert(sum(vi_omega_delta_j) == 0);

      // The constraint Newton update

      omega_j += step_size * vi_omega_delta_j;

      if (verbose > 1){
        cout << " tau-sum: " << sum(omega_j);
      }
      // assert(sum(omega_j.col(j)) == 1);

      /////////////////////////////////////////////////////////////////////////
      // Newton update for a_j
      /////////////////////////////////////////////////////////////////////////
      tmp = (
        alpha_h
        + collection_size(j)
        - a_j * omega_j
        - ( gamma_h * collection_size(j) * omega_j )
      );

      grad_lb_a_j = sum (
        ( omega_j % trigamma_vec(a_j * omega_j) - Rf_trigamma(a_j) ) % tmp
      );
      grad_lb_a_j += ( (num_topics - 1) * gamma_h * collection_size(j) / (a_j * a_j) );

      hess_lb_a_j = sum (
        (
            omega_j % omega_j % tetragamma_vec(a_j * omega_j)
            - Rf_tetragamma(a_j)
        ) % tmp
      );
      hess_lb_a_j -= sum (
        omega_j % ( omega_j % trigamma_vec(a_j * omega_j) - Rf_trigamma(a_j) )
      );
      hess_lb_a_j -= ( 2. * (num_topics - 1) * gamma_h * collection_size(j) / (a_j * a_j * a_j) );

      vi_a_delta_j = - grad_lb_a_j / hess_lb_a_j;
      assert(is_finite(vi_a_delta_j));

      // The Newton update

      a_j += step_size * vi_a_delta_j;

      /////////////////////////////////////////////////////////////////////////
      // Computes the variational lowerbound using the new a_j and \omega_j
      /////////////////////////////////////////////////////////////////////////
      vi_lb_tau_j = sum(
        (alpha_h - a_j * omega_j) % ( digamma_vec(a_j * omega_j) - Rf_digamma(a_j) )
        + log_gamma_vec(a_j * omega_j)
      );
      vi_lb_tau_j -= lgamma(a_j);
      vi_lb_tau_j -= collection_size(j) * gamma_h * (num_topics - 1.) / a_j;
      vi_lb_tau_j -= (gamma_h - num_topics) * (
        collection_size(j) * ( log(a_j) - Rf_digamma(a_j) )
        + sum(vi_tau_rho_ss.col(j))
      );
      vi_lb_tau_j -= sum(
        collection_size(j) * log_gamma_vec(gamma_h * omega_j)
        + (1. - gamma_h * omega_j) % (
            collection_size(j) * ( log(a_j * omega_j) - digamma_vec(a_j * omega_j) )
            + vi_tau_rho_ss.col(j)
        )
      );
      if (verbose > 1){
        cout << " a-j: " << a_j;
        cout << " tau-lb-j: " << vi_lb_tau_j;
      }

      if (opt_iter > 5){ // minimum number of iterations
        conv_ratio = (vi_lb_tau_j_old - vi_lb_tau_j) / vi_lb_tau_j_old;
        if (verbose > 1){ cout << " tau-conv-ratio: " << conv_ratio; }
      }

      if (verbose > 1){ cout << endl; }

    } // for each Newton step

    vi_tau.col(j) = a_j * omega_j;


    if (verbose > 1){ cout << endl; }

  } // for each collection


  return vi_tau;
}
