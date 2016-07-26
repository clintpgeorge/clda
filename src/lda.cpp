// Compound Latent Dirichlet Allocation Model
//
// Copyright (c) 2016  Clint P. George
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

# include "lda.h"

///////////////////////////////////////////////////////////////////////////////
// LDA: Helper functions
///////////////////////////////////////////////////////////////////////////////

/***
* Computes the log posterior probability of the LDA model up to a
* multiplicative constant.
*
*/
double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    double alpha_h,
    double eta_h){

  double lp = 0.0;
  unsigned int d, i, num_docs, num_topics, vocab_size;
  arma::mat log_theta = log(theta_samples);
  arma::mat log_beta = log(beta_samples);
  num_topics = beta_samples.n_rows;
  vocab_size = beta_samples.n_cols;
  num_docs = theta_samples.n_cols;
  arma::vec n_dj(num_topics);
  arma::mat m_djt(num_topics, vocab_size);

  for (d = 0; d < num_docs; d++){ // for each document

    vector < unsigned int > word_ids = doc_word_ids[d];
    vector < unsigned int > word_zids = doc_word_zids[d];

    n_dj.fill(0.);
    m_djt.fill(0.);
    for (i = 0; i < doc_word_counts(d); i++){
      n_dj(word_zids[i]) += 1;
      m_djt(word_zids[i], word_ids[i]) += 1;
    }

    lp += arma::accu(m_djt % log_beta);
    lp += arma::accu((n_dj + alpha_h - 1.0) % log_theta.col(d));

  }

  lp += arma::accu((eta_h - 1.0) * log_beta);

  return lp;
}

double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    vector < vector < unsigned int > > doc_word_class,
    double alpha_h,
    double eta_h){

  double lp = 0.0;
  unsigned int d, i, num_docs, num_topics, vocab_size;
  arma::mat log_theta = log(theta_samples);
  arma::mat log_beta = log(beta_samples);
  num_topics = beta_samples.n_rows;
  vocab_size = beta_samples.n_cols;
  num_docs = theta_samples.n_cols;
  arma::vec n_dj(num_topics);
  arma::mat m_djt(num_topics, vocab_size);

  for (d = 0; d < num_docs; d++){ // for each document

    vector < unsigned int > word_ids = doc_word_ids[d];
    vector < unsigned int > word_zids = doc_word_zids[d];
    vector < unsigned int > word_class = doc_word_class[d];

    n_dj.fill(0.);
    m_djt.fill(0.);
    for (i = 0; i < doc_word_counts(d); i++){
      if (word_class[i]) continue; // ignores test words from sampling
      n_dj(word_zids[i]) += 1;
      m_djt(word_zids[i], word_ids[i]) += 1;
    }

    lp += arma::accu(m_djt % log_beta);
    lp += arma::accu((n_dj + alpha_h - 1.0) % log_theta.col(d));

  }

  lp += arma::accu((eta_h - 1.0) * log_beta);

  return lp;
}

