// /////////////////////////////////////////////////////////////////////////////
//
// cLDA: Variational Expectation Maximization
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
# include "optimize_h.h"
# include "optimize_tau.h"

//' cLDA: Variational Expectation Maximization
//'
//' This implements the Variational Expectation Maximization (EM) algorithm for
//' the compound latent Dirichlet allocation (cLDA) model.
//'
//' To compute perplexity, we first partition words in a corpus into two sets:
//' (a) a test set (held-out set), which is selected from the set of words in
//' the test (held-out) documents (identified via \code{test_doc_share} and
//' \code{test_word_share}) and (b) a training set, i.e., the remaining words in
//' the corpus. We then run the variational EM algorithm based on the training
//' set. Finally, we compute per-word perplexity based on the held-out set.
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size Vocabulary size
//' @param docs_cid Documents collection IDs (ID indices starts 0)
//' @param docs_tf A list of corpus documents read from the Blei corpus using
//'             \code{\link{read_docs}} (term indices starts with 0)
//' @param alpha_h Hyperparameter for collection-level Dirichlets \eqn{\pi}
//' @param gamma_h Hyperparameter for document-level Dirichlets \eqn{\theta}
//' @param eta_h Hyperparameter for corpus level topic Dirichlets \eqn{\beta}
//' @param vi_max_iter Maximum number of iterations for variational inference
//' @param em_max_iter Maximum number of iterations for variational EM
//' @param vi_conv_thresh Convergence threshold for the document variational inference loop
//' @param em_conv_thresh Convergence threshold for the variational EM loop
//' @param tau_max_iter Maximum number of iterations for the constraint Newton updates of \eqn{\tau}
//' @param tau_step_size the step size for the constraint Newton updates of \eqn{\tau}
//' @param estimate_alpha If true, run hyperparameter \eqn{\alpha} optimization
//' @param estimate_gamma dummy parameter [not implemented]
//' @param estimate_eta If true, run hyperparameter \eqn{\eta} optimization
//' @param verbose from {0, 1, 2, 3}
//' @param init_pi the initial configuration for the collection level topic
//'                mixtures, i.e., \eqn{\pi} samples
//' @param test_doc_share proportion of the test documents in the corpus.
//'   Must be from [0., 1.)
//' @param test_word_share proportion of the test words in each test document.
//'   Must be from [0., 1.)
//'
//' @return A list of variational parameters
//'
//' @export
//'
//' @family Variational Inference
//'
//' @note Created on May 13, 2016
//'
// [[Rcpp::export]]
List clda_vem(
    unsigned int num_topics,
    unsigned int vocab_size,
    NumericVector docs_cid,
    List docs_tf,
    double alpha_h,
    double gamma_h,
    double eta_h,
    unsigned int vi_max_iter,
    unsigned int em_max_iter,
    double vi_conv_thresh,
    double em_conv_thresh,
    unsigned int tau_max_iter,
    double tau_step_size,
    bool estimate_alpha,
    bool estimate_gamma,
    bool estimate_eta,
    int verbose,
    arma::mat init_pi, // TODO: to be deleted in the final version
    double test_doc_share = 0.,
    double test_word_share = 0.
) {

  cout << endl << endl;
  if (verbose > 0){
    cout << "clda-vem (c++): Initializes variables and count statistics....";
  }

  assert(tau_step_size > 0);
  assert(test_word_share < 1.);
  assert(test_word_share >= 0.);
  assert(test_doc_share < 1.);
  assert(test_doc_share >= 0.);

  NumericVector cids = unique(docs_cid);
  unsigned int num_collections = cids.size(); // number of collections
  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int d, em_iter, vi_iter, c, word_id, word_count, k, j, i, n_jd;
  unsigned int num_words = 0; // number of words in the corpus
  unsigned int num_uwords;
  unsigned int num_test_words = 0;
  unsigned int num_train_docs = (double)num_docs * (1. - test_doc_share);

  double alpha_ss = 0;
  double eta_ss = 0;
  double doc_vi_lb;
  double doc_vi_cr;
  double doc_vi_lb_old;
  double em_conv_ratio;
  double em_lb_old;
  double em_lb_current;

  vector < vector < unsigned int > > doc_uword_ids; // unique word ids
  vector < vector < unsigned int > > doc_uword_counts; // unique word counts
  vector < vector < unsigned int > > doc_uword_test_counts; // unique test word counts
  vector < double > vi_lb; // variational lower-bound

  vec phi_ui;
  vec vi_rho_ss;
  vec vi_rho_new;
  vec perplexity;
  vec pred_likelihood;
  vec collection_size = zeros<vec>(num_collections); // number of documents in
                                                     // each collection
  rowvec collection_word_counts = zeros<rowvec>(num_collections);
  uvec doc_word_counts = zeros<uvec>(num_docs); // doc lengths
  uvec doc_train_word_counts = zeros<uvec>(num_docs); // doc lengths
  mat vi_lambda = zeros<mat>(num_topics, vocab_size); // K x V matrix
  mat vi_lambda_ss = zeros<mat>(num_topics, vocab_size); // K x V matrix
  mat vi_rho = zeros<mat>(num_topics, num_docs); // K x D matrix
  mat vi_tau = zeros<mat>(num_topics, num_collections); // K x J matrix
  mat vi_tau_ss = zeros<mat>(num_topics, num_collections); // K x J matrix
  mat vi_tau_rho_ss; // K x J matrix

  cube vi_tau_t = zeros<cube>(num_topics, num_collections, em_max_iter);

  //////////////////////////////////////////////////////////////////////////////
  // Reads corpus statistics
  //////////////////////////////////////////////////////////////////////////////

  for (d = 0; d < num_docs; d++){ // for each document
    j = docs_cid(d); // document d's collection id
    umat document = as<umat>(docs_tf(d));
    vector < unsigned int > uword_ids;
    vector < unsigned int > uword_counts;
    vector < unsigned int > uword_test_counts;
    vector < unsigned int > uword_idx;

    for (c = 0; c < document.n_cols; c++){
      word_id = document(0,c);
      word_count = document(1,c);
      num_words += word_count; // increments number of words in the corpus
      doc_word_counts(d) += word_count; // increments doc word counts
      uword_ids.push_back(word_id); // saves unique word ids
      uword_counts.push_back(word_count); // saves unique word counts
      uword_test_counts.push_back(0); // initializes to zero
      for (i = 0; i < word_count; i++){
        uword_idx.push_back(c); // the word's index in uword_test_counts
      }
    }

    // random selection of test words
    if (test_word_share > 0 && d > num_train_docs){
      n_jd = doc_word_counts(d); // the document length
      uvec rp_d = randperm(n_jd); // gets random permutations
      unsigned int num_train_words = (1. - test_word_share) * n_jd;
      for (i = num_train_words; i < n_jd; i++){
        uword_test_counts[uword_idx[rp_d(i)]] += 1; //
        num_test_words++;
      }
      doc_train_word_counts(d) = num_train_words;
    }
    else {
      doc_train_word_counts(d) = doc_word_counts(d);
    }

    doc_uword_ids.push_back(uword_ids);
    doc_uword_counts.push_back(uword_counts);
    doc_uword_test_counts.push_back(uword_test_counts);
    collection_size(j) += 1;

    collection_word_counts(j) += doc_train_word_counts(d);
  } // for each document

  if (test_word_share > 0 && test_doc_share > 0){
    perplexity = arma::zeros<arma::vec>(em_max_iter);
    pred_likelihood = arma::zeros<arma::vec>(num_test_words);
  }

  if (verbose > 0){ cout << "DONE" << endl; }

  //////////////////////////////////////////////////////////////////////////////

  if (verbose > 1){
    cout << "clda-vem (c++): Number of docs: " << num_docs << endl;
    cout << "clda-vem (c++): Number of total words: " << num_words << endl;
    cout << "clda-vem (c++): Number of test words: " << num_test_words << endl;
    cout << "clda-vem (c++): Number of topics: " << num_topics << endl;
    cout << "clda-vem (c++): Vocabulary size: " << vocab_size << endl;
    cout << "clda-vem (c++): alpha_h: " << alpha_h << endl;
    cout << "clda-vem (c++): gamma_h: " << gamma_h << endl;
    cout << "clda-vem (c++): eta_h: " << eta_h << endl;
  }

  if (verbose > 0){
    cout << "clda-vem (c++): Variational-EM..." << endl << endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Initializes variational-EM parameters
  //////////////////////////////////////////////////////////////////////////////

  // vi_lambda - the variational Dirichlet parameter for \beta
  // random initiliazation, adapted from
  //     void random_initialize_ss(lda_suffstats* ss, lda_model* model)
  // given in the LDA-C package
  vi_lambda.randu();
  vi_lambda += (1.0 / (double) vocab_size);

  // initializes the \lambda statistics
  for (k = 0; k < num_topics; k++) { // for each topic
    vi_lambda_ss.row(k) = digamma_rowvec(vi_lambda.row(k)) - Rf_digamma(sum(vi_lambda.row(k)));
  }
  assert(vi_lambda_ss.is_finite());

  // vi_tau - the variational Dirichlet parameter for \pi
//   if (init_pi.n_cols > 0){
//     vi_tau = mat(init_pi);
//   } else {
  vi_tau.fill(alpha_h);
  vi_tau.each_row() += ( collection_word_counts / (double)num_topics );
  vi_tau.each_row() /= sum(vi_tau, 0); // makes each column sums to 1
  // }


  em_iter = 0;
  em_conv_ratio = 1.0;
  em_lb_old = 0;
  cout.precision(10);
  // while (((em_conv_ratio < 0) || (em_conv_ratio > em_conv_thresh) || (em_iter <= 2)) && (em_iter <= em_max_iter)) {
  while (((em_conv_ratio < 0) || (em_conv_ratio > em_conv_thresh)) && (em_iter < em_max_iter)) {

    vi_tau_t.slice(em_iter) = vi_tau;
    em_iter++;

    ////////////////////////////////////////////////////////////////////////////
    // E Step
    ////////////////////////////////////////////////////////////////////////////

    if (verbose > 0) {
      cout << "clda-vem (c++): em_iter #" << em_iter;
      cout << " alpha: " << alpha_h << " gamma: " << gamma_h << " eta: " << eta_h << endl << endl;
    }

    // resets variational Dirichlets in each EM iteration
    // vi_rho - the variational Dirichlet parameter for \theta
    // vi_lambda - the variational Dirichlet parameter for \beta

    rowvec vi_tau_colsums = sum(vi_tau, 0); // sums elements in each column
    for (d = 0; d < num_docs; d++){
      j = docs_cid(d); // document d's collection id
      vi_rho.col(d) = (
        ((gamma_h / vi_tau_colsums(j)) * vi_tau.col(j))
        + ((double) doc_train_word_counts(d) / (double) num_topics) // + ((double) doc_word_counts(d) / (double) num_topics)
      );
    }

    mat vi_lambda_new = zeros<mat>(num_topics, vocab_size); // K x V matrix
    vi_lambda_new.fill(eta_h); // fills with the current \eta

    vi_tau_rho_ss = zeros<mat>(num_topics, num_collections);
    em_lb_current = 0.0; // corpus variational lowerbound

    for (d = 0; d < num_docs; d++) { // for each document

      j = docs_cid(d); // document d's collection id
      vector < unsigned int > uword_ids = doc_uword_ids[d];
      vector < unsigned int > uword_counts = doc_uword_counts[d];
      vector < unsigned int > uword_test_counts = doc_uword_test_counts[d];
      num_uwords = uword_ids.size();

      // variational inference for document d

      mat phi_d = zeros<mat>(num_topics, num_uwords);
      doc_vi_cr = 1.0;
      doc_vi_lb = 0;
      doc_vi_lb_old = 0;
      vi_iter = 0;

      while ((doc_vi_cr > vi_conv_thresh) && (vi_iter < vi_max_iter)) { // for each doc VI iteration

        vi_iter++;
        if (verbose > 2) {
          cout << "clda-vem (c++): doc #" << (d + 1) << " vi_iter # " << vi_iter;
        } else if (verbose > 1) {
          cout << ".";
        }

        // vi updates for \rho and \phi

        vi_rho_ss = digamma_vec(vi_rho.col(d)) - Rf_digamma(sum(vi_rho.col(d)));
        vi_rho_new = (gamma_h / sum(vi_tau.col(j))) * vi_tau.col(j);

        for (c = 0; c < num_uwords; c++){ // for each unique word
          word_id = uword_ids[c];
          word_count = uword_counts[c] - uword_test_counts[c];
          if (word_count > 0) {
            phi_ui = exp(vi_lambda_ss.col(word_id) + vi_rho_ss);
            phi_ui.elem( find(phi_ui < 1e-15) ).fill(1e-15); // this is a HACK!!!!
            phi_ui /= sum(phi_ui); // normalize to sum to 1.
            phi_d.col(c) = phi_ui;

            vi_rho_new += word_count * phi_ui;
          }
        } // for each unique word

        vi_rho.col(d) = vi_rho_new;


        // computes log document d's variational lower-bound

        vi_rho_ss = digamma_vec(vi_rho_new) - Rf_digamma(sum(vi_rho_new)); // K x 1 vector

        doc_vi_lb = lgamma(gamma_h) - ( gamma_h / vi_tau_colsums(j) ) * (num_topics - 1.) ;
        doc_vi_lb -= ( (gamma_h - num_topics) * (log(vi_tau_colsums(j)) - Rf_digamma(vi_tau_colsums(j)) + Rf_digamma(sum(vi_rho_new))) );
        doc_vi_lb -= sum(
          log_gamma_vec( (gamma_h / vi_tau_colsums(j))  * vi_tau.col(j) )
          + ( 1. - ( (gamma_h / vi_tau_colsums(j))  * vi_tau.col(j) ) )
          % ( log(vi_tau.col(j)) - digamma_vec(vi_tau.col(j)) + digamma_vec(vi_rho_new) )
        );

        doc_vi_lb += (
          sum(log_gamma_vec(vi_rho_new))
          - lgamma(sum(vi_rho_new))
          - sum( (vi_rho_new - 1.) % vi_rho_ss )
        ); //

        for (c = 0; c < num_uwords; c++) { // for each unique word
          word_id = uword_ids[c];
          word_count = uword_counts[c] - uword_test_counts[c];
          if (word_count > 0) {
            phi_ui = phi_d.col(c);
            doc_vi_lb += (sum(phi_ui % vi_rho_ss) * word_count); //
            doc_vi_lb += (sum(phi_ui % vi_lambda_ss.col(word_id)) * word_count); //
            doc_vi_lb -= (sum(phi_ui % log(phi_ui)) * word_count); //
          }
        }
        assert(!isnan(doc_vi_lb));

        // checks for convergence

        if (verbose > 2){ cout << " vi-lb: " << doc_vi_lb; }
        if (vi_iter > 1){
          doc_vi_cr = (doc_vi_lb_old - doc_vi_lb) / doc_vi_lb_old;
          if (verbose > 2){ cout << " conv-ratio: " << doc_vi_cr; }
        }
        doc_vi_lb_old = doc_vi_lb;


        if (verbose > 2){ cout << endl; }


      } // for each doc VI iteration


      em_lb_current += doc_vi_lb;

      // doc statistics to update \tau_j's
      vi_tau_rho_ss.col(j) += digamma_vec(vi_rho.col(d));

      // vi updates for \lambda

      for (c = 0; c < num_uwords; c++) {
        word_count = uword_counts[c] - uword_test_counts[c];
        if (word_count > 0) {
          vi_lambda_new.col(uword_ids[c]) += word_count * phi_d.col(c);
        }
      }

      if (verbose > 1){ cout << endl; }

    } // for each document


    if (verbose > 1){ cout << "doc-vi-lb: " << em_lb_current << endl; }

    // updates variational Dirichlet parameter \lambda, for topics (i.e. \beta_k's)

    vi_lambda = vi_lambda_new;


    ////////////////////////////////////////////////////////////////////////////
    // M Step
    ////////////////////////////////////////////////////////////////////////////

    // vi updates for each collection - \tau_j's.
    // a constraint Newton method...

    vi_tau = optimize_tau(vi_tau, vi_tau_rho_ss, collection_size,
                          num_collections, num_topics, alpha_h, gamma_h,
                          tau_max_iter, tau_step_size, verbose);


    // computes \eta sufficient statistics

    eta_ss = 0;
    for (k = 0; k < num_topics; k++) { // for each topic
      vi_lambda_ss.row(k) = digamma_rowvec(vi_lambda.row(k)) - Rf_digamma(sum(vi_lambda.row(k))); // 1 x K vector
      eta_ss += sum( vi_lambda_ss.row(k) );
    }

    // computes \alpha sufficient statistics

    alpha_ss = 0;
    for (j = 0; j < num_collections; j++){ // for each collection
      vi_tau_ss.col(j) = digamma_vec(vi_tau.col(j)) - Rf_digamma(sum(vi_tau.col(j)));
      alpha_ss += sum(vi_tau_ss.col(j));
    }


    // hyperparameter optimization for \alpha and \eta

    if(estimate_alpha){
      // alpha_h = opt_hp(100, alpha_ss, num_collections, num_topics); // resets \alpha each iteration
      alpha_h = opt_hp(alpha_h, alpha_ss, num_collections, num_topics); // initializes with the prior \alpha
      if (verbose > 1){ cout << endl; }
    }

    if(estimate_eta){
      // eta_h = opt_hp(100, eta_ss, num_topics, vocab_size); // resets \eta each iteration
      eta_h = opt_hp(eta_h, eta_ss, num_topics, vocab_size); // initializes with the prior \eta
      if (verbose > 1){ cout << endl; }
    }

    // computes the em-lower-bound based on the updates \alpha and \eta

    for (k = 0; k < num_topics; k++) { // for each topic
      em_lb_current += (
        lgamma(vocab_size * eta_h) -
          vocab_size * lgamma(eta_h) +
          sum( (eta_h - 1.) * vi_lambda_ss.row(k) )
      );
      em_lb_current += (
        sum( log_gamma_rowvec( vi_lambda.row(k) ) ) -
          lgamma( sum( vi_lambda.row(k) ) ) -
          sum( (vi_lambda.row(k) - 1.) % vi_lambda_ss.row(k) )
      );
    }

    for (j = 0; j < num_collections; j++){ // for each collection
      em_lb_current += (
        lgamma(num_topics * alpha_h) -
          num_topics * lgamma(alpha_h) +
          sum( (alpha_h - 1.) * vi_tau_ss.col(j) )
      );
      em_lb_current += (
        sum( log_gamma_vec( vi_tau.col(j) ) ) -
          lgamma( sum( vi_tau.col(j) ) ) -
          sum( (vi_tau.col(j) - 1.) % vi_tau_ss.col(j) )
      );
    }

    if (verbose > 1){
      cout << "clda-vem (c++): em_iter #" << em_iter;
      cout << " vi-lb: " << em_lb_current;
      cout << " alpha_ss: " << alpha_ss;
      cout << " eta_ss: " << eta_ss;
      cout << " opt-alpha: " << alpha_h;
      cout << " opt-gamma: " << gamma_h;
      cout << " opt-eta: " << eta_h;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Estimates of the predictive likelihood
    // p(w^{\text{test}}_{jdi} | \bw^{\text{train}})
    // for each test word w^{\text{test}}_{jdi} in the corpus
    ////////////////////////////////////////////////////////////////////////////
    if (test_word_share > 0 && test_doc_share > 0){ // perplexity
      mat beta_hat = vi_lambda;
      beta_hat.each_col() /= sum(beta_hat, 1); // normalizes the rows to sum to 1
      assert(beta_hat.is_finite());

      unsigned int test_word_idx = 0;
      for (d = 0; d < num_docs; d++) { // for each document

        j = docs_cid(d); // document d's collection id
        vector < unsigned int > uword_ids = doc_uword_ids[d];
        vector < unsigned int > uword_test_counts = doc_uword_test_counts[d];
        vec theta_d_hat = vi_rho.col(d);
        theta_d_hat /= sum(theta_d_hat); // normalizes the column to sum to 1
        assert(theta_d_hat.is_finite());

        for (c = 0; c < uword_ids.size(); c++) { // for each unique word
          word_id = uword_ids[c];
          for (unsigned int twc = 0; twc < uword_test_counts[c]; twc++){ // only for test words
            double pll_t = sum(theta_d_hat % beta_hat.col(word_id));
            pred_likelihood(test_word_idx) = (
              ( (em_iter - 1) * pred_likelihood(test_word_idx) + pll_t )
              / em_iter
            ); // calculates the predictive likelihood via online average
            test_word_idx++;
          }
        }

      }

      // perplexity of the test set for the current iteration
      perplexity(em_iter - 1) = exp(-mean(log(pred_likelihood)));
      if (verbose > 1){
        cout << " perp: " << perplexity(em_iter - 1);
      }
    } // perplexity


    ////////////////////////////////////////////////////////////////////////////
    // EM: check for convergence
    ////////////////////////////////////////////////////////////////////////////




    if (em_iter > 1){
      em_conv_ratio = (em_lb_old - em_lb_current) / em_lb_old;
      if (verbose > 1){ cout << " em-conv-ratio: " << em_conv_ratio; }
      // if (em_conv_ratio < 0.0) { vi_max_iter *= 2; }
    }
    if (verbose > 1){  cout << endl << endl; }

    em_lb_old = em_lb_current;
    vi_lb.push_back(em_lb_current);


  } // for each EM iteration



  if (verbose > 0){
    cout << "clda-vem (c++): Completed Variational-EM." << endl << endl;
  }

  return List::create(
    Named("vi_lambda") = wrap(vi_lambda),
    Named("vi_tau") = wrap(vi_tau),
    Named("vi_tau_t") = wrap(vi_tau_t),
    Named("vi_rho") = wrap(vi_rho),
    Named("alpha_h") = wrap(alpha_h),
    Named("gamma_h") = wrap(gamma_h),
    Named("eta_h") = wrap(eta_h),
    Named("perplexity") = wrap(perplexity)
  );

}

