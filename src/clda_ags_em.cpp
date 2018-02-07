// /////////////////////////////////////////////////////////////////////////////
//
// cLDA: Auxiliary Variable Update within Collpased Gibbs Sampler with
// hyperparameter alpha sampling and EM updates for hyperparameters eta and
// gamma
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

//' cLDA: Auxiliary Variable Update within Collpased Gibbs Sampler with
//' hyperparameter \eqn{\alpha} sampling  and EM updates for hyperparameters
//' \eqn{\eta} and \eqn{\gamma}
//'
//' This implements a Markov chain on \eqn{(z, \pi)} via the collapsed Gibbs
//' sampling with auxiliary variable updates for the compound latent Dirichlet
//' allocation (cLDA) model.
//'
//' To compute perplexity, we first partition words in a corpus into two sets:
//' (a) a test set (held-out set), which is selected from the set of words in
//' the test (held-out) documents (identified via \code{test_doc_share} and
//' \code{test_word_share}) and (b) a training set, i.e., the remaining words in
//' the corpus. We then run the variational EM algorithm based on the training
//' set. Finally, we compute per-word perplexity based on the held-out set.
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_cid Collection ID for each document in the corpus (indices starts 0)
//' @param docs_tf Corpus documents read from the Blei corpus format, e.g., via \code{\link{read_docs}} (indices starts with 0)
//' @param alpha_h Hyperparameter for \eqn{\pi}. When \code{sample_alpha_h} is \code{true} this variable is used to initialize hyperparameter \eqn{\alpha}
//' @param gamma_h Hyperparameter for \eqn{\theta}
//' @param eta_h Hyperparameter for \eqn{\beta}
//' @param em_max_iter Maximum number of EM iterations to be performed
//' @param gibbs_max_iter Maximum number of Gibbs iterations to be performed
//' @param burn_in Burn-in-period for the Gibbs sampler
//' @param spacing Spacing between the stored samples (to reduce correlation)
//' @param save_pi if 0 the function does not save \eqn{\pi} samples
//' @param save_beta if 0 the function does not save \eqn{\beta} samples
//' @param save_theta if 0 the function does not save \eqn{\theta} samples
//' @param save_lp if 0 the function does not save computed log posterior for iterations
//' @param verbose from {0, 1, 2}
//' @param init_pi the initial configuration for the collection level topic mixtures, i.e., \eqn{\pi} samples
//' @param test_doc_share proportion of the test documents in the corpus. Must be from [0., 1.)
//' @param test_word_share proportion of the test words in each test document. Must be from [0., 1.)
//' @param burn_in_pi burn in iterations until pi sampling
//' @param sample_alpha_h sample hyperparameter \eqn{\alpha} (true) or not (false)
//' @param gamma_shape hyperparameter \code{shape} for the Gamma prior on \eqn{\alpha}. Default is 1.
//' @param gamma_rate hyperparameter \code{rate} for the Gamma prior on \eqn{\alpha}. Default is 1.
//'
//' @return A list of
//'   \item{corpus_topic_counts}{corpus-level topic counts from last iteration of the Markov chain}
//'   \item{pi_counts}{collection-level topic counts from the last iteration of the Markov chain}
//'   \item{theta_counts}{document-level topic counts from last iteration of the Markov chain}
//'   \item{beta_counts}{topic word counts from last iteration of the Markov chain}
//'   \item{pi_samples}{\eqn{\pi} samples after the burn in period, if \code{save_pi} is set}
//'   \item{theta_samples}{\eqn{\theta} samples after the burn in period, if \code{save_theta} is set}
//'   \item{beta_samples}{\eqn{\beta} samples after the burn in period, if \code{save_beta} is set}
//'   \item{log_posterior}{the log posterior (upto a constant multiplier) of the hidden variable \eqn{\psi = (\beta, \pi, \theta, z)} in the LDA model, if \code{save_lp} is set}
//'   \item{log_posterior_pi_z}{the log posterior (upto a constant multiplier) of the hidden variables \eqn{(\pi, z)} in the LDA model, if \code{save_lp} is set}
//'   \item{perplexity}{perplexity of the set of held-out words}
//'   \item{alpha_h_samples}{\eqn{\alpha} samples if \code{sample_alpha_h} is \code{true}}
//'   \item{gamma_h_estimates}{\eqn{\gamma} estimates from each EM iteration}
//'   \item{eta_h_estimates}{\eqn{\eta} estimates from each EM iteration}
//'
//'
//' @export
//'
//' @family MCMC
//'
//' @note
//'
//' Updated on: December 18, 2017 -- Added hyperparameter alpha sampling and AGS EM updates
//'
//' Updated on: June 02, 2016
//'
//' Created on: May 18, 2016
//'
//' Created by: Clint P. George
//'
// [[Rcpp::export]]
List clda_ags_em(
    unsigned int num_topics,
    unsigned int vocab_size,
    NumericVector docs_cid,
    List docs_tf,
    double alpha_h,
    double gamma_h,
    double eta_h,
    unsigned int em_max_iter,
    unsigned int gibbs_max_iter,
    unsigned int burn_in,
    unsigned int spacing,
    bool save_pi,
    bool save_theta,
    bool save_beta,
    bool save_lp,
    int verbose,
    arma::mat init_pi, // TODO: to be deleted in the final version
    double test_doc_share = 0.,
    double test_word_share = 0.,
    unsigned int burn_in_pi = 10,
    bool sample_alpha_h = false,
    double gamma_shape = 1.,
    double gamma_rate = 1.) {

  assert(test_word_share < 1.);
  assert(test_word_share >= 0.);
  assert(test_doc_share < 1.);
  assert(test_doc_share >= 0.);
  assert(gamma_shape > 0.);
  assert(gamma_rate > 0.);

  burn_in_pi = (burn_in_pi >= gibbs_max_iter) ? 0 : burn_in_pi;

  NumericVector cids = unique(docs_cid);
  unsigned int num_collections = cids.size(); // number of collections
  unsigned int num_docs = docs_cid.size(); // number of documents
  unsigned int valid_samples = ceil((gibbs_max_iter - burn_in) / (double) spacing);
  unsigned int d, i, j, k, iter, c, word_id, word_count, topic_id, new_topic_id, v;
  unsigned int n_jd; // number of words in document d_j
  unsigned int num_words = 0; // number of words in the corpus
  unsigned int num_train_docs = (double)num_docs * (1. - test_doc_share);

  double doc_denom;

  vector < vector < unsigned int > > doc_word_ids;
  vector < vector < unsigned int > > doc_word_zids;
  vector < vector < unsigned int > > doc_word_class;

  arma::vec prob;
  arma::uvec num_accept = arma::zeros<arma::uvec>(num_collections);
  arma::uvec doc_word_counts = arma::zeros<arma::uvec>(num_docs); // doc lengths
  uvec doc_train_word_counts = zeros<uvec>(num_docs); // doc lengths, for training
  arma::mat beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  arma::vec corpus_topic_counts = arma::zeros<arma::vec>(num_topics); // corpus-level topic counts
  arma::mat pi_counts = arma::zeros<arma::mat>(num_topics, num_collections); // K x J matrix
  arma::mat pi_t = arma::zeros<arma::mat>(num_topics, num_collections); // K x J matrix
  arma::mat sterling_counts; // K x D matrix

  arma::vec alpha_h_samples;
  arma::cube pi_samples;
  arma::cube theta_samples;
  arma::cube beta_samples;
  arma::vec log_posterior;
  arma::vec log_posterior_pi_z;
  arma::vec perplexity;
  arma::vec perplexity2;
  arma::vec pred_likelihood;
  arma::vec pred_likelihood2;
  vector < double > gamma_h_estimates;
  vector < double > eta_h_estimates;

  arma::mat PI_PRIOR(num_topics, num_collections); // K x J matrix
  PI_PRIOR.fill(alpha_h - 1.);

  if (sample_alpha_h) {
    alpha_h_samples = arma::zeros<arma::vec>(valid_samples);
  }
  if (save_pi) {
    pi_samples = arma::cube(num_topics, num_collections, valid_samples);
  }
  else {
    pi_samples = arma::cube(num_topics, num_collections, 1); // to save the last sample
  }
  if (save_theta) {
    theta_samples = arma::cube(num_topics, num_docs, valid_samples);
  }
  else {
    theta_samples = arma::cube(num_topics, num_docs, 1); // to save the last sample
  }
  if (save_beta) {
    beta_samples = arma::cube(num_topics, vocab_size, valid_samples);
  }
  else {
    beta_samples = arma::cube(num_topics, vocab_size, 1); // to save the last sample
  }
  if (save_lp) {
    log_posterior = arma::zeros<arma::vec>(valid_samples);
    log_posterior_pi_z = arma::zeros<arma::vec>(valid_samples);
  }


  cout << endl << endl;
  if (verbose > 1){
    cout << "clda-ags-em (c++): Number of saved samples - "
         << valid_samples << endl;
  }
  // Calculates the document word indices

  if (verbose > 0){
    cout << "clda-ags-em (c++): Initializes variables and count statistics....";
  }

  // corpus level topic mixture is used to initialize count statistics
  arma::vec alpha_vec = arma::zeros<arma::vec>(num_topics);
  alpha_vec.fill(alpha_h);

  //   for (j = 0; j < num_collections; j++){
  //     pi_t.col(j) = sample_dirichlet(num_topics, alpha_vec);
  //   }
  pi_t = arma::mat(init_pi);


  for (d = 0; d < num_docs; d++){

    j = docs_cid(d); // document d's collection id
    assert(j < num_collections);
    arma::umat document = as<arma::umat>(docs_tf(d));
    vector < unsigned int > word_ids;
    vector < unsigned int > word_zids;
    vector < unsigned int > word_class;

    for (c = 0; c < document.n_cols; c++){
      word_id = document(0,c);
      word_count = document(1,c);

      for (i = 0; i < word_count; i++){
        // samples z for each word
        arma::vec theta_c = sample_dirichlet(num_topics, alpha_vec);
        topic_id = sample_multinomial(theta_c);

        word_zids.push_back(topic_id);
        word_ids.push_back(word_id);
        word_class.push_back(0); // train word
        num_words++; // increments number of words in the corpus
      }

      doc_word_counts(d) += word_count; // increments doc word counts
    }

    // random selection of test words
    if (test_word_share > 0 && d > num_train_docs){
      n_jd = doc_word_counts(d); // the document length
      arma::uvec rp_d = randperm(n_jd); // gets random permutations
      unsigned int num_train_words = (1. - test_word_share) * n_jd;
      for (i = num_train_words; i < n_jd; i++){
        word_class[rp_d(i)] = 1; // test word
      }
      doc_train_word_counts(d) = num_train_words;
    } else {
      doc_train_word_counts(d) = doc_word_counts(d);
    }


    // doc_word_indices.push_back(word_indices);
    doc_word_ids.push_back(word_ids);
    doc_word_zids.push_back(word_zids);
    doc_word_class.push_back(word_class);

  }

  //////////////////////////////////////////////////////////////////////////////
  // updates count statistics for training words
  //////////////////////////////////////////////////////////////////////////////
  unsigned int num_test_words = 0;
  for (d = 0; d < num_docs; d++){

    j = docs_cid(d); // document d's collection id

    for (i = 0; i < doc_word_counts(d); i++){
      topic_id = doc_word_zids[d][i];
      word_id = doc_word_ids[d][i];
      if (doc_word_class[d][i] == 0){ // train words
        pi_counts(topic_id, j) += 1.;
        corpus_topic_counts(topic_id) += 1.;
        theta_counts(topic_id, d) += 1.;
        beta_counts(topic_id, word_id) += 1.;
      }
      else { // test words
        num_test_words++;
      }
    }

  }
  //////////////////////////////////////////////////////////////////////////////

  if (test_word_share > 0 && test_doc_share > 0){
    perplexity = arma::zeros<arma::vec>(valid_samples);
    perplexity2 = arma::zeros<arma::vec>(valid_samples);
    pred_likelihood = arma::zeros<arma::vec>(num_test_words);
    pred_likelihood2 = arma::zeros<arma::vec>(num_test_words);
  }

  if (verbose > 0){
    cout << "DONE" << endl;
  }

  if (verbose > 1){
    cout << "clda-ags-em (c++): Number of collections: " << num_collections << endl;
    cout << "clda-ags-em (c++): Number of docs: " << num_docs << endl;
    cout << "clda-ags-em (c++): Number of total words: " << num_words << endl;
    cout << "clda-ags-em (c++): Number of test words: " << num_test_words << endl;
    cout << "clda-ags-em (c++): Number of topics: " << num_topics << endl;
    cout << "clda-ags-em (c++): Vocabulary size: " << vocab_size << endl;
    cout << "clda-ags-em (c++): alpha_h: " << alpha_h << endl;
    cout << "clda-ags-em (c++): gamma_h: " << gamma_h << endl;
    cout << "clda-ags-em (c++): eta_h: " << eta_h << endl;
    cout << "clda-ags-em (c++): burn_in_pi: " << burn_in_pi << endl;
    cout << "clda-ags-em (c++): gibbs_max_iter: " << gibbs_max_iter << endl;
    cout << "clda-ags-em (c++): valid_samples: " << valid_samples << endl;
  }


  //////////////////////////////////////////////////////////////////////////////
  // EM iterations
  //////////////////////////////////////////////////////////////////////////////

  unsigned int em_iter = 0;
  double gamma_num_avg = 0;
  double gamma_den_avg = 0;
  double eta_num_avg = 0;
  double eta_den_avg = 0;
  double eta_num_sum = 0;
  double eta_den_sum = 0;
  double gamma_num_sum = 0;
  double gamma_den_sum = 0;
  double gamma_h_old;
  double eta_h_old;
  cout.precision(10);

  while ((em_iter <= 2) || (em_iter <= em_max_iter)) { // for EM iteration

    em_iter++;

    ////////////////////////////////////////////////////////////////////////////
    // E Step: Gibbs sampling
    ////////////////////////////////////////////////////////////////////////////


    if (verbose > 0){
      cout << "clda-ags-em (c++): EM iter #" << em_iter << endl << endl;
      cout << "clda-ags-em (c++): Augmented Gibbs sampling..." << endl << endl;
    }
    unsigned int ss_idx = 0;


    for (iter = 0; iter < gibbs_max_iter; iter++) { // for each Gibbs iteration

      if (verbose > 1) { cout << "clda-ags-em (c++): AGS iter# " << iter + 1; }
      else if (verbose > 0) { cout << "."; }

      bool save_flag = (iter >= burn_in) && (iter % spacing == 0);

      if (save_flag && save_pi) {
        pi_samples.slice(ss_idx) = pi_t;
      }
      else if (save_flag && ((iter + 1) == gibbs_max_iter)) {
        pi_samples.slice(0) = pi_t;
      }

      // samples \beta
      if (save_flag && save_beta) {
        for(k = 0; k < num_topics; k++)
          beta_samples.slice(ss_idx).row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);
      }
      else if (save_flag && ((iter + 1) == gibbs_max_iter)) {
        for(k = 0; k < num_topics; k++)
          beta_samples.slice(0).row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);
      }


      ////////////////////////////////////////////////////////////////////////////
      // updates z's
      ////////////////////////////////////////////////////////////////////////////

      gamma_num_sum = 0; // resets the sum in every Gibbs iteration
      gamma_den_sum = 0; // resets the sum in every Gibbs iteration

      for (d = 0; d < num_docs; d++) { // for each document

        j = docs_cid(d); // document d's collection id
        // assert(j < num_collections); // check
        n_jd = doc_word_counts(d); // number of words in document d
        doc_denom = doc_train_word_counts(d) - 1. + gamma_h; // it's a constant for a term
        vector < unsigned int > word_ids = doc_word_ids[d];
        vector < unsigned int > word_zids = doc_word_zids[d];
        vector < unsigned int > word_class = doc_word_class[d];
        vec pi_j = pi_t.col(j);

        // samples \theta
        if (save_flag && save_theta) {
          theta_samples.slice(ss_idx).col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + gamma_h * pi_j);
        }
        else if (save_flag && ((iter + 1) == gibbs_max_iter)) {
          theta_samples.slice(0).col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + gamma_h * pi_j);
        }


        for (i = 0; i < n_jd; i++) { // for each word

          if (word_class[i]) continue; // ignores test words from sampling

          topic_id = word_zids[i];
          word_id = word_ids[i];
          prob = arma::zeros <arma::vec> (num_topics); // initialize with zero

          // decrements the counts by one, to ignore the current sampling word
          theta_counts(topic_id, d) -= 1.;
          beta_counts(topic_id, word_id) -= 1.;
          pi_counts(topic_id, j) -= 1.;
          corpus_topic_counts(topic_id) -= 1.;

          // samples z's
          // TODO: need to check whether we can vectorize this for loop.

          for (k = 0; k < num_topics; k++){ // for each topic
            // compute p(z_{jdi} == j | \bz^{(-jdi)}, \bw, \pi)
            prob(k) = (((theta_counts(k, d) + gamma_h * pi_j(k)) / doc_denom)
                         * ((beta_counts(k, word_id) + eta_h)
                              / (corpus_topic_counts(k) + eta_h * vocab_size)));
          }
          new_topic_id = sample_multinomial(prob); // the ***new*** topic

          // assert(new_topic_id < num_topics); // check

          // increments the counts by one
          theta_counts(new_topic_id, d) += 1.;
          beta_counts(new_topic_id, word_id) += 1.;
          pi_counts(new_topic_id, j) += 1.;
          corpus_topic_counts(new_topic_id) += 1.;

          // updates newly generated topic to the database
          word_zids[i] = new_topic_id;

        } // for each word

        doc_word_zids[d] = word_zids; // updates global variable


        // Compute hyperparameter \gamma statistics

        gamma_num_sum += sum(pi_j % (digamma_vec(theta_counts.col(d) + gamma_h * pi_j) - digamma_vec(gamma_h * pi_j)));
        gamma_den_sum += Rf_digamma(doc_word_counts(d) + gamma_h) - Rf_digamma(gamma_h);

      } // for each document

      //////////////////////////////////////////////////////////////////////////
      // Computes \gamma and \eta statistics
      //////////////////////////////////////////////////////////////////////////
      if (save_flag){
        gamma_num_avg = ((ss_idx * gamma_num_avg + (gamma_num_sum / (double) num_docs)) / (ss_idx + 1.));
        gamma_den_avg = ((ss_idx * gamma_den_avg + (gamma_den_sum / (double) num_docs)) / (ss_idx + 1.));

        eta_num_sum = 0;
        eta_den_sum = 0;
        for(k = 0; k < num_topics; k++){
          eta_num_sum += sum(digamma_rowvec(beta_counts.row(k) + eta_h) - Rf_digamma(eta_h)); // sum_{v = 1}^V - numerator
          eta_den_sum += vocab_size * (Rf_digamma(corpus_topic_counts(k) + vocab_size * eta_h) - Rf_digamma(vocab_size * eta_h)); // denominator sum
        }
        eta_num_avg = ((ss_idx * eta_num_avg + (eta_num_sum / (double) num_topics)) / (ss_idx + 1.));
        eta_den_avg = ((ss_idx * eta_den_avg + (eta_den_sum / (double) num_topics)) / (ss_idx + 1.));
      }



      ////////////////////////////////////////////////////////////////////////////
      // updates \pi's
      //
      // auxiliary variable update
      ////////////////////////////////////////////////////////////////////////////

      if (iter > burn_in_pi) { // updates \pi's

        // Computes the Stirling number of the first kind via sampling from the
        // Antoniak distribution
        sterling_counts = arma::zeros<arma::mat>(num_topics, num_collections); // K x J matrix
        for (d = 0; d < num_docs; d++) { // for each document
          j = docs_cid(d); // document d's collection id
          for (k = 0; k < num_topics; k++){
            sterling_counts(k, j) += sample_antoniak(
              theta_counts(k, d), // the number of samples
              gamma_h * pi_t(k, j) // the strength parameter
            );
          }
        }

        // Sample \pi_j's from Diriclets
        for (j = 0; j < num_collections; j++){
          pi_t.col(j) = sample_dirichlet(num_topics, sterling_counts.col(j) + alpha_h);
        }

      } // updates \pi's


      ////////////////////////////////////////////////////////////////////////////
      // samples alpha_h from the Gamma posterior
      //
      // Added on: December 17, 2017
      ////////////////////////////////////////////////////////////////////////////

      if (sample_alpha_h) {
        double post_gamma_rate = gamma_rate - arma::accu(arma::log(pi_t));
        alpha_h = rgamma(1, gamma_shape, 1. / post_gamma_rate)(0); // rgamma(n, shape, scale)
        if (save_flag) {
          alpha_h_samples(ss_idx) = alpha_h;
        }
        if (verbose > 1){
          cout << " sum_log_pi: " << arma::accu(arma::log(pi_t));
          cout << " post_gamma_rate: " << post_gamma_rate;
          cout << " alpha_h: " << alpha_h;
        }
      }


      ////////////////////////////////////////////////////////////////////////////
      // Computes the estimate of the predictive likelihood
      // p(w^{\text{test}}_{jdi} | \bw^{\text{train}})
      // for each test word w^{\text{test}}_{jdi} in the corpus
      ////////////////////////////////////////////////////////////////////////////
      if (save_flag && (test_word_share > 0 && test_doc_share > 0)){

        arma::mat beta_hat = beta_counts + eta_h;
        beta_hat.each_col() /= (corpus_topic_counts + eta_h * vocab_size); // arma::sum(beta_hat, 1); // 1 x K vector

        unsigned int test_word_idx = 0;
        for (d = 0; d < num_docs; d++) { // for each document

          j = docs_cid(d); // document d's collection id
          n_jd = doc_word_counts(d); // number of words in document d
          vector < unsigned int > word_ids = doc_word_ids[d];
          vector < unsigned int > word_class = doc_word_class[d];

          arma::vec theta_d_hat = theta_counts.col(d) + gamma_h * pi_t.col(j);
          theta_d_hat /= (doc_train_word_counts(d) + gamma_h);

          for (i = 0; i < n_jd; i++) { // for each word
            if (word_class[i] == 0) continue; // only for test words

            // Calculates the predictive likelihood via online average
            // Method-1
            double pll_t = arma::sum(theta_d_hat % beta_hat.col(word_ids[i]));
            pred_likelihood(test_word_idx) = ((ss_idx * pred_likelihood(test_word_idx) + pll_t) / (ss_idx + 1.));
            // Method-2
            double pll_t2 = arma::sum(pi_t.col(j) % beta_hat.col(word_ids[i]));
            pred_likelihood2(test_word_idx) = ((ss_idx * pred_likelihood2(test_word_idx) + pll_t2) / (ss_idx + 1.));

            test_word_idx++;
          }

        }

        // perplexity of the test set for the current iteration
        perplexity(ss_idx) = exp(-arma::mean(arma::log(pred_likelihood)));
        perplexity2(ss_idx) = exp(-arma::mean(arma::log(pred_likelihood2)));
        if (verbose > 1){
          cout << " perp: " << perplexity(ss_idx);
          cout << " perp2: " << perplexity2(ss_idx);
        }

      }

      ////////////////////////////////////////////////////////////////////////////
      // Computes the log posterior of (\beta, \pi, \theta, z) in the c-LDA model
      // up to a normalizing constant, based on the updated z's and \pi's
      ////////////////////////////////////////////////////////////////////////////
      if (save_flag && save_beta && save_theta && save_lp){

        double lp = arma::accu((beta_counts + eta_h - 1.0) % log(beta_samples.slice(ss_idx)));
        lp += arma::accu((alpha_h - 1.0) * log(pi_t));
        for (d = 0; d < num_docs; d++) { // \sum_{j=1}^J \sum_{k=1}^K
          j = docs_cid(d); // document d's collection id
          // Note: theta_samples can be zero some times and taking log will give
          // NAN or INF. To avoid that case adding a small constant 1e-24
          arma::vec log_theta_d = log(theta_samples.slice(ss_idx).col(d) + 1e-24);
          lp += arma::accu(
            (theta_counts.col(d) + gamma_h * pi_t.col(j)  - 1.0) % log_theta_d -
              log_gamma_vec(gamma_h * pi_t.col(j))
          );
        }

        log_posterior(ss_idx) = lp;

        if (verbose > 1){ cout << " lp: " << lp; }
      }

      ////////////////////////////////////////////////////////////////////////////
      // Computes the log posterior of (\pi, z) in the c-LDA model up to a
      // normalizing constant, based on the updated z's and \pi's
      ////////////////////////////////////////////////////////////////////////////
      if (save_flag && save_lp){
        double lppi = arma::accu(PI_PRIOR % log(pi_t));
        for (d = 0; d < num_docs; d++) { // \sum_{j=1}^J \sum_{k=1}^K
          j = docs_cid(d); // document d's collection id
          lppi += (lgamma(gamma_h) - lgamma(gamma_h + doc_train_word_counts(d)));
          lppi += arma::sum(
            log_gamma_vec(theta_counts.col(d) + gamma_h * pi_t.col(j)) -
              log_gamma_vec(gamma_h * pi_t.col(j))
          ); // sum_{k=1}^K
        }
        lppi += num_topics * (lgamma(vocab_size * eta_h) - vocab_size * lgamma(eta_h));
        for (k = 0; k < num_topics; k++){
          for (v = 0; v < vocab_size; v++){
            lppi += lgamma(beta_counts(k, v) + eta_h);
          }
          lppi -= lgamma(corpus_topic_counts(k) + vocab_size * eta_h);
        }
        log_posterior_pi_z(ss_idx) = lppi;

        if (verbose > 1){ cout << " lp-pi-z: " << lppi; }
      }

      ////////////////////////////////////////////////////////////////////////////


      if (verbose > 1) { cout << endl; }

      if (save_flag) { ss_idx += 1; }

    } // for each Gibbs iteration

    if (verbose > 0){
      cout << endl;
      cout << "clda-ags-em (c++): Completed sampling." << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    // M Step: Hyperparameter Optimization
    ////////////////////////////////////////////////////////////////////////////

    // new alpha
    gamma_h_old = gamma_h;
    gamma_h = gamma_h * gamma_num_avg / gamma_den_avg;

    // new eta
    eta_h_old = eta_h;
    eta_h = eta_h * eta_num_avg / eta_den_avg;

    gamma_h_estimates.push_back(gamma_h);
    eta_h_estimates.push_back(eta_h);

    cout << "clda-ags-em (c++): new gamma: " << gamma_h << " old-alpha: " << gamma_h_old;
    cout << " gamma_num_avg: " << gamma_num_avg << " gamma_den_avg: " << gamma_den_avg << endl;
    cout << "clda-ags-em (c++): new eta: " << eta_h << " old-eta: " << eta_h_old;
    cout << " eta_num_avg: " << eta_num_avg << " eta_den_avg: " << eta_den_avg << endl;
    cout << endl;

  } // for EM iteration


  return List::create(
    Named("corpus_topic_counts") = wrap(corpus_topic_counts),
    Named("pi_counts") = wrap(pi_counts),
    Named("theta_counts") = wrap(theta_counts),
    Named("beta_counts") = wrap(beta_counts),
    Named("num_accept") = wrap(num_accept),
    Named("beta_samples") = wrap(beta_samples),
    Named("pi_samples") = wrap(pi_samples),
    Named("theta_samples") = wrap(theta_samples),
    Named("log_posterior") = wrap(log_posterior),
    Named("log_posterior_pi_z") = wrap(log_posterior_pi_z),
    Named("perplexity") = wrap(perplexity),
    Named("perplexity2") = wrap(perplexity2), // not used.
    Named("alpha_h_samples") = wrap(alpha_h_samples),
    Named("gamma_h_estimates") = wrap(gamma_h_estimates),
    Named("eta_h_estimates") = wrap(eta_h_estimates)
  );

}

