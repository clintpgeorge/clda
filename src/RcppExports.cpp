// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// clda_ags
List clda_ags(unsigned int num_topics, unsigned int vocab_size, NumericVector docs_cid, List docs_tf, double alpha_h, double gamma_h, double eta_h, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_pi, bool save_theta, bool save_beta, bool save_lp, int verbose, arma::mat init_pi, double test_doc_share, double test_word_share, unsigned int burn_in_pi);
RcppExport SEXP clda_clda_ags(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_cidSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP gamma_hSEXP, SEXP eta_hSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_piSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP init_piSEXP, SEXP test_doc_shareSEXP, SEXP test_word_shareSEXP, SEXP burn_in_piSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type docs_cid(docs_cidSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_h(gamma_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_pi(save_piSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init_pi(init_piSEXP);
    Rcpp::traits::input_parameter< double >::type test_doc_share(test_doc_shareSEXP);
    Rcpp::traits::input_parameter< double >::type test_word_share(test_word_shareSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in_pi(burn_in_piSEXP);
    __result = Rcpp::wrap(clda_ags(num_topics, vocab_size, docs_cid, docs_tf, alpha_h, gamma_h, eta_h, max_iter, burn_in, spacing, save_pi, save_theta, save_beta, save_lp, verbose, init_pi, test_doc_share, test_word_share, burn_in_pi));
    return __result;
END_RCPP
}
// clda_mgs
List clda_mgs(unsigned int num_topics, unsigned int vocab_size, NumericVector docs_cid, List docs_tf, double alpha_h, double gamma_h, double eta_h, double step_size, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_pi, bool save_theta, bool save_beta, bool save_lp, int verbose, arma::mat init_pi, double test_doc_share, double test_word_share, unsigned int burn_in_pi);
RcppExport SEXP clda_clda_mgs(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_cidSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP gamma_hSEXP, SEXP eta_hSEXP, SEXP step_sizeSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_piSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP init_piSEXP, SEXP test_doc_shareSEXP, SEXP test_word_shareSEXP, SEXP burn_in_piSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type docs_cid(docs_cidSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_h(gamma_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_pi(save_piSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init_pi(init_piSEXP);
    Rcpp::traits::input_parameter< double >::type test_doc_share(test_doc_shareSEXP);
    Rcpp::traits::input_parameter< double >::type test_word_share(test_word_shareSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in_pi(burn_in_piSEXP);
    __result = Rcpp::wrap(clda_mgs(num_topics, vocab_size, docs_cid, docs_tf, alpha_h, gamma_h, eta_h, step_size, max_iter, burn_in, spacing, save_pi, save_theta, save_beta, save_lp, verbose, init_pi, test_doc_share, test_word_share, burn_in_pi));
    return __result;
END_RCPP
}
// clda_vem
List clda_vem(unsigned int num_topics, unsigned int vocab_size, NumericVector docs_cid, List docs_tf, double alpha_h, double gamma_h, double eta_h, unsigned int vi_max_iter, unsigned int em_max_iter, double vi_conv_thresh, double em_conv_thresh, unsigned int tau_max_iter, double tau_step_size, bool estimate_alpha, bool estimate_gamma, bool estimate_eta, int verbose, arma::mat init_pi, double test_doc_share, double test_word_share);
RcppExport SEXP clda_clda_vem(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_cidSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP gamma_hSEXP, SEXP eta_hSEXP, SEXP vi_max_iterSEXP, SEXP em_max_iterSEXP, SEXP vi_conv_threshSEXP, SEXP em_conv_threshSEXP, SEXP tau_max_iterSEXP, SEXP tau_step_sizeSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP estimate_etaSEXP, SEXP verboseSEXP, SEXP init_piSEXP, SEXP test_doc_shareSEXP, SEXP test_word_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type docs_cid(docs_cidSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_h(gamma_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vi_max_iter(vi_max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type em_max_iter(em_max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type vi_conv_thresh(vi_conv_threshSEXP);
    Rcpp::traits::input_parameter< double >::type em_conv_thresh(em_conv_threshSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type tau_max_iter(tau_max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tau_step_size(tau_step_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_eta(estimate_etaSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init_pi(init_piSEXP);
    Rcpp::traits::input_parameter< double >::type test_doc_share(test_doc_shareSEXP);
    Rcpp::traits::input_parameter< double >::type test_word_share(test_word_shareSEXP);
    __result = Rcpp::wrap(clda_vem(num_topics, vocab_size, docs_cid, docs_tf, alpha_h, gamma_h, eta_h, vi_max_iter, em_max_iter, vi_conv_thresh, em_conv_thresh, tau_max_iter, tau_step_size, estimate_alpha, estimate_gamma, estimate_eta, verbose, init_pi, test_doc_share, test_word_share));
    return __result;
END_RCPP
}
// lda_cgs
List lda_cgs(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_theta, bool save_beta, bool save_lp, int verbose, double test_doc_share, double test_word_share);
RcppExport SEXP clda_lda_cgs(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP test_doc_shareSEXP, SEXP test_word_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type test_doc_share(test_doc_shareSEXP);
    Rcpp::traits::input_parameter< double >::type test_word_share(test_word_shareSEXP);
    __result = Rcpp::wrap(lda_cgs(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in, spacing, save_theta, save_beta, save_lp, verbose, test_doc_share, test_word_share));
    return __result;
END_RCPP
}
// sample_antoniak
double sample_antoniak(unsigned int N, double alpha);
RcppExport SEXP clda_sample_antoniak(SEXP NSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(sample_antoniak(N, alpha));
    return __result;
END_RCPP
}
// sample_multinomial
unsigned int sample_multinomial(arma::vec theta);
RcppExport SEXP clda_sample_multinomial(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    __result = Rcpp::wrap(sample_multinomial(theta));
    return __result;
END_RCPP
}
// sample_dirichlet
arma::vec sample_dirichlet(unsigned int num_elements, arma::vec alpha);
RcppExport SEXP clda_sample_dirichlet(SEXP num_elementsSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type num_elements(num_elementsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(sample_dirichlet(num_elements, alpha));
    return __result;
END_RCPP
}