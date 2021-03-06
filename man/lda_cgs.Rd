% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lda_cgs}
\alias{lda_cgs}
\title{LDA: Collapsed Gibbs Sampler with Perplexity Computation}
\usage{
lda_cgs(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in,
  spacing, save_theta, save_beta, save_lp, verbose, test_doc_share = 0,
  test_word_share = 0)
}
\arguments{
\item{num_topics}{Number of topics in the corpus}

\item{vocab_size}{Vocabulary size}

\item{docs_tf}{A list of corpus documents read from the Blei corpus using
\code{\link{read_docs}} (term indices starts with 0)}

\item{alpha_h}{Hyperparameter for \eqn{\theta} sampling}

\item{eta_h}{Smoothing parameter for the \eqn{\beta} matrix}

\item{max_iter}{Maximum number of Gibbs iterations to be performed}

\item{burn_in}{Burn-in-period for the Gibbs sampler}

\item{spacing}{Spacing between the stored samples (to reduce correlation)}

\item{save_theta}{if 0 the function does not save \eqn{\theta} samples}

\item{save_beta}{if 0 the function does not save \eqn{\beta} samples}

\item{save_lp}{if 0 the function does not save computed log posterior for
iterations}

\item{verbose}{from {0, 1, 2}}

\item{test_doc_share}{proportion of the test documents in the corpus.
Must be from [0., 1.)}

\item{test_word_share}{proportion of the test words in each test document.
Must be from [0., 1.)}
}
\value{
The Markov chain output as a list of
  \item{corpus_topic_counts}{corpus-level topic counts from last iteration
  of the Markov chain}
  \item{theta_counts}{document-level topic counts from last iteration
  of the Markov chain}
  \item{beta_counts}{topic word counts from last iteration of the Markov chain}
  \item{theta_samples}{\eqn{\theta} samples after the burn in period, if
  \code{save_theta} is set}
  \item{beta_samples}{\eqn{\beta} samples after the burn in period, if
  \code{save_beta} is set}
  \item{log_posterior}{the log posterior (upto a constant multiplier) of
  the hidden variable \eqn{\psi = (\beta, \theta, z)} in the LDA model,
  if \code{save_lp} is set}
  \item{perplexity}{perplexity of the held-out words' set}
}
\description{
This implements of the collapsed Gibbs sampler for the LDA model---a Markov
chain on \eqn{z}.
}
\details{
To compute perplexity, we first partition words in a corpus into two sets:
(a) a test set (held-out set), which is selected from the set of words in
the test (held-out) documents (identified via \code{test_doc_share} and
\code{test_word_share}) and (b) a training set, i.e., the remaining words in
the corpus. We then run the variational EM algorithm based on the training
set. Finally, we compute per-word perplexity based on the held-out set.
}
\seealso{
Other MCMC: \code{\link{clda_ags_em}},
  \code{\link{clda_ags_sample_alpha}},
  \code{\link{clda_ags}}, \code{\link{clda_mgs}}
}
