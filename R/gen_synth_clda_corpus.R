################################################################################
#
# Generates a Synthetic c-LDA Corpus: This file is part of clda
#
# Copyright (c) 2016  Clint P. George
#
# clda is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# clda is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


#' Generates a Synthetic c-LDA Corpus
#'
#' Generates documents using the c-LDA generative process based on a set of
#' predefined values.
#'
#' @param K number of topics
#' @param V vocabulary size
#' @param J number of collections
#' @param collection.size number of documents in each collection (a list of J
#' elements)
#' @param doc.size number of words in each document
#' @param alpha.h hyperparameter for collection-level Dirichlet sampling
#' @param gamma.h hyperparameter for document-level Dirichlet sampling
#' @param eta.h hyperparameter for topic Diriclet sampling
#'
#' @return a list of generated corpus and their statistics
#'
#' @export
#'
#' @family corpus
#'
#' @details Last modified on: March 04, 2016
#'
#' @examples
#
#' ## Generates documents with given parameters
#'
#' J                  <- 2
#' K                  <- 4
#' V                  <- 20
#' alpha.h            <- 2
#' gamma.h            <- .2
#' eta.h              <- .25
#' doc.size           <- 80
#' collection.size    <- c(40, 40) # number of documents in each collection
#' ds.name            <- paste("synth-J", J, "-K", K, "-D", D, "-V", V, sep = "")
#'
#' ds                 <- gen_synth_clda_corpus(K, V, J, collection.size, doc.size, alpha.h, gamma.h, eta.h)
#'
gen_synth_clda_corpus <- function(K, V, J, collection.size, doc.size, alpha.h, gamma.h, eta.h)
{

  calc_topic_counts <- function(Z, K)
  {
    Nt <- array(0, c(1, K));
    for (k in 1:K) Nt[k] <- sum(Z == k);

    return(Nt);
  }


  num_docs <- sum(collection.size) # number of documents
  theta.counts <- matrix(0, nrow = K, ncol = num_docs) # document topic word counts
  beta.counts <- matrix(0, nrow = K, ncol = V) # topic word counts
  pi.counts <- matrix(0, nrow = K, ncol = J)

  pi.samples <- matrix(0, nrow = K, ncol = J)
  theta.samples <- matrix(0, nrow = K, ncol = num_docs)
  beta.samples <- matrix(1e-2, nrow = K, ncol = V)

  alpha.v <- array(alpha.h/K, c(K, 1))
  gamma.v <- array(gamma.h, c(K, 1))
  eta.v <- array(eta.h, c(1, V))

  cids <- rep(0, num_docs)
  did <- c()
  wid <- c()
  zid <- c()
  doc.N <- array(doc.size, dim = c(num_docs, 1))
  doc.idx <- vector("list", num_docs)
  docs <- vector("list", num_docs)
  word.idx <- 1; # initialize the corpus word index

  ptm <- proc.time()

  # Topic Dirichlet sampling
  for (k in 1:K) {
    beta.samples[k,]  <- sample_dirichlet(V, eta.v)
  }


  # collection sampling
  #
  d.index <- 1
  for (j in 1:J) {

    pi.samples[, j] <- sample_dirichlet(K, alpha.v)

    # Document sampling
    for (d in 1:collection.size[j]) {

      theta.samples[, d.index] <- sample_dirichlet(K, gamma.v * pi.samples[, j]);

      did <- cbind(did, array(1, c(1, doc.N[d.index])) * d.index); # document instances
      z_d <- c();
      indices <- c();

      # Word sampling
      word_ids <- rep(0, doc.N[d.index])
      for (i in 1:doc.N[d.index]) {

        # samples topic
        z_dn <- which(rmultinom(1, size = 1, prob = theta.samples[, d.index]) == 1)
        z_d <- cbind(z_d, z_dn)

        # samples word
        w_dn <- which(rmultinom(1, size = 1, beta.samples[z_dn,]) == 1)
        word_ids[i] <- w_dn - 1
        wid <- cbind(wid, w_dn)

        indices <- cbind(indices, word.idx)
        word.idx <- word.idx + 1
        pi.counts[z_dn,j] <- pi.counts[z_dn,j] + 1

      }

      cids[d.index] <- (j - 1) # collection id
      doc <- as.data.frame(table(word_ids))
      doc <- rbind(as.integer(levels(doc$word_ids)), doc$Freq)
      docs[[d.index]] <- doc
      doc.idx[[d.index]] <- as.integer(indices) # document word indices

      theta.counts[, d.index] <- calc_topic_counts(z_d, K) # document topic counts
      zid <- cbind(zid, z_d)

      d.index <- d.index + 1
    }

  }

  total.N <- sum(doc.N);
  for (i in 1:total.N) {
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1;
  }

  ptm <- proc.time() - ptm
  cat("Corpus generation time: ", ptm[3], ", number of total words: ",
      total.N, "\n", sep = "")

  # returns a list
  list(
    docs = docs,
    cids = cids,
    did = as.vector(did),
    wid = as.vector(wid),
    zid = as.vector(zid),
    pi.counts = pi.counts,
    theta.counts = theta.counts,
    beta.counts = beta.counts,
    pi.samples = pi.samples,
    theta.samples = theta.samples,
    beta.samples = beta.samples,
    total.N = total.N,
    doc.N = doc.N,
    doc.idx = doc.idx
  )

}


#' Generates a Synthetic c-LDA Corpus (given the \eqn{\pi)} matrix)
#'
#' Generates documents using the c-LDA generative process based on a set of
#' predefined values.
#'
#' @param K number of topics
#' @param V vocabulary size
#' @param J number of collections
#' @param collection.size number of documents in each collection (a list of J
#' elements)
#' @param doc.size number of words in each document
#' @param pi.prior the \eqn{\pi)} matrix
#' @param gamma.h hyperparameter for document-level Dirichlet sampling
#' @param eta.h hyperparameter for topic Diriclet sampling
#'
#' @return a list of generated corpus and their statistics
#'
#' @export
#'
#' @family corpus
#'
#' @details Last modified on: March 04, 2016
#'
#' @examples
#
#' ## Generates documents with given parameters
#'
#' J                  <- 2
#' K                  <- 4
#' V                  <- 20
#' alpha.h            <- 2
#' gamma.h            <- .2
#' eta.h              <- .25
#' doc.size           <- 80
#' collection.size    <- c(40, 40) # number of documents in each collection
#' ds.name            <- paste("synth-J", J, "-K", K, "-D", D, "-V", V, sep = "")
#'
#' ds                 <- gen_synth_clda_corpus(K, V, J, collection.size, doc.size, alpha.h, gamma.h, eta.h)
#'
gen_synth_clda_corpus_pi <- function(K, V, J, collection.size, doc.size, pi.prior, gamma.h, eta.h)
{

  calc_topic_counts <- function(Z, K)
  {
    Nt <- array(0, c(1, K));
    for (k in 1:K) Nt[k] <- sum(Z == k);

    return(Nt);
  }


  num_docs <- sum(collection.size) # number of documents
  theta.counts <- matrix(0, nrow = K, ncol = num_docs) # document topic word counts
  beta.counts <- matrix(0, nrow = K, ncol = V) # topic word counts
  pi.counts <- matrix(0, nrow = K, ncol = J)

  theta.samples <- matrix(0, nrow = K, ncol = num_docs)
  beta.samples <- matrix(1e-2, nrow = K, ncol = V)

  gamma.v <- array(gamma.h, c(K, 1))
  eta.v <- array(eta.h, c(1, V))

  cids <- rep(0, num_docs)
  did <- c()
  wid <- c()
  zid <- c()
  doc.N <- array(doc.size, dim = c(num_docs, 1))
  doc.idx <- vector("list", num_docs)
  docs <- vector("list", num_docs)
  word.idx <- 1; # initialize the corpus word index

  ptm <- proc.time()

  # Topic Dirichlet sampling
  for (k in 1:K) {
    beta.samples[k,]  <- sample_dirichlet(V, eta.v)
  }


  # collection sampling
  #
  d.index <- 1
  for (j in 1:J) {

    # Document sampling
    for (d in 1:collection.size[j]) {

      theta.samples[, d.index] <- sample_dirichlet(K, gamma.v * pi.prior[, j]);

      did <- cbind(did, array(1, c(1, doc.N[d.index])) * d.index); # document instances
      z_d <- c();
      indices <- c();

      # Word sampling
      word_ids <- rep(0, doc.N[d.index])
      for (i in 1:doc.N[d.index]) {

        # samples topic
        z_dn <- which(rmultinom(1, size = 1, prob = theta.samples[, d.index]) == 1)
        z_d <- cbind(z_d, z_dn)

        # samples word
        w_dn <- which(rmultinom(1, size = 1, beta.samples[z_dn,]) == 1)
        word_ids[i] <- w_dn - 1
        wid <- cbind(wid, w_dn)

        indices <- cbind(indices, word.idx)
        word.idx <- word.idx + 1
        pi.counts[z_dn,j] <- pi.counts[z_dn,j] + 1

      }

      cids[d.index] <- (j - 1) # collection id
      doc <- as.data.frame(table(word_ids))
      doc <- rbind(as.integer(levels(doc$word_ids)), doc$Freq)
      docs[[d.index]] <- doc
      doc.idx[[d.index]] <- as.integer(indices) # document word indices

      theta.counts[, d.index] <- calc_topic_counts(z_d, K) # document topic counts
      zid <- cbind(zid, z_d)

      d.index <- d.index + 1
    }

  }

  total.N <- sum(doc.N);
  for (i in 1:total.N) {
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1;
  }

  ptm <- proc.time() - ptm
  cat("Corpus generation time: ", ptm[3], ", number of total words: ",
      total.N, "\n", sep = "")

  # returns a list
  list(
    docs = docs,
    cids = cids,
    did = as.vector(did),
    wid = as.vector(wid),
    zid = as.vector(zid),
    pi.counts = pi.counts,
    theta.counts = theta.counts,
    beta.counts = beta.counts,
    theta.samples = theta.samples,
    beta.samples = beta.samples,
    total.N = total.N,
    doc.N = doc.N,
    doc.idx = doc.idx
  )

}
