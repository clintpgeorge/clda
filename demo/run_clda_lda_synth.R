#' #############################################################################
#' This script runs the cLDA algorithms Variational Expectation Maximization
#' (VEM), Metropolis Adjusted Langevin Algorithm within Gibbs Sampler (MGS), and
#' Auxiliary Variable Update within Gibbs Sampler (AGS), and the LDA Collapsed
#' Gibbs Samper (CGS), using a synthetic dataset.
#'
#' Created on:
#' April 01, 2016
#'
#' Last modified on:
#' June 02, 2016
#'
#' Created by:
#' Clint P. George
#'
#' Example run:
#' Rscript run_clda_lda_synth.R 1983 .1 1 .25 3 2 40 100 200
#' #############################################################################

rm(list = ls()); # Removes all objects in the current R state

library(clda)

setwd('~')


# Handles commandline arguments

args               <- commandArgs(TRUE)
SEED               <- as.numeric(args[1]) # seed
alpha.h            <- as.numeric(args[2]) # alpha
gamma.h            <- as.numeric(args[3]) # gamma
eta.h              <- as.numeric(args[4]) # eta
K                  <- as.numeric(args[5]) # number of topics
J                  <- as.numeric(args[6]) # number of collections
V                  <- as.numeric(args[7]) # vocab size
coll.size          <- as.numeric(args[8]) # number of documents in a collection
doc.size           <- as.numeric(args[9]) # number of words in a document

test.doc.share     <- .2   # percentage of corpus documents used for test set
test.word.share    <- .2   # percentage of words used in each test document
step.size          <- 1e-1 # step size for MGS
max.iter           <- 1000 # the maximum number of Gibbs iterations
burn.in            <- 0    # burn in period
spacing            <- 1    # thinning
store.pi           <- 1    # store pi samples ?
store.beta         <- 0    # store beta samples ?
store.theta        <- 0    # store theta samples ?
store.lp           <- 0    # store log posterior for each iteration
verbose            <- 2
collection.size    <- rep(coll.size, J)  # number of documents in each collection

vi.max.iter        <- 50   # maximum number of variational inference iterations
em.max.iter        <- 100  # maximum number of VEM iterations
vi.conv.thresh     <- 1e-6 # threshold of convergence for variational inference
em.conv.thresh     <- 1e-4 # threshold of convergence for VEM
tau.max.iter       <- 20   # maximum number of tau update iterations
tau.step.size      <- 1    # step size for tau updates
estimate.alpha     <- 0    # estimate alpha ?
estimate.gamma     <- 0    # estimate gamma ?
estimate.eta       <- 0    # estimate eta ?

burn.in.pi         <- 50   # burn in period for \pi

fn.prefix          <-
  paste(
    "clda-lda-J", J, "-K", K, "-D", sum(collection.size), "-V", V, "-cs", coll.size,
    "-ds", doc.size, "-seed", SEED, "-a", alpha.h, "-g", gamma.h, "-e", eta.h,
    "-", format(Sys.time(), "%Y%b%d%H%M%S"), sep = ""
  )
fn.prefix          <- gsub("\\.", "d", fn.prefix)

# Generates data -----------------------------------------------------------

set.seed(SEED)
pi.prior           <- gen_synth_pi(K, J, alpha.h)
if (dim(pi.prior)[1] == J) {
  pi.prior <- t(pi.prior)
}

set.seed(SEED) # data generation seed
ds                 <-
  gen_synth_clda_corpus_pi(
    K, V, J, collection.size, doc.size, pi.prior, gamma.h, eta.h
  )

set.seed(SEED)
init.pi            <- matrix(0, nrow = K, ncol = J)
for (j in 1:J) {
  init.pi[, j]     <- sample_dirichlet(K, array(alpha.h, c(K, 1)))
}

# Variational inference for the CLDA model ------------------------------------

ptm                <- proc.time()
set.seed(SEED)
clda.vem           <-
  clda_vem(
    K, V, ds$cids, ds$docs, alpha.h, gamma.h, eta.h, vi.max.iter, em.max.iter,
    vi.conv.thresh, em.conv.thresh, tau.max.iter, tau.step.size, estimate.alpha,
    estimate.gamma, estimate.eta, verbose, init.pi, test.doc.share,
    test.word.share
  )
clda.vem.ptm       <- proc.time() - ptm


# Gibbs sampling for the CLDA model -----------------------------------------

init.pi            <- clda.vem$vi_tau_t[,,1] # to keep the same initial point for all algorithms


ptm                <- proc.time();
set.seed(SEED)
clda.aux           <-
  clda_ags(
    K, V, ds$cids, ds$docs, alpha.h, gamma.h, eta.h, max.iter, burn.in, spacing,
    store.pi, store.beta, store.theta, store.lp, verbose, init.pi,
    test.doc.share, test.word.share, 0
  )
clda.aux.ptm       <- proc.time() - ptm;


ptm                <- proc.time();
set.seed(SEED)
clda.mmala         <-
  clda_mgs(
    K, V, ds$cids, ds$docs, alpha.h, gamma.h, eta.h, step.size, max.iter,
    burn.in, spacing, store.pi, store.beta, store.theta, store.lp, verbose,
    init.pi, test.doc.share, test.word.share, 0
  )
clda.mmala.ptm     <- proc.time() - ptm;



# Gibbs sampling for the CLDA model (Initializes with VEM estimate) -----------

vem.init.pi <- t(apply(clda.vem$vi_tau, 1, function(x)(x / colSums(clda.vem$vi_tau))))

ptm                <- proc.time();
set.seed(SEED)
clda.aux1          <-
  clda_ags(
    K, V, ds$cids, ds$docs, alpha.h, gamma.h, eta.h, max.iter, burn.in, spacing,
    store.pi, store.beta, store.theta, store.lp, verbose, vem.init.pi,
    test.doc.share, test.word.share, burn.in.pi
  )
clda.aux1.ptm      <- proc.time() - ptm;


ptm                <- proc.time();
set.seed(SEED)
clda.mmala1        <-
  clda_mgs(
    K, V, ds$cids, ds$docs, alpha.h, gamma.h, eta.h, step.size, max.iter,
    burn.in, spacing, store.pi, store.beta, store.theta, store.lp, verbose,
    vem.init.pi, test.doc.share, test.word.share, burn.in.pi
  )
clda.mmala1.ptm    <- proc.time() - ptm;


# Gibbs sampling for the LDA model -----------------------------------------

ptm                <- proc.time();
set.seed(SEED)
lda                <-
  lda_cgs(
    K, V, ds$docs, alpha.h, eta.h, max.iter, burn.in, spacing, store.beta,
    store.theta, store.lp, verbose, test.doc.share, test.word.share
  )
lda.ptm            <- proc.time() - ptm;

cat("Execution time (clda-vem)       = ", clda.vem.ptm[3], "\n")
cat("Execution time (clda-cgs-mmala) = ", clda.mmala.ptm[3], "\n");
cat("Execution time (clda-cgs-aux)   = ", clda.aux.ptm[3], "\n");
cat("Execution time (clda-cgs-mmala1)= ", clda.mmala1.ptm[3], "\n");
cat("Execution time (clda-cgs-aux1)  = ", clda.aux1.ptm[3], "\n");
cat("Execution time (lda-cgs)        = ", lda.ptm[3], "\n");


# Saves every object into a file -------------------------------------------


save.image(paste(fn.prefix, ".RData", sep = ""))

num.samples <- dim(clda.mmala$pi_samples)[3]

cat("True pi\n")
pi.prior

cat("\nVariational tau\n")
clda.vem$vi_tau

cat("\nVEM Estimate of pi\n")
t(apply(clda.vem$vi_tau, 1, function(x) (x / colSums(clda.vem$vi_tau))))

cat("\nMMALA CGS Estimate of pi\n")
clda.mmala$pi_samples[,,num.samples]

cat("\nAUX CGS Estimate of pi\n")
clda.aux$pi_samples[,,num.samples]

cat("\nMMALA CGS (Initializes with VEM) Estimate of pi\n")
clda.mmala1$pi_samples[,,num.samples]

cat("\nAUX CGS (Initializes with VEM) Estimate of pi\n")
clda.aux1$pi_samples[,,num.samples]
