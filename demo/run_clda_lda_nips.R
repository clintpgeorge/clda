#' #############################################################################
#' This runs the following cLDA and LDA algorithms on the 16newsgroups dataset.
#'  * cLDA MMALA-Gibbs with J = 1
#'  * cLDA AUX-Gibbs with J = 1
#'  * cLDA MMALA-Gibbs with J = 4
#'  * cLDA AUX-Gibbs with J = 4
#'  * cLDA VEM with J = 4
#'  * LDA CGS
#'
#'
#' Versions:
#'  Februray 18, 2016 - Created
#'  July 26, 2016     - Added to clda
#'
#' Created by:
#'  Clint P. George
#'
#' Example run:
#'  Rscript run_clda_lda_nips.R 1983 .25 .5 1 90
#' #############################################################################

rm(list = ls()); # Removes all objects in the current R state

library(clda)
data("nips") # Loads data
setwd("~")

args               <- commandArgs(TRUE)
SEED               <- as.numeric(args[1])  # seed
eta.h              <- as.numeric(args[2])  # hyperparameter eta
alpha.h            <- as.numeric(args[3])  # hyperparamter alpha
gamma.h            <- as.numeric(args[4])  # hyperparameter gamma
K                  <- as.numeric(args[5])  # number of topics
J                  <- length(unique(cids)) # number of collections
step.size          <- 1e-1   # step size for MGS
max.iter           <- 10   # the maximum number of Gibbs iterations
burn.in            <- 0      # burn in period
spacing            <- 1      # thinning
store.pi           <- 0      # store pi samples ?
store.beta         <- 0      # store beta samples ?
store.theta        <- 0      # store theta samples ?
store.lp           <- 0      # store log posterior for each iteration
verbose            <- 1      # verbose

vi.max.iter        <- 50     # maximum number of variational inference iterations
em.max.iter        <- 100    # maximum number of VEM iterations
vi.conv.thresh     <- 1e-5   # threshold of convergence for variational inference
em.conv.thresh     <- 1e-4   # threshold of convergence for VEM
tau.max.iter       <- 20     # maximum number of tau update iterations
tau.step.size      <- 1      # step size for tau updates
estimate.alpha     <- 0      # estimate alpha ?
estimate.gamma     <- 0      # estimate gamma ?
estimate.eta       <- 0      # estimate eta ?

test.doc.share     <- .2     # percentage of corpus documents used for test set
test.word.share    <- .2     # percentage of words used in each test document

fn.prefix          <-
  paste(
    "clda-lda-", ds.name, "-J", J, "-K", K, "-D", nrow(docs.metadata), "-V", V,
    "-seed", SEED, "-e", eta.h, "-a", alpha.h, "-g", gamma.h, "-",
    format(Sys.time(), "%Y%b%d%H%M%S"), sep = ""
  )
fn.prefix          <- gsub("\\.", "d", fn.prefix)



# Gibbs sampling for the CLDA model ----------------------------------------

set.seed(SEED)
init.pi0           <- matrix(0, nrow = K, ncol = 1)
init.pi0[, 1]      <- sample_dirichlet(K, array(alpha.h, c(K, 1)))

cids0              <- rep(0, length(cids)) # one collection

ptm                <- proc.time();
set.seed(SEED)
clda.aux0          <-
  clda_ags(
    K, V, cids0, docs, alpha.h, gamma.h, eta.h, max.iter, burn.in, spacing,
    store.pi, store.beta, store.theta, store.lp, verbose, init.pi0,
    test.doc.share, test.word.share
  )
clda.aux0.ptm     <- proc.time() - ptm;

ptm                <- proc.time();
set.seed(SEED)
clda.mmala0        <-
  clda_mgs(
    K, V, cids0, docs, alpha.h, gamma.h, eta.h, step.size, max.iter,
    burn.in, spacing, store.pi, store.beta, store.theta, store.lp, verbose,
    init.pi0, test.doc.share, test.word.share
  )
clda.mmala0.ptm    <- proc.time() - ptm;





# Variational inference for the CLDA model ------------------------------------

# Note: This is a dummy initialization

set.seed(SEED)
init.pi            <- matrix(0, nrow = K, ncol = J)
for (j in 1:J) {
  init.pi[, j]     <- sample_dirichlet(K, array(alpha.h, c(K, 1)))
}

ptm                <- proc.time()
set.seed(SEED)
clda.vem           <-
  clda_vem(
    K, V, cids, docs, alpha.h, gamma.h, eta.h, vi.max.iter, em.max.iter,
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
    K, V, cids, docs, alpha.h, gamma.h, eta.h, max.iter, burn.in, spacing,
    store.pi, store.beta, store.theta, store.lp, verbose, init.pi,
    test.doc.share, test.word.share
  )
clda.aux.ptm       <- proc.time() - ptm;



ptm                <- proc.time();
set.seed(SEED)
clda.mmala         <-
  clda_mgs(
    K, V, cids, docs, alpha.h, gamma.h, eta.h, step.size, max.iter,
    burn.in, spacing, store.pi, store.beta, store.theta, store.lp, verbose,
    init.pi, test.doc.share, test.word.share
  )
clda.mmala.ptm     <- proc.time() - ptm;




# Gibbs sampling for the LDA model -----------------------------------------

ptm                <- proc.time();
set.seed(SEED)
lda                <-
  lda_cgs(
    K, V, docs, alpha.h, eta.h, max.iter, burn.in, spacing, store.beta,
    store.theta, store.lp, verbose, test.doc.share, test.word.share
  )
lda.ptm            <- proc.time() - ptm;


cat("Time elapsed (c-lda-mmala, J = 1) = ", clda.mmala0.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");
cat("Time elapsed (c-lda-mmala, J = ", J, ") = ", clda.mmala.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");
cat("Time elapsed (c-lda-aux, J = 1) = ", clda.aux0.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");
cat("Time elapsed (c-lda-aux, J = ", J, ") = ", clda.aux.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");
cat("Time elapsed (c-lda-vem, J = ", J, ") = ", clda.vem.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");
cat("Time elapsed (lda) = ", lda.ptm[3]/(60 * (max.iter - burn.in) / spacing), " min\n", sep = "");


# Saves every object into a file -------------------------------------------


save.image(paste(fn.prefix, ".RData", sep = ""))
