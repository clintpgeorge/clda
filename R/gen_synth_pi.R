#' @export
#'
gen_synth_pi <- function(K, J, alpha.h, num.samples = 100){
  pi.prior <- array(0, dim = c(num.samples, K))  # num.samples x K

  for (s in 1:num.samples) {
    pi.prior[s, ] <- sample_dirichlet(K, rep(alpha.h, K))
  }

  dm <- as.matrix(dist(pi.prior))
  x <- rep(0, num.samples*num.samples)
  y <- rep(0, num.samples*num.samples)
  d <- rep(0, num.samples*num.samples)
  cnt <- 1
  for (i in 1:num.samples) {
    for (j in 1:num.samples) {
      x[cnt] <- i
      y[cnt] <- j
      d[cnt] <- dm[i,j]
      cnt <- cnt + 1
    }
  }

  fth.q <- quantile(d, probs = .5)
  sth.q <- quantile(d, probs = .75)
  selected <- c()
  for (cnt in 1:(num.samples*num.samples)) {

    if ((d[cnt] >= fth.q) && (d[cnt] < sth.q)) {
      xc <- x[cnt]
      yc <- y[cnt]

      if (!is.element(xc, selected)) {
        selected <- rbind(selected, xc)
        if (length(selected) == J) {
          break
        }
      }

      if (!is.element(yc, selected)) {
        selected <- rbind(selected, yc)
        if (length(selected) == J) {
          break
        }
      }
    }


  }

  sel.pi.prior <- pi.prior[selected,] # J x K
  t(sel.pi.prior)
}

