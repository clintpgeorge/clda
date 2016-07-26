################################################################################
#
# Distance Functions: This file is part of clda
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



#' Entropy
#'
#' @param v a distribution as a vector
#'
#' @family utils
#'
#' @export
#'
H <- function(v) {
  v <- v[v > 0]
  return(sum(-v * log(v)))
}

#' Jensen-Shannon Divergence (for more than two distributions)
#'
#' For more than two distributions, we need a function to compute the Entropy
#' H(v)
#'
#' @param w A vector of weights which should sum up to 1
#' @param M a matrix with the input distributions as columns
#'
#' @note Adapted from the discussion
#' http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
#'
#' @family utils
#'
#' @export
#'
JSD.matrix <- function(w, M) {
  return(H(M %*% w) - apply(M, 2, H) %*% w)
}

#' Jensen-Shannon Divergence (pairwise)
#'
#' @param p a distribution as a vector
#' @param q another distribution as a vector
#'
#' @note Adapted from the discussion
#' http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
#'
#' @family utils
#'
#' @export
#'
JSD <- function(p, q){
  m <- 0.5 * (p + q)
  JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  return(JS)
}

#' Kullback-Leibler Divergence
#'
#'
#' @param p a distribution as a vector
#' @param q another distribution as a vector
#'
#' @note Adapted from the discussion
#' http://enterotype.embl.de/enterotypes.html
#'
#' @family utils
#'
#' @export
#'
KLD <- function(p, q) sum(p * log(p / q))
