################################################################################
#
# This file is part of clda
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

#' Align Topics
#'
#' Align the rows of the \code{dest} matrix in terms of the \code{src} matrix
#'
#' @param src a K x V matrix
#' @param dest a K x V matrix
#' @param dist.fn the distance function that is used to find distance between rows
#'
#' @return  \code{align.indices} vector
#'
#' @note Created on April 14, 2016
#'
#' @export
#'
align_beta <- function(src, dest, dist.fn=KLD){
  align.indices <- rep(0, K)
  is.selected <- rep(0, K)
  align.dist <- rep(0, K)
  for (i in 1:K) {

    min.dist <- Inf
    for (j in 1:K) {
      dd <- dist.fn(src[i,], dest[j,])
      # cat(i, j, dd, "-> ")
      if (dd < min.dist && is.selected[j] == 0) {
        # cat(align.indices[i], is.selected[align.indices[i]])
        # first reverts the previous selection, if any
        if (align.indices[i] > 0) {
          is.selected[align.indices[i]] <- 0
        }

        # selects
        align.indices[i] <- j
        is.selected[j] <- 1
        align.dist[i] <- dd
        min.dist <- dd
      }
      # cat("\n")
    }
  }
  align.indices
}
