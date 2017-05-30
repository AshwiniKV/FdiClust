#' Compute the dissimilarity function
#'
#' This function calculates the dissimilarity between CDFs for junctions 
#' @param i Junction i 
#' @param j Junction j
#' @param cdf.junction List of cumulative distribution functions for junctions
#' @export
#' @examples
#' compute_badness()

compute_badness<-function(i, j, cdf.junction) {
  cdf1<-cdf.junction[[i]]
  cdf2<-cdf.junction[[j]]
  badness <- mean(colMeans(abs(cdf1 - cdf2), na.rm = TRUE),na.rm= TRUE)
}