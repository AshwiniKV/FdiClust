#' Compute the shape penalty function 
#'
#' This function calculates the dissimilarity between CDFS for junctions 
#' @param i Junction i 
#' @param j Junction j
#' @param clustering List of cumulative distribution functions
#' @export
#' @examples
#' shape_penalty()

shape_penalty <- function(i, j, clustering) {
  clust1 <- which(clustering==i)
  clust2 <- which(clustering==j)
  -closedness(c(clust1,clust2))+.5*closedness(clust1)+.5*closedness(clust2)
}