#' Compute the closedness
#'
#' This function calculates the dissimilarity between CDFs for junctions
#' @param cluster cluster
#' @export
#' @examples
#' closedness()

closedness <- function(cluster) {
  neighbours <- c(adjacency[adjacency[,1]%in%cluster,2],adjacency[adjacency[,2]%in%cluster,1])
  result <- mean(neighbours%in%cluster)
  #  attr(result,"neighbours") <- length(neighbours)
  result
}
