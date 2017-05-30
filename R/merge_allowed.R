
#' Check if the two junctions within a cluster are adjacent 
#'
#' This function allows only adjacenct junctions to merge to ensure spatially contiguous clusters 
#' @param i Junction i 
#' @param j Junction j
#' @param clustering Set of clusters 
#' @param adjacency Matrix to represent links from Junction j to Junction i
#' @export
#' @examples
#' merge_allowed()

merge_allowed <- function(i, j, clustering, adjacency) {
  junc.i <- which(clustering==i)  
  junc.j <- which(clustering==j)
  any(adjacency[,1]%in%junc.i & adjacency[,2]%in%junc.j) | any(adjacency[,1]%in%junc.j & adjacency[,2]%in%junc.i)
}
