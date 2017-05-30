
#' Clustering algorithm
#'
#' This function clusters the network using CDFs calculated for individual junctions aggregated over relevant incoming links.
#' @param clustering Set of clusters 
#' @param adjacency Matrix to represent links from junction j to junction i
#' @param cdf.junction Matrix to represent CDF for each junction
#' @export
#' @examples
#' junction_link()

junction_link<-function(adjacency, clustering, cdf.junction)
{
  history <- list(clustering)
  counter <- max(clustering)
  position <- character(length(clustering))
  cumulative.badness <- 0
  while (length(unique(clustering))>1) {
    clusters <- unique(clustering)
    #newcdf<-cdf.junction
    best.badness <- Inf
    best.idx <- NULL
    for (i in 1:(length(clusters)-1))
      for (j in (i+1):length(clusters)) {
        clust.i <- clusters[i]
        clust.j <- clusters[j]
        if (merge_allowed(clust.i, clust.j, clustering, adjacency)) {
          clusters.new <- clustering          
          clusters.new[clusters.new%in%c(clust.i, clust.j)] <- max(clusters.new)+1
          attributes(clusters.new) <- NULL
          current.badness <- compute_badness(clust.i, clust.j, cdf.junction) + 1.5*
            shape_penalty(clust.i, clust.j, clustering) + 2
          cat(paste(clust.i,"+",clust.j,"=",current.badness,"\n"))
          attr(current.badness, "clusters.new") <- clusters.new
          if (!is.nan(current.badness) && ((current.badness < best.badness-1e-8) ||
                                           ((current.badness < best.badness+1e-8) && (runif(1)<0.5)))) {
            #Keep least worst merge
            best.idx <- c(clust.i, clust.j)
            best.badness <- current.badness
          }
        }
      }
    counter <- counter + 1
    cat(".")
    select <- clustering==best.idx[1]
    if (is.null(best.idx) || mean(which(clustering==best.idx[1]))<=mean(which(clustering==best.idx[2]))) {
      c1 <- "l"
      c2 <- "r"
    } else {
      c1 <- "r"
      c2 <- "l"
    }
    position[select] <- paste0(c1,position[select])
    select <- clustering==best.idx[2]
    position[select] <- paste0(c2,position[select])   
    clustering[clustering%in%best.idx] <- counter
    cdf.junction[[counter]]<-(cdf.junction[[best.idx[1]]]+ cdf.junction[[best.idx[2]]])/2
    cumulative.badness <- cumulative.badness + best.badness
    attr(clustering,"badness") <- cumulative.badness
    history[[length(history)+1]] <- clustering
    if (is.null(best.idx)) break
    
  }
  attr(history, "position") <- position
  history
}