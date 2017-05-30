#' Cumulative distribution function for each junction 
#' @param dataset Data for incoming links over the entire junction
#' @param sd.x Vector of standard deviations used for calculating the distribution function 
#' @param sd.t Vector of standard deviations used for calculating the density function over differences in time   
#' @param times Number of observations
#' @param t0 Time points recorded corresponding to relevant sampling rate
#' @param x0 Vector of quantiles
#' @param type CDF calculated with functional only (type = 2), distributional only (type = 1) and both functional and distributional (type = 3)
#' @export
#' @examples
#' cdf_junction()
#' 
cdf_junction<-function(dataset, sd.x, sd.t, times, t0, x0, type){
if(type == 1 | type == 2| type == 3){
  if(type == 1){
  dataset<-matrix(rep(colMeans(dataset), each = nrow(dataset)), nrow = nrow(dataset))

  # Create a list with # entries = # junctions
  cdf <- matrix(nrow=length(t0), ncol=length(x0))

  cdfj_mean_vec<-vector("list", 400)
  for(j in unique(adjacency[,1]))
  {
  sel <- which(adjacency[,1] == j)
  time.weights <- matrix(dnorm(outer(times, t0, "-"),sd=sd.t), nrow=length(times))
  time.weights.sum <- colSums(time.weights)
  for (i in 1:length(x0)) {
    cdfs <- rowMeans(matrix(pnorm(x0[i], mean=as.numeric(as.matrix(dataset[, sel])), sd=sd.x), ncol = length(sel)))
    cdf[,i] <-  colSums(cdfs*time.weights)/time.weights.sum
  }
  cat(".")
  cdfj_mean_vec[[j]]<-cdf
  }
return(cdfj_mean_vec)
}
  
if(type == 2){
  cdf <- matrix(nrow=length(t0), ncol=length(x0))
  
  sd.x <- 2.5
  sd.t <- 10
  cdfj_fun<-vector("list", 400)
  for(j in unique(adjacency[,1]))
  {
    sel <- which(adjacency[,1] == j)
    for (i in 1:length(x0)){
      if(length(sel)>1){cdfs<-rowMeans(dataset[,sel])}else{
        cdfs<-dataset[,sel]
      }
      cdf[,i] <-  cdfs
    }
    cat(".")
    cdfj_fun[[j]]<-cdf
  }
  return(cdfj_fun)
}
  
if(type == 3){
  cdf <- matrix(nrow=length(t0), ncol=length(x0))
  
  cdfj_funmean<-vector("list", 400)
  for(j in unique(adjacency[,1]))
  {
    sel <- which(adjacency[,1] == j)
    time.weights <- matrix(dnorm(outer(times, t0, "-"),sd=sd.t), nrow=length(times))
    time.weights.sum <- colSums(time.weights)
    for (i in 1:length(x0)) {
      cdfs <- rowMeans(matrix(pnorm(x0[i], mean=as.numeric(as.matrix(dataset[,sel])), sd=sd.x), ncol=length(sel)))
      std.cdfs <-  round(colSums(cdfs*time.weights)/time.weights.sum,3)
      #std.cdfs[which(is.nan(std.cdfs))] <- 0
      cdf[,i] <-std.cdfs
    }
    cat(".")
    cdfj_funmean[[j]]<-cdf
  }
  return(cdfj_funmean)
}
  }
    }