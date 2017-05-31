#setwd("./SF_Network/Data/")
occupancy<-read.csv("NewClusterData/Occupancy2.csv", header = FALSE)

# Origin - Destination
origin.idx<-as.numeric(occupancy[3,-1])
origin.idy<-as.numeric(occupancy[4,-1])
dest.idx<-as.numeric(occupancy[5,-1])
dest.idy<-as.numeric(occupancy[6,-1])
times<-as.numeric(as.character(occupancy[-(1:12),1]))

cluster<-rep(NA, length(origin.idx))
market.st <- 3940
cluster[origin.idy > market.st]<-1
cluster[origin.idy < market.st]<-2
cluster[(origin.idx>= ((origin.idy - 2381.194)/0.16)) & origin.idx >= 3300 |((origin.idx >= (origin.idy + 1439.754)/1.5) & origin.idy >=4293.81) |(origin.idx >= 3656.891 & origin.idy <= 3936.467) |(origin.idx >= (origin.idy + 13300.21)/3.323) |(origin.idy <= -503.8062 + 1.21 * origin.idx)& (origin.idx >= 36720.891 & origin.idx <= 3988.736)]<-3
cluster[origin.idy<=-470.8062+1.21*origin.idx & origin.idy>3800]<-3
cluster[origin.idy <= -2745.561 +1.87*origin.idx & origin.idy > 4300]<-3
cluster[origin.idx>3400 & origin.idy<3100]<-3

clustering<-1:158
nodes.1<-as.numeric(occupancy[11,-1])
nodes.2<-as.numeric(occupancy[12,-1])
# Time Periods
adjacency<-matrix(cbind(t(nodes.1), t(nodes.2)), ncol = 2)

adj.mat<-matrix(0, ncol = 158, nrow = 158)
for (i in 1:nrow(adjacency))
  adj.mat[adjacency[i,1], adjacency[i,2]]<-1

#adj.mat<-adj.mat[adjacency[,1],adjacency[,2]]
adj.mat <- adj.mat + t(adj.mat)
adj.mat[adj.mat>1] <- 1

A<-adj.mat

index.clust1<-adjacency[which(cluster==1),1]
index.clust2<-adjacency[which(cluster==2),1]
index.clust3<-adjacency[which(cluster==3),1]
A.1<-A[index.clust1, index.clust1] #obs1
A.2<-A[index.clust2, index.clust2] #obs2
A.3<-A[index.clust3, index.clust3] #obs3

##############################################################################
########################################################################


set.seed(1)
alpha <- 0.01
beta <- 0.5
nobs <- 360
nvars <- dim(A.1)[1]

states <- matrix(sample(0:1, size=nvars, replace=TRUE), ncol=nvars, nrow=nobs)
for (i in 2:nrow(states))
  for (j in 1:ncol(states)) {
    m <- beta*states[i-1,j] + (1-beta)*mean(states[i,which(A.1[j,]==1)])
    states[i,j] <- runif(1)< alpha+(1-2*alpha)*m
  }
image(states)

rho <- 0.8
tau<-0.2
phi<-0.5
spatial.prec1 <- -rho*A.1+diag(rowSums(A.1)*rho+1-rho)
spat.1<-solve(spatial.prec1)* tau^2
times<-seq(1, 21600, length  =  nobs)
temp.1 <- zapsmall(.5^(abs(outer(times,times,"-")/phi)))

y1<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
y2<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
L_S<-t(chol(spat.1))
L_T<-t(chol(temp.1))
#L<-(L_T%x%L_S)
sigma1<-40*L_T%*%y1%*%t(L_S)
sigma2<-40*L_T%*%y2%*%t(L_S)
#image(sigma)
mu.1<-40+1:nobs/nobs*15
mu.2<-90-1:nobs/nobs*10

obs1<-matrix(0, nrow = nrow(states), ncol = ncol(states))
for(i in 1:nrow(states)){
  for(j in 1:ncol(states)){
    if(states[i, j] == 1){obs1[i, j] <- (mu.1+sigma1)[i, j]}else{
      obs1[i, j] <- (mu.2+sigma2)[i, j]
    }
  }
}

for(i in 1:ncol(obs1)){
  obs1[which(obs1[,i] > 100),i] <- 100
  obs1[which(obs1[,i] < 0), i] <- 0
}

#############################################################################

#set.seed(105)
alpha <- 0.01
beta <- 0.5
nobs <- 360
nvars <- dim(A.2)[1]

states <- matrix(sample(0:1, size=nvars, replace=TRUE), ncol=nvars, nrow=nobs)
for (i in 2:nrow(states))
  for (j in 1:ncol(states)) {
    m <- beta*states[i-1,j] + (1-beta)*mean(states[i,which(A.2[j,]==1)])
    states[i,j] <- runif(1)< alpha+(1-2*alpha)*m
  }
image(states)

rho <- 0.8
tau<-0.2
phi<-0.5
spatial.prec1 <- -rho*A.2+diag(rowSums(A.2)*rho+1-rho)
spat.1<-solve(spatial.prec1)* tau^2
times<-seq(1, 21600, length  =  nobs)
temp.1 <- zapsmall(.5^(abs(outer(times,times,"-")/phi)))

y1<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
y2<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
L_S<-t(chol(spat.1))
L_T<-t(chol(temp.1))
#L<-(L_T%x%L_S)
sigma1<-40*L_T%*%y1%*%t(L_S)
sigma2<-35*L_T%*%y2%*%t(L_S)
#image(sigma)
mu.1<-50-1:nobs/nobs*15
mu.2<-25+1:nobs/nobs*15

obs2<-matrix(0, nrow = nrow(states), ncol = ncol(states))
for(i in 1:nrow(states)){
  for(j in 1:ncol(states)){
    if(states[i, j] == 1){obs2[i, j] <- (mu.1+sigma1)[i, j]}else{
      obs2[i, j] <- (mu.2+sigma2)[i, j]
    }
  }
}

for(i in 1:ncol(obs2)){
  obs2[which(obs2[,i] > 100),i] <- 100
  obs2[which(obs2[,i] < 0), i] <- 0
}

##################################################################################

alpha <- 0.01
beta <- 0.5
nobs <- 360
nvars <- dim(A.3)[1]

#A.3[85,85]<-A.3[86, 86]<-1
states <- matrix(sample(0:1, size=nvars, replace=TRUE), ncol=nvars, nrow=nobs)
for (i in 2:nrow(states))
  for (j in 1:ncol(states)) {
    m <- beta*states[i-1,j] + (1-beta)*mean(states[i,which(A.3[j,]==1)], na.rm = TRUE)
    states[i,j] <- runif(1)< alpha+(1-2*alpha)*m
  }
image(states)

rho <- 0.8
tau<-0.2
phi<-0.5
spatial.prec1 <- -rho*A.3+diag(rowSums(A.3)*rho+1-rho)
spat.1<-solve(spatial.prec1)* tau^2
times<-seq(1, 21600, length  =  nobs)
temp.1 <- zapsmall(.5^(abs(outer(times,times,"-")/phi)))

y1<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
y2<-matrix(rnorm(nobs*nvars), nrow=nobs, ncol=nvars)
L_S<-t(chol(spat.1))
L_T<-t(chol(temp.1))
#L<-(L_T%x%L_S)
sigma1<-45*L_T%*%y1%*%t(L_S)
sigma2<-40*L_T%*%y2%*%t(L_S)
#image(sigma)
mu.1<-85-1:nobs/nobs*10
mu.2<-75-1:nobs/nobs*25

obs3<-matrix(0, nrow = nrow(states), ncol = ncol(states))
for(i in 1:nrow(states)){
  for(j in 1:ncol(states)){
    if(states[i, j] == 1){obs3[i, j] <- (mu.1+sigma1)[i, j]}else{
      obs3[i, j] <- (mu.2+sigma2)[i, j]
    }
  }
}

for(i in 1:ncol(obs3)){
  obs3[which(obs3[,i] > 100),i] <- 100
  obs3[which(obs3[,i] < 0), i] <- 0
}


#obs3<-cbind(obs3, obs3)
##########################################################

#obs3<-cbind(obs3, obs3)
##########################################################

dataset<-matrix(0, nrow = 360, ncol = nrow(adjacency))

obs1.l<-which(cluster == 1)
obs2.l<-which(cluster == 2)
obs3.l<-which(cluster == 3)

dataset[,obs1.l]<-obs1
dataset[,obs2.l]<-obs2
dataset[,obs3.l]<-obs3
#dataset[,unlist(obs4.ind)]<-obs4

#######################################################

times <- seq(1, nrow(dataset))
t0 <- seq(1, 21600, len=360)
x0 <- seq(0, 100, len=1000)
sd.x <- 2.5
sd.t <- 10

library(devtools)
install_github("AshwiniKV/fdclust")
library(fdclust)

cdf.junction<-cdf_junction(dataset, sd.x, sd.t, times, t0, x0, 3)
jl<-junction_link(adjacency, clustering, cdf.junction)
junctions<-1:158
