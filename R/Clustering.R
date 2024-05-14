#' K-means Clustering
#'
#' @param RawData A data frame containing the raw data as read in by [pooledpeaks::LoadData]
#' @param K An integer specifying the number of clusters.
#'
#' @importFrom stats kmeans
#' @importFrom utils capture.output
#'
#' @return A list containing the results of the K-means cluster analysis, including cluster assignments and original data.
#' @export
#'
#' @examples

cluster <- function(RawData=data.frame,K=2) {
  d <- RawData
  N <- TypedLoci(d)
  Ns <- 0
  for (i in 1:ncol(N)) Ns[i] <- N[i,i]
  w <- which(Ns==max(N))
  w <- w+2
  d <- d[,c(1:2,w)]
  dat <- d
  names <- colnames(d)
  names <- names[-(1:2)]
  PopLabel <- substr(names,1,1)
  PopFactor <- as.factor(PopLabel)
  d <- subset(d,(d[,2] != 'n')&(d[,3]!='NA'))
  d <- d[,-c(1:2)]
  fit <- stats::kmeans(t(d),centers=K,iter.max=20,nstart=50)
  utils::capture.output(fit, file = "clusterOutput.txt")

  fit <- fit$cluster
  fit <- list(PopFactor=PopFactor,clust=fit,dat=dat)

  return(fit)
}


#' Sample Of Loci
#'
#' An internal function that supports ClusterFromSamples. Sample loci from a dataset based on the number of loci specified.
#'
#' @param aaax A data frame containing loci information.
#' @param NLoci An integer specifying the number of loci to sample.
#'
#' @return A data frame containing the sampled loci.
#'
SampleOfLoci <- function(aaax=data.frame,NLoci=max(aaax[,1])) {
  x <- sample(1:max(aaax[,1]),NLoci,replace = TRUE)
  z <- aaax[1,]
  for (i in x) {
    b <- subset(aaax,aaax[,1]==i)
    z <- rbind(z,b)
  }
  aaax <- z[-1,]
  return(aaax)
}



#' Cluster From Samples
#'
#' Perform clustering on samples of loci from a data frame and calculate statistics.
#'
#' @param datafile A data frame containing the data.
#' @param numloci An integer specifying the number of loci to sample.
#' @param reps An integer specifying the number of repetitions.
#'
#' @return A matrix containing statistics calculated from the clustering results.
#' @export
#'
#' @examples

ClusterFromSamples <- function(datafile=data.frame,numloci=5,reps=100) {
  for (j in 1:reps) {
    A <- SampleOfLoci(datafile,numloci)
    names <- colnames(datafile)
    names <- names[-(1:2)]
    PopLabel <- substr(names,1,1)
    PopFactor <- as.factor(PopLabel)
    H <- cluster(A,K=length(levels(PopFactor)))
    T <- table(H$clust,H$PopFactor)
    S <- colSums(T)
    m <-0
    for (i in 1:ncol(T)) {
      m[i] <- max(T[,i]/S[i])
    }

    if (j==1)
    {M <- m} else
    {M <- rbind(M,m)}
  }

  M <- round(M,4)
  return(M)
}



#' Multi Dimensional Scaling (MDS) Plot
#'
#' Generate a multidimensional scaling (MDS) plot from genetic distance data.
#'
#' @param distance A matrix containing the genetic distance data.
#' @param pcs A numeric vector specifying the principal coordinates to plot.
#' @param PF A factor vector specifying population labels.
#' @param y A character vector specifying colors for population labels.
#'
#' @importFrom stats cmdscale
#' @importFrom graphics par
#' @importFrom graphics title
#'
#' @return The ouput is the MDS plot for the samples for the specified principal coordinates.
#'
#' @export
#'
#' @examples

MDSplot<- function(distance=matrix,pcs=c(1,2),PF=NULL, y= c('dodgerblue','red','turquoise3','purple','olivedrab3') ) {
  graphics::par(mar=c(5, 4, 4, 8), xpd=TRUE)
  K <- nrow(distance)
  if (K < 11){
    K <- K - 1
  } else {
    K <- 10
  }
  names <- colnames(distance)
  PopLabel <- substr(names,1,1)

  if (is.null(PF)) PF <- as.factor(PopLabel)

  if (is.null(y)){ y <- c('dodgerblue','red','turquoise3','purple','olivedrab3')
  } else {y<- y}


  z <- rep(1,length(PF))

  for (i in 1:length(PF)) z[i] <- y[PF[i]]

  E <- stats::cmdscale(distance,eig=TRUE, k=K)   # multidimensional scaling from R

  T <- 0
  for (i in 1:length(E$eig)) if (E$eig[i]>0) T <- T + E$eig[i]

  percent <- round(E$eig/T,3)*100
  LtextX <-  paste('Principal Coordinate (',pcs[1],') :  ',percent[pcs[1]],'% Dispersion')
  LtextY <-paste('Principal Coordinate (',pcs[2],') :  ',percent[pcs[2]],'% Dispersion')

  plot(E$points[,pcs[1]],E$points[,pcs[2]],pch=21,col='gray25',bg=z,cex=1,
       xlab=NA,ylab=NA,las=1)
  legend('topright', legend = c(substr(levels(PF), 1,6)), xjust=1, yjust=1,
         inset=c(-0.2, 0), fill = y, bty = "n", cex = 1)
  graphics::title(xlab=LtextX)
  graphics::title(ylab=LtextY)
  graphics::title(main='Principal Coordinates of Genetic Distances',font.main=1,cex.main=1.25)
}
