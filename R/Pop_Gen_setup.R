#' Load Genetic Data
#'
#' This function imports data for genetic analysis.
#'
#' @param datafile The path to your datafile. The format of your data should be .txt or .csv.
#' @importFrom utils read.table
#'
#' @return A data frame containing the imported data formatted in the way necessary for downstream population genetic functions.
#' @export
#'
#' @examples
LoadData <- function(datafile=NULL)
{
  if (is.null(datafile))
    datafile <- readline(" Please enter the path to your data file :    ")
  N1453 <- utils::read.table(datafile,header=T)

  Loci <- 0
  Locus <- 1

  for (i in 1:(nrow(N1453)-1)) {
    if (N1453[i+1,1]=='n') {
      Loci <- Loci+1
    }
    Locus[i] <- Loci
  }
  Locus[nrow(N1453)]<- Loci

  N1453 <- cbind(Locus,N1453)

  return(N1453)
}


#' Typed Loci
#'
#' This function calculates the number of loci successfully genotyped by each individual included in our data set
#'
#' @param datafile A data frame containing the input data must be in LoadData style [pooledpeaks::LoadData].
#'
#' @return A matrix representing processed data.
#' @export
#'
#' @examples

TypedLoci <- function(datafile=data.frame) {

  b <- subset(datafile,datafile[,3]!='NA')

  for (i in 1:nrow(b))
    for (j in 3:ncol(b))
      if (b[i,j]>1) b[i,j]<-1

  b <- subset(b,b[,2]=='n')
  b <- b[,-(1:2)]
  X <- as.matrix(b)
  N <- t(X)%*%X

}

#' Gene Identity Matrix
#'
#' Using the number of typed loci, this function calculates the gene identity between all possible pairwise combinations between individuals for all markers creating a matrix.
#'
#' @param RawData A data frame containing the input data must be in LoadData style [pooledpeaks::LoadData].
#' @param LociGenotyped The Output from the TypedLoci function
#'
#' @return The Gene Identity Matrix
#' @export
#'
#' @examples

GeneIdentityMatrix <- function(RawData=data.frame,LociGenotyped=matrix) {
  d <- subset(RawData,(RawData[,2] != 'n')&(RawData[,3]!='NA'))
  d <- d[,-(1:2)]
  X <- as.matrix(d)
  J <- t(X)%*%X
  J <- J/LociGenotyped
  return(J)
}


#' Genetic Distance Matrix
#'
#' This function calculates the genetic distance matrix from a given gene identity matrix.
#'
#' @param J The Gene Identity Matrix created using [pooledpeaks::GeneIdentityMatrix]
#'
#' @return The Genetic Distance Matrix
#' @export
#'
#' @examples

GeneticDistanceMatrix <- function(J=matrix) {
  D <- rep(0,nrow(J)*ncol(J))
  D <- array(D,dim=c(nrow(J),ncol(J)))
  rownames(D) <- rownames(J)
  colnames(D) <- colnames(J)

  for (i in 1:nrow(J))
    for (j in 1:ncol(J))
      D[i,j] <- (J[i,i]+J[j,j]-2*J[i,j])/2
  return(D)
}

#' Random Walk Covariance Distance Matrix
#'
#' This function calculates the RWC (Random Walk Covariance) distance matrix from a given matrix of genetic distances.
#'
#' @param J The Genetic Distance Matrix calculated using [pooledpeaks::GeneticDistanceMatrix]
#'
#' @return A matrix representing the distance matrix calculated using the Random Walk Covariance method.
#' @export
#'
#' @examples

RWCDistanceMatrix <- function(J=matrix){
  D <- rep(0,nrow(J)*ncol(J))
  D <- array(D,dim=c(nrow(J),ncol(J)))
  rownames(D) <- rownames(J)
  colnames(D) <- colnames(J)

  for (i in 1:nrow(J))
    for (j in 1:ncol(J)) {
      D[i,j] <- (J[i,i]+J[j,j]-2*J[i,j])/2
      D[i,j] <- D[i,j]/(1-J[i,j])
    }
  return(D)
}
