#' Distance Correlation
#'
#' Calculate the correlation between expected and realized genetic distances and plot them.
#'
#' @param GD A matrix containing the genetic distance data.
#'
#' @importFrom stats cophenetic
#' @importFrom stats cor
#' @importFrom ape nj
#' @importFrom graphics par
#' @importFrom graphics abline
#' @importFrom graphics title
#'
#' @return A plot showing the Expected Genetic Distance versus Realized Genetic Distance
#' @export
#'
#' @examples
DistCor <- function(GD=matrix) {
  T <- ape::nj(GD)
  ED <- stats::cophenetic(T)
  graphics::par(mar=c(5,5,5,5))
  plot(GD,ED)
  graphics::abline(0,1,col=2)

  k<-0
  x<-0
  y<-0
  for (i in 1:(nrow(GD)-1))
    for (j in 2:ncol(GD))
    {k<-k+1
    x[k]<- GD[i,j]
    y[k]<- ED[i,j] }
  print(stats::cor(x,y),las=1)
  graphics::title(main='Expected Genetic Distance \n versus Realized Genetic Distance',
        font.main=1,cex.main=1.25)
}




#' Calculate Empirical Standard Error
#'
#' This function calculates the empirical standard error based on repeated sampling.
#'
#' @param datafile A data frame containing genetic data from the [pooledpeaks::LoadData]
#' @param NLoci Number of loci to sample in each iteration.
#'
#' @return A numeric vector containing the empirical standard error estimates.
#' @export
#'
#' @examples

EmpiricalSE <- function(datafile=data.frame,NLoci=10) {
  G<-0
  for ( i in 1:100)
  { if(i %% 10==0) print(i)
    A <- SampleOfLoci(datafile,NLoci)
    names <- colnames(A)
    names <- names[-(1:2)]

    PopLabel <- substr(names,1,1)
    PopFactor <- as.factor(PopLabel)

    Nx <- TypedLoci(A)
    Jx <- GeneIdentityMatrix(A,Nx)
    G[i] <- GST(Jx)
  }
  return(G)
}




#' Calculate Pre-Jost's D
#'
#' This function calculates the pre-Jost's D measure from a genetic distance matrix.
#'
#' @param G A square matrix representing a genetic distance matrix.
#'
#' @return The Jost's D value.

preJostD <- function(G=matrix) {

  JS <- 0
  for (i in 1:nrow(G)) JS <- JS+G[i,i]
  JS <- JS/nrow(G)

  JT <-0
  for (i in 1:nrow(G))
    for (j in 1:nrow(G))
      JT <- JT+G[i,j]
  JT <- JT/nrow(G)/nrow(G)
  D <- (JT/JS -1)/(1/nrow(G)-1)
  return(D)
}




#' Calculate Jost's D
#'
#' This function calculates Jost's D measure from a genetic distance matrix.
#'
#' @param J A genetic distance matrix.
#' @param pairwise Logical indicating whether to calculate pairwise Jost's D.
#'
#' @return If pairwise = TRUE, returns a matrix of pairwise Jost's D values.
#' If pairwise = FALSE, returns the overall Jost's D value.
#'
#' @export
#'
#' @examples

JostD <- function(J=matrix, pairwise=TRUE) {
  if (pairwise==TRUE){
    N <- nrow(J)
    PJostD <- array(dim=c(N,N))
    for (i in 1:N)
      for (j in 1:N) {
        G <- J[c(i,j),c(i,j)]
        PJostD[i,j] <- preJostD(G)
      }
    return(PJostD)}

  if (pairwise==FALSE)
  {preJostD()}

}


#' Pre GST Calculation
#'
#' This function calculates the GST from a genetic distance matrix.
#'
#' @param G The genetic distance matrix
#'
#' @return The GST value.

preGST <- function(G=matrix) {
  JS <- 0
  for (i in 1:nrow(G)) JS <- JS+G[i,i]
  JS <- JS/nrow(G)

  JT <-0
  for (i in 1:nrow(G))
    for (j in 1:nrow(G))
      JT <- JT+G[i,j]
  JT <- JT/nrow(G)/nrow(G)

  GST <- (JS-JT)/(1-JT)
  return(GST)
}

#' Nei's GST
#'
#' This function calculates GST (Nei's standard genetic distance) measure from a genetic distance matrix.
#'
#' @param J A square matrix representing a genetic distance matrix.
#' @param pairwise Logical indicating whether to calculate pairwise GST.
#'
#' @return If pairwise = TRUE, returns a matrix of pairwise GST values.
#' If pairwise = FALSE, returns the overall GST value.
#' @export
#'
#' @examples

GST <- function(J=matrix, pairwise=TRUE) {
  if (pairwise==TRUE){
    N <- nrow(J)
    PGST <- array(dim=c(N,N))
    for (i in 1:N)
      for (j in 1:N) {
        G <- J[c(i,j),c(i,j)]
        PGST[i,j] <- preGST(G)
      }
    return(PGST)}

  if (pairwise==FALSE)
  {preGST()}

}

#' Calculate Two-Level GST
#'
#' This function calculates two-level GST (Nei's standard genetic distance) measure from a genetic distance matrix.
#'
#' @param G A square matrix representing a genetic distance matrix.
#'
#' @return A list containing the components of two-level GST including within-group gene identity, between-group gene identity, and GST values.
#' @export
#'
#' @examples

TwoLevelGST <- function(G=matrix) {
  names <- colnames(G)
  names <- names[-(1:2)]
  PopLabel <- substr(names,1,1)
  PF <- as.factor(PopLabel)
  z <- levels(PF)
  Gsc <- rep(0,length(z))
  Jc <- rep(0,length(z))
  n <- 0
  Jx <- 0  # Gene Identity Matrix for Component
  Js <- rep(0,length(z)) # Gene Identity Within

  for (k in 1:length(z)) {
    group <- which(PF==z[k])
    Jx <- G[group,group]
    n[k] <- nrow(Jx)
    Js[k]<-0
    for (i in 1:length(group))
      Js[k] <- Js[k]+Jx[i,i]
    Js[k] <- Js[k]/length(group)

    for (i in 1:length(group))
      for (j in 1:length(group))
        Jc[k] <- Jc[k]+Jx[i,j]

    Jc[k] <- Jc[k]/length(group)/length(group)

    Gsc[k] <- (Js[k]-Jc[k])/(1-Jc[k])
  }

  JS <- 0
  for (i in 1:nrow(G))
    JS <- JS+G[i,i]
  JS <- JS/nrow(G)

  JC <-0
  for (k in 1:length(z))
    JC <- JC + n[k]*Jc[k]
  JC <- JC/sum(n)

  JT <-0
  for (i in 1:nrow(G))
    for (j in 1:nrow(G))
      JT <- JT + G[i,j]
  JT <- JT/nrow(G)/nrow(G)

  GG <- list(Js=Js,Jc=Jc,Gsc=Gsc,JS=JS,JC=JC,JT=JT,
             GSC=(JS-JC)/(1-JC),GCT=(JC-JT)/(1-JT),GST=(JS-JT)/(1-JT))
  return(GG)
}


#' Calculate Allelic Richness
#'
#' This function calculates allelic richness based on the provided genetic data.
#'
#' @param datafile A data frame containing the data as read in by [pooledpeaks::LoadData]
#' @param n Vector specifying the number of alleles per locus.
#'
#' @return A vector containing the allelic richness for each locus.
#' @export
#'
#' @examples

AlRich <- function(datafile=data.frame,n=c()) {

  loci <- rep(0,nrow(n))
  loci <- diag(n)

  g <- subset(datafile,(datafile[,2] != 'n')&(datafile[,3]!='NA'))
  g <- g[,-(1:2)]
  g <- as.matrix(g)

  allnum <- rep(0,ncol(g))

  for (i in 1:ncol(g)) {
    h <- g[,i]
    h <- subset(h,h>0)
    allnum[i] <- length(h)
  }

  allnum <- allnum/loci

  return(allnum)
}

#' Perform Bootstrap Analysis
#'
#' This function performs bootstrap analysis on genetic data.
#'
#' @param A A data frame containing the data as read in by [pooledpeaks::LoadData]
#' @param Rep Number of bootstrap replicates.
#' @param Stat Type of statistic to compute (1 for AlRich, 2 for TwoLevelGST)
#'
#' @return Either a matrix of AlRich statistics or a list containing various statistics computed using TwoLevelGST.
#' @export
#'
#' @examples

BootStrap3 <- function(A=data.frame,Rep=20,Stat=1) {
  A <- A
  Rep <- 5
  Stat=2

  Loci <- max(A[,1])

  names <- colnames(A)
  names <- names[-(1:2)]

  PopLabel <- substr(names,1,1)
  PopFactor <- as.factor(PopLabel)
  K <- length(levels(PopFactor))

  gst <- 0

  if (Stat==1) w <- rep(0,length(A[1,])-2)

  if (Stat==2) {
    js <- rep(0,Rep*K)
    js <- array(js,dim=c(Rep,K))
    jc <- rep(0,Rep*K)
    jc <- array(jc,dim=c(Rep,K))
    gsc <- rep(0,Rep*K)
    gsc <- array(gsc,dim=c(Rep,K))
    JS <- rep(0,Rep)
    JC <- 0
    JT <- 0
    GSC <-0
    GCT <- 0
    GST <- 0
  }

  for (j in 1:Rep) {
    if (j%%10==0) print(j)
    z <-A[1,]
    x <- sample(1:Loci,Loci,replace=TRUE)
    for (i in x) {
      b <- subset(A,A[,1]==i)
      z <- rbind(z,b)
    }
    z <- z[-1,]


    if (Stat==1) {
      N <- TypedLoci(z)
      v <- AlRich(z,N)
      if (j==1) w <- v
      if (j>1)  w <- rbind(w,v)
    }  # End of Stat = 1 loop

    if (Stat==2) {
      nn <- TypedLoci(z)
      Jx<-GeneIdentityMatrix(z,nn)
      G <- TwoLevelGST(Jx)
      for (k in 1:length(G$Js)) js[j,k] <- G$Js[k]
      for (k in 1:length(G$Jc)) jc[j,k] <- G$Jc[k]
      for (k in 1:length(G$Gsc)) gsc[j,k] <- G$Gsc[k]
      JS[j] <- G$JS
      JC[j] <- G$JC
      JT[j] <- G$JT
      GST[j] <-G$GST
      GSC[j] <- (JS[j]-JC[j])/(1-JC[j])
      GCT[j] <- (JC[j] -JT[j] )/(1-JT[j])

    }  # end of stat=2 loop
  }   # end of Replicates loop

  if (Stat==1) send_out <- w
  if (Stat == 2) send_out <- list(Js=js,Jc=jc,Gsc=gsc,JS=JS,JC=JC,JT=JT,
                                  GSC=GSC,GCT=GCT,GST=GST)

  return(send_out)

}  # End of Function Bootstrap3
