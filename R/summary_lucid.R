#' Plot Sankey diagram for integrative clustering
#'
#'\code{summary_lucid} generates a summary for the results of integrative clustering based on an \code{IntClust} object.
#' @param x An \code{IntClust} class object
#' @param ... Additional parameters
#' @return \code{summary_lucid} returns a list containing important outputs from an \code{IntClust} object.
#' \item{Beta}{Estimates of genetic effects, matrix}
#' \item{Mu}{Estimates of cluster-specific biomarker means, matrix}
#' \item{Gamma}{Estimates of cluster-specific disease risk, vector}
#' \item{select_G}{A logical vector indicates non-zero genetic features}
#' \item{select_Z}{A logical vector indicates non-zero bio-features}
#' \item{No0G}{A total # of non-zero genetic features}
#' \item{No0Z}{A total # of non-zero bio-features}
#' \item{BIC}{Model BIC}
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Yang, Z., Conti, D.V.

summary_lucid <- function(x, ...) {
  K <- x$K
  family <- x$family

  Beta <- x$beta
  Mu <- x$mu
  Gamma <- x$gamma[1:K]
  Pred <- x$pred

  colnames(Beta) <- c("Int", x$Gnames)
  rownames(Beta) <- paste0("Cluster", 1:K)
  colnames(Mu) <- x$Znames
  rownames(Mu) <- paste0("Cluster", 1:K)
  names(Gamma) <- paste0("Cluster", 1:K)

  G_diff <- apply(apply(x$beta,2,range),2,function(x){x[2]-x[1]})[-1]
  select_G <- G_diff != 0
  InvSigmaMu <- solve(x$sigma[[1]])%*%x$mu[1,]
  for(i in 2:K){
    InvSigmaMu <- cbind(InvSigmaMu, solve(x$sigma[[i]])%*%x$mu[i,])
  }
  Z_diff <- apply(apply(InvSigmaMu,1,range),2,function(x){x[2]-x[1]})
  select_Z <- Z_diff != 0
  No0G <- sum(select_G)
  No0Z <- sum(select_Z)

  model_LL <- sum(log(rowSums(x$Likelihood)))
  Nparm <- (No0G+1)*(K-1) + (No0Z*K + No0Z*(No0Z+1)/2*K) + K*2
  BIC <- -2*model_LL + Nparm*log(nrow(x$pred))

  SumResults <- list(Beta, Mu, Gamma, select_G, select_Z, No0G, No0Z, BIC)

  if(family == "normal"){
    names(SumResults) <- c("Beta", "Mu", "Gamma", "select_G", "select_Z", "No0G", "No0Z", "BIC")
  }
  if(family == "binary"){
    names(SumResults) <- c("Beta", "Mu", "OR_Gamma", "select_G", "select_Z", "No0G", "No0Z", "BIC")
  }
  return(SumResults)
}
