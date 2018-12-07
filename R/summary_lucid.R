#' Summarize results for integrative clustering
#'
#'\code{summary_lucid} generates a summary for the results of integrative clustering based on an \code{IntClust} object.
#' @param x An \code{IntClust} class object
#' @param switch An indicator to do label switching or not, the default is FALSE
#' @param order A customized order for label switching, a vector with a length of K; the default is NULL, which is a descending order in gamma
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
#' @examples
#' # For a testing dataset with 10 genetic features (5 causal) and 4 biomarkers (2 causal)
#'
#' # Integrative clustering without feature selection
#' set.seed(10)
#' IntClusFit <- est_lucid(G=G1,Z=Z1,Y=Y1,K=2,family="binary",Pred=TRUE)
#'
#' # Check important model outputs
#' summary_lucid(IntClusFit)

summary_lucid <- function(x, switch=FALSE, order=NULL) {
  K <- x$K
  family <- x$family

  Beta <- x$beta
  Mu <- x$mu
  Gamma <- x$gamma[1:K]
  Pred <- x$pred

  if(switch){
    if(is.null(order)){
      Beta <- Beta[order(Gamma),]
      Base_Beta <- Beta[1,]
      Beta <- t(apply(Beta, 1, function(m) return(m-Base_Beta)))
      Mu <- Mu[order(Gamma),]
      Pred <- Pred[,order(Gamma)]
      Gamma <- Gamma[order(Gamma)]
    }
    else{
      Beta <- Beta[order,]
      Base_Beta <- Beta[1,]
      Beta <- t(apply(Beta, 1, function(m) return(m-Base_Beta)))
      Mu <- Mu[order,]
      Pred <- Pred[,order]
      Gamma <- Gamma[order]
    }
  }

  colnames(Beta) <- c("Int", x$Gnames)
  rownames(Beta) <- paste0("Cluster", 1:K)
  colnames(Mu) <- x$Znames
  rownames(Mu) <- paste0("Cluster", 1:K)
  names(Gamma) <- paste0("Cluster", 1:K)
  colnames(Pred) <- paste0("Cluster", 1:K)

  if(family == "binary"){
    if(is.null(x$YFIT)){
      OR_Gamma <- c(1,exp(Gamma[2:K]-Gamma[1]))
    }else{
      OR_Gamma <- c(1,exp(coef(x$YFIT)[2:K]))
    }
    names(OR_Gamma) <- paste0("Cluster", 1:K)
  }

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

  lambda <- apply(x$pred, 2, mean)
  model_LL <- sum(log(rowSums(lambda*x$Likelihood)))
  Nparm <- (No0G+1)*(K-1) + (No0Z*K + No0Z*(No0Z+1)/2*K) + K*2
  # Nparm <- ifelse(is.null(G),0,(ncol(G)+1)*(K-1)) + ifelse(is.null(Z),0,ncol(Z)*K + ncol(Z)*(ncol(Z)+1)/2*K) + K*2
  BIC <- -2*model_LL + Nparm*log(nrow(x$pred))
  #Other types of GIC
  GIC1 <- -2*model_LL + Nparm*log(log(nrow(x$pred)))*log(nrow(x$pred))
  GIC2 <- -2*model_LL + Nparm*log(log(nrow(x$pred)))*log(ncol(x$beta)-1+ncol(x$mu))

  if(family == "normal"){
    SumResults <- list(Beta, Mu, Gamma, select_G, select_Z, No0G, No0Z, BIC, GIC1, GIC2)
    names(SumResults) <- c("Beta", "Mu", "Gamma", "select_G", "select_Z", "No0G", "No0Z", "BIC", "GIC1", "GIC2")
  }
  if(family == "binary"){
    SumResults <- list(Beta, Mu, OR_Gamma, select_G, select_Z, No0G, No0Z, BIC, GIC1, GIC2)
    names(SumResults) <- c("Beta", "Mu", "OR_Gamma", "select_G", "select_Z", "No0G", "No0Z", "BIC", "GIC1", "GIC2")
  }
  return(SumResults)
}
