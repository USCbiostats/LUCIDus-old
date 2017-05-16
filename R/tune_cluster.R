#' Parallel Grid Search for Tuning Parameters in Latent Cluster Analysis
#'
#' \code{tune_cluster} estimates latent clusters of genetic effects using complete biomarker data with/without disease outcomes. Options to produce sparse solutions for cluster-specific parameter estimates under a circumstance of analyzing high-dimensional data are also provided.
#' @param G Genetic effects, a matrix
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param K Pre-specified # of latent clusters, default is 2
#' @param family "binary" or "normal" for Y
#' @param useY Using Y or not, default is TRUE
#' @param Rho_G Penalty for selection on genetic data, numeric, default is -9 using a sequence of penalties
#' @param Rho_Z_InvCov Penalty for the inverse of covariance of biomarkers, numeric, default is 0
#' @param Rho_Z_CovMu Penalty for the product of covariance and mean of biomarkers, numeric, default is 0
#' @param MAX_ITR Maximum number of iterations, integer, default is 100
#' @param MAX_TOT_ITR Maximum number of total iterations, integer, default is 10000
#' @param reltol Convergence cut-off using a relative tolerance, default is 1e-8
#' @param Pred Flag to compute predicted disease probability with fitted model, boolean, default is FALSE
#' @keywords Tunning Parameter, Grid-search
#' @return \code{tune_cluster} returns an object of list containing Modelfits, Results, and Optimal:
#' \item{Modelfits}{Latent cluster model fits for a combination of given tuning parameters}
#' \item{Results}{Summary results of grid-search}
#' \item{Optimal}{Fatures of the optimal model with minimun BIC in the grid-search results}
#' @importFrom mvtnorm dmvnorm
#' @importFrom nnet multinom
#' @importFrom glmnet glmnet
#' @importFrom glasso glasso
#' @importFrom lbfgs lbfgs
#' @importFrom Matrix bdiag
#' @importFrom stats kmeans
#' @importFrom stats runif
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats sd
#' @importFrom stats dnorm
#' @import doParallel
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Yang, Z., Conti, D.V.
#' @examples
#' # For a testing dataset with 10 genetic features (5 causal) and 4 biomarkers (2 causal)
#' est_cluster(G=G1,Z=Z1,Y=Y1,K=2,
#'             init_b = NULL, init_m = NULL, init_s = NULL, init_g = NULL,
#'             family="binary",Pred=TRUE,Select_G=TRUE,Select_Z=TRUE,
#'             Rho_G=0.02,Rho_Z_InvCov=0.1,Rho_Z_CovMu=93,
#'             tol_m = 1e-8,tol_b=1e-8,tol_s=1e-8,tol_g=1e-8,MAX_ITR = 800,MAX_TOT_ITR=800)
library(doParallel)

tune_cluster <- function(G, Z, Y, K, Family,
                         LRho_g, URho_g, NoRho_g,
                         LRho_z_invcov, URho_z_invcov, NoRho_z_invcov,
                         LRho_z_covmu, URho_z_covmu, NoRho_z_covmu){

  parallel_cluster <- function(Lrho_g, Urho_g, Norho_g, Lrho_z_invcov, Urho_z_invcov, Norho_z_invcov, Lrho_z_covmu, Urho_z_covmu, Norho_z_covmu){
    foreach(rho_g = seq(Lrho_g, Urho_g, length.out=Norho_g)) %:%
      foreach(rho_z_invcov = seq(Lrho_z_invcov, Urho_z_invcov, length.out=Norho_z_invcov)) %:%
        foreach(rho_z_covmu = seq(Lrho_z_covmu, Urho_z_covmu, length.out=Norho_z_covmu), .combine = list, .multicombine = TRUE,
                .export=c("est_cluster", "G", "Z", "Y", "K", "Family", "start_b", "start_m", "start_s", "start_g"),
                .packages = c("glmnet", "glasso", "mvtnorm", "nnet", "lbfgs", "stats", "Matrix"))  %dopar%{
                  est_cluster(G=G,Z=Z,Y=Y,K=K,
                              init_b = start_b, init_m = start_m, init_s = start_s, init_g = start_g,
                              family=Family,Pred=TRUE,Select_G=TRUE,Select_Z=TRUE,
                              Rho_G=rho_g,Rho_Z_InvCov=rho_z_invcov,Rho_Z_CovMu=rho_z_covmu,
                              tol_m = 1e-8,tol_b=1e-8,tol_s=1e-8,tol_g=1e-8,MAX_ITR = 800,MAX_TOT_ITR=800)
                }
  }

  parallel_results <- function(Norho_g, Norho_z_invcov, Norho_z_covmu){
    foreach(e=1:Norho_g, .combine = 'rbind') %:%
      foreach(f=1:Norho_z_invcov, .combine = 'rbind') %:%
        foreach(g=1:Norho_z_covmu, .combine = 'rbind', .export=c("modelfits", "K", "G", "Z", "Y")) %dopar%{
          G_diff <- apply(apply(modelfits[[e]][[f]][[g]]$beta,2,range),2,function(x){x[2]-x[1]})[-1]
          select_G <- G_diff != 0

          InvSigmaMu <- solve(modelfits[[e]][[f]][[g]]$sigma[[1]])%*%modelfits[[e]][[f]][[g]]$mu[1,]
          for(j in 2:K){
            InvSigmaMu <- cbind(InvSigmaMu, solve(modelfits[[e]][[f]][[g]]$sigma[[j]])%*%modelfits[[e]][[f]][[g]]$mu[j,])
          }

          Z_diff <- apply(apply(InvSigmaMu,1,range),2,function(x){x[2]-x[1]})
          select_Z <- Z_diff != 0

          G_select <- NULL
          Z_select <- NULL

          if(!all(select_G==FALSE)){
            G_select <- matrix(G[,select_G],ncol=sum(select_G))
          }
          if(!all(select_Z==FALSE)){
            Z_select <- matrix(Z[,select_Z],ncol=sum(select_Z))
          }

          Non0g <- sum(select_G)
          Non0z <- sum(select_Z)
          model_LL <- sum(log(rowSums(modelfits[[e]][[f]][[g]]$Likelihood)))
          Nparm <- ifelse(is.null(G_select),0,(ncol(G_select)+1)*(K-1)) + ifelse(is.null(Z_select),0,ncol(Z_select)*K + ncol(Z_select)*(ncol(Z_select)+1)/2*K) + K
          bic <- -2*model_LL + Nparm*log(length(Y))
          data.frame(Rho_G=modelfits[[e]][[f]][[g]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]][[g]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]][[g]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic)
        }
  }


  modelfits <- parallel_cluster(Lrho_g = LRho_g, Urho_g = URho_g, Norho_g = NoRho_g,
                                Lrho_z_invcov = LRho_z_invcov, Urho_z_invcov = URho_z_invcov, Norho_z_invcov = NoRho_z_invcov,
                                Lrho_z_covmu = LRho_z_covmu, Urho_z_covmu = URho_z_covmu, Norho_z_covmu = NoRho_z_covmu)


  results <- parallel_results(Norho_g = NoRho_g, Norho_z_invcov = NoRho_z_invcov, Norho_z_covmu = NoRho_z_covmu)

  #Show the BEST Model with minBIC
  optimal <- results[which.min(results$BIC), ]

  return(list(Modelfits = modelfits, Results = results, Optimal = optimal))
}
