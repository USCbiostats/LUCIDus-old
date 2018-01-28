#' Parallel Grid Search for Tuning Parameters in Latent Cluster Analysis
#'
#' \code{tune_cluster} fits regularized latent cluster models with various combinations of three tuning parameters based on joint inference across data types to perform a grid-search helping determine an optimal choice of three tuning parameters with minimum model BIC.
#' @param G Genetic effects, a matrix
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param K Pre-specified # of latent clusters
#' @param Family "binary" or "normal" for Y
#' @param USEY Using Y or not, default is TRUE
#' @param initial A list of initial model parameters will be returned for integrative clustering
#' @param LRho_g Lower limit of the penalty for selection on genetic data
#' @param URho_g Upper limit of the penalty for selection on genetic data
#' @param NoRho_g Number of \code{Rho_g} for grid-search
#' @param LRho_z_invcov Lower limit of the penalty for the inverse of covariance of biomarkers
#' @param URho_z_invcov Upper limit of the penalty for the inverse of covariance of biomarkers
#' @param NoRho_z_invcov Number of \code{Rho_z_invcov} for grid-search
#' @param LRho_z_covmu Lower limit of the penalty for the product of covariance and mean of biomarkers
#' @param URho_z_covmu Upper limit of the penalty for the product of covariance and mean of biomarkers
#' @param NoRho_z_covmu Number of \code{Rho_z_covmu} for grid-search
#' @param NoCores Number of CPU cores for parallel grid-search, default is total number of cores minus 1
#' @keywords Tunning Parameter, Grid-search
#' @return \code{tune_cluster} returns an object of list containing Modelfits, Results, and Optimal:
#' \item{Modelfits}{Latent cluster model fits for a combination of given tuning parameters}
#' \item{Results}{Summary results of grid-search}
#' \item{Optimal}{Features of the optimal model with minimun BIC in the grid-search summary}
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
#' @import parallel
#' @import foreach
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Yang, Z., Conti, D.V.
#' @examples
#' # For a testing dataset with 10 genetic features (5 causal) and 4 biomarkers (2 causal)
#' # Grid-search with 30 combinations of tuning parameters
#' GridSearch <- tune_cluster(G=G1, Z=Z1, Y=Y1, K=2, Family="binary", USEY = TRUE,
#'                            LRho_g = 0.008, URho_g = 0.01, NoRho_g = 3,
#'                            LRho_z_invcov = 0.1, URho_z_invcov = 0.2, NoRho_z_invcov = 2,
#'                            LRho_z_covmu = 86, URho_z_covmu = 90, NoRho_z_covmu = 5)
#' GridSearch$Results
#' GridSearch$Optimal

tune_cluster <- function(G = NULL, Z = NULL, Y, K, Family, USEY = TRUE,
                         initial = def_initial(),
                         LRho_g, URho_g, NoRho_g,
                         LRho_z_invcov, URho_z_invcov, NoRho_z_invcov,
                         LRho_z_covmu, URho_z_covmu, NoRho_z_covmu,
                         NoCores = detectCores()-1){

  M <- dim(G)[2]
  Q <- dim(Z)[2]

  parallel_cluster <- function(Lrho_g, Urho_g, Norho_g, Lrho_z_invcov, Urho_z_invcov, Norho_z_invcov, Lrho_z_covmu, Urho_z_covmu, Norho_z_covmu){
    foreach(rho_g = seq(Lrho_g, Urho_g, length.out=Norho_g)) %:%
      foreach(rho_z_invcov = seq(Lrho_z_invcov, Urho_z_invcov, length.out=Norho_z_invcov)) %:%
        foreach(rho_z_covmu = seq(Lrho_z_covmu, Urho_z_covmu, length.out=Norho_z_covmu),
                .combine = list, .multicombine = TRUE, .maxcombine = 2000, .errorhandling = 'pass',
                .export=c("est_cluster", "def_initial", "def_tune", "def_tol", "G", "Z", "Y", "K", "Family", "USEY", "initial"),
                .packages = c("glmnet", "glasso", "mvtnorm", "nnet", "lbfgs", "stats", "Matrix"))  %dopar%{
                  est_cluster(G=G,Z=Z,Y=Y,K=K,useY=USEY,family=Family,Pred=TRUE,
                              initial = initial, tunepar = def_tune(Select_G=T,Select_Z=T,Rho_G=rho_g,Rho_Z_InvCov=rho_z_invcov,Rho_Z_CovMu=rho_z_covmu),
                              def_tol(MAX_ITR = 500,MAX_TOT_ITR=1000))
                }
  }

  parallel_results <- function(Norho_g, Norho_z_invcov, Norho_z_covmu){
    foreach(e=1:Norho_g, .combine = 'rbind') %:%
      foreach(f=1:Norho_z_invcov, .combine = 'rbind') %:%
        foreach(g=1:Norho_z_covmu, .combine = 'rbind', .errorhandling = 'remove', .export=c("modelfits", "K", "G", "Z", "Y")) %dopar%{
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
            G_select <- G[,select_G]
          }
          if(!all(select_Z==FALSE)){
            Z_select <- Z[,select_Z]
          }

          Non0g <- sum(select_G)
          Non0z <- sum(select_Z)

          model_LL <- sum(log(rowSums(modelfits[[e]][[f]][[g]]$Likelihood)))
          #Nparm <- ifelse(is.null(G_select),0,(ncol(G_select)+1)*(K-1)) + ifelse(is.null(Z_select),0,ncol(Z_select)*K + ncol(Z_select)*(ncol(Z_select)+1)/2*K) + K*2
          Nparm <- ifelse(is.null(G),0,(ncol(G)+1)*(K-1)) + ifelse(is.null(Z),0,ncol(Z)*K + ncol(Z)*(ncol(Z)+1)/2*K) + K*2
          bic <- -2*model_LL + Nparm*log(length(Y))
          #Other types of GIC
          gic1 <- -2*model_LL + Nparm*log(log(length(Y)))*log(length(Y))
          gic2 <- -2*model_LL + Nparm*log(log(length(Y)))*log(ncol(G)+ncol(Z))

          data.frame(Rho_G=modelfits[[e]][[f]][[g]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]][[g]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]][[g]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic, GIC1=gic1, GIC2=gic2)
        }
  }

  # Initiate cluster
  cl <- makeCluster(NoCores)
  #Start parallel computing
  registerDoParallel(cl)

  modelfits <- parallel_cluster(Lrho_g = LRho_g, Urho_g = URho_g, Norho_g = NoRho_g,
                                Lrho_z_invcov = LRho_z_invcov, Urho_z_invcov = URho_z_invcov, Norho_z_invcov = NoRho_z_invcov,
                                Lrho_z_covmu = LRho_z_covmu, Urho_z_covmu = URho_z_covmu, Norho_z_covmu = NoRho_z_covmu)


  results <- parallel_results(Norho_g = NoRho_g, Norho_z_invcov = NoRho_z_invcov, Norho_z_covmu = NoRho_z_covmu)

  stopCluster(cl)

  #Show the BEST Model with minBIC
  optimal <- results[which.min(results$BIC), ]

  return(list(Modelfits = modelfits, Results = results, Optimal = optimal))
}
