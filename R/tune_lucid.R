#' Parallel Grid Search for Tuning Parameters in Latent Cluster Analysis
#'
#' \code{tune_lucid} fits regularized latent cluster models with various combinations of three tuning parameters based on joint inference across data types to perform a grid-search helping determine an optimal choice of three tuning parameters with minimum model BIC.
#' @param G Genetic effects, a matrix
#' @param CoG Covariates to be added in G->X path
#' @param Z Biomarker data, a matrix
#' @param CoY Covariates to be added in X->Y path
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
#' @keywords Tuning Parameter, Grid-search
#' @return \code{tune_lucid} returns an object of list containing Modelfits, Results, and Optimal:
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
#' # Parallel grid-search with 100 combinations of tuning parameters
#' GridSearch <- tune_lucid(G=G1, Z=Z1, Y=Y1, K=2, Family="binary", USEY = TRUE, NoCores = 2,
#'                          LRho_g = 0.008, URho_g = 0.012, NoRho_g = 3,
#'                          LRho_z_invcov = 0.04, URho_z_invcov = 0.06, NoRho_z_invcov = 3,
#'                          LRho_z_covmu = 90, URho_z_covmu = 100, NoRho_z_covmu = 3)
#' GridSearch$Results
#' # Determine the best tuning parameters
#' GridSearch$Optimal

tune_lucid <- function(G = NULL, CoG = NULL, Z = NULL, CoY = NULL, Y, K, Family, USEY = TRUE,
                         initial = def_initial(),
                         LRho_g, URho_g, NoRho_g,
                         LRho_z_invcov, URho_z_invcov, NoRho_z_invcov,
                         LRho_z_covmu, URho_z_covmu, NoRho_z_covmu,
                         NoCores = detectCores()-1){

  if(getRversion() >= "2.15.1")  utils::globalVariables(c("e","f","g","rho_g","rho_z_covmu","rho_z_invcov"))

  M <- dim(G)[2]
  Q <- dim(Z)[2]

  parallel_cluster <- function(Lrho_g, Urho_g, Norho_g, Lrho_z_invcov, Urho_z_invcov, Norho_z_invcov, Lrho_z_covmu, Urho_z_covmu, Norho_z_covmu){
    foreach(rho_g = seq(Lrho_g, Urho_g, length.out=Norho_g)) %:%
      foreach(rho_z_invcov = seq(Lrho_z_invcov, Urho_z_invcov, length.out=Norho_z_invcov)) %:%
        foreach(rho_z_covmu = seq(Lrho_z_covmu, Urho_z_covmu, length.out=Norho_z_covmu),
                .combine = list, .multicombine = TRUE, .maxcombine = 2000, .errorhandling = 'pass',
                .export=c("G", "CoG", "Z", "CoY", "Y", "K", "Family", "USEY", "initial"),
                .packages = c("glmnet", "glasso", "mvtnorm", "nnet", "lbfgs", "stats", "Matrix", "LUCid"))  %dopar%{
                  set.seed(rho_g*rho_z_invcov*rho_z_covmu)
                  est_lucid(G=G,CoG=CoG,Z=Z,CoY=CoY,Y=Y,K=K,useY=USEY,family=Family,Pred=TRUE,
                              initial = initial, tunepar = def_tune(Select_G=T,Select_Z=T,Rho_G=rho_g,Rho_Z_InvCov=rho_z_invcov,Rho_Z_CovMu=rho_z_covmu),
                              def_tol(MAX_ITR = 500,MAX_TOT_ITR=1000))
                }
  }

  parallel_results <- function(Norho_g, Norho_z_invcov, Norho_z_covmu){
    if(Norho_z_covmu==1){
      foreach(e=1:Norho_g, .combine = 'rbind') %:%
        foreach(f=1:Norho_z_invcov, .combine = 'rbind', .errorhandling = 'remove',
                .export=c("modelfits", "K"), .packages=c("stats", "LUCid")) %do%{

                  Non0g <- summary_lucid(modelfits[[e]][[f]])$No0G
                  Non0z <- summary_lucid(modelfits[[e]][[f]])$No0Z

                  bic <- summary_lucid(modelfits[[e]][[f]])$BIC
                  #Other types of GIC
                  gic1 <- summary_lucid(modelfits[[e]][[f]])$GIC1
                  gic2 <- summary_lucid(modelfits[[e]][[f]])$GIC2

                  data.frame(Rho_G=modelfits[[e]][[f]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic, GIC1=gic1, GIC2=gic2)
                }
    }else{
      foreach(e=1:Norho_g, .combine = 'rbind') %:%
        foreach(f=1:Norho_z_invcov, .combine = 'rbind') %:%
          foreach(g=1:Norho_z_covmu, .combine = 'rbind', .errorhandling = 'remove',
                  .export=c("modelfits", "K"), .packages=c("stats", "LUCid")) %do%{

                    Non0g <- summary_lucid(modelfits[[e]][[f]][[g]])$No0G
                    Non0z <- summary_lucid(modelfits[[e]][[f]][[g]])$No0Z

                    bic <- summary_lucid(modelfits[[e]][[f]][[g]])$BIC
                    #Other types of GIC
                    gic1 <- summary_lucid(modelfits[[e]][[f]][[g]])$GIC1
                    gic2 <- summary_lucid(modelfits[[e]][[f]][[g]])$GIC2

                    data.frame(Rho_G=modelfits[[e]][[f]][[g]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]][[g]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]][[g]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic, GIC1=gic1, GIC2=gic2)
                  }
    }
  }

  # Initiate cluster
  cl <- makeCluster(NoCores)
  #Start parallel computing
  registerDoParallel(cl)

  modelfits <- parallel_cluster(Lrho_g = LRho_g, Urho_g = URho_g, Norho_g = NoRho_g,
                                Lrho_z_invcov = LRho_z_invcov, Urho_z_invcov = URho_z_invcov, Norho_z_invcov = NoRho_z_invcov,
                                Lrho_z_covmu = LRho_z_covmu, Urho_z_covmu = URho_z_covmu, Norho_z_covmu = NoRho_z_covmu)

  stopCluster(cl)

  results <- parallel_results(Norho_g = NoRho_g, Norho_z_invcov = NoRho_z_invcov, Norho_z_covmu = NoRho_z_covmu)

  #Show the BEST Model with minBIC
  optimal <- results[which.min(results$BIC), ]

  return(list(Modelfits = modelfits, Results = results, Optimal = optimal))
}
