#' Bootstrap method to estimate variability of latent clusters
#'
#' \code{boot_lucid} provides SEs of parameter estimates from a LUCID model through bootstrapping.
#' @param G Genetic features, a matrix
#' @param CoG Covariates to be included in the G->X path
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param CoY Covariates to be included in the X->Y path
#' @param family "binary" or "normal" for Y
#' @param useY Using Y or not, default is TRUE
#' @param K Pre-specified # of latent clusters, default is 2
#' @param initial A list of initial model parameters will be returned for integrative clustering
#' @param itr_tol A list of tolerance settings will be returned for integrative clustering
#' @param tunepar A list of tuning parameters and settings will be returned for integrative clustering
#' @param Pred Flag to compute posterior probability of latent cluster with fitted model, default is TRUE
#' @param R The number of bootstrap replicates, default is 100
#' @param DeltaE Flag to return the difference in parameter estimate across latent clusters, default is TRUE
#' @param NCPUs The number of processes to be used in parallel computing, default is total number of cores minus 1
#' @keywords bootstrap, latent cluster
#' @return \code{boot_lucid} returns an object of list containing a "boot" class object of LUCID fit and a summary of bootstrap results:
#' \item{Bootstrap}{an object of "boot" class after bootstrapping a LUCID model}
#' \item{Results}{A summary of bootstrap includes original estimate, a bias of the bootstrap estimate, standard error of the bootstrap estimate, and three types of bootstrap confidence intervals based on normal approximation, basic, percentile bootstrap methods}
#' @import boot
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
#'
#' Efron, B., Tibshirani, R. Bootstrap Methods for Standard Errors, Confidence Intervals, and Other Measures of Statistical Accuracy. Statist. Sci. 1 (1986), no. 1, 54--75. doi:10.1214/ss/1177013815.
#' @examples
#' \dontrun{
#' boot_lucid(G = G1, CoG = CoG, Z = Z1, Y = Y1, CoY = CoY, useY = TRUE, family = "binary", K = 2, R=500)
#' }

boot_lucid <- function(G = NULL, CoG = NULL, Z = NULL, Y, CoY = NULL, useY = TRUE, family = "binary", K = 2, Pred = TRUE,
                       initial = def_initial(), itr_tol = def_tol(), tunepar = def_tune(),
                       R = 100, DeltaE = TRUE, NCPUs = detectCores()-1){

  M <- dim(G)[2]
  Q <- dim(Z)[2]

  cM <- dim(CoG)[2]
  cY <- dim(CoY)[2]

  Overalldata <- cbind(Y,G,Z,CoG,CoY)

  lucid_par <- function(data, indices) {
    d <- data[indices,] # allows boot to select sample

    if(is.null(CoG)&&is.null(CoY)){
      try_lucid <- try(try_lucid <- est_lucid(G=as.matrix(d[,2:(M+1)]), Z=as.matrix(d[,(M+2):(M+Q+1)]), Y=as.matrix(d[,1]),
                                              K=K, useY=useY, family=family, Pred=Pred,
                                              initial=initial, tunepar=tunepar, itr_tol = itr_tol))

      par_lucid <- c(try_lucid$beta[2,],as.vector(t(try_lucid$mu)),try_lucid$gamma[1:2]) # original parameters

      if(DeltaE == TRUE){
        par_lucid <- c(abs(try_lucid$beta[2,]),abs(diff(try_lucid$mu)),abs(diff(try_lucid$gamma[1:2]))) # delta Z&Y
      }
    }

    if(is.null(CoG)&&!is.null(CoY)){
      try_lucid <- try(try_lucid <- est_lucid(G=as.matrix(d[,2:(M+1)]), Z=as.matrix(d[,(M+2):(M+Q+1)]), Y=as.matrix(d[,1]), CoY=as.matrix(d[,(M+Q+2):(cY+M+Q+1)]),
                                              K=K, useY=useY, family=family, Pred=Pred,
                                              initial=initial, tunepar=tunepar, itr_tol = itr_tol))

      par_lucid <- c(try_lucid$beta[2,],as.vector(t(try_lucid$mu)),try_lucid$gamma[1:2]) # original parameters

      if(DeltaE == TRUE){
        par_lucid <- c(abs(try_lucid$beta[2,]),abs(diff(try_lucid$mu)),abs(diff(try_lucid$gamma[1:2]))) # delta Z&Y
      }
    }

    if(!is.null(CoG)&&is.null(CoY)){
      try_lucid <- try(try_lucid <- est_lucid(G=as.matrix(d[,2:(M+1)]), Z=as.matrix(d[,(M+2):(M+Q+1)]), Y=as.matrix(d[,1]), CoG=as.matrix(d[,(M+Q+2):(cM+M+Q+1)]),
                                              K=K, useY=useY, family=family, Pred=Pred,
                                              initial=initial, tunepar=tunepar, itr_tol = itr_tol))

      par_lucid <- c(try_lucid$beta[2,],as.vector(t(try_lucid$mu)),try_lucid$gamma[1:2]) # original parameters

      if(DeltaE == TRUE){
        par_lucid <- c(abs(try_lucid$beta[2,]),abs(diff(try_lucid$mu)),abs(diff(try_lucid$gamma[1:2]))) # delta Z&Y
      }
    }

    if(!is.null(CoG)&&!is.null(CoY)){
      try_lucid <- try(try_lucid <- est_lucid(G=as.matrix(d[,2:(M+1)]), Z=as.matrix(d[,(M+2):(M+Q+1)]), Y=as.matrix(d[,1]),
                                              CoG=as.matrix(d[,(M+Q+2):(cM+M+Q+1)]), CoY=as.matrix(d[,(cM+M+Q+2):(cY+cM+M+Q+1)]),
                                              K=K, useY=useY, family=family, Pred=Pred,
                                              initial=initial, tunepar=tunepar, itr_tol = itr_tol))

      par_lucid <- c(try_lucid$beta[2,],as.vector(t(try_lucid$mu)),try_lucid$gamma[1:2]) # original parameters

      if(DeltaE == TRUE){
        par_lucid <- c(abs(try_lucid$beta[2,]),abs(diff(try_lucid$mu)),abs(diff(try_lucid$gamma[1:2]))) # delta Z&Y
      }
    }

    return(par_lucid)
  }

  # Parallel Realization - Bootstrap

  # # Calculate the number of cores
  # no_cores <- detectCores()
  # # Initiate cluster
  # cl <- makeCluster(no_cores)
  #
  # #Start parallel computing
  # registerDoParallel(cl)
  # getDoParWorkers()

  #Set # of bootstrap replicates
  BootRun <- boot::boot(data=Overalldata, statistic=lucid_par, R=R, parallel="multicore", ncpus=NCPUs)

  # stopCluster(cl)

  if(is.null(CoG)){
    BootResults <- as.matrix(cbind(as.vector(BootRun$t0)[2:length(BootRun$t0)],
                                   as.vector(colMeans(BootRun$t)-BootRun$t0)[2:length(BootRun$t0)],
                                   as.vector(apply(BootRun$t,2,sd)[2:length(BootRun$t0)]),
                                   t(sapply(lapply(2:length(BootRun$t0), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$normal[2:3]), cbind)),
                                   t(sapply(lapply(2:length(BootRun$t0), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$basic[4:5]), cbind)),
                                   t(sapply(lapply(2:length(BootRun$t0), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$perc[4:5]), cbind))))

    colnames(BootResults) <- c("original", "bias", "std. error", "lower.norm", "upper.norm", "lower.basic", "upper.basic", "lower.perc", "upper.perc")
    rownames(BootResults) <- c(colnames(G), colnames(Z), "Y")

  }else{
    BootResults <- as.matrix(cbind(as.vector(BootRun$t0)[c(2:(M+1), (cM+M+2):length(BootRun$t0))],
                                   as.vector(colMeans(BootRun$t)-BootRun$t0)[c(2:(M+1), (cM+M+2):length(BootRun$t0))],
                                   as.vector(apply(BootRun$t,2,sd)[c(2:(M+1), (cM+M+2):length(BootRun$t0))]),
                                   t(sapply(lapply(c(2:(M+1), (cM+M+2):length(BootRun$t0)), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$normal[2:3]), cbind)),
                                   t(sapply(lapply(c(2:(M+1), (cM+M+2):length(BootRun$t0)), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$basic[4:5]), cbind)),
                                   t(sapply(lapply(c(2:(M+1), (cM+M+2):length(BootRun$t0)), function(i) (boot::boot.ci(BootRun, type = c("norm", "basic", "perc"), index = i))$perc[4:5]), cbind))))

    colnames(BootResults) <- c("original", "bias", "std. error", "lower.norm", "upper.norm", "lower.basic", "upper.basic", "lower.perc", "upper.perc")

    rownames(BootResults) <- c(colnames(G), colnames(Z), "Y")

  }



  return(list(Results=BootResults, Bootstrap=BootRun))

}
