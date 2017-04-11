#' Estimating latent clusters with multi-omics data
#'
#' \code{est_cluster} estimates latent clusters of genetic effects using complete biomarker data with/without disease outcomes. Options to produce sparse solutions for cluster-specific parameter estimates under a circumstance of analyzing high-dimensional data are also provided.
#' @param G Genetic effects, a matrix
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param family "binary" or "normal" for Y
#' @param useY Using Y or not, default is TRUE
#' @param K Pre-specified # of latent clusters, default is 2
#' @param Select_G Flag to do model selection on genetic data, default is FALSE
#' @param Select_Z Flag to do model selection on biomarker data, default is FALSE
#' @param Rho_G Penalty for selection on genetic data, numeric, default is -9 using a sequence of penalties
#' @param Rho_Z_InvCov Penalty for the inverse of covariance of biomarkers, numeric, default is 0
#' @param Rho_Z_CovMu Penalty for the product of covariance and mean of biomarkers, numeric, default is 0
#' @param MAX_ITR Maximum number of iterations, integer, default is 100
#' @param MAX_TOT_ITR Maximum number of total iterations, integer, default is 10000
#' @param reltol Convergence cut-off using a relative tolerance, default is 1e-8
#' @param Pred Flag to compute predicted disease probability with fitted model, boolean, default is FALSE
#' @keywords latent cluster
#' @return \code{est_cluster} returns an object of list containing parameters estimates, predicted probability of latent clusters, and other features:
#' \item{beta}{Estimates of genetic effects, matrix}
#' \item{mu}{Estimates of cluster-specific biomarker means, matrix}
#' \item{sigma}{Estimates of cluster-specific biomarker covariance matrix, list}
#' \item{gamma}{Estimates of cluster-specific disease risk, vector}
#' \item{pcluster}{Probability of cluster, when G is null}
#' \item{pred}{Predicted probability of belonging to each latent cluster}
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

est_cluster <- function(G=NULL, Z=NULL, Y,
                        family="binary", useY = TRUE,
                        init_b = NULL, init_m = NULL, init_s = NULL, init_g = NULL,init_pcluster=NULL,
                        K = 2, MAX_ITR = 100, MAX_TOT_ITR = 10000, reltol=1e-8,
                        tol_b = 1e-4, tol_m = 1e-4, tol_s = 1e-4, tol_g = 1e-4, tol_p = 1e-4, tol_sem = 1e-4,
                        Pred = FALSE, Select_G = FALSE, Select_Z = FALSE, Rho_G = -9, Rho_Z_InvCov = 0, Rho_Z_CovMu = 0){

  # check input
  if(family != "binary" && family!= "normal"){
    print("family can only be 'binary' or 'linear'...")
    return (list(err = -99))
  }

  #Initialize E-M algorithm
  if(is.null(G)&&is.null(Z)){
    if(useY == FALSE){
      cat("ERROR: Y mush be used when G and Z are empty!","\n")
      return (list(err = -99))
    }
  }

  if(is.null(Y)){
    cat("ERROR: Y cannot be empty!","\n")
    return (list(err = -99))
  }

  N <- length(Y) #sample size
  if(!is.null(G)){
    M <- ncol(G) #dimension of G
    P <- 0 #indicator for whether include pcluster as parameters
    new_beta <- mat.or.vec(K,M+1) #store estimated beta in M-step
    new_pcluster <- NULL
  }else{
    M <- -1
    P <- 1 #indicator for whether include pcluster as parameters
    new_beta <- NULL
    new_pcluster <- rep(0,K)
  }

  if(useY){
    if(family=="normal"){
      new_gamma <- array(0,2*K) #store estimated gamma in M-step
    }
    if(family=="binary"){
      new_gamma <- array(0,K)
    }
  }else{
    new_gamma <- NULL
  }

  if(!is.null(Z)){
    Q <- ncol(Z) #dimension of biomarkers
    new_mu <- mat.or.vec(K,Q) #store estimated mu in M-step
    new_sigma <- replicate(K,diag(0,nrow=Q,ncol=Q),simplify=F) #store estiamted sigma in M-step
  }else{
    Q <- 0
    new_mu <- NULL
    new_sigma <- NULL
  }

  r <- mat.or.vec(N,K) #prob of cluster given other vars
  total_itr <- 0 #total number of iterations

  success <- FALSE #flag for successful convergence

  likelihood <- function(Beta=NULL,Mu=NULL,Sigma=NULL,Gamma=NULL,Family){

    if(!is.null(Gamma)){
      if(!is.null(Beta)){
        pAgG <- exp(as.matrix(cbind(intercept=rep(1,N),G))%*%t(Beta))
        pAgG <- data.frame(t(apply(pAgG,1,function(x)return(x/sum(x)))))
      }

      if(!is.null(Mu)){
        pZgA <- mat.or.vec(N,K)
        for(k in 1:K){
          pZgA[,k] <- apply(Z,1,function(x)return(dmvnorm(x,Mu[k,],Sigma[[k]])))
        }
      }

      pYgA <- mat.or.vec(N,K)

      if(Family == "binary"){
        for(k in 1:K){
          pYgA[,k] <- ((exp(Gamma[k])/(1+exp(Gamma[k])))^Y)*(1/(1+exp(Gamma[k])))^(1-Y)
        }
      }

      if(Family == "normal"){
        for(k in 1:K){
          pYgA[,k] <- dnorm(Y,mean = Gamma[k],sd = sqrt(Gamma[k+K]))
        }
      }

      if(!is.null(Mu)&&is.null(Beta)){
        return (pZgA*pYgA)
      }
      if(is.null(Mu)&&!is.null(Beta)){
        return (pAgG*pYgA)
      }
      if(is.null(Mu)&&is.null(Beta)){
        return (pYgA)
      }
      return (pAgG*pZgA*pYgA)
    }else{
      if(!is.null(Beta)){
        pAgG <- exp(as.matrix(cbind(intercept=rep(1,N),G))%*%t(Beta))
        pAgG <- data.frame(t(apply(pAgG,1,function(x)return(x/sum(x)))))
      }

      if(!is.null(Mu)){
        pZgA <- mat.or.vec(N,K)
        for(k in 1:K){
          pZgA[,k] <- apply(Z,1,function(x)return(dmvnorm(x,Mu[k,],Sigma[[k]])))
        }
      }

      if(!is.null(Mu)&&is.null(Beta)){
        return (pZgA)
      }
      if(is.null(Mu)&&!is.null(Beta)){
        return (pAgG)
      }

      return (pAgG*pZgA)
    }
  }

  #Choose starting point for mu
  if(!is.null(Z)){
    if(is.null(init_m)){

      mu_start <- kmeans(Z,K)$centers

    }
    else{
      mu_start <- as.matrix(init_m)
    }
  }

  while(!success && total_itr < MAX_TOT_ITR){
    #choose starting point
    mu <- NULL
    sigma <- NULL
    gamma <- NULL


    if(!is.null(Z)){
      mu <- mu_start

      if(is.null(init_s)){
        sigma <- replicate(K,diag(runif(1)*2,nrow=Q,ncol=Q),simplify=F) # use diag for now; use more general form in the future
      }
      else{
        sigma <- init_s
      }
    }


    if(useY){
      if(is.null(init_g)){
        if(family == "binary"){
          gamma <- array(runif(K)*2)
        }

        if(family == "normal"){
          gamma <- array(runif(2*K)*2)
          gamma[(K+1):(2*K)] <- abs(gamma[(K+1):(2*K)])
        }
      }
      else{
        gamma <- init_g
      }
    }

    beta <- NULL
    pcluster <- NULL
    if(!is.null(G)){
      if(is.null(init_b)){
        beta <- array(runif((K)*(M+1))*2,dim=c(K,M+1))
      }
      else{
        beta <- init_b
      }
    }else{
      if(is.null(init_pcluster)){
        pcluster <- rep(0,K)
        pcluster[1] <- runif(1)
        if(K>2){
          for(k in 2:(K-1)){
            pcluster[k] <- runif(1)*(1-sum(pcluster[1:(k-1)]))
          }
        }
        pcluster[K] <- 1-sum(pcluster[1:(K-1)])
      }else{
        pcluster <- init_pcluster
      }
    }

    #flag of breakdown
    breakdown <- FALSE

    #flag of convergence
    convergence <- FALSE

    #number of iterations
    itr <- 0

    while(!convergence && !breakdown && itr < MAX_ITR){
      total_itr <- total_itr + 1
      itr <- itr + 1

      #------E-step------#
      r <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=gamma, Family=family)

      if(is.null(G)){
        r <- t(t(r)*pcluster)
      }
      r <- r/rowSums(r)

      #------check r------#
      valid_r <- all(is.finite(as.matrix(r)))

      if(!valid_r){
        breakdown <- TRUE
        print(paste(itr,": invalid r"))
      }
      else{
        print(paste(itr,": E-step finished"))

        #------M-step------#

        # estimate new mu and sigma
        if(!is.null(Z)){
          if(Select_Z){
            k <- 1

            while(k <= K && !breakdown){
              #estimate E(S_k) to be used by glasso
              Z_mu <- t(t(Z)-mu[k,])
              E_S <- (matrix(colSums(r[,k]*t(apply(Z_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(r[,k])

              #use glasso and E(S_k) to estimate new_sigma and new_sigma_inv
              try_glasso <- try(l_cov <- glasso(E_S,Rho_Z_InvCov))

              if("try-error" %in% class(try_glasso)){
                breakdown <- TRUE
                print(paste(itr,": glasso failed"))
              }
              else{
                new_sigma[[k]] <- l_cov$w
                new_sigma_inv <- l_cov$wi
                new_sigma_est <- l_cov$w

                #use lbfgs to estimate mu with L1 penalty
                fn <- function(a){
                  mat <- Z
                  cov_inv <- new_sigma_inv
                  cov <- new_sigma_est
                  Mu <- cov%*%a
                  return(sum(r[,k]*apply(mat,1,function(v) return(t(v-Mu)%*%cov_inv%*%(v-Mu)))))
                }

                gr <- function(a){
                  mat <- Z
                  cov <- new_sigma_est
                  Mu <- cov%*%a
                  return(2.0*apply(r[,k]*t(apply(mat,1,function(v) return(Mu-v))),2,sum))
                }

                try_optim_mu <- try(lbfgs(call_eval=fn,call_grad = gr, vars = rep(0,Q), invisible=1, orthantwise_c = Rho_Z_CovMu))

                if("try-error" %in% class(try_optim_mu)){
                  breakdown <- true
                }
                else{
                  new_mu[k,] <- new_sigma[[k]]%*%(try_optim_mu$par)
                }
              }
              k <- k + 1
            }
          }else{

            new_mu <- matrix(t(apply(r,2,function(x) return(colSums(x*Z)/sum(x)))),ncol = ncol(Z))

            if(ncol(Z) > 1){
              for(k in 1:K){
                Z_mu <- t(t(Z)-new_mu[k,])
                new_sigma[[k]] <- (matrix(colSums(r[,k]*t(apply(Z_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(r[,k])
              }
            }else{
              for(k in 1:K){
                Z_mu <- t(t(Z)-new_mu[k,])
                new_sigma[[k]] <- (matrix(sum(r[,k]*t(apply(Z_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(r[,k])
              }
            }

          }
        }


        if(!is.null(G)){
          #estimate new beta
          if(Select_G){
            if(Rho_G == -9){
              tryLasso <- try(glmnet(as.matrix(G),as.matrix(r),family="multinomial"))
            }
            else{
              tryLasso <- try(glmnet(as.matrix(G),as.matrix(r),family="multinomial",lambda=Rho_G))
            }

            if("try-error" %in% class(tryLasso)){
              breakdown <- TRUE
              print(paste(itr,": lasso failed"))
            }
            else{
              if(Rho_G == -9){
                new_beta[,1] <- tryLasso$a0[,length(tryLasso$lambda)]
                new_beta[,-1] <- t(matrix(unlist(lapply(tryLasso$beta,function(x) return(x[,ncol(x)]))),ncol = K))
              }
              else{
                new_beta[,1] <- tryLasso$a0
                new_beta[,-1] <- t(matrix(unlist(lapply(tryLasso$beta,function(x) return(x[,1]))),ncol=K))
              }
            }

          }else{
            fit_MultiLogit <- multinom(as.matrix(r)~as.matrix(G))
            new_beta[1,] <- 0
            new_beta[-1,] <- coef(fit_MultiLogit)
          }
        }else{
          new_pcluster <- colSums(r)/sum(r)
        }


        if(useY){
          #estimate new gamma
          if(family == "binary"){
            new_gamma <- apply(r,2,function(x) return(log(sum(x*Y)/(sum(x)-sum(x*Y)))))
          }
          if(family == "normal"){
            new_gamma[1:K] <- apply(r,2,function(x) return(sum(x*Y)/sum(x)))
            new_gamma[(K+1):(2*K)] <- colSums(r*apply(matrix(new_gamma[1:K]),1,function(x){(x-Y)^2}))/colSums(r)
          }
        }

        #stop criteria
        if(useY){
          if(!is.null(Z)&&is.null(G)){
            checkvalue <- all(is.finite(new_mu),is.finite(new_gamma),is.finite(unlist(new_sigma)),is.finite(new_pcluster))
            checksingular <- try(lapply(new_sigma,solve),silent=T)
            is_singular <- "try-error" %in% class(checksingular)
          }
          if(is.null(Z)&&!is.null(G)){
            checkvalue <- all(is.finite(new_beta),is.finite(new_gamma))
            is_singular <- FALSE
          }
          if(is.null(Z)&&is.null(G)){
            checkvalue <- all(is.finite(new_gamma),is.finite(new_pcluster))
            is_singular <- FALSE
          }
          if(!is.null(Z)&&!is.null(G)){
            checkvalue <- all(is.finite(new_beta),is.finite(new_mu),is.finite(new_gamma),is.finite(unlist(new_sigma)))
            checksingular <- try(lapply(new_sigma,solve),silent=T)
            is_singular <- "try-error" %in% class(checksingular)
          }

        }else{
          if(!is.null(Z)&&is.null(G)){
            checkvalue <- all(is.finite(new_mu),is.finite(unlist(new_sigma)),is.finite(new_pcluster))
            checksingular <- try(lapply(new_sigma,solve),silent=T)
            is_singular <- "try-error" %in% class(checksingular)
          }
          if(is.null(Z)&&!is.null(G)){
            checkvalue <- all(is.finite(new_beta))
            is_singular <- FALSE
          }
          if(!is.null(Z)&&!is.null(G)){
            checkvalue <- all(is.finite(new_beta),is.finite(new_mu),is.finite(unlist(new_sigma)))
            checksingular <- try(lapply(new_sigma,solve),silent=T)
            is_singular <- "try-error" %in% class(checksingular)
          }

        }

        if(!checkvalue || is_singular){
          breakdown <- TRUE
          print(paste(itr,": Invalid estimates"))
        }
        else{
          # Termination criteria
          if(!is.null(Z)){
            diffm <- matrix(mu-new_mu,nrow=1)
            dm <- sum(diffm^2)

            diffs <- unlist(sigma)-unlist(new_sigma)
            ds <- sum(diffs^2)
          }else{
            dm <- 0
            ds <- 0
          }

          if(useY){
            diffg <- matrix(gamma-new_gamma,nrow=1)
            dg <- sum(diffg^2)
          }else{
            dg <- 0
          }


          if(!is.null(G)){
            diffb <- matrix(beta-new_beta,nrow=1)
            db <- sum(diffb^2)
            dp <- 0
          }else{
            db <- 0
            diffp <- matrix(pcluster - new_pcluster,nrow=1)
            dp <- sum(diffp^2)
          }


          # Termination criteria

          if(dm<tol_m&&dg<tol_g&&db<tol_b&&ds<tol_s&&dp<tol_p){
            convergence <- 1
            mu <- new_mu
            sigma <- new_sigma
            gamma <- new_gamma
            beta <- new_beta
            pcluster <- new_pcluster
          }else{
            mu <- new_mu
            sigma <- new_sigma
            gamma <- new_gamma
            beta <- new_beta
            pcluster <- new_pcluster
          }

          print(paste(itr,": Updated parameters"))
        }
      }
    }

    if(convergence){
      success <- TRUE
      print(paste(itr,": Converged!"))
    }

  }

  if(!success){
    print("Fitting failed: exceed MAX_TOT_ITR...")
    err <- -99
    return (list(err = err))
  }
  else{
    if(useY){
      jointP <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=gamma, Family=family)

      if(!Pred){
        return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster,
                    Likelihood = jointP))
      }
      else{

        if(is.null(G)){
          preR <- t(t(jointP)*pcluster)
        }
        preR <- jointP/rowSums(jointP)
        return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster, pred = preR,
                    Likelihood = jointP))
      }

    }else{
      # Y not used in EM algorithm to estimate beta, mu and sigma

      estX <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=gamma, Family=family)

      if(is.null(G)){
        estX <- t(t(estX)*pcluster)
      }
      estX <- estX/rowSums(estX)

      #estimate gamma by regressing Y on estX
      if(family == "normal"){

        fit_y_model <- glm(Y~as.matrix(estX[,-1]),family=gaussian)
        gamma <- mat.or.vec(1,2*K)
        gamma[1:K] <- coef(fit_y_model)
        gamma_for_likelihood <- gamma
        gamma_for_likelihood[2:K] <- gamma[2:K]+gamma[1]
        gamma_for_likelihood[(K+1):(2*K)] <- (sd(fit_y_model$res))^2
        se_gamma <- summary(fit_y_model)$coef[,2]
      }
      if(family == "binary"){

        fit_y_model <- glm(Y~as.matrix(estX[,-1]),family="binomial")
        gamma <- mat.or.vec(1,K)
        gamma <- coef(fit_y_model)
        gamma_for_likelihood <- gamma
        gamma_for_likelihood[2:K] <- gamma[2:K]+gamma[1]
        se_gamma <- summary(fit_y_model)$coef[,2]
      }

      jointP <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=gamma_for_likelihood, Family=family)

      # gamma are deleted in the following section of estimating SE by SEM
      if(!Pred){
        return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster,
                    Likelihood = jointP))
      }
      else{

        if(is.null(G)){
          preR <- t(t(jointP)*pcluster)
        }
        preR <- jointP/rowSums(jointP)
        return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster, pred = preR,
                    Likelihood = jointP))
      }
    }

  }

}
