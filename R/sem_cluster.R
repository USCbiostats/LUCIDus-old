#' SEM for latent cluster estimation
#'
#' \code{sem_cluster} provides standard errors (SE) of parameter estimates when performing latent cluster analysis with multi-omics data. SEs are obtained through supplemented EM-algorithm (SEM).
#' @param G Genetic effects, a matrix
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param family "binary" or "normal" for Y
#' @param useY Using Y or not, default is TRUE
#' @param K Pre-specified # of latent clusters, default is 2
#' @param Get_SE Flag to perform SEM to get SEs of parameter estimates, default is TRUE
#' @param Ad_Hoc_SE Flag to fit ad hoc regression models to get SEs of parameter estimates, default is FALSE
#' @param Pred Flag to compute predicted disease probability with fitted model, boolean, default is FALSE
#' @param initial A list of initial model parameters will be returned for integrative clustering
#' @param itr_tol A list of tolerance settings will be returned for integrative clustering
#' @keywords SEM, latent cluster
#' @return \code{sem_cluster} returns an object of list containing parameters estimates, their corresponding standard errors, and other features:
#' \item{beta}{Estimates of genetic effects, matrix}
#' \item{se_beta}{SEM standard errors of Beta}
#' \item{se_ah_beta}{Ad hoc standard errors of Beta}
#' \item{mu}{Estimates of cluster-specific biomarker means, matrix}
#' \item{se_mu}{SEM standard errors of Mu}
#' \item{se_ah_mu}{Ad hoc standard errors of Mu}
#' \item{sigma}{Estimates of cluster-specific biomarker covariance matrix, list}
#' \item{gamma}{Estimates of cluster-specific disease risk, vector}
#' \item{se_gamma}{SEM standard errors of Gamma}
#' \item{se_ah_gamma}{Ad hoc standard errors of Gamma}
#' \item{pcluster}{Probability of cluster, when G is null}
#' \item{pred}{Predicted probability of belonging to each latent cluster}
#' @importFrom mvtnorm dmvnorm
#' @importFrom nnet multinom
#' @importFrom glmnet glmnet
#' @importFrom glasso glasso
#' @importFrom lbfgs lbfgs
#' @import Matrix
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Meng, X., & Rubin, D. B. (1991). Using EM to Obtain Asymptotic Matrices : The SEM Algorithm.
#' Journal of the American Statistical Association, 86(416), 899-909. http://doi.org/10.2307/2290503
#' @examples
#' sem_cluster(G=G2,Z=Z2,Y=Y2,useY=TRUE,K=2,Pred=TRUE,family="normal",Get_SE=TRUE,
#'             def_initial(),def_tol(MAX_ITR=1000,MAX_TOT_ITR=3000))


sem_cluster <- function(G=NULL, Z=NULL, Y, family="binary", useY = TRUE, K = 2,
                        initial = def_initial(), itr_tol = def_tol(),
                        Pred = FALSE, Get_SE = TRUE, Ad_Hoc_SE = FALSE){

  init_b <- initial$init_b; init_m <- initial$init_m; init_s <- initial$init_s; init_g <- initial$init_g; init_pcluster <- initial$init_pcluster
  MAX_ITR <- itr_tol$MAX_ITR; MAX_TOT_ITR <- itr_tol$MAX_TOT_ITR; reltol <- itr_tol$reltol
  tol_b <- itr_tol$tol_b; tol_m <- itr_tol$tol_m; tol_s <- itr_tol$tol_s; tol_g <- itr_tol$tol_g; tol_p <- itr_tol$tol_p; tol_sem <- itr_tol$tol_sem

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

  #define likelihood function
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

    if(Get_SE == TRUE){
      # create lists to store parameter values at each iteration

      if(family=="binary"){
        theta_itr <- mat.or.vec(MAX_ITR, (K-1)*(M+1) + useY*K + Q*K + Q*(Q+1)*0.5*K + P*K)
      }else{
        theta_itr <- mat.or.vec(MAX_ITR, (K-1)*(M+1) + useY*2*K + Q*K + Q*(Q+1)*0.5*K + P*K)
      }
    }

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

        if(!is.null(Z)){
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

        if(!is.null(G)){
          fit_MultiLogit <- multinom(as.matrix(r)~as.matrix(G))
          new_beta[1,] <- 0
          new_beta[-1,] <- coef(fit_MultiLogit)
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

          if(Get_SE){

            # store parameter values

            if(!is.null(G)){
              itr_beta <- array(t(beta[-1,]))
            }else{
              itr_pcluster <- pcluster
            }
            if(!is.null(Z)){
              itr_mu <- array(t(mu))
              itr_sigma <- sigma[[1]][lower.tri(sigma[[1]],diag=T)]

              if(K > 1){
                k = 2
                while(k <= K){
                  itr_sigma <- c(itr_sigma, sigma[[k]][lower.tri(sigma[[k]],diag=T)])
                  k <- k + 1
                }
              }
            }

            if(useY){
              itr_gamma <- gamma
              if(!is.null(Z)&&is.null(G)){
                theta_itr[itr,] <- c(itr_gamma, itr_mu, itr_sigma, itr_pcluster)
              }
              if(is.null(Z)&&!is.null(G)){
                theta_itr[itr,] <- c(itr_beta, itr_gamma)
              }
              if(is.null(Z)&&is.null(G)){
                theta_itr[itr,] <- c(itr_gamma,itr_pcluster)
              }
              if(!is.null(Z)&&!is.null(G)){
                theta_itr[itr,] <- c(itr_beta, itr_gamma, itr_mu, itr_sigma)
              }
            }else{
              if(!is.null(Z)&&is.null(G)){
                theta_itr[itr,] <- c(itr_mu, itr_sigma, itr_pcluster)
              }
              if(is.null(Z)&&!is.null(G)){
                theta_itr[itr,] <- c(itr_beta)
              }
              if(is.null(Z)&&is.null(G)){
                theta_itr[itr,] <- c(itr_pcluster)
              }
              if(!is.null(Z)&&!is.null(G)){
                theta_itr[itr,] <- c(itr_beta, itr_mu, itr_sigma)
              }
            }

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

      if(Get_SE || Ad_Hoc_SE){

        se_beta <- NULL
        se_gamma <- NULL
        se_mu <- NULL

        se_msg <- NULL
        se_ah_beta <- NULL
        se_ah_gamma <- NULL
        se_ah_mu <- NULL

        ## compute the inverse of Conditional Expectation of Complete data Observed Information Matrix, (I_oc)^-1
        pred_r <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=gamma, Family=family)

        if(is.null(G)){
          pred_r <- t(t(pred_r)*pcluster)
        }
        pred_r <- pred_r/rowSums(pred_r)
        # beta
        if(!is.null(G)){
          fit_beta <- multinom(as.matrix(pred_r)~as.matrix(G), Hess = TRUE)
          # to begin with, this version only works for 2 clusters
          Hess_beta <- fit_beta$Hess
          Var_beta <- solve(Hess_beta)
        }else{
          #use zeros for pcluster
          Var_pcluster <- diag(0,K)
        }

        # gamma
        if(family=="binary"){
          cross_product_gamma <- unlist(lapply(gamma, function(x){exp(x)/((1+exp(x))^2)}))
          Var_gamma <- 1/rowSums(apply(pred_r,1,function(x) return(x*cross_product_gamma)))
        }
        if(family=="normal"){
          Var_gamma <- rep(0,2*K)
          Var_gamma[1:K] <- gamma[(K+1):(2*K)]/colSums(pred_r)
        }

        if(!is.null(Z)){
          # mu
          Var_mu <- replicate(K,mat.or.vec(Q,Q),simplify=F)
          weight <- colSums(pred_r)
          for(k in 1:K){
            Var_mu[[k]] <- sigma[[k]]/weight[k]
          }

          # use zeros for sigma
          Var_sigma <- diag(0, nrow = Q*(Q+1)/2*K, ncol=Q*(Q+1)/2*K)
        }

        # write them in a single matrix
        if(!is.null(Z)&&is.null(G)){
          inv_I_oc <- bdiag(diag(Var_gamma),bdiag(Var_mu),Var_sigma,Var_pcluster)
        }
        if(is.null(Z)&&!is.null(G)){
          inv_I_oc <- bdiag(Var_beta, diag(Var_gamma))
        }
        if(is.null(Z)&&is.null(G)){
          inv_I_oc <- bdiag(diag(Var_gamma),Var_pcluster)
        }
        if(!is.null(Z)&&!is.null(G)){
          inv_I_oc <- bdiag(Var_beta, diag(Var_gamma),bdiag(Var_mu),Var_sigma)
        }


        # ad hoc variance
        if(Ad_Hoc_SE){
          V_ad_hoc <- inv_I_oc

          se_ah <- sqrt(diag(V_ad_hoc))

          if(!is.null(G)){
            se_ah_beta <- se_ah[1:((K-1)*(M+1))]
          }

          if(family=="binary"){
            se_ah_gamma <- se_ah[((K-1)*(M+1)+1):((K-1)*(M+1)+K)]
            if(!is.null(Z)){
              se_ah_mu <- se_ah[((K-1)*(M+1)+K+1):((K-1)*(M+1)+K+K*Q)]
              se_ah_mu <- t(matrix(se_ah_mu, ncol = K))
            }
          }

          if(family=="normal"){
            se_ah_gamma <- se_ah[((K-1)*(M+1)+1):((K-1)*(M+1)+K)]
            if(!is.null(Z)){
              se_ah_mu <- se_ah[((K-1)*(M+1)+2*K+1):((K-1)*(M+1)+2*K+K*Q)]
              se_ah_mu <- t(matrix(se_ah_mu, ncol = K))
            }
          }

        }

        if(Get_SE){
          ## Compute the DM matrix
          # store parameters estimated from EM algorithm in theta_star
          if(!is.null(G)){
            array_beta <- beta[-1,]
          }else{
            array_pcluster <- pcluster
          }

          array_gamma <- gamma

          if(!is.null(Z)){
            array_mu <- array(t(mu))
            array_sigma <- sigma[[1]][lower.tri(sigma[[1]],diag=T)]
            if(K > 1){
              k = 2
              while(k <= K){
                array_sigma <- c(array_sigma, sigma[[k]][lower.tri(sigma[[k]],diag=T)])
                k <- k + 1
              }
            }
          }



          if(!is.null(Z)&&is.null(G)){
            theta_star <- c(array_gamma, array_mu, array_sigma, array_pcluster)
          }
          if(is.null(Z)&&!is.null(G)){
            theta_star <- c(array_beta, array_gamma)
          }
          if(is.null(Z)&&is.null(G)){
            theta_star <- c(array_gamma, array_pcluster)
          }
          if(!is.null(Z)&&!is.null(G)){
            theta_star <- c(array_beta, array_gamma, array_mu, array_sigma)
          }

          # use saved parameters from EM iterations

          itr_sem <- 1

          n_parm <- length(theta_star) # total number of parameters
          sem_converged <- matrix(FALSE, ncol=n_parm, nrow=n_parm) # initialize indicator
          sem_success <- 0 # initialize indicator
          # initialize
          sem_mu <- NULL
          sem_sigma <- NULL
          sem_gamma <- NULL
          sem_beta <- NULL

          new_sem_mu <- NULL
          new_sem_sigma <- NULL
          new_sem_gamma <- NULL
          new_sem_beta <- NULL


          # SEM to compute the DM matrix
          if(!is.null(Z)){
            sem_mu <- mat.or.vec(K,Q) #store estimated mu SEM
            sem_sigma <- replicate(K,diag(0,nrow=Q,ncol=Q),simplify=F) #store estiamted sigma SEM
            new_sem_mu <- mat.or.vec(K,Q) #store estimated mu SEM
            new_sem_sigma <- replicate(K,diag(0,nrow=Q,ncol=Q),simplify=F) #store estiamted sigma SEM
          }

          if(!is.null(G)){
            sem_beta <- mat.or.vec(K,M+1) #store estimated beta in SEM
            new_sem_beta <- mat.or.vec(K,M+1) #store estimated beta in SEM
          }else{
            sem_pcluster <- mat.or.vec(1,K)
            new_sem_pcluster <- mat.or.vec(1,K)
          }


          if(family=="binary"){
            sem_gamma <- array(0,K) #store estimated gamma in SEM
            new_sem_gamma <- array(0,K) #store estimated gamma in SEM
          }
          if(family=="normal"){
            sem_gamma <- array(0,2*K) #store estimated gamma in SEM
            new_sem_gamma <- array(0,2*K) #store estimated gamma in SEM
          }

          DM <- mat.or.vec(n_parm, n_parm)
          new_DM <- mat.or.vec(n_parm, n_parm)

          while(itr_sem <= itr && any(sem_converged == FALSE)){
            for(i in 1:n_parm){
              if(any(sem_converged[i,] == FALSE)){

                curr_theta <- theta_star
                curr_theta[i] <- theta_itr[itr_sem,i]

                if(!is.null(G)){
                  sem_beta[-1,] <- t(array(curr_theta[1:((M+1)*(K-1))],dim=c((M+1),(K-1))))
                }else{
                  sem_pcluster <- curr_theta[(length(curr_theta)-K+1) : length(curr_theta)]
                }

                if(family=="binary"){
                  sem_gamma <- curr_theta[((M+1)*(K-1)+1):((M+1)*(K-1)+K)]
                  if(!is.null(Z)){
                    sem_mu <- t(matrix(curr_theta[((M+1)*(K-1)+K+1):((M+1)*(K-1)+K+Q*K)],ncol=K))
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+K+Q*K+1):(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2)]
                    sem_sigma[[1]] <- t(sem_sigma[[1]])
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+K+Q*K+1):(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2)]

                    if(K>1){
                      k=2
                      while(k <= K){
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2*k)]
                        sem_sigma[[k]] <- t(sem_sigma[[k]])
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+K+Q*K)+Q*(Q+1)/2*k)]
                        k <- k + 1
                      }
                    }
                  }
                }
                if(family=="normal"){
                  sem_gamma <- curr_theta[((M+1)*(K-1)+1):((M+1)*(K-1)+2*K)]

                  if(!is.null(Z)){
                    sem_mu <- t(matrix(curr_theta[((M+1)*(K-1)+2*K+1):((M+1)*(K-1)+2*K+Q*K)],ncol=K))
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+2*K+Q*K+1):(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2)]
                    sem_sigma[[1]] <- t(sem_sigma[[1]])
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+2*K+Q*K+1):(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2)]

                    if(K>1){
                      k=2
                      while(k <= K){
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2*k)]
                        sem_sigma[[k]] <- t(sem_sigma[[k]])
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+2*K+Q*K)+Q*(Q+1)/2*k)]
                        k <- k + 1
                      }
                    }
                  }

                }




                #browser()
                #------E-step------#
                sem_r <- likelihood(Beta=sem_beta, Mu=sem_mu, Sigma=sem_sigma, Gamma=sem_gamma, Family=family)
                if(is.null(G)){
                  sem_r <- t(t(sem_r)*sem_pcluster)
                }
                sem_r <- sem_r/rowSums(sem_r)

                #------check r------#
                valid_sem_r <- all(is.finite(as.matrix(sem_r)))

                if(!valid_sem_r){
                  sem_breakdown <- TRUE
                  print(paste(itr,": invalid sem_r"))
                }
                else{
                  print(paste(itr,": E-step finished for SEM"))
                }


                #------M-STEP------#
                # mu and sigma

                if(!is.null(Z)){
                  new_sem_mu <- matrix(t(apply(sem_r,2,function(x) return(colSums(x*Z)/sum(x)))),ncol=ncol(Z))
                  if(ncol(Z)>1){
                    for(k in 1:K){
                      Z_sem_mu <- t(t(Z)-new_sem_mu[k,])
                      new_sem_sigma[[k]] <- (matrix(colSums(sem_r[,k]*t(apply(Z_sem_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(sem_r[,k])
                    }
                  }else{
                    for(k in 1:K){
                      Z_sem_mu <- t(t(Z)-new_sem_mu[k,])
                      new_sem_sigma[[k]] <- (matrix(sum(sem_r[,k]*t(apply(Z_sem_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(sem_r[,k])
                    }
                  }
                }





                # beta
                if(!is.null(G)){
                  fit_sem_MultiLogit <- multinom(as.matrix(sem_r)~as.matrix(G))
                  new_sem_beta[1,] <- 0
                  new_sem_beta[-1,] <- coef(fit_sem_MultiLogit)
                }else{
                  new_sem_pcluster <- colSums(sem_r)/sum(sem_r)
                }


                # gamma
                if(family=="binary"){
                  new_sem_gamma <- apply(sem_r,2,function(x) return(log(sum(x*Y)/(sum(x)-sum(x*Y)))))
                }
                if(family=="normal"){
                  new_sem_gamma[1:K] <- apply(sem_r,2,function(x) return(sum(x*Y)/sum(x)))
                  new_sem_gamma[(K+1):(2*K)] <- colSums(sem_r*apply(matrix(new_sem_gamma[1:K]),1,function(x){(x-Y)^2}))/colSums(sem_r)
                }

                # write updated parameter into a vector
                if(!is.null(G)){
                  #sem_itr_beta <- new_sem_beta[2,]
                  sem_itr_beta <- array(t(new_sem_beta[-1,]))
                }else{
                  sem_itr_pcluster <- new_sem_pcluster
                }

                sem_itr_gamma <- new_sem_gamma

                if(!is.null(Z)){
                  sem_itr_mu <- array(t(new_sem_mu))
                  sem_itr_sigma <- new_sem_sigma[[1]][lower.tri(new_sem_sigma[[1]],diag=T)]
                  if(K > 1){
                    k = 2
                    while(k <= K){
                      sem_itr_sigma <- c(sem_itr_sigma, new_sem_sigma[[k]][lower.tri(new_sem_sigma[[k]],diag=T)])
                      k <- k + 1
                    }
                  }
                }

                if(!is.null(Z)&&is.null(G)){
                  curr_theta_tilta <- c(sem_itr_gamma, sem_itr_mu, sem_itr_sigma, sem_itr_pcluster)
                }
                if(is.null(Z)&&!is.null(G)){
                  curr_theta_tilta <- c(sem_itr_beta, sem_itr_gamma)
                }
                if(is.null(Z)&&is.null(G)){
                  curr_theta_tilta <- c(sem_itr_gamma, sem_itr_pcluster)
                }
                if(!is.null(Z)&&!is.null(G)){
                  curr_theta_tilta <- c(sem_itr_beta, sem_itr_gamma, sem_itr_mu, sem_itr_sigma)
                }




                # compute DM_ij
                denom <- theta_itr[itr_sem,i]-theta_star[i]
                num <- curr_theta_tilta - theta_star

                new_DM[i,] <- num/denom
                #browser()
              }
            }

            d_DM <- new_DM - DM
            sem_converged <- ifelse(d_DM^2<tol_sem, TRUE, FALSE)
            DM <- new_DM
            itr_sem <- itr_sem + 1
          }


          ## Compute Observed Variance-Covariance matrix, V = (I_oc)^-1 + (I_oc)^-1 %*% DM %*% (I-DM)^-1

          V = inv_I_oc + inv_I_oc %*% DM %*% solve((diag(n_parm)-DM))

          se <- sqrt(diag(V))

          if(all(sem_converged == TRUE) && sum(is.na(se))==0){
            sem_success = 1
            print("SEM succeed!")
          }
          else{
            sem_success = 0
            print("SEM failed!")
          }

          if(!is.null(G)){
            se_beta <- se[1:((K-1)*(M+1))]
          }

          se_gamma <- se[((K-1)*(M+1)+1):((K-1)*(M+1)+K)]

          if(!is.null(Z)){
            if(family=="binary"){
              se_mu <- se[((K-1)*(M+1)+K+1):((K-1)*(M+1)+K+K*Q)]
            }
            if(family=="normal"){
              se_mu <- se[((K-1)*(M+1)+2*K+1):((K-1)*(M+1)+2*K+K*Q)]
            }

            se_mu <- t(matrix(se_mu, ncol = K))
          }

        }


        if(!Pred){
          return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster,
                      se_beta = se_beta, se_gamma = se_gamma, se_mu = se_mu, se_msg = sem_success,
                      se_ah_beta = se_ah_beta, se_ah_gamma = se_ah_gamma, se_ah_mu = se_ah_mu,
                      Likelihood = jointP))
        }
        else{

          if(is.null(G)){
            preR <- t(t(jointP)*pcluster)
          }
          preR <- jointP/rowSums(jointP)
          return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma, pcluster = pcluster,
                      se_beta = se_beta, se_gamma = se_gamma, se_mu = se_mu, se_msg = sem_success,
                      se_ah_beta = se_ah_beta, se_ah_gamma = se_ah_gamma, se_ah_mu = se_ah_mu,
                      pred = preR, Likelihood = jointP))
        }
      }

    }else{
      # Y not used in EM algorithm to estimate beta, mu and sigma

      estX <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=NULL, Family=family)

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
      if(Get_SE || Ad_Hoc_SE){

        se_beta <- NULL
        se_mu <- NULL

        se_msg <- NULL
        se_ah_beta <- NULL
        se_ah_mu <- NULL

        ## compute the inverse of Conditional Expectation of Complete data Observed Information Matrix, (I_oc)^-1
        pred_r <- likelihood(Beta=beta, Mu=mu, Sigma=sigma, Gamma=NULL, Family=family)

        if(is.null(G)){
          pred_r <- t(t(pred_r)*pcluster)
        }
        pred_r <- pred_r/rowSums(pred_r)
        # beta
        if(!is.null(G)){
          fit_beta <- multinom(as.matrix(pred_r)~as.matrix(G), Hess = TRUE)
          # to begin with, this version only works for 2 clusters
          Hess_beta <- fit_beta$Hess
          Var_beta <- solve(Hess_beta)
        }else{
          #use zeros for pcluster
          Var_pcluster <- diag(0,K)
        }

        if(!is.null(Z)){
          # mu
          Var_mu <- replicate(K,mat.or.vec(Q,Q),simplify=F)
          weight <- colSums(pred_r)
          for(k in 1:K){
            Var_mu[[k]] <- sigma[[k]]/weight[k]
          }

          # use zeros for sigma
          Var_sigma <- diag(0, nrow = Q*(Q+1)/2*K, ncol=Q*(Q+1)/2*K)
        }

        # write them in a single matrix
        if(!is.null(Z)&&is.null(G)){
          inv_I_oc <- bdiag(bdiag(Var_mu),Var_sigma,Var_pcluster)
        }
        if(is.null(Z)&&!is.null(G)){
          inv_I_oc <- bdiag(Var_beta)
        }
        if(is.null(Z)&&is.null(G)){
          inv_I_oc <- bdiag(Var_pcluster)
        }
        if(!is.null(Z)&&!is.null(G)){
          inv_I_oc <- bdiag(Var_beta,bdiag(Var_mu),Var_sigma)
        }


        # ad hoc variance
        if(Ad_Hoc_SE){
          V_ad_hoc <- inv_I_oc

          se_ah <- sqrt(diag(V_ad_hoc))

          if(!is.null(G)){
            se_ah_beta <- se_ah[1:((K-1)*(M+1))]
          }

          if(family=="binary"){

            if(!is.null(Z)){
              se_ah_mu <- se_ah[((K-1)*(M+1)+1):((K-1)*(M+1)+K*Q)]
              se_ah_mu <- t(matrix(se_ah_mu, ncol = K))
            }
          }

          if(family=="normal"){

            if(!is.null(Z)){
              se_ah_mu <- se_ah[((K-1)*(M+1)+1):((K-1)*(M+1)+K*Q)]
              se_ah_mu <- t(matrix(se_ah_mu, ncol = K))
            }
          }

        }

        if(Get_SE){
          ## Compute the DM matrix
          # store parameters estimated from EM algorithm in theta_star
          if(!is.null(G)){
            array_beta <- beta[-1,]
          }else{
            array_pcluster <- pcluster
          }

          if(!is.null(Z)){
            array_mu <- array(t(mu))
            array_sigma <- sigma[[1]][lower.tri(sigma[[1]],diag=T)]
            if(K > 1){
              k = 2
              while(k <= K){
                array_sigma <- c(array_sigma, sigma[[k]][lower.tri(sigma[[k]],diag=T)])
                k <- k + 1
              }
            }
          }


          if(!is.null(Z)&&is.null(G)){
            theta_star <- c(array_mu, array_sigma, array_pcluster)
          }
          if(is.null(Z)&&!is.null(G)){
            theta_star <- c(array_beta)
          }
          if(is.null(Z)&&is.null(G)){
            theta_star <- c(array_pcluster)
          }
          if(!is.null(Z)&&!is.null(G)){
            theta_star <- c(array_beta,array_mu, array_sigma)
          }

          # use saved parameters from EM iterations

          itr_sem <- 1

          n_parm <- length(theta_star) # total number of parameters
          sem_converged <- matrix(FALSE, ncol=n_parm, nrow=n_parm) # initialize indicator
          sem_success <- 0 # initialize indicator
          # initialize
          sem_mu <- NULL
          sem_sigma <- NULL
          sem_beta <- NULL

          new_sem_mu <- NULL
          new_sem_sigma <- NULL
          new_sem_beta <- NULL


          # SEM to compute the DM matrix
          if(!is.null(Z)){
            sem_mu <- mat.or.vec(K,Q) #store estimated mu SEM
            sem_sigma <- replicate(K,diag(0,nrow=Q,ncol=Q),simplify=F) #store estiamted sigma SEM
            new_sem_mu <- mat.or.vec(K,Q) #store estimated mu SEM
            new_sem_sigma <- replicate(K,diag(0,nrow=Q,ncol=Q),simplify=F) #store estiamted sigma SEM
          }

          if(!is.null(G)){
            sem_beta <- mat.or.vec(K,M+1) #store estimated beta in SEM
            new_sem_beta <- mat.or.vec(K,M+1) #store estimated beta in SEM
          }else{
            sem_pcluster <- mat.or.vec(1,K)
            new_sem_pcluster <- mat.or.vec(1,K)
          }

          DM <- mat.or.vec(n_parm, n_parm)
          new_DM <- mat.or.vec(n_parm, n_parm)

          while(itr_sem <= itr && any(sem_converged == FALSE)){
            for(i in 1:n_parm){
              if(any(sem_converged[i,] == FALSE)){

                curr_theta <- theta_star
                curr_theta[i] <- theta_itr[itr_sem,i]

                if(!is.null(G)){
                  sem_beta[-1,] <- t(array(curr_theta[1:((M+1)*(K-1))],dim=c((M+1),(K-1))))
                }else{
                  sem_pcluster <- curr_theta[(length(curr_theta)-K+1) : length(curr_theta)]
                }

                if(family=="binary"){
                  if(!is.null(Z)){
                    sem_mu <- t(matrix(curr_theta[((M+1)*(K-1)+1):((M+1)*(K-1)+Q*K)],ncol=K))
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+Q*K+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2)]
                    sem_sigma[[1]] <- t(sem_sigma[[1]])
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+Q*K+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2)]

                    if(K>1){
                      k=2
                      while(k <= K){
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*k)]
                        sem_sigma[[k]] <- t(sem_sigma[[k]])
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*k)]
                        k <- k + 1
                      }
                    }
                  }
                }
                if(family=="normal"){

                  if(!is.null(Z)){
                    sem_mu <- t(matrix(curr_theta[((M+1)*(K-1)+1):((M+1)*(K-1)+Q*K)],ncol=K))
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+Q*K+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2)]
                    sem_sigma[[1]] <- t(sem_sigma[[1]])
                    sem_sigma[[1]][lower.tri(sem_sigma[[1]],diag=T)] <- curr_theta[((M+1)*(K-1)+Q*K+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2)]

                    if(K>1){
                      k=2
                      while(k <= K){
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*k)]
                        sem_sigma[[k]] <- t(sem_sigma[[k]])
                        sem_sigma[[k]][lower.tri(sem_sigma[[k]],diag=T)] <- curr_theta[(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*(k-1)+1):(((M+1)*(K-1)+Q*K)+Q*(Q+1)/2*k)]
                        k <- k + 1
                      }
                    }
                  }

                }



                #------E-step------#
                sem_r <- likelihood(Beta=sem_beta, Mu=sem_mu, Sigma=sem_sigma, Gamma=NULL, Family=family)
                if(is.null(G)){
                  sem_r <- t(t(sem_r)*sem_pcluster)
                }
                sem_r <- sem_r/rowSums(sem_r)

                #------check r------#
                valid_sem_r <- all(is.finite(as.matrix(sem_r)))

                if(!valid_sem_r){
                  sem_breakdown <- TRUE
                  print(paste(itr,": invalid sem_r"))
                }
                else{
                  print(paste(itr,": E-step finished for SEM"))
                }


                #------M-STEP------#
                # mu and sigma

                if(!is.null(Z)){
                  new_sem_mu <- matrix(t(apply(sem_r,2,function(x) return(colSums(x*Z)/sum(x)))),ncol=ncol(Z))
                  if(ncol(Z)>1){
                    for(k in 1:K){
                      Z_sem_mu <- t(t(Z)-new_sem_mu[k,])
                      new_sem_sigma[[k]] <- (matrix(colSums(sem_r[,k]*t(apply(Z_sem_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(sem_r[,k])
                    }
                  }else{
                    for(k in 1:K){
                      Z_sem_mu <- t(t(Z)-new_sem_mu[k,])
                      new_sem_sigma[[k]] <- (matrix(sum(sem_r[,k]*t(apply(Z_sem_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(sem_r[,k])
                    }
                  }
                }





                # beta
                if(!is.null(G)){
                  fit_sem_MultiLogit <- multinom(as.matrix(sem_r)~as.matrix(G))
                  new_sem_beta[1,] <- 0
                  new_sem_beta[-1,] <- coef(fit_sem_MultiLogit)
                }else{
                  new_sem_pcluster <- colSums(sem_r)/sum(sem_r)
                }

                # write updated parameter into a vector
                if(!is.null(G)){
                  sem_itr_beta <- array(t(new_sem_beta[-1,]))
                }else{
                  sem_itr_pcluster <- new_sem_pcluster
                }

                if(!is.null(Z)){
                  sem_itr_mu <- array(t(new_sem_mu))
                  sem_itr_sigma <- new_sem_sigma[[1]][lower.tri(new_sem_sigma[[1]],diag=T)]
                  if(K > 1){
                    k = 2
                    while(k <= K){
                      sem_itr_sigma <- c(sem_itr_sigma, new_sem_sigma[[k]][lower.tri(new_sem_sigma[[k]],diag=T)])
                      k <- k + 1
                    }
                  }
                }

                if(!is.null(Z)&&is.null(G)){
                  curr_theta_tilta <- c(sem_itr_mu, sem_itr_sigma, sem_itr_pcluster)
                }
                if(is.null(Z)&&!is.null(G)){
                  curr_theta_tilta <- c(sem_itr_beta)
                }
                if(is.null(Z)&&is.null(G)){
                  curr_theta_tilta <- c(sem_itr_pcluster)
                }
                if(!is.null(Z)&&!is.null(G)){
                  curr_theta_tilta <- c(sem_itr_beta, sem_itr_mu, sem_itr_sigma)
                }




                # compute DM_ij
                denom <- theta_itr[itr_sem,i]-theta_star[i]
                num <- curr_theta_tilta - theta_star

                new_DM[i,] <- num/denom

              }
            }


            d_DM <- new_DM - DM
            sem_converged <- ifelse(d_DM^2<tol_sem, TRUE, FALSE)
            DM <- new_DM
            itr_sem <- itr_sem + 1
          }

          ## Compute Observed Variance-Covariance matrix, V = (I_oc)^-1 + (I_oc)^-1 %*% DM %*% (I-DM)^-1

          V = inv_I_oc + inv_I_oc %*% DM %*% solve((diag(n_parm)-DM))

          se <- sqrt(diag(V))

          if(all(sem_converged == TRUE) && sum(is.na(se))==0){
            sem_success = 1
            print("SEM succeed!")
          }
          else{
            sem_success = 0
            print("SEM failed!")
          }

          if(!is.null(G)){
            se_beta <- se[1:((K-1)*(M+1))]
          }

          if(!is.null(Z)){
            if(family=="binary"){
              se_mu <- se[((K-1)*(M+1)+1):((K-1)*(M+1)+K*Q)]
            }
            if(family=="normal"){
              se_mu <- se[((K-1)*(M+1)+1):((K-1)*(M+1)+K*Q)]
            }

            se_mu <- t(matrix(se_mu, ncol = K))
          }

        }


        if(!Pred){
          return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma_for_likelihood, pcluster = pcluster,
                      se_beta = se_beta, se_gamma = se_gamma, se_mu = se_mu, se_msg = sem_success,
                      se_ah_beta = se_ah_beta, se_ah_gamma = NULL, se_ah_mu = se_ah_mu,
                      Likelihood = jointP))
        }
        else{

          if(is.null(G)){
            preR <- t(t(jointP)*pcluster)
          }
          preR <- jointP/rowSums(jointP)
          return(list(beta = beta, mu = mu, sigma = sigma, gamma = gamma_for_likelihood, pcluster = pcluster,
                      se_beta = se_beta, se_gamma = se_gamma, se_mu = se_mu, se_msg = sem_success,
                      se_ah_beta = se_ah_beta, se_ah_gamma = NULL, se_ah_mu = se_ah_mu,
                      pred = preR, Likelihood = jointP))
        }
      }

    }

  }

}

