#' Model Predictions for LUCid
#'
#'\code{pred_lucid} produces predicted values for latent clusers and outcome with an \code{IntClust} object and new data.
#' @param Fit An \code{IntClust} class object
#' @param G Genetic effects, a matrix
#' @param CoG Covariates to be included in the G->X path
#' @param Z Biomarker data, a matrix
#' @param Y Disease outcome, a vector
#' @param CoY Covariates to be included in the X->Y path
#' @return \code{pred_lucid} returns a list containing predicted values.
#' \item{pred_cluster}{predicted probabilities for latent clusters}
#' \item{pred_outcome}{predicted values for outcome}
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Yang, Z., Conti, D.V.
#' @examples
#' set.seed(10)
#' IntClusFit <- est_lucid(G=G1,Z=Z1,Y=Y1,K=2,family="binary",Pred=TRUE)
#' PRED <- pred_lucid(Fit = IntClusFit, G=G1, CoG = NULL, Z=Z1, Y=Y1, CoY = NULL)

pred_lucid <- function(Fit=NULL, G=NULL, CoG=NULL, Z=NULL, Y, CoY=NULL){

  if(!is.null(CoG)){
    G <- cbind(G, CoG)
  }

  Beta <- Fit$beta
  Mu <- Fit$mu
  Sigma <- Fit$sigma
  Gamma <- Fit$gamma
  family <- Fit$family
  K <- Fit$K

  N <- nrow(G)

  likelihood <- function(Beta=NULL,Mu=NULL,Sigma=NULL,Gamma=NULL,Family,total_ITR=2){

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

      if(total_ITR==1){
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
      }
      else{
        if(Family == "binary"){
          for(k in 1:K){
            if(is.null(CoY)){
              SetK <- as.data.frame(mat.or.vec(N,K-1))
              colnames(SetK) <- names(coef(Fit$YFIT))[-1]
            }
            else{
              SetK <- as.data.frame(cbind(mat.or.vec(N,K-1), CoY))
              colnames(SetK) <- names(coef(Fit$YFIT))[-1]
            }
            if(k>1){
              SetK[,k-1] <- 1
            }
            pYgA[,k] <- dbinom(Y, 1, predict(Fit$YFIT, newdata = SetK, type = "response"))
          }
        }

        if(Family == "normal"){
          for(k in 1:K){
            pYgA[,k] <- dnorm(Y,mean = Gamma[k],sd = sigma(Fit$YFIT))
          }
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

  PYgX <- function(X=NULL, Gamma=NULL){
    lincom <- as.matrix(cbind(rep(1,N),X[,-1]))%*%as.matrix(Gamma)
    # lincom <- as.matrix(X)%*%as.matrix(Gamma)
    pi <- exp(lincom)/(1+exp(lincom))
  }

  #Posterior Probabilities
  r <- likelihood(Beta=Beta,Mu=Mu,Sigma=Sigma,Gamma=Gamma,Family=family)
  pred_cluster <- r/rowSums(r)

  #Predicted Outcome
  if(family == "binary"){
    SetF <- as.data.frame(cbind(Y, pred_cluster[,-1], CoY))
    pred_outcome <- predict(Fit$YFIT, newdata = SetF, type = "response")
  }

  if(family == "normal"){
    SetF <- as.data.frame(cbind(Y, pred_cluster, CoY))
    pred_outcome <- predict(Fit$YFIT, newdata = SetF)
  }

  Preds <- list(pred_cluster,pred_outcome)
  names(Preds) <- c("pred_cluster", "pred_outcome")
  return(Preds)
}