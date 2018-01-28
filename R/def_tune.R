#' Define tuning parameters for regularization during integrative clustering
#'
#' @description Defines selection options and tuning parameters in \code{\link{est_cluster}}, \code{\link{sem_cluster}} fitting.
#' @param Select_G Flag to do model selection on genetic data, default is FALSE
#' @param Select_Z Flag to do model selection on biomarker data, default is FALSE
#' @param Rho_G Penalty for selection on genetic data, numeric, default is -9 using a sequence of penalties
#' @param Rho_Z_InvCov Penalty for the inverse of covariance of biomarkers, numeric, default is 0
#' @param Rho_Z_CovMu Penalty for the product of covariance and mean of biomarkers, numeric, default is 0
#' @return A list of tolerance settings will be returned for integrative clustering.
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti

def_tol <- function(Rho_G = -9,
                    Rho_Z_InvCov = 0,
                    Rho_Z_CovMu = 0,
                    Select_G = FALSE,
                    Select_Z = FALSE) {

  if(!Select_G){
    Rho_G == 0
  }else{
    Rho_G == Rho_G
  }

  if(!Select_Z){
    Rho_Z_InvCov == 0
    Rho_Z_CovMu == 0
  }else{
    Rho_Z_InvCov = Rho_Z_InvCov
    Rho_Z_CovMu = Rho_Z_CovMu
  }

  structure(list(Select_G = Select_G,
                 Select_Z = Select_Z,
                 Rho_G = Rho_G,
                 Rho_Z_InvCov = Rho_Z_InvCov,
                 Rho_Z_CovMu = Rho_Z_CovMu)
  )
}
