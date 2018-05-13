#' Define initial values of parameters for clustering
#'
#' @description Defines initial values of model parameters in \code{\link{est_lucid}}, \code{\link{sem_lucid}} , & \code{\link{tune_lucid}} fitting.
#'
#' @param init_b Initial model parameters of \eqn{\beta}, genetic effects parameter
#' @param init_m Initial model parameters of \eqn{\mu}, biomarker mean effects parameters
#' @param init_s Initial model parameters of \eqn{\Sigma}, biomarker covariance matrix
#' @param init_g Initial model parameters of \eqn{\gamma}, outcome effects parameter
#' @param init_pcluster Initial probabilites of latent clusters
#' @return A list of initial model parameters will be returned for integrative clustering.
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti

def_initial <- function(init_b = NULL,
                        init_m = NULL,
                        init_s = NULL,
                        init_g = NULL,
                        init_pcluster = NULL) {

  structure(list(init_b = init_b,
                 init_m = init_m,
                 init_s = init_s,
                 init_g = init_g,
                 init_pcluster = init_pcluster)
  )
}
