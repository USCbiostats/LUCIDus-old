#' Define initial values of parameters for clustering
#'
#' @description Defines initial values of model parameters in \code{\link{est_cluster}}, \code{\link{sem_cluster}} , & \code{\link{tune_cluster}} fitting.
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

  if (!is.null(init_b)){
    if (ncol(init_b) != ncol(G)+1 | nrow(init_b) != K ) {
      stop("Invalid initial beta; reset to NULL")
      init_b = NULL
    } else {
      init_b <- init_b
    }
  }

  if (!is.null(init_m)){
    if (ncol(init_m) != ncol(Z) | nrow(init_m) != K ) {
      stop("Invalid initial mu; reset to NULL")
      init_m = NULL
    } else {
      init_m <- init_m
    }
  }

  if (!is.null(init_s)){
    if (length(init_s) != K && !is.list(init_s)) {
      stop("Invalid initial sigma; reset to NULL")
      init_s = NULL
    } else {
      init_s <- init_s
    }
  }

  if (!is.null(init_g)){
    if (!is.vector(init_g)) {
      stop("Invalid initial gamma; reset to NULL")
      init_g = NULL
    } else {
      init_g <- init_g
    }
  }

  if (!is.null(init_pcluster)){
    if (length(init_pcluster) != K | !is.vector(init_g)) {
      stop("Invalid initial gamma; reset to NULL")
      init_pcluster = NULL
    } else {
      init_pcluster <- init_pcluster
    }
  }

  structure(list(init_b = init_b,
                 init_m = init_m,
                 init_s = init_s,
                 init_g = init_g,
                 init_pcluster = init_pcluster)
  )
}
