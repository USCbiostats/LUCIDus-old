#' Define maximum number of iteration and convergence
#'
#' @description Defines tolerance settings in \code{\link{est_cluster}}, \code{\link{sem_cluster}} , & \code{\link{tune_cluster}} fitting.
#' @param MAX_ITR Maximum number of iterations, integer, default is 100
#' @param MAX_TOT_ITR Maximum number of total iterations, integer, default is 10000
#' @param reltol Convergence cut-off using a relative tolerance, default is 1e-8
#' @param tol_b Convergence criteria of \eqn{\beta}, genetic effects parameter, default is 1e-4
#' @param tol_m Convergence criteria of \eqn{\mu}, biomarker mean effects parameters, default is 1e-4
#' @param tol_s Convergence criteria of \eqn{\Sigma}, biomarker covariance matrix, default is 1e-4
#' @param tol_g Convergence criteria of \eqn{\gamma}, outcome effects parameter, default is 1e-4
#' @param tol_p Convergence criteria of the probability of latent clusters, default is 1e-4
#' @param tol_sem Convergence criteria of SEM, default is 1e-3
#' @return A list of tolerance settings will be returned for integrative clustering.
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti

def_tol <- function(MAX_ITR = 100,
                    MAX_TOT_ITR = 10000,
                    reltol = 1e-8,
                    tol_b = 1e-4,
                    tol_m = 1e-4,
                    tol_s = 1e-4,
                    tol_g = 1e-4,
                    tol_p = 1e-4,
                    tol_sem = 1e-3) {

  structure(list(MAX_ITR = MAX_ITR,
                 MAX_TOT_ITR = MAX_TOT_ITR,
                 reltol = reltol,
                 tol_b = tol_b,
                 tol_m = tol_m,
                 tol_s = tol_s,
                 tol_g = tol_g,
                 tol_p = tol_p,
                 tol_sem = tol_sem)
  )
}
