#' Genetic Features Set 1
#'
#' A simulated dataset containing one of the components to run est_lucid,
#' plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 10 variables:
#' \describe{
#'   \item{CG1 - CG5}{Causal SNPs}
#'   \item{NG1 - NG5}{Null SNPs}
#' }
"G1"

#' Biomarker Set 1
#'
#' A simulated dataset containing one of the components to run est_lucid,
#' plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 4 variables:
#' \describe{
#'   \item{CZ1, CZ2}{Causal biomarkers}
#'   \item{NZ1, NZ2}{Null biomarkers}
#' }
"Z1"

#' Outcome Set 1
#'
#' A simulated dataset containing one of the components to run est_lucid,
#' plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 1 variable:
#' \describe{
#'   \item{Y1}{A binary outcome}
#' }
"Y1"

#' Covariate Set in the G->X path
#'
#' A simulated dataset containing one of the optional components to run est_lucid,
#' plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 5 variables:
#' \describe{
#'   \item{C1 - C3}{Three continuous covariates}
#'   \item{C4, C5}{Two binary covariates}
#' }
"CoG"

#' Genetic Features Set 2
#'
#' A simulated dataset containing one of the components to run sem_lucid.
#' The variables are as follows:
#'
#' @format A set with 2000 rows and 10 variables:
#' \describe{
#'   \item{CG1 - CG5}{Causal SNPs}
#'   \item{NG1 - NG5}{Null SNPs}
#' }
"G2"

#' Biomarker Set 2
#'
#' A simulated dataset containing one of the components to run sem_lucid.
#' The variables are as follows:
#'
#' @format A set with 2000 rows and 4 variables:
#' \describe{
#'   \item{CZ1, CZ2}{Causal biomarkers}
#'   \item{NZ1, NZ2}{Null biomarkers}
#' }
"Z2"

#' Outcome Set 2
#'
#' A simulated dataset containing one of the components to run sem_lucid.
#' The variables are as follows:
#'
#' @format A set with 2000 rows and 1 variable:
#' \describe{
#'   \item{Y2}{A continuous outcome}
#' }
"Y2"
