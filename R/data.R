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

#' Incomplete Biomarker Set 1
#'
#' A simulated dataset containing one of the components to run est_lucid incomplete
#' option, plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 4 variables, 500 rows are NAs:
#' \describe{
#'   \item{CZ1, CZ2}{Causal biomarkers}
#'   \item{NZ1, NZ2}{Null biomarkers}
#' }
"Z1_Incomp"

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
#'   \item{GC1 - GC3}{Three continuous covariates}
#'   \item{GC4, GC5}{Two binary covariates}
#' }
"CoG"

#' Covariate Set in the X->Y path
#'
#' A simulated dataset containing one of the optional components to run est_lucid,
#' plot_lucid, and tune_lucid. The variables are as follows:
#'
#' @format A set with 2000 rows and 5 variables:
#' \describe{
#'   \item{YC1 - YC3}{Three continuous covariates}
#'   \item{YC4, YC5}{Two binary covariates}
#' }
"CoY"

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
