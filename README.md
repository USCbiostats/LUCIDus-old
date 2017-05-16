# LCid

A R package allows users to achieve joint estimation of latent clusters using multi-omics data with/without outcome of interest.

## Getting Started

This supplementary package is based on the research paper "Integrated Analysis of Germline, Omic and Disease Data" (Under development, will replace this field with a citation when we have one). 

### Prerequisites

R (>= 3.1.0).

```
Package Dependencies: "mvtnorm", "nnet", "glmnet", "glasso", "lbfgs".
```

### Installing

We have a plan to submit the package to CRAN/Bioconductor when the development cycle is done.

For now, it can be installed from GitHub using the following codes:

```
install.packages("devtools")
devtools::install_github("USCbiostats/LCid")
```
Otherwise, one can download the package from GitHub, and run the following codes from the parent working directory that contains the LCid folder:

```
install.packages("devtools")
setwd("..")
devtools::install("LCid")
```

## Fitting the latent cluster models

Three functions, including *est_cluster*, *sem_cluster*, & *tune_cluster*, are currently available for model fitting. 

### *est_cluster*

Estimating latent clusters with multi-omics data

#### Example

```
est_cluster(G=G1,Z=Z1,Y=Y1,K=2,
            init_b = NULL, init_m = NULL, init_s = NULL, init_g = NULL,
            family="binary",Pred=TRUE,Select_G=TRUE,Select_Z=TRUE,
            Rho_G=0.02,Rho_Z_InvCov=0.1,Rho_Z_CovMu=93,
            tol_m = 1e-8,tol_b=1e-8,tol_s=1e-8,tol_g=1e-8,MAX_ITR = 800,MAX_TOT_ITR=800)
```

### *sem_cluster*

Supplemented EM-algorithm for latent cluster estimation

#### Example

```
sem_cluster(G=G2,Z=Z2,Y=Y2,useY=TRUE,K=2,
            init_b = NULL, init_m = NULL, init_s = NULL, init_g = NULL,
            Pred=TRUE,family="normal",Get_SE=TRUE,Ad_Hoc_SE=TRUE,
            tol_m = 1e-6,tol_b=1e-6,tol_s=1e-6,tol_g=1e-6,tol_sem=1e-3,
            MAX_ITR=3000,MAX_TOT_ITR=3000)
```

### *tune_cluster*

Grid search for tuning parameters using parallel computing

#### Example

```
tune_cluster(G=G1, Z=Z1, Y=Y1, K=2, Family="binary", USEY = TRUE,
             LRho_g = 0.005, URho_g = 0.015, NoRho_g = 3,
             LRho_z_invcov = 0.1, URho_z_invcov = 0.2, NoRho_z_invcov = 2,
             LRho_z_covmu = 81, URho_z_covmu = 85, NoRho_z_covmu = 5)
```
For more details, see documentations for each function in the R package.


## Built With

* [devtools](https://cran.r-project.org/web/packages/devtools/index.html) - Tools to Make Developing R Packages Easier
* [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html) - In-Line Documentation for R

## Versioning

The current version is 0.3.0.

For the versions available, see the [Version on this repository](https://github.com/your/project/Version) (not available yet). 

## Authors

* **Cheng Peng**

## License

This project is licensed under the GPL-2 License.

## Acknowledgments

* Dr. Zhao Yang
* Dr. David V. Conti
* etc.
