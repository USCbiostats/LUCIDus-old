LUCid: Latent or Unobserved Clustering with integrated data
================

Introduction
============

An R package allows users to achieve a joint estimation of latent or unobserved clusters using multi-omics data with/without the outcome of interest.

Getting Started
---------------

This supplementary package is based on the research paper "Integrative Latent Cluster Assignment Using Multi-Omics Data with Phenotypic Traits" (Under development, will replace this field with a citation when we have one).

### Prerequisites

R (&gt;= 3.1.0).

    Package Dependencies: "mvtnorm", "nnet", "glmnet", "glasso", "lbfgs".

### Installing

We have a plan to submit the package to CRAN/Bioconductor when the development cycle is done.

For now, it can be installed from GitHub using the following codes:

    install.packages("devtools")
    devtools::install_github("USCbiostats/LUCid")

Otherwise, one can download the package from GitHub, and run the following codes from the parent working directory that contains the LUCid folder:

    install.packages("devtools")
    setwd("..")
    devtools::install("LUCid")

``` r
library(LUCid)
```

Fitting the latent cluster models
---------------------------------

Three functions, including *est\_cluster*, *sem\_cluster*, & *tune\_cluster*, are currently available for model fitting and selection.

### *est\_cluster*

Estimating latent clusters with multi-omics data

#### Example

``` r
# For a testing dataset with 10 genetic features (5 causal) and 4 biomarkers (2 causal)
ncluster <- 2
IntClusFit <- est_cluster(G=G1,Z=Z1,Y=Y1,K=ncluster,family="binary",Pred=TRUE,
                          initial=def_initial(), itr_tol=def_tol(),
                          tunepar = def_tune(Select_G=TRUE,Select_Z=TRUE,
                                             Rho_G=0.02,Rho_Z_InvCov=0.1,Rho_Z_CovMu=93))

# Identify selected features
G_diff <- apply(apply(IntClusFit$beta,2,range),2,function(x){x[2]-x[1]})[-1]
select_G <- G_diff != 0
InvSigmaMu <- solve(IntClusFit$sigma[[1]])%*%IntClusFit$mu[1,]
for(i in 2:ncluster){
    InvSigmaMu <- cbind(InvSigmaMu, solve(IntClusFit$sigma[[i]])%*%IntClusFit$mu[i,])
}
Z_diff <- apply(apply(InvSigmaMu,1,range),2,function(x){x[2]-x[1]})
select_Z <- Z_diff != 0
sum(select_G); sum(select_Z)
colnames(G1)[select_G]; colnames(Z1)[select_Z]

# Select the features
if(!all(select_G==FALSE)){
    G_select <- G1[,select_G]
}
if(!all(select_Z==FALSE)){
    Z_select <- Z1[,select_Z]
}

# Re-fit with selected features
set.seed(10)
IntClusFitFinal <- est_cluster(G=G_select,Z=Z_select,Y=Y1,K=ncluster,family="binary",Pred=TRUE,
                               initial=def_initial(), itr_tol=def_tol(), tunepar = def_tune())
```

#### *Visualize the results using Sankey diagram*

``` r
plot_cluster(IntClusFitFinal)
```

![](README_files/figure-markdown_github/Sankey1.png)

### *sem\_cluster*

Supplemented EM-algorithm for latent cluster estimation

#### Example

    sem_cluster(G=G2,Z=Z2,Y=Y2,useY=TRUE,K=2,Pred=TRUE,family="normal",Get_SE=TRUE,
                def_initial(),def_tol(MAX_ITR=1000,MAX_TOT_ITR=3000))

### *tune\_cluster*

Grid search for tuning parameters using parallel computing

#### Example

    # Better be run on a server or HPC
    GridSearch <- tune_cluster(G=G1, Z=Z1, Y=Y1, K=2, Family="binary", USEY = TRUE,
                               LRho_g = 0.001, URho_g = 0.1, NoRho_g = 10,
                               LRho_z_invcov = 0.05, URho_z_invcov = 0.3, NoRho_z_invcov = 6,
                               LRho_z_covmu = 75, URho_z_covmu = 100, NoRho_z_covmu = 6)
    GridSearch$Results
    GridSearch$Optimal

For more details, see documentations for each function in the R package.

Built With
----------

-   [devtools](https://cran.r-project.org/web/packages/devtools/index.html) - Tools to Make Developing R Packages Easier
-   [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html) - In-Line Documentation for R

Versioning
----------

The current version is 0.3.0.

For the versions available, see the [Version on this repository](https://github.com/your/project/Version) (not available yet).

Authors
-------

-   **Cheng Peng**

License
-------

This project is licensed under the GPL-2 License.

Acknowledgments
---------------

-   Dr. David V. Conti
-   Dr. Zhao Yang
-   USC IMAGE P01 Group.
