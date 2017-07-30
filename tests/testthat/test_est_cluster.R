context("LCid Testing")

test_that("Test for est_cluster", {
  fit_selection <- est_cluster(G=G1,Z=Z1,Y=Y1,K=2,
                               init_b = NULL, init_m = NULL, init_s = NULL, init_g = NULL,
                               family="binary",Pred=TRUE,Select_G=TRUE,Rho_G=0.02,
                               Select_Z=TRUE,Rho_Z_InvCov=0.1,Rho_Z_CovMu=93,
                               tol_m = 1e-8,tol_b=1e-8,tol_s=1e-8,tol_g=1e-8,MAX_ITR = 800,MAX_TOT_ITR=800)

  beta0 <- 5
  beta1 <- abs(fit_selection$mu[2,1])
  expect_length(fit_selection, 10)
  expect_equal(beta0, beta1, tolerance = 0.05)
})
