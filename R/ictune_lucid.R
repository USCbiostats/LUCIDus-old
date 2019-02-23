#Grid-search ictune_lucid
ictune_lucid <- function(G = NULL, CoG = NULL, Z = NULL, CoY = NULL, Y, IG=NULL, IY=NULL, K, Family, USEY = TRUE,
                         initial = def_initial(),
                         LRho_g, URho_g, NoRho_g,
                         LRho_z_invcov, URho_z_invcov, NoRho_z_invcov,
                         LRho_z_covmu, URho_z_covmu, NoRho_z_covmu,
                         NoCores = detectCores()-1){
  
  e=f=g=rho_g=rho_z_covmu=rho_z_invcov=NULL
  
  M <- dim(G)[2]
  Q <- dim(Z)[2]
  
  parallel_cluster <- function(Lrho_g, Urho_g, Norho_g, Lrho_z_invcov, Urho_z_invcov, Norho_z_invcov, Lrho_z_covmu, Urho_z_covmu, Norho_z_covmu){
    foreach(rho_g = seq(Lrho_g, Urho_g, length.out=Norho_g)) %:%
      foreach(rho_z_invcov = seq(Lrho_z_invcov, Urho_z_invcov, length.out=Norho_z_invcov)) %:%
      foreach(rho_z_covmu = seq(Lrho_z_covmu, Urho_z_covmu, length.out=Norho_z_covmu),
              .combine = list, .multicombine = TRUE, .maxcombine = 2000, .errorhandling = 'pass',
              .export=c("G", "CoG", "Z", "CoY", "Y", "IG", "IY", "K", "Family", "USEY", "initial"),
              .packages = c("glmnet", "glasso", "mvtnorm", "nnet", "lbfgs", "stats", "Matrix", "LUCIDus"))  %dopar%{
                set.seed(rho_g*rho_z_invcov*rho_z_covmu)
                icest_lucid(G=G,CoG=CoG,Z=Z,CoY=CoY,Y=Y,IG=IG,IY=IY,K=K,useY=USEY,family=Family,Pred=TRUE,
                            initial = initial, tunepar = def_tune(Select_G=T,Select_Z=T,Rho_G=rho_g,Rho_Z_InvCov=rho_z_invcov,Rho_Z_CovMu=rho_z_covmu),
                            def_tol(MAX_ITR = 500,MAX_TOT_ITR=1000),itr_tol = def_tol(tol_b = 1e-03, tol_m = 1e-03, tol_s = 1e-03, tol_g = 1e-03))
              }
  }
  
  parallel_results <- function(Norho_g, Norho_z_invcov, Norho_z_covmu){
    if(Norho_z_covmu==1){
      foreach(e=1:Norho_g, .combine = 'rbind') %:%
        foreach(f=1:Norho_z_invcov, .combine = 'rbind', .errorhandling = 'remove',
                .export=c("modelfits", "K"), .packages=c("stats", "LUCIDus")) %do%{
                  
                  Non0g <- summary_lucid(modelfits[[e]][[f]])$No0G
                  Non0z <- summary_lucid(modelfits[[e]][[f]])$No0Z
                  
                  bic <- summary_lucid(modelfits[[e]][[f]])$BIC
                  #Other types of GIC
                  gic1 <- summary_lucid(modelfits[[e]][[f]])$GIC1
                  gic2 <- summary_lucid(modelfits[[e]][[f]])$GIC2
                  
                  data.frame(Rho_G=modelfits[[e]][[f]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic, GIC1=gic1, GIC2=gic2)
                }
    }else{
      foreach(e=1:Norho_g, .combine = 'rbind') %:%
        foreach(f=1:Norho_z_invcov, .combine = 'rbind') %:%
        foreach(g=1:Norho_z_covmu, .combine = 'rbind', .errorhandling = 'remove',
                .export=c("modelfits", "K"), .packages=c("stats", "LUCIDus")) %do%{
                  
                  Non0g <- summary_lucid(modelfits[[e]][[f]][[g]])$No0G
                  Non0z <- summary_lucid(modelfits[[e]][[f]][[g]])$No0Z
                  
                  bic <- summary_lucid(modelfits[[e]][[f]][[g]])$BIC
                  #Other types of GIC
                  gic1 <- summary_lucid(modelfits[[e]][[f]][[g]])$GIC1
                  gic2 <- summary_lucid(modelfits[[e]][[f]][[g]])$GIC2
                  
                  data.frame(Rho_G=modelfits[[e]][[f]][[g]]$rho_g, Rho_Z_InvCov=modelfits[[e]][[f]][[g]]$rho_z_InvCov, Rho_Z_CovMu=modelfits[[e]][[f]][[g]]$rho_z_CovMu, Non0G=Non0g, Non0Z=Non0z, BIC=bic, GIC1=gic1, GIC2=gic2)
                }
    }
  }
  
  # Initiate cluster
  cl <- makeCluster(NoCores)
  #Start parallel computing
  registerDoParallel(cl)
  
  modelfits <- parallel_cluster(Lrho_g = LRho_g, Urho_g = URho_g, Norho_g = NoRho_g,
                                Lrho_z_invcov = LRho_z_invcov, Urho_z_invcov = URho_z_invcov, Norho_z_invcov = NoRho_z_invcov,
                                Lrho_z_covmu = LRho_z_covmu, Urho_z_covmu = URho_z_covmu, Norho_z_covmu = NoRho_z_covmu)
  
  stopCluster(cl)
  
  results <- parallel_results(Norho_g = NoRho_g, Norho_z_invcov = NoRho_z_invcov, Norho_z_covmu = NoRho_z_covmu)
  
  #Show the BEST Model with minBIC
  optimal <- results[which.min(results$BIC), ]
  
  return(list(Modelfits = modelfits, Results = results, Optimal = optimal))
}