library(boot)
library(MASS)
library(nnet)
library(pROC)
library(methods)
library(LUCIDus)

##### Simulate Data with ContY #####
Npop <- 200000
Nsample <- 2000

Nrep <- 100

M <- 5 #number of causal genes
NM <- 5 #number of non-causal genes
Q <- 2 #number of causal biomarkers
NQ <- 2 #number of non-causal biomarkers

K <- 2 #number of latent clusters

# MAF of genes
MAF_M <- rep(0.1,M)
MAF_NM <- rep(0.1,NM)

# Beta: A|G
Beta <- mat.or.vec(K,M+1)
OR_G <- 1.9
Beta[2,] <- c(-0.0,rep(log(OR_G),M))

# W: Z|A
W <- mat.or.vec(K,Q)
W[1,] <- rep(-0.4,Q)
W[2,] <- W[1,]+0.8

NW <- mat.or.vec(K,NQ)

W <- cbind(W,NW)

# Sigma: Z|A
Sigma <- replicate(K,diag(Q+NQ),simplify=F)

# Gamma: Y|A
Gamma <- c(-0.1,0.1,1,1)

# Storage Objectives
sample_w_hat <- replicate(Nrep,mat.or.vec(K,Q+NQ),simplify=F)
sample_g_hat <- replicate(Nrep,mat.or.vec(1,2*K),simplify=F)
sample_b_hat <- replicate(Nrep,mat.or.vec(K,M+NM+1),simplify=F)
sample_s_hat <- replicate(Nrep,replicate(K,diag(0,nrow=Q+NQ,ncol=Q+NQ),simplify=F),simplify=F)

sample_auc_no_selection_se <- mat.or.vec(Nrep,1)

# use pop parms as start point
start_b <- mat.or.vec(K,M+NM+1)
start_b[,1:(M+1)] <- Beta
start_m <- W
start_g <- Gamma
start_s <- Sigma

set.seed(15)

# generate population
# generate genotype
c_gene1 <- array(0,dim=c(Npop,M))
c_gene2 <- array(0,dim=c(Npop,M))
c_gene1 <- t(apply(c_gene1,1,function(x)return(ifelse(runif(length(x))<MAF_M,1,0))))
c_gene2 <- t(apply(c_gene2,1,function(x)return(ifelse(runif(length(x))<MAF_M,1,0))))

c_gene <- data.frame(c_gene1+c_gene2)
names(c_gene) <- c(paste("CG",1:M,sep=""))

n_gene1 <- array(0,dim=c(Npop,NM))
n_gene2 <- array(0,dim=c(Npop,NM))
n_gene1 <- t(apply(n_gene1,1,function(x)return(ifelse(runif(length(x))<MAF_NM,1,0))))
n_gene2 <- t(apply(n_gene2,1,function(x)return(ifelse(runif(length(x))<MAF_NM,1,0))))

n_gene <- data.frame(n_gene1+n_gene2)
names(n_gene) <- c(paste("NG",1:NM,sep=""))

gene <- cbind(c_gene,n_gene)
apply(gene,2,table)/Npop #check the genotype distribution of simulated population

e_cgene <- apply(c_gene,2,mean)

pop <- gene

# generate latent cluster

prA <-exp(as.matrix(cbind(intercept=rep(1,Npop),t(t(c_gene)-e_cgene)))%*%t(Beta))
prA <-data.frame(t(apply(prA,1,function(x)return(x/sum(x)))))
names(prA) <- c(paste("P",1:K,sep=""))

cluster <- apply(prA,1,function(x){sample(c(1:K),size=1,replace=F,prob=x)})

e_cluster <- table(cluster)/length(cluster)

# generate biomarkers
gen_biomarker <- function(x){
  return(mvrnorm(1,mu=W[x,],Sigma=Sigma[[x]]))
}

biomarker <- data.frame(t(apply(matrix(cluster),1,gen_biomarker)))

names(biomarker) <- c(paste("CZ",1:Q,sep=""),paste("NZ",1:NQ,sep=""))

apply(biomarker[cluster==1,],2,mean) #check population mean
apply(biomarker[cluster==2,],2,mean) #check population mean

cov(biomarker[cluster==1,]) #check population cov
cov(biomarker[cluster==2,]) #check population cov

# combine gene, cluster and biomarker
pop <- cbind(gene,cluster,biomarker)

# generate continuous outcome Y

Y <- data.frame(apply(matrix(cluster),1,function(x){rnorm(1,mean=Gamma[x],sd=sqrt(Gamma[x+K]))}))

names(Y) <- "Y"

pop <- cbind(gene,cluster,biomarker,Y)

mean(pop[which(pop$cluster==1), "Y"])
mean(pop[which(pop$cluster==2), "Y"])

#Bootstapping Estimate and CI

lucid_par <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample 
  try_no_selection <- try(t_no_selection_se <- est_lucid(G=d[,1:10],Z=d[,12:15],Y=as.matrix(d[,16]),K=2,useY=TRUE,family="normal",Pred=TRUE,
                                                         initial = def_initial(init_b=start_b,init_m=start_m,init_s=start_s,init_g=start_g),
                                                         tunepar = def_tune(),itr_tol = def_tol(MAX_ITR=2000)))
  # par_lucid <- c(try_no_selection$beta[2,],as.vector(t(try_no_selection$mu)),try_no_selection$gamma[1:2]) # original parameters
  par_lucid <- c(try_no_selection$beta[2,],diff(try_no_selection$mu),diff(try_no_selection$gamma[1:2])) # delta Z&Y
  return(par_lucid) 
}

# Parallel Realization
library(doParallel)

# Calculate the number of cores
no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)

#Start parallel computing
registerDoParallel(cl)
getDoParWorkers()

#Set # of Replication
BootCI <- foreach(icount(Nrep), 
                  .combine = list,
                  .multicombine = TRUE,
                  .maxcombine = 2000,
                  .packages = c("parallel", "glmnet", "glasso", "mvtnorm", "nnet", "lbfgs", "stats", "Matrix", "boot", "LUCIDus"))  %dopar%{
                    sampledata <- pop[sample(nrow(pop),Nsample,replace=F),]
                    bootcpu <- detectCores()
                    results <- boot(data=sampledata, statistic=lucid_par, R=50, parallel="multicore", ncpus=bootcpu)
                    return(list(Sampledata = sampledata, Boot = results))
                  }

stopCluster(cl)

#Save the results to local
setwd("/auto/pmd-01/chengpen/LUCid/Results")
save(BootCI, file = "Bootstrapping_Results_Rep100.RData")

# #Summarize the bootstrapping results
# sapply(2:6, function(i) plot(results, index = i))
# mean(as.vector(results$t[,2:6])); sd(as.vector(results$t[,2:6]))
# hist(as.vector(results$t[,2:6]), freq = F)
# abline(v=mean(as.vector(results$t[,2:6])), col=2, lty=2, lwd=2)
# abline(v=log(1.9), col=3, lty=1, lwd=2)
# 
# sapply(7:11, function(i) plot(results, index = i))
# mean(as.vector(results$t[,7:11])); sd(as.vector(results$t[,7:11]))
# hist(as.vector(results$t[,7:11]), freq = F)
# abline(v=mean(as.vector(results$t[,7:11])), col=2, lty=2, lwd=2)
# abline(v=0, col=3, lty=1, lwd=2)
# 
# sapply(12:13, function(i) plot(results, index = i))
# mean(as.vector(results$t[,12:13])); sd(as.vector(results$t[,12:13]))
# hist(as.vector(results$t[,12:13]), freq = F)
# abline(v=mean(as.vector(results$t[,12:13])), col=2, lty=2, lwd=2)
# abline(v=-0.4, col=3, lty=1, lwd=2)
# 
# sapply(16:17, function(i) plot(results, index = i))
# mean(as.vector(results$t[,16:17])); sd(as.vector(results$t[,16:17]))
# hist(as.vector(results$t[,16:17]), freq = F)
# abline(v=mean(as.vector(results$t[,16:17])), col=2, lty=2, lwd=2)
# abline(v=0.4, col=3, lty=1, lwd=2)
# 
# sapply(20:21, function(i) plot(results, index = i))
# mean(as.vector(results$t[,20])); sd(as.vector(results$t[,20]))
# mean(as.vector(results$t[,21])); sd(as.vector(results$t[,21]))
# 
# 
# boot.ci(results, type = c("norm","basic", "perc"), index = 21)
# 
# sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$normal[2:3]), cbind)
# sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$basic[4:5]), cbind)
# sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$perc[4:5]), cbind)
# 
# cbind(t(sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$normal[2:3]), cbind)),
#       t(sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$basic[4:5]), cbind)),
#       t(sapply(lapply(2:21, function(i) (boot.ci(results, type = c("norm","basic", "perc"), index = i))$perc[4:5]), cbind)))
