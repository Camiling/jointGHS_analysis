rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")
run_SSJGL = TRUE

source("Monocytes/SSJGL/R/SSJGL.R")
source("Monocytes/SSJGL/R/JGL.R")
source("Monocytes/SSJGL/R/admm.iters.R")
source("Monocytes/SSJGL/R/eval.R")
source("Monocytes/SSJGL/R/gete.R")

# Illustrate that the SSJGL fails


# Use jointGHS on data ---------------------------------------------------------

if(run_SSJGL){
  lambda1 <- 1 
  lambda2 <- 1 
  lambda.eff2 <- lambda1 + 1 + c(0:20)*10
  v0s <- lambda1/lambda.eff2
  res.ssjgl = SSJGL(Y=expr_all,penalty='fused',lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = 1, v0s = v0s, tol.em=1e-3, a=1, b=1, doubly=FALSE, c = 0.01)
  save(res.ssjgl,file='Monocytes/data/monocytes_SSJGL.RData')
}

load('Monocytes/data/monocytes_SSJGL.RData')

