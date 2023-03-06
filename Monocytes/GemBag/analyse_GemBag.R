rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5_notOnGit.RData")
source('Monocytes/GemBag/GemBag_implementation/BIC_GemBag.R')
Rcpp::sourceCpp('Monocytes/GemBag/GemBag_implementation/GemBag-algo.cpp')


run_GemBag = TRUE

# Analyse data and do some descriptive analyses


# Use jointGHS on data ---------------------------------------------------------

if(run_GemBag){
  # Perform GemBag
  S_l = lapply(expr_all, cov)
  n.vals = unlist(lapply(expr_all, nrow))
  p = ncol(expr_all[[1]])
  # Set of hyperparameters
  v0_l <- c(0.25, 0.5, 0.75, 1) * sqrt(1/n.vals[1]/log(p))
  v1_l <- c(2.5, 5, 7.5, 10) * sqrt(1/n.vals[1]/log(p))
  p1 <- 0.5
  # Tuning by BIC
  hyper <- Tune_GemBag(v0_l, v1_l,S_l, n.vals, maxiter=20, p1, p_2=0.5)
  v0 <- hyper$v0
  v1 <- hyper$v1
  # Final estimate
  res.gembag <- GemBag(S_l=S_l, n= n.vals, 
                       v_0=v0, v_1=v1, tau=v0, 
                       p_1=p1, p_2=1,
                       maxiter=20)
  names(res.gembag) <- c('Theta', 'P', 'W')
  save(res.gembag,file='Monocytes/data/monocytes_GemBag.RData')
}

load('Monocytes/data/monocytes_GemBag.RData.RData')
