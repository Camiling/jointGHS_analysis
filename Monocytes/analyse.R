rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")
run_jointGHS = FALSE
run_boot = TRUE

# Use jointGHS on data ---------------------------------------------------------

if(run_jointGHS){
  res.joint = jointGHS::jointGHS(expr_all, scale=T, AIC_selection = T, AIC_eps = 1,eps=1e-3) # AIC_eps=0.1 for try 2
  save(res.joint,file='data/expression_data_monocytes_LYZ_region_FDR_5_jointGHS.RData')
}
#load('data/expression_data_monocytes_LYZ_region_FDR_5_jointGHS.RData')





# Investigate properties of results --------------------------------------------
 
condition.names = c('IFN-gamma', 'LPS-2h', 'LPS-24h', 'Unstim')
# Get precision matrices from joint estimates
thetas.est.mono = lapply(res.joint$theta, cov2cor)
for(k in 1:length(thetas.est.mono)){
  thetas.est.mono[[k]][which(abs(thetas.est.mono[[k]]) < 1e-5, arr.ind = T)] = 0
}
# Get precision matrices of single network estimates
thetas.est.mono.single = lapply(res.joint$theta_single, cov2cor)
for(k in 1:length(thetas.est.mono.single)){
  thetas.est.mono.single[[k]][which(abs(thetas.est.mono.single[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Sparsities
unlist(lapply(thetas.est.mono, FUN = function(s) tailoredGlasso::sparsity(s!=0)))

# Sparsities of single network estimates
unlist(lapply(thetas.est.mono.single, FUN = function(s) tailoredGlasso::sparsity(s!=0)))


# Confusion matrices of joint estimates

# 1 and 2 (IFN-gamma and LPS-2h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[2]]=!0)

# 1 and 3 (IFN-gamma and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[3]]=!0)

# 1 and 4 (IFN-gamma and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[4]]=!0)

# 2 and 3 (LPS-2h and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[3]]=!0)

# 2 and 4 (LPS-2h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[4]]=!0)

# 3 and 4 (LPS-24h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,thetas.est.mono[[4]]=!0)


# How much did the networks change when using a joint approach instead of a single?

# IFN-gamma
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono.single[[1]]=!0)

# LPS-2h
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono.single[[2]]=!0)

# LPS-24h
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,thetas.est.mono.single[[3]]=!0)

# Unstim
tailoredGlasso::confusion.matrix(thetas.est.mono[[4]]!=0,thetas.est.mono.single[[4]]=!0)

# Plot MCC of results -----------------------------------------------------------------

# MCC of joint estimates
MCC_vals_all = sapply(thetas.est.mono, function(x) sapply(thetas.est.mono, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
rownames(MCC_vals_all)=colnames(MCC_vals_all) = condition.names

# MCC of single estimates
MCC_vals_all_single = sapply(thetas.est.mono.single, function(x) sapply(thetas.est.mono.single, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
rownames(MCC_vals_all_single)=colnames(MCC_vals_all_single) = condition.names


png("monocytes_MCC_jointGHS.png", width=1200, height =1200)
superheat(MCC_vals_all,title = 'jointGHS')
dev.off()

png("monocytes_MCC_fastGHS.png", width=1200, height =1200)
superheat(MCC_vals_all_single,title = 'fastGHS')
dev.off()

# Perform clustering

png("monocytes_MCC_jointGHS_clustered.png", width=1200, height =1200)
superheat(MCC_vals_all,title = 'jointGHS', row.dendrogram = T, col.dendrogram = T)
dev.off()

png("monocytes_MCC_fastGHS_clustered.png", width=1200, height =1200)
superheat(MCC_vals_all_single,title = 'fastGHS', row.dendrogram = T, col.dendrogram = T)
dev.off()

# Plot networks ----------------------------------------------------------------





# Investigate what is unique to the different networks -------------------------






# Use Bayesian bootstrap -------------------------------------------------------

nCores=56
if(run_boot){
  res.boot = jointGHS::jointGHS(expr_all, scale=T, AIC_selection = T, AIC_eps = 1,eps=1e-3, boot_check = T, boot_lambda = F,B=100, nCores = nCores) # AIC_eps=0.1 for try 2
  save(res.boot,file='data/expression_data_monocytes_LYZ_region_FDR_5_jointGHS_boot.RData')
}

plot(res.boot,k=1)
plot(res.boot,k=2)
plot(res.boot,k=3)
plot(res.boot,k=4)

print(res.boot,k=1)
print(res.boot,k=2)
print(res.boot,k=3)
print(res.boot,k=4)


