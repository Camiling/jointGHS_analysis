rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")
run_jointGHS = FALSE
run_boot = FALSE

# Analyse data and do some descriptive analyses


# Use jointGHS on data ---------------------------------------------------------

if(run_jointGHS){
  res.joint = jointGHS::jointGHS(expr_all, scale=T, AIC_selection = T, AIC_eps = 5,eps=1e-3, stop_overflow = TRUE) 
  save(res.joint,file='Monocytes/data/monocytes_jointGHS_smalleps_largeAIC.RData')
}

load('Monocytes/data/monocytes_jointGHS_smalleps_largeAIC.RData')






# Investigate properties of results --------------------------------------------

p = dim(expr_all[[1]])[2]
n.vals = unlist(lapply(expr_all, nrow))
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
as.data.frame(rbind(condition.names, unlist(lapply(thetas.est.mono, FUN = function(s) tailoredGlasso::sparsity(s!=0)))))
#           IFN-gamma             LPS-2h            LPS-24h             Unstim
#   0.0177372565271446 0.0202928581295759 0.0208868628263572 0.016369664318276

# Sparsities of single network estimates
as.data.frame(rbind(condition.names, unlist(lapply(thetas.est.mono.single, FUN = function(s) tailoredGlasso::sparsity(s!=0)))))
#         IFN-gamma             LPS-2h            LPS-24h             Unstim
#   0.0305843348528802 0.0394805912418842 0.037007874015748 0.0260809504075149

# Number of edges in each network


# Confusion matrices of joint estimates

# 1 and 2 (IFN-gamma and LPS-2h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[2]]!=0)
#     [,1]  [,2]
#[1,]  587   882
#[2,]  697 70224

# 1 and 3 (IFN-gamma and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[3]]!=0)
#     [,1]  [,2]
#[1,]  606   906
#[2,]  678 70200

# 1 and 4 (IFN-gamma and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[4]]!=0)
#     [,1]  [,2]
#[1,]  608   577
#[2,]  676 70529

# 2 and 3 (LPS-2h and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[3]]!=0)
#     [,1]  [,2]
#[1,]  645   867
#[2,]  824 70054

# 2 and 4 (LPS-2h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[4]]!=0)
#    [,1]  [,2]
#[1,]  576   609
#[2,]  893 70312

# 3 and 4 (LPS-24h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,thetas.est.mono[[4]]!=0)
#     [,1]  [,2]
#[1,]  579   606
#[2,]  933 70272

# How much did the networks change when using a joint approach instead of a single?

# IFN-gamma
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono.single[[1]]!=0)
#     [,1]  [,2]
#[1,] 1024  1190
#[2,]  260 69916

# LPS-2h
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono.single[[2]]!=0)
#     [,1]  [,2]
#[1,] 1217  1641
#[2,]  252 69280

# LPS-24h
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,thetas.est.mono.single[[3]]!=0)
#     [,1]  [,2]
#[1,] 1240  1439
#[2,]  272 69439

# Unstim
tailoredGlasso::confusion.matrix(thetas.est.mono[[4]]!=0,thetas.est.mono.single[[4]]!=0)
#     [,1]  [,2]
#[1,]  917   971
#[2,]  268 70234


# Did the loglikelihood improve compared to single network version? -------------------------------------------------------------

sample.covs = lapply(expr_all, FUN = function(x) cov(scale(x)))

loglik.ifn = gaussianloglik_util(sample.covs[[1]], res.joint$theta[[1]], n.vals[1])
loglik.ifn 
# -60516.67
loglik.lps2 = gaussianloglik_util(sample.covs[[2]], res.joint$theta[[2]], n.vals[2])
loglik.lps2
# -77967.44
loglik.lps24 = gaussianloglik_util(sample.covs[[3]], res.joint$theta[[3]], n.vals[3])
loglik.lps24
# -91633.02
loglik.unstim = gaussianloglik_util(sample.covs[[4]], res.joint$theta[[4]], n.vals[4])
loglik.unstim 
# -105617.3

res.joint$tau_sq
# 8.001 6.401 9.001 3.201

# Get single networks to same sparsity
set.seed(123)
res.single.1 = fastGHS::fastGHS(scale(expr_all[[1]]), tau_sq=0.01397, fix_tau = T, epsilon = 1e-3) 
theta.single.1 = cov2cor(res.single.1$theta)
theta.single.1[which(abs(theta.single.1)<1e-5,arr.ind=T)]=0
tailoredGlasso::sparsity(theta.single.1)
# 0.01769581
loglik.ifn.single = gaussianloglik_util(sample.covs[[1]], res.single.1$theta, n.vals[1])
loglik.ifn.single
# -63271.96

set.seed(123)
res.single.2 = fastGHS::fastGHS(scale(expr_all[[2]]), tau_sq=0.01674, fix_tau = T, epsilon = 1e-3) # 0.01673001045 0.01673001047
theta.single.2 = cov2cor(res.single.2$theta)
theta.single.2[which(abs(theta.single.2)<1e-5,arr.ind=T)]=0
tailoredGlasso::sparsity(theta.single.2)
# 0.02030667
loglik.lps2.single =gaussianloglik_util(sample.covs[[2]], res.single.2$theta, n.vals[2])
loglik.lps2.single
# -78843.15

set.seed(123)
res.single.3 = fastGHS::fastGHS(scale(expr_all[[3]]), tau_sq=0.015042, fix_tau = T, epsilon = 1e-3) # 0.1504 # 0.15044
theta.single.3 = cov2cor(res.single.3$theta)
theta.single.3[which(abs(theta.single.3)<1e-5,arr.ind=T)]=0
tailoredGlasso::sparsity(theta.single.3)
# 0.02088686
loglik.lps24.single = gaussianloglik_util(sample.covs[[3]], res.single.3$theta, n.vals[3])
loglik.lps24.single 
# -91731.75

set.seed(123)
res.single.4 = fastGHS::fastGHS(scale(expr_all[[4]]), tau_sq=0.01105, fix_tau = T, epsilon = 1e-3) #0.01105 #0.0111
theta.single.4 = cov2cor(res.single.4$theta)
theta.single.4[which(abs(theta.single.4)<1e-5,arr.ind=T)]=0
tailoredGlasso::sparsity(theta.single.4)
# 0.01635585
loglik.unstim.single = gaussianloglik_util(sample.covs[[4]], res.single.4$theta, n.vals[4])
loglik.unstim.single
# -105582.7

# Save single network estimates 
res.single = list(res.single.1, res.single.2, res.single.3, res.single.4)
save(res.single, file='Monocytes/data/res_joint_single_samespars.RData')

# Plot for better visualization
logliks.joint = c(loglik.ifn, loglik.lps2, loglik.lps24, loglik.unstim)
logliks.single = c(loglik.ifn.single, loglik.lps2.single, loglik.lps24.single, loglik.unstim.single)
df.loglik = data.frame(loglikelihood = c(logliks.joint, logliks.single), method=c(rep('joint',length(logliks.joint)), rep('single',length(logliks.single))),
                       condition = rep(condition.names,2))

pdf("Monocytes/plots/loglik_jointGHS.pdf")
ggplot2::ggplot(df.loglik, aes(y=loglikelihood, x='', group=method))+geom_point(aes(color=method))+facet_wrap(~ condition, ncol=4)+
  labs(y= "Loglikelihood", x = "")
dev.off()

df.loglik

# Plot MCC of results -----------------------------------------------------------------

# MCC of joint estimates
MCC_vals_all = sapply(thetas.est.mono, function(x) sapply(thetas.est.mono, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
diag(MCC_vals_all)=NA
rownames(MCC_vals_all)=colnames(MCC_vals_all) = condition.names

# MCC of single estimates
MCC_vals_all_single = sapply(thetas.est.mono.single, function(x) sapply(thetas.est.mono.single, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
diag(MCC_vals_all_single)=NA
rownames(MCC_vals_all_single)=colnames(MCC_vals_all_single) = condition.names


png("Monocytes/plots/monocytes_MCC_jointGHS.png", width=1200, height =1200)
superheat(MCC_vals_all,heat.lim = c(min(MCC_vals_all,na.rm = T)-0.00001, max(MCC_vals_all,na.rm = T)+0.018), heat.col.scheme = "red")
dev.off()

png("Monocytes/plots/monocytes_MCC_fastGHS.png", width=1200, height =1200)
superheat(MCC_vals_all_single,heat.lim = c(min(MCC_vals_all_single,na.rm = T), max(MCC_vals_all_single,na.rm = T)), heat.col.scheme = "red")
dev.off()

# Perform clustering

png("Monocytes/plots/monocytes_MCC_jointGHS_clustered.png", width=1200, height =1200)
superheat(MCC_vals_all, row.dendrogram = T, col.dendrogram = T,heat.lim = c(min(MCC_vals_all,na.rm = T)-0.00001, max(MCC_vals_all,na.rm = T)+0.018), heat.col.scheme = "red")
dev.off()

png("Monocytes/plots/monocytes_MCC_fastGHS_clustered.png", width=1200, height =1200)
superheat(MCC_vals_all_single, row.dendrogram = T, col.dendrogram = T,heat.lim = c(min(MCC_vals_all_single,na.rm = T), max(MCC_vals_all_single,na.rm = T)), heat.col.scheme = "red")
dev.off()

# Plot networks ----------------------------------------------------------------

# Simple plot

nets = list()
for(i in 1:length(thetas.est.mono)){
  net = network::network(thetas.est.mono[[i]]!=0,directed=F)
  nets[[i]] = GGally::ggnet2(net,node.size = 2, edge.size = 0.3,alpha=0.9,mode = "fruchtermanreingold",color = 'dodgerblue')+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

pdf("Monocytes/plots/monocytes_net_jointGHS.pdf", 12, 12)
gridExtra::grid.arrange(grobs=nets)
dev.off()

# Plot with individual edges marked

unique.list = list()
theta.joint.array = array(unlist(thetas.est.mono), dim=c(ncol(thetas.est.mono[[1]]),ncol(thetas.est.mono[[1]]),length(thetas.est.mono)))
theta.joint.array = theta.joint.array!=0
for(k in 1:length(thetas.est.mono)){
  unique.list[[k]] = ifelse(theta.joint.array[,,k] !=0, 1, 0) # 1 means edge
  unique.list[[k]][which((theta.joint.array[,,k] !=0 ) & (apply(theta.joint.array[,,-k],c(1,2),sum)==0))] = 2 # only present in network k
  unique.list[[k]][which(apply(theta.joint.array[,,],c(1,2),sum)==length(thetas.est.mono))] = 3 # present in all
}
nets = list()
for(i in 1:length(thetas.est.mono)){
  net = network::network(unique.list[[i]],directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[(net %e% "weights")+1])
  nets[[i]] = GGally::ggnet2(net,node.size = 2, edge.size = 0.3,alpha=0.9,mode = "fruchtermanreingold",color = 'dodgerblue', edge.color = 'color')+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}
pdf("Monocytes/plots/monocytes_net_unique_jointGHS.pdf", 12, 12)
gridExtra::grid.arrange(grobs=nets)
dev.off()

# Same layout for all?

# Get layout for Unstim
net = network::network(unique.list[[4]],directed=F, ignore.eval=F,names.eval='weights')
x = sna::gplot.layout.fruchtermanreingold(net, NULL)

nets = list()
for(i in 1:length(thetas.est.mono)){
  net = network::network(unique.list[[i]],directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[(net %e% "weights")+1])
  net %v% "x" = x[, 1]
  net %v% "y" = x[, 2]
  nets[[i]] = GGally::ggnet2(net,node.size = 2, edge.size = 0.3,alpha=0.9,mode = c('x','y'),color = 'dodgerblue', edge.color = 'color')+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

pdf("Monocytes/plots/monocytes_net_unique_samelayout_jointGHS.pdf", 12, 12)
gridExtra::grid.arrange(grobs=nets)
dev.off()



# Circular plot

nets = list()
for(i in 1:length(thetas.est.mono)){
  net = network::network(unique.list[[i]],directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[(net %e% "weights")+1])
  nets[[i]] = GGally::ggnet2(net,node.size = 2, edge.size = 0.3,alpha=0.9,mode = "circle",color = 'dodgerblue', edge.color = 'color')+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}
pdf("Monocytes/plots/monocytes_net_circle_jointGHS.pdf", 12, 12)
gridExtra::grid.arrange(grobs=nets)
dev.off()

# Look at selected parameters --------------------------------------------------

# Tau_sq
res.joint$tau_sq
# 8.001 6.401 9.001 3.201

hist(res.joint$E_NuInv[which(abs(res.joint$E_NuInv)>1e-7)],breaks=100)


thetas.offdiag = list()
Lambdas.offdiag = list()
for(k in 1:length(res.joint$theta)){
  theta.tmp = cov2cor(res.joint$theta[[k]])
  theta.tmp[!upper.tri(theta.tmp)] = NA
  thetas.offdiag[[k]] = c(theta.tmp)
  Lambda.tmp = res.joint$Lambda_sq[[k]]
  Lambda.tmp[!upper.tri(Lambda.tmp)] = NA
  Lambdas.offdiag[[k]] = c(Lambda.tmp)
}
# NuInv is the same for all
Nu.tmp = res.joint$E_NuInv
Nu.tmp[!upper.tri(Nu.tmp)] = NA
NuInv.offdiag = c(Nu.tmp)

df.all = data.frame(theta=unlist(thetas.offdiag), Lambda_sq=unlist(Lambdas.offdiag), NuInv = NuInv.offdiag, graph=factor(rep(1:4,each=p)))


# Plot theta_ijk's as functions of nu_ij's

p.all <- ggplot2::ggplot(na.omit(df.all),  aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.6))+
  geom_point(aes(colour=graph, shape=graph)) + geom_hline(yintercept=0, linetype='dashed', color='darkgrey') + geom_vline(xintercept=0,linetype='dashed', color='darkgrey')           

pdf('Monocytes/plots/monocytes_theta_vs_NuInv.pdf')
p.all
dev.off()





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





