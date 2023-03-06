rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_and_bcells_LYZ_region_notOnGit.RData")
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5_notOnGit.RData")
load('Monocytes/data/monocytes_GemBag.RData')


# Investigate properties of results --------------------------------------------

p = dim(expr_all[[1]])[2]
n.vals = unlist(lapply(expr_all, nrow))
condition.names = c('IFN-gamma', 'LPS-2h', 'LPS-24h', 'Unstim')
# Get precision matrices from joint estimates
thetas.est.mono = lapply(res.gembag$Theta, cov2cor)

# Sparsities
as.data.frame(rbind(condition.names, unlist(lapply(thetas.est.mono, FUN = function(s) tailoredGlasso::sparsity(s!=0)))))
#V1                  V2                  V3                  V4
#condition.names           IFN-gamma              LPS-2h             LPS-24h              Unstim
#                0.00513883133029424 0.00542892664732698 0.00519408758115762 0.00520790164387346


# Confusion matrices of joint estimates

# 1 and 2 (IFN-gamma and LPS-2h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[2]]!=0)
#[,1]  [,2]
#[1,]  313    80
#[2,]   59 71938


# 1 and 3 (IFN-gamma and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[3]]!=0)
#[,1]  [,2]
#[1,]  294    82
#[2,]   78 71936

# 1 and 4 (IFN-gamma and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,thetas.est.mono[[4]]!=0)
#    [,1]  [,2]
#[1,]  314    63
#[2,]   58 71955

# 2 and 3 (LPS-2h and LPS-24h)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[3]]!=0)
#[,1]  [,2]
#[1,]  333    43
#[2,]   60 71954

# 2 and 4 (LPS-2h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,thetas.est.mono[[4]]!=0)
#[,1]  [,2]
#[1,]  326    51
#[2,]   67 71946

# 3 and 4 (LPS-24h and Unstim)
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,thetas.est.mono[[4]]!=0)
#[,1]  [,2]
#[1,]  308    69
#[2,]   68 71945


sample.covs = lapply(expr_all, FUN = function(x) cov(scale(x)))

loglik.ifn = gaussianloglik_util(sample.covs[[1]], thetas.est.mono[[1]], n.vals[1])
loglik.ifn 
#-185518.7

loglik.lps2 = gaussianloglik_util(sample.covs[[2]], thetas.est.mono[[2]], n.vals[2])
loglik.lps2
# -133104.6

loglik.lps24 = gaussianloglik_util(sample.covs[[3]], thetas.est.mono[[3]], n.vals[3])
loglik.lps24
# -164492.9

loglik.unstim = gaussianloglik_util(sample.covs[[4]], thetas.est.mono[[4]], n.vals[4])
loglik.unstim 
# -210026.2

# Very small loglik...



# How many edges in common with jointGHS estimate?
load('Monocytes/data/monocytes_jointGHS_smalleps_largeAIC.RData')
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0,res.joint$theta[[1]] !=0)
#[,1]  [,2]
#[1,]  372 72018
#[2,]    0     0
tailoredGlasso::confusion.matrix(thetas.est.mono[[2]]!=0,res.joint$theta[[2]] !=0)
#[,1]  [,2]
#[1,]  393 71997
#[2,]    0     0
tailoredGlasso::confusion.matrix(thetas.est.mono[[3]]!=0,res.joint$theta[[3]]!=0)
#[,1]  [,2]
#[1,]  376 72014
#[2,]    0     0
tailoredGlasso::confusion.matrix(thetas.est.mono[[4]]!=0,res.joint$theta[[4]]!=0)
#[,1]  [,2]
#[1,]  377 72013
#[2,]    0     0
tailoredGlasso::confusion.matrix(thetas.est.mono[[1]]!=0 & thetas.est.mono[[2]]!=0 & thetas.est.mono[[3]]!=0 & thetas.est.mono[[4]]!=0,
                                 res.joint$theta[[1]]!=0 & res.joint$theta[[2]]!=0 & res.joint$theta[[3]]!=0 & res.joint$theta[[4]]!=0)

#[,1]  [,2]
#[1,]  254 72136
#[2,]    0     0


# Add information about whether a gene is cis or trans to top hotspot ---------------------------------

top_hotspot = 69757429 # On chromosome 12 (hotspot rs6581889)

# Find coordinates for all 
coords.all = annot_expr[,c('Gene_start','Gene_end')]
# Shortest dist to top hotspot region
dists.all = matrix(c(abs(top_hotspot-coords.all[,2]),abs(top_hotspot-coords.all[,1])),ncol=2, byrow = F)
dists.all = apply(dists.all, 1, min)

cis.which = dists.all<1e6 & annot_expr[,'Chr']==12
sum(cis.which) # 6 genes are CIS. (not all were used in analysis, as we used FDR 0.05. 5 were used.)
rownames(annot_expr)[which(cis.which)]
cis.genes.names = rownames(annot_expr)[which(cis.which)]
annot_expr[which(cis.which),]
annot_expr_full = annot_expr
# Check that annot table is right
annot_expr[which(annot_expr$cis_to_LYZ=='yes'),]

# Add information about whether a gene is controlled by the top hotspot --------------------------------

# FDR<0.05 only
df_hits_reduced = df_hits[df_hits$Gene_symbol %in% genes_id,]
sort(table(df_hits_reduced$SNP))

top_hotspot_name = 'rs6581889'
genes.controlled = df_hits_reduced[df_hits_reduced$SNP==top_hotspot_name,]
length(unique(genes.controlled$Gene_symbol)) # 299 (out of 381)

# Stratify by condition
genes.flagged = list()
for(i in 1:length(condition.names)){
  genes.flagged[[i]] = genes.controlled[genes.controlled$Condition==condition.names[i],'Gene_symbol']
  cat(condition.names[i], ': ', length(genes.flagged[[i]]), '\n')
}
# IFN-gamma :  294 
# LPS-2h :  88 
# LPS-24h :  16 
# Unstim :  215 

# Not as many when we stratify

# Which cis genes are controlled in each group?
cis.genes.names[which(cis.genes.names %in% genes.flagged[[1]])]
cis.genes.names[which(cis.genes.names %in% genes.flagged[[2]])]
cis.genes.names[which(cis.genes.names %in% genes.flagged[[3]])]
cis.genes.names[which(cis.genes.names %in% genes.flagged[[4]])]
# Only LYZ and YEATS4 in all. So these are the ones we will use.

# Only marking the cis genes controlled by the top hotspot
cis.genes.names = c('LYZ', 'YEATS4')



# What are the neighbours of the cis genes? ------------------------------------

# Look at edges from LYZ only
mean(genes_id == colnames(expr_all[[1]])) # Same order

ind.lyz = which(colnames(expr_all[[1]])=='LYZ')
# In IFN-GAMMA
edges.lyz.ifn = which(thetas.est.mono[[1]][ind.lyz,]!=0)
edges.lyz.ifn = edges.lyz.ifn[edges.lyz.ifn!=ind.lyz]
genes_id[edges.lyz.ifn]
annot_expr[genes_id[edges.lyz.ifn],]
# lps-2h
edges.lyz.lps2 = which(thetas.est.mono[[2]][ind.lyz,]!=0)
edges.lyz.lps2 = edges.lyz.lps2[edges.lyz.lps2!=ind.lyz]
genes_id[edges.lyz.lps2]
annot_expr[genes_id[edges.lyz.lps2],]
# lps-24h
edges.lyz.lps24 = which(thetas.est.mono[[3]][ind.lyz,]!=0)
edges.lyz.lps24 = edges.lyz.lps24[edges.lyz.lps24!=ind.lyz]
genes_id[edges.lyz.lps24]
annot_expr[genes_id[edges.lyz.lps24],]
# Unstim
edges.lyz.unstim = which(thetas.est.mono[[4]][ind.lyz,]!=0)
edges.lyz.unstim = edges.lyz.unstim[edges.lyz.unstim!=ind.lyz]
genes_id[edges.lyz.unstim]
annot_expr[genes_id[edges.lyz.unstim],]

# Very few..

# What neighbours are common across networks?
common.neighs = annot_expr[genes_id[which(genes_id %in% genes_id[edges.lyz.lps2] & genes_id %in% genes_id[edges.lyz.lps24] & genes_id %in% genes_id[edges.lyz.ifn] & genes_id %in% genes_id[edges.lyz.unstim])],]
common.neighs
#          Probe_Id Gene_start Gene_end Chr Symbol cis_to_LYZ
#AFMID ILMN_2095653   76203637 76203686  17  AFMID         no

# Far fewer than in jointGHS estimate

# Look at degree distribution ----------------------------------------------------


# For all cis genes
cis.names = cis.genes.names[which(cis.genes.names %in% genes_id)] # Only those used in sim
cis.names = c(cis.names, 'CREB1') # Also looking at CREB1
for(i in 1:length(cis.names)){
  g.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected', diag=F)
  df.deg.ifn = data.frame(degree=igraph::degree(g.ifn))
  rownames(df.deg.ifn) = genes_id
  deg.cis = df.deg.ifn[cis.names[i],]
  g.deg.ifn.noline = ggplot2::ggplot(df.deg.ifn, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.ifn = ggplot2::ggplot(df.deg.ifn, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected', diag=F)
  df.deg.lps2 = data.frame(degree=igraph::degree(g.lps2))
  rownames(df.deg.lps2) = genes_id
  deg.cis = df.deg.lps2[cis.names[i],]
  g.deg.lps2.noline =ggplot2::ggplot(df.deg.lps2, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.lps2 =ggplot2::ggplot(df.deg.lps2, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected', diag=F)
  df.deg.lps24 = data.frame(degree=igraph::degree(g.lps24))
  rownames(df.deg.lps24) = genes_id
  deg.cis = df.deg.lps24[cis.names[i],]
  g.deg.lps24.noline = ggplot2::ggplot(df.deg.lps24, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.lps24 = ggplot2::ggplot(df.deg.lps24, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected', diag=F)
  df.deg.unstim = data.frame(degree=igraph::degree(g.unstim ))
  rownames(df.deg.unstim) = genes_id
  deg.cis = df.deg.unstim[cis.names[i],]
  g.deg.unstim.noline =ggplot2::ggplot(df.deg.unstim , aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.unstim =ggplot2::ggplot(df.deg.unstim , aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  #gridExtra::grid.arrange(g.deg.ifn, g.deg.lps2, g.deg.lps24, g.deg.unstim)
}

# Also make a combined plot
degs.cis.ifn = df.deg.ifn[cis.names,] # LYZ, YEATS4, CREB1
degs.cis.lps2 = df.deg.lps2[cis.names,] # LYZ, YEATS4, CREB1
degs.cis.lps24 = df.deg.lps24[cis.names,] # LYZ, YEATS4, CREB1
degs.cis.unstim = df.deg.unstim[cis.names,] # LYZ, YEATS4, CREB1

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

df.test = data.frame(x=1:6,y=1:6,gene=rep(cis.names,2))
legend=get_legend(ggplot(df.test, aes(x=x,y=y,group=gene),palette=c('green','red','blue'))+geom_line(aes(color=gene))+theme(legend.position="right"))

pdf(paste0("Monocytes/GemBag/plots/degreedistr_hist_allcis_GemBag.pdf"), 8, 8)
gridExtra::grid.arrange(g.deg.ifn.noline+geom_vline(xintercept = degs.cis.ifn[1]+0.2, color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.ifn[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.ifn[3], color='red', alpha=0.7), 
                        g.deg.lps2.noline+geom_vline(xintercept = degs.cis.lps2[1]+0.2, color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.lps2[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.lps2[3], color='red', alpha=0.7),
                        legend,
                        g.deg.lps24.noline+geom_vline(xintercept = degs.cis.lps24[1], color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.lps24[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.lps24[3], color='red', alpha=0.7),
                        g.deg.unstim.noline+geom_vline(xintercept = degs.cis.unstim[1], color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.unstim[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.unstim[3], color='red', alpha=0.7),
                        ncol=3,widths=c(1,1,0.3))
dev.off()


# Also make density plot
g.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected', diag=F)
g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected', diag=F)
g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected', diag=F)
g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected', diag=F)

df.degree = data.frame(degree=c(igraph::degree(g.ifn), igraph::degree(g.lps2), igraph::degree(g.lps24), igraph::degree(g.unstim)), 
                       condition=factor(c(rep(condition.names[1], p), rep(condition.names[2], p),rep(condition.names[3], p),rep(condition.names[4], p))))

gg.dens = ggplot2::ggplot(df.degree, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))  +scale_x_continuous(limits = c(0,70))

pdf(paste0("Monocytes/GemBag/plots/degreedistr_density_GemBag.pdf"), 8, 8)
gg.dens+theme_bw()
dev.off()



# Not with CREB1, final plot, divided up

gg.dens.split = ggplot2::ggplot(df.degree, aes(degree))+geom_density(color='darkgray', fill='darkgray', alpha=0.2)+scale_fill_manual(values='darkgray')+
  ggplot2::theme(legend.position = 'none')  + scale_x_continuous(limits = c(0,70))+facet_wrap(~ condition)+
  geom_vline(data=dplyr::filter(df.degree, condition=="IFN-gamma"), aes(xintercept=degs.cis.ifn[1]), colour="darkorange",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="IFN-gamma"), aes(xintercept=degs.cis.ifn[2]), linetype='dashed', colour="darkolivegreen3",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-2h"), aes(xintercept=degs.cis.lps2[1]-0.2), colour="darkorange",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-2h"), aes(xintercept=degs.cis.lps2[2]), linetype='dashed', colour="darkolivegreen3",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-24h"), aes(xintercept=degs.cis.lps24[1]), colour="darkorange",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-24h"), aes(xintercept=degs.cis.lps24[2]), linetype='dashed', colour="darkolivegreen3",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="Unstim"), aes(xintercept=degs.cis.unstim[1]), colour="darkorange",alpha=0.7, size=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="Unstim"), aes(xintercept=degs.cis.unstim[2]), linetype='dashed', colour="darkolivegreen3",alpha=0.7, size=0.7)

df.test = data.frame(x=1:4,y=1:4,gene=rep(cis.names[1:2],2))
legend=get_legend(ggplot(df.test, aes(x=x,y=y,group=gene))+geom_line(aes(color=gene,linetype=gene), alpha=0.7, size=0.7)+theme(legend.position="right")+
                    scale_colour_manual(values=c("darkorange",'darkolivegreen3'))+ theme(text = element_text(size = 15)))


pdf("Monocytes/GemBag/plots/degreedistr_density_cisgenes_GemBag.pdf", 10, 8)
gridExtra::grid.arrange(gg.dens.split+theme_bw()+ggplot2::theme(text = element_text(size = 15)), legend ,ncol=2,widths=c(1,0.15))
dev.off()


# What are the hubs? ------------------------------------------------------------

# Defining hubs as having top 10% highest degree

df.degree.ifn = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.ifn)))], rev(sort(igraph::degree(g.ifn)))))
names(df.degree.ifn)=c('Gene', 'Degree')
top.ifn = as.numeric(df.degree.ifn$Degree)>=quantile(as.numeric(df.degree.ifn$Degree), 0.9)
df.degree.ifn[top.ifn,]
length(df.degree.ifn[top.ifn,1]) # 42
df.degree.lps2 = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.lps2)))], rev(sort(igraph::degree(g.lps2)))))
names(df.degree.lps2)=c('Gene', 'Degree')
top.lps2 = as.numeric(df.degree.lps2$Degree)>=quantile(as.numeric(df.degree.lps2$Degree), 0.9)
df.degree.lps2[top.lps2,]
length(df.degree.lps2[top.lps2,1]) # 45
df.degree.lps24 = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.lps24)))], rev(sort(igraph::degree(g.lps24)))))
names(df.degree.lps24)=c('Gene', 'Degree')
top.lps24 = as.numeric(df.degree.lps24$Degree)>=quantile(as.numeric(df.degree.lps24$Degree), 0.9)
df.degree.lps24[top.lps24,]
length(df.degree.lps24[top.lps24,1]) # 43
df.degree.unstim = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.unstim)))], rev(sort(igraph::degree(g.unstim)))))
names(df.degree.unstim) =c('Gene', 'Degree')
top.unstim = as.numeric(df.degree.unstim$Degree)>=quantile(as.numeric(df.degree.unstim$Degree), 0.9)
df.degree.unstim[top.unstim,]
length(df.degree.unstim[top.unstim,1]) # 42

# How much do the hubs agree?
sum(df.degree.ifn[top.ifn,1]%in% df.degree.lps2[top.ifn,1])
# 38
sum(df.degree.ifn[top.lps2,1]%in% df.degree.lps24[top.lps24,1])
# 39
sum(df.degree.ifn[top.lps2,1]%in% df.degree.unstim[top.unstim,1])
# 40
sum(df.degree.lps2[top.ifn,1]%in% df.degree.lps24[top.lps24,1])
# 37
sum(df.degree.lps2[top.ifn,1]%in% df.degree.unstim[top.unstim,1])
# 38
sum(df.degree.lps24[top.lps24,1]%in% df.degree.unstim[top.unstim,1])
# 37

# Almost agree on all hubs

# Is YEATS4 a hub in all? 
'YEATS4' %in% df.degree.ifn[top.ifn,1]
'YEATS4' %in% df.degree.lps2[top.lps2,1]
'YEATS4' %in% df.degree.lps24[top.lps24,1]
'YEATS4' %in% df.degree.unstim[top.unstim,1]
# None!

# Is LYZ a hub in all? 
'LYZ' %in% df.degree.ifn[top.ifn,1]
'LYZ' %in% df.degree.lps2[top.lps2,1]
'LYZ' %in% df.degree.lps24[top.lps24,1]
'LYZ' %in% df.degree.unstim[top.unstim,1]
# None!

# Their degrees
cbind(condition.names, c(sum(thetas.est.mono[[1]][which(colnames(expr_all[[1]])=='LYZ'),]!=0)-1,sum(thetas.est.mono[[2]][which(colnames(expr_all[[2]])=='LYZ'),]!=0)-1,
                         sum(thetas.est.mono[[3]][which(colnames(expr_all[[3]])=='LYZ'),]!=0)-1, sum(thetas.est.mono[[4]][which(colnames(expr_all[[4]])=='LYZ'),]!=0)-1))
#condition.names    
#[1,] "IFN-gamma"     "2"
#[2,] "LPS-2h"        "2"
#[3,] "LPS-24h"       "2"
#[4,] "Unstim"        "1"

cbind(condition.names, c(sum(thetas.est.mono[[1]][which(colnames(expr_all[[1]])=='YEATS4'),]!=0)-1,sum(thetas.est.mono[[2]][which(colnames(expr_all[[2]])=='YEATS4'),]!=0)-1,
                         sum(thetas.est.mono[[3]][which(colnames(expr_all[[3]])=='YEATS4'),]!=0)-1, sum(thetas.est.mono[[4]][which(colnames(expr_all[[4]])=='YEATS4'),]!=0)-1))
#     condition.names    
#[1,] "IFN-gamma"     "0"
#[2,] "LPS-2h"        "0"
#[3,] "LPS-24h"       "0"
#[4,] "Unstim"        "0"

# Are any genes a hub in all conditions?
hub.in.all = genes_id[which(genes_id %in% df.degree.ifn[top.ifn,1] &  genes_id %in% df.degree.lps2[top.lps2,1] &  genes_id %in%  df.degree.lps24[top.lps24,1] &  genes_id %in%  df.degree.unstim[top.unstim,1])]
hub.in.all
#[1] "AFMID"        "BMS1P5"       "C3orf34"      "C8orf45"      "CATSPER2"     "CCBE1"        "CDKN2AIPNL"   "CHRNA5"       "CREB1"        "DDX51"        "DUSP19"      
#[12] "EID2B"        "FAM119A"      "FAM175A"      "GGA1"         "GSDM1"        "HCG2P7"       "HSPC268"      "KCNH6"        "LOC100128288" "LOC100132391" "LOC100190938"
#[23] "LOC441087"    "LOC729090"    "LOC730313"    "LRAP"         "PRIM2"        "QRFPR"        "RRP7B"        "SNHG10"       "TP53BP2"      "TRIM16L"      "ZMAT3"       
#[34] "ZNF394"   
# Almost all!

# Save as RData object 
top.hubs.ifn = df.degree.ifn[top.ifn,1]
top.hubs.lps2 = df.degree.lps2[top.lps2,1]
top.hubs.lps24 = df.degree.lps24[top.lps24,1]
top.hubs.unstim = df.degree.unstim[top.unstim,1]
save(top.hubs.ifn, top.hubs.lps2, top.hubs.lps24, top.hubs.unstim, file='Monocytes/GemBag/data/top_hubs_GemBag.RData')



# What are their neighbours? Looking at most influential hub

# IFN-GAMMA
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[1]),]!=0)],]
# LYZ a neighbour

# lps2
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[1]),]!=0)],]

# lps24
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[1]),]!=0)],]

# Unstim
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[1]),]!=0)],]


# Second most influential

# IFN-GAMMA
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[2]),]!=0)],]

# lps2
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[2]),]!=0)],]

# lps24
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[2]),]!=0)],]

# Unstim
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[2]),]!=0)],]

# Third most influential

# IFN-GAMMA
df.degree.ifn$Gene[3] # AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[3]),]!=0)],]

# lps2
df.degree.lps2$Gene[3]  #AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[3]),]!=0)],]

# lps24
df.degree.lps24$Gene[3] # AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[3]),]!=0)],]

# Unstim 
df.degree.unstim$Gene[3] # PHRF1
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[3]),]!=0)],]

# So LYZ is only a neighbour of one of the top three hubs in IFNg...


# Are all the cis genes a neighbour to the top hubs?


for(i in 1:length(cis.names)){
  neigh.ifn = genes_id[which(thetas.est.mono[[1]][genes_id==cis.names[i],which(genes_id %in% df.degree.ifn[top.ifn,1])]!=0)]
  cat('IFN: ', cis.names[i], ': ', neigh.ifn, '\n')
  neigh.lps2 = genes_id[which(thetas.est.mono[[1]][genes_id==cis.names[i],which(genes_id %in% df.degree.lps2[top.lps2,1])]!=0)]
  cat('LPS2h: ', cis.names[i], ': ', neigh.lps2, '\n')
  neigh.lps24 = genes_id[which(thetas.est.mono[[1]][genes_id==cis.names[i],which(genes_id %in% df.degree.lps24[top.lps24,1])]!=0)]
  cat('LPS24h: ', cis.names[i], ': ', neigh.lps24, '\n')
  neigh.unstim = genes_id[which(thetas.est.mono[[1]][genes_id==cis.names[i],which(genes_id %in% df.degree.unstim[top.unstim,1])]!=0)]
  cat('LPS2h: ', cis.names[i], ': ', neigh.unstim, '\n')
}

# LYZ is, but not YEATS4 


# What are their degrees?
degree.list.ifn = list()
degree.list.lps2 = list()
degree.list.lps24 = list()
degree.list.unstim = list()
for(i in 1:(length(cis.names)-1)){
  neighs.all.ifn = genes_id[which(thetas.est.mono[[1]][genes_id==cis.names[i],]!=0)]
  neighs.all.ifn = neighs.all.ifn[-which(neighs.all.ifn==cis.names[i])]
  degs = df.deg.ifn[neighs.all.ifn,]
  df.save.ifn = data.frame(gene=neighs.all.ifn[rev(order(degs))], degree= df.deg.ifn[neighs.all.ifn,][rev(order(degs))])
  neighs.all.lps2 = genes_id[which(thetas.est.mono[[2]][genes_id==cis.names[i],]!=0)]
  neighs.all.lps2 = neighs.all.lps2[-which(neighs.all.lps2==cis.names[i])]
  degs = df.deg.lps2[neighs.all.lps2,]
  df.save.lps2 = data.frame(gene = neighs.all.lps2[rev(order(degs))], degree = df.deg.lps2[neighs.all.lps2,][rev(order(degs))])
  neighs.all.lps24 = genes_id[which(thetas.est.mono[[3]][genes_id==cis.names[i],]!=0)]
  neighs.all.lps24 = neighs.all.lps24[-which(neighs.all.lps24==cis.names[i])]
  degs = df.deg.lps24[neighs.all.lps24,]
  df.save.lps24 = data.frame(gene = neighs.all.lps24[rev(order(degs))], degree = df.deg.lps24[neighs.all.lps24,][rev(order(degs))])
  neighs.all.unstim = genes_id[which(thetas.est.mono[[4]][genes_id==cis.names[i],]!=0)]
  neighs.all.unstim = neighs.all.unstim[-which(neighs.all.unstim==cis.names[i])]
  degs = df.deg.unstim[neighs.all.unstim,]
  df.save.unstim = data.frame(gene = neighs.all.unstim[rev(order(degs))], degree = df.deg.unstim[neighs.all.unstim,][rev(order(degs))])
  degree.list.ifn[[i]] = df.save.ifn
  degree.list.lps2[[i]] = df.save.lps2
  degree.list.lps24[[i]] = df.save.lps24
  degree.list.unstim[[i]] = df.save.unstim
}


# LYZ (degree of hub neigbours)
dim.neighs = max(unlist(lapply(list(degree.list.ifn[[1]],degree.list.lps2[[1]],degree.list.lps24[[1]],degree.list.unstim[[1]]), FUN = function(s) nrow(s))))
for(i in 1:dim.neighs){
  if(is.na(degree.list.ifn[[1]][i,1])) cat(' & &&')
  else cat(degree.list.ifn[[1]][i,1], ' & ', degree.list.ifn[[1]][i,2], '&&')
  if(is.na(degree.list.lps2[[1]][i,1])) cat(' & &&')
  else cat(degree.list.lps2[[1]][i,1], ' & ', degree.list.lps2[[1]][i,2], '&&')
  if(is.na(degree.list.lps24[[1]][i,1])) cat(' & &&')
  else cat(degree.list.lps24[[1]][i,1], ' & ', degree.list.lps24[[1]][i,2], '&&')
  if(is.na(degree.list.unstim[[1]][i,1])) cat(' & \\\\ \n')
  else cat(degree.list.unstim[[1]][i,1], ' & ', degree.list.unstim[[1]][i,2], ' \\\\ \n')
}
#HCG2P7  &  60 &&AFMID  &  12 &&AFMID  &  15 &&AFMID  &  12  \\ 
#AFMID  &  12 &&TP53BP2  &  9 &&CCL20  &  8 && & \\ 

# YEATS4
dim.neighs = max(unlist(lapply(list(degree.list.ifn[[2]],degree.list.lps2[[2]],degree.list.lps24[[2]],degree.list.unstim[[2]]), FUN = function(s) nrow(s))))
for(i in 1:dim.neighs){
  if(is.na(degree.list.ifn[[2]][i,1])) cat(' & &&')
  else cat(degree.list.ifn[[2]][i,1], ' & ', degree.list.ifn[[2]][i,2], '&&')
  if(is.na(degree.list.lps2[[2]][i,1])) cat(' & &&')
  else cat(degree.list.lps2[[2]][i,1], ' & ', degree.list.lps2[[2]][i,2], '&&')
  if(is.na(degree.list.lps24[[2]][i,1])) cat(' & &&')
  else cat(degree.list.lps24[[2]][i,1], ' & ', degree.list.lps24[[2]][i,2], '&&')
  if(is.na(degree.list.unstim[[2]][i,1])) cat(' & \\\\ \n')
  else cat(degree.list.unstim[[2]][i,1], ' & ', degree.list.unstim[[2]][i,2], '\\\\\n')
}
# None


# Which are controlled by the top hotspot? 

# LYZ
degree.list.ifn[[1]]$gene[which(degree.list.ifn[[1]]$gene %in% genes.flagged[[1]])]
degree.list.lps2[[1]]$gene[which(degree.list.lps2[[1]]$gene %in% genes.flagged[[2]])]
degree.list.lps24[[1]]$gene[which(degree.list.lps24[[1]]$gene %in% genes.flagged[[3]])]
degree.list.unstim[[1]]$gene[which(degree.list.unstim[[1]]$gene %in% genes.flagged[[4]])]
# AFMID in all conditions, plus one more for ifn and lps2

# YEATS4
degree.list.ifn[[2]]$gene[which(degree.list.ifn[[2]]$gene %in% genes.flagged[[1]])]
degree.list.lps2[[2]]$gene[which(degree.list.lps2[[2]]$gene %in% genes.flagged[[2]])]
degree.list.lps24[[2]]$gene[which(degree.list.lps24[[2]]$gene %in% genes.flagged[[3]])]
degree.list.unstim[[2]]$gene[which(degree.list.unstim[[2]]$gene %in% genes.flagged[[4]])]
# none

# Now the same, but for all nodes with degree above xx

top.frac=0.9
df.deg.ifn.ordered = data.frame(genes = rownames(df.deg.ifn)[rev(order(df.deg.ifn))], degree = df.deg.ifn$degree[rev(order(df.deg.ifn))])
df.deg.ifn.ordered.top = df.deg.ifn.ordered[which(df.deg.ifn.ordered$degree>quantile(df.deg.ifn.ordered$degree,top.frac)),]
df.deg.lps2.ordered = data.frame(genes = rownames(df.deg.lps2)[rev(order(df.deg.lps2))], degree = df.deg.lps2$degree[rev(order(df.deg.lps2))])
df.deg.lps2.ordered.top = df.deg.lps2.ordered[which(df.deg.lps2.ordered$degree>quantile(df.deg.lps2.ordered$degree,top.frac)),]
df.deg.lps24.ordered = data.frame(genes = rownames(df.deg.lps24)[rev(order(df.deg.lps24))], degree = df.deg.lps24$degree[rev(order(df.deg.lps24))])
df.deg.lps24.ordered.top = df.deg.lps24.ordered[which(df.deg.lps24.ordered$degree>quantile(df.deg.lps24.ordered$degree,top.frac)),]
df.deg.unstim.ordered = data.frame(genes = rownames(df.deg.unstim)[rev(order(df.deg.unstim))], degree = df.deg.unstim$degree[rev(order(df.deg.unstim))])
df.deg.unstim.ordered.top = df.deg.unstim.ordered[which(df.deg.unstim.ordered$degree>quantile(df.deg.unstim.ordered$degree,top.frac)),]


# Print for table (THIS WE MIGHT USE!) Table of top hubs
dim.neighs.top = max(c(length(df.deg.ifn.ordered.top$degree), length(df.deg.lps2.ordered.top$degree), length(df.deg.lps24.ordered.top$degree), length(df.deg.unstim.ordered.top$degree)))
for(i in 1:dim.neighs.top){
  if(is.na(df.deg.ifn.ordered.top[i,1])) cat(' & &&')
  else cat(df.deg.ifn.ordered.top[i,1], ' & ', df.deg.ifn.ordered.top[i,2], '&&')
  if(is.na(df.deg.lps2.ordered.top[i,1])) cat(' & &&')
  else cat(df.deg.lps2.ordered.top[i,1], ' & ', df.deg.lps2.ordered.top[i,2], '&&')
  if(is.na(df.deg.lps24.ordered.top[i,1])) cat(' & &&')
  else cat(df.deg.lps24.ordered.top[i,1], ' & ', df.deg.lps24.ordered.top[i,2], '&&')
  if(is.na(df.deg.unstim.ordered.top[i,1])) cat(' & \\\\ \n')
  else cat(df.deg.unstim.ordered.top[i,1], ' & ', df.deg.unstim.ordered.top[i,2], ' \\\\ \n')
}

# Which are common to all?
genes_id[which(genes_id %in% df.deg.ifn.ordered.top$genes &  genes_id %in% df.deg.lps2.ordered.top$genes & genes_id %in% df.deg.lps24.ordered.top$genes & genes_id %in% df.deg.unstim.ordered.top$genes)]
# unique to IFN
genes_id[which(genes_id %in% df.deg.ifn.ordered.top$genes &  ! genes_id %in% df.deg.lps2.ordered.top$genes & ! genes_id %in% df.deg.lps24.ordered.top$genes & ! genes_id %in% df.deg.unstim.ordered.top$genes)]
# unique to LPS2
genes_id[which(! genes_id %in% df.deg.ifn.ordered.top$genes &  genes_id %in% df.deg.lps2.ordered.top$genes & ! genes_id %in% df.deg.lps24.ordered.top$genes & !genes_id %in% df.deg.unstim.ordered.top$genes)]
# unique to LPS24
genes_id[which(! genes_id %in% df.deg.ifn.ordered.top$genes &  !genes_id %in% df.deg.lps2.ordered.top$genes &  genes_id %in% df.deg.lps24.ordered.top$genes & !genes_id %in% df.deg.unstim.ordered.top$genes)]
# unique to unstim
genes_id[which(! genes_id %in% df.deg.ifn.ordered.top$genes & ! genes_id %in% df.deg.lps2.ordered.top$genes & ! genes_id %in% df.deg.lps24.ordered.top$genes & genes_id %in% df.deg.unstim.ordered.top$genes)]
# Most hubs are in common. How many?
hubs.common=genes_id[which(genes_id %in% df.deg.ifn.ordered.top$genes &  genes_id %in% df.deg.lps2.ordered.top$genes & genes_id %in% df.deg.lps24.ordered.top$genes & genes_id %in% df.deg.unstim.ordered.top$genes)]
length(hubs.common) # 24
# Out of how many hubs?
length(df.deg.ifn.ordered.top$degree) # 34
length(df.deg.lps2.ordered.top$degree) # 36
length(df.deg.lps24.ordered.top$degree) # 34
length(df.deg.unstim.ordered.top$degree) # 29
# So most...


# What are their hotspot statuses?
df.deg.ifn.ordered.top$genes[which(df.deg.ifn.ordered.top$genes %in% genes.flagged[[1]])]
df.deg.lps2.ordered.top$genes[which(df.deg.lps2.ordered.top$genes %in% genes.flagged[[2]])]
df.deg.lps24.ordered.top$genes[which(df.deg.lps24.ordered.top$genes %in% genes.flagged[[3]])]
df.deg.unstim.ordered.top$genes[which(df.deg.unstim.ordered.top$genes %in% genes.flagged[[4]])]


# Evaluate significance of density of top hotspot mediated subnetworks ------------------------

# IFN
spars.ifn = tailoredGlasso::sparsity(thetas.est.mono[[1]][which(genes_id %in% genes.flagged[[1]]),which(genes_id %in% genes.flagged[[1]])])
spars.ifn 
# 0.008567249
# LPS2
spars.lps2 = tailoredGlasso::sparsity(thetas.est.mono[[2]][which(genes_id %in% genes.flagged[[2]]),which(genes_id %in% genes.flagged[[2]])])
spars.lps2
# 0.0514629
# LPS24
spars.lps24 = tailoredGlasso::sparsity(thetas.est.mono[[3]][which(genes_id %in% genes.flagged[[3]]),which(genes_id %in% genes.flagged[[3]])])
spars.lps24
# 0.2083333
# Unstim
spars.unstim = tailoredGlasso::sparsity(thetas.est.mono[[4]][which(genes_id %in% genes.flagged[[4]]),which(genes_id %in% genes.flagged[[4]])])
spars.unstim
# 0.0154749

# Sampling to assess significance
n.sample = 10000
spars.vals.ifn = c()
spars.vals.lps2 = c()
spars.vals.lps24 = c()
spars.vals.unstim = c()
set.seed(123)
for(i in 1:n.sample){
  # Same as many genes as are flagged
  genes.sampled.ifn = sample(1:p, length(genes.flagged[[1]]))
  spars.vals.ifn[i] = tailoredGlasso::sparsity(thetas.est.mono[[1]][genes.sampled.ifn, genes.sampled.ifn])
  genes.sampled.lps2 = sample(1:p, length(genes.flagged[[2]]))
  spars.vals.lps2[i] = tailoredGlasso::sparsity(thetas.est.mono[[2]][genes.sampled.lps2, genes.sampled.lps2])
  genes.sampled.lps24 = sample(1:p, length(genes.flagged[[3]]))
  spars.vals.lps24[i] = tailoredGlasso::sparsity(thetas.est.mono[[3]][genes.sampled.lps24, genes.sampled.lps24])
  genes.sampled.unstim = sample(1:p, length(genes.flagged[[4]]))
  spars.vals.unstim[i] = tailoredGlasso::sparsity(thetas.est.mono[[4]][genes.sampled.unstim, genes.sampled.unstim])
}



(sum(spars.vals.ifn>spars.ifn)+1)/(n.sample+1)
# 9.999e-05
(sum(spars.vals.lps2>spars.lps2)+1)/(n.sample+1)
# 9.999e-05
(sum(spars.vals.lps24>spars.lps24)+1)/(n.sample+1)
# 9.999e-05
(sum(spars.vals.unstim>spars.unstim)+1)/(n.sample+1)
# 9.999e-05


# Make upset plots


dat.edges = data.frame(IFN =thetas.est.mono[[1]][upper.tri(thetas.est.mono[[1]])]!=0, LPS2 =thetas.est.mono[[2]][upper.tri(thetas.est.mono[[2]])]!=0, 
                       LPS24 =thetas.est.mono[[3]][upper.tri(thetas.est.mono[[3]])]!=0, Unstim =thetas.est.mono[[4]][upper.tri(thetas.est.mono[[4]])]!=0)
names(dat.edges) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

plot.upset = upset(dat.edges, names(dat.edges)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=50000, sort_sets=F)
plot.upset.wrapped = upset(dat.edges, names(dat.edges)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=50000, sort_sets=F, wrap=T, 
                           themes=upset_default_themes(text=element_text(size=15)))


# Change order of intersections:

inter.order.grouped = list()
inter.order.grouped[[1]] = c(combn(condition.names,4))
comb.3 = combn(condition.names,3)
order.3 = c(4,2,1,3) # Order manually according to size
for(k in 1:ncol(comb.3)){
  #inter.order.grouped[[k+1]] = c(comb.3[,k])
  inter.order.grouped[[k+1]] = c(comb.3[,order.3[k]])
}
comb.2 = combn(condition.names,2)
#order.2 = c(4,3,1, 2, 6,5)
order.2 = c(3,4,2,1,6,5)
for(k in 1:ncol(comb.2)){
  inter.order.grouped[[order.2[k]+ncol(comb.3)+1]] = c(comb.2[,k])
}
order.1 = c(1, 3, 2, 4)
for(k in 1:4){
  inter.order.grouped[[order.1[k]+ncol(comb.3)+ncol(comb.2)+1]] = condition.names[k]
}

#inter.order.grouped = list(c("IFN-gamma","LPS-2h", "LPS-24h", "Unstim"),  c("LPS-2h", "LPS-24h", "Unstim" ), c("IFN-gamma", "LPS-2h", "Unstim"),
 #                          c("IFN-gamma", "LPS-2h", "LPS-24h"), 
#                           
#                           c( "LPS-2h", "LPS-24h"), c("IFN-gamma", "Unstim"), c("IFN-gamma"), c("LPS-24h" ),
#                           c("LPS-2h"), c("Unstim"), 
#                           )

plot.upset.wrapped.ordered = upset(dat.edges, names(dat.edges)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=50000, sort_intersections=F,sort_sets=F, wrap=T, 
                                   themes=upset_default_themes(text=element_text(size=15)), intersections=inter.order.grouped)


pdf("Monocytes/GemBag/plots/intersection_GemBag_ordered.pdf", 14, 7)
plot.upset.wrapped.ordered
dev.off()



