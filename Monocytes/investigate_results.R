rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_and_bcells_LYZ_region.RData")
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")
load('Monocytes/data/monocytes_jointGHS_smalleps_largeAIC.RData')
load('Monocytes/data/res_joint_single_samespars.RData')

p = dim(expr_all[[1]])[2]
n.vals = unlist(lapply(expr_all, nrow))
condition.names = c('IFN-gamma', 'LPS-2h', 'LPS-24h', 'Unstim')
# Get precision matrices from joint estimates
thetas.est.mono = lapply(res.joint$theta, cov2cor)
for(k in 1:length(thetas.est.mono)){
  thetas.est.mono[[k]][which(abs(thetas.est.mono[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Get precision matrices from single estimates
thetas.est.mono.single = lapply(res.single, FUN = function(s) cov2cor(s$theta))
for(k in 1:length(thetas.est.mono.single)){
  thetas.est.mono.single[[k]][which(abs(thetas.est.mono.single[[k]]) < 1e-5, arr.ind = T)] = 0
}

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

save(genes.flagged, file='Monocytes/data/genes_flagged.RData')


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

# Only in lps-24h is another cis gene a neighbour to LYZ (YEATS4)

# What neighbours are common across networks?
common.neighs = annot_expr[genes_id[which(genes_id %in% genes_id[edges.lyz.lps2] & genes_id %in% genes_id[edges.lyz.lps24] & genes_id %in% genes_id[edges.lyz.ifn] & genes_id %in% genes_id[edges.lyz.unstim])],]
common.neighs
#                 Probe_Id Gene_start Gene_end Chr       Symbol cis_to_LYZ
#AFMID        ILMN_2095653   76203637 76203686  17        AFMID         no
#KLHL28       ILMN_3251605   46326235 46326284  14       KLHL28         no
#LOC100128098 ILMN_3261439   17451419 17451468  10 LOC100128098         no
#MAFF         ILMN_2322375   38612067 38612116  22         MAFF         no
#SNRNP48      ILMN_3237516    7611650  7611699   6      SNRNP48         no


# Look at degree distribution ----------------------------------------------------


# For all cis genes
cis.names = cis.genes.names[which(cis.genes.names %in% genes_id)] # Only those used in sim
cis.names = c(cis.names, 'CREB1') # Also looking at CREB1
for(i in 1:length(cis.names)){
  g.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected')
  df.deg.ifn = data.frame(degree=igraph::degree(g.ifn))
  rownames(df.deg.ifn) = genes_id
  deg.cis = df.deg.ifn[cis.names[i],]
  g.deg.ifn = ggplot2::ggplot(df.deg.ifn, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected')
  df.deg.lps2 = data.frame(degree=igraph::degree(g.lps2))
  rownames(df.deg.lps2) = genes_id
  deg.cis = df.deg.lps2[cis.names[i],]
  g.deg.lps2 =ggplot2::ggplot(df.deg.lps2, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected')
  df.deg.lps24 = data.frame(degree=igraph::degree(g.lps24))
  rownames(df.deg.lps24) = genes_id
  deg.cis = df.deg.lps24[cis.names[i],]
  g.deg.lps24 = ggplot2::ggplot(df.deg.lps24, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected')
  df.deg.unstim = data.frame(degree=igraph::degree(g.unstim ))
  rownames(df.deg.unstim) = genes_id
  deg.cis = df.deg.unstim[cis.names[i],]
  g.deg.unstim =ggplot2::ggplot(df.deg.unstim , aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_", cis.names[i],"_combined_jointGHS.pdf"), 8, 8)
  gridExtra::grid.arrange(g.deg.ifn, g.deg.lps2, g.deg.lps24, g.deg.unstim)
  dev.off()
}

# Also make density plot
g.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected')
g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected')
g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected')
g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected')

df.degree = data.frame(degree=c(igraph::degree(g.ifn), igraph::degree(g.lps2), igraph::degree(g.lps24), igraph::degree(g.unstim)), 
                       condition=factor(c(rep(condition.names[1], p), rep(condition.names[2], p),rep(condition.names[3], p),rep(condition.names[4], p))))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_jointGHS.pdf"), 8, 8)
ggplot2::ggplot(df.degree, aes(degree, group=condition, colour=condition))+geom_density()+
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) 
dev.off()

# Only considering degree of genes controlled by top hotspot

# Get indices
ind.flagget.ifn = which(genes_id %in% genes.flagged[[1]])
ind.flagget.lps2 = which(genes_id %in% genes.flagged[[2]])
ind.flagget.lps24 = which(genes_id %in% genes.flagged[[3]])
ind.flagget.unstim = which(genes_id %in% genes.flagged[[4]])

df.degree = data.frame(degree=c(igraph::degree(g.ifn)[ind.flagget.ifn], igraph::degree(g.lps2)[ind.flagget.lps2], igraph::degree(g.lps24)[ind.flagget.lps24], igraph::degree(g.unstim)[ind.flagget.unstim]), 
                       condition=factor(c(rep(condition.names[1], length(ind.flagget.ifn)), rep(condition.names[2], length(ind.flagget.lps2)),rep(condition.names[3], length(ind.flagget.lps24)),rep(condition.names[4], length(ind.flagget.unstim)))))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_onlyflagged_jointGHS.pdf"), 8, 8)
ggplot2::ggplot(df.degree, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) 
dev.off()



# Look at degree distribution of single network with same sparsity ----------------------------------------------------


# For all cis genes
for(i in 1:length(cis.names)){
  g.ifn.single = igraph::graph.adjacency(thetas.est.mono.single[[1]]!=0, mode='undirected')
  df.deg.ifn.single = data.frame(degree=igraph::degree(g.ifn.single))
  rownames(df.deg.ifn.single) = genes_id
  deg.cis.single = df.deg.ifn.single[cis.names[i],]
  g.deg.ifn.single = ggplot2::ggplot(df.deg.ifn.single, aes(x=degree))+geom_histogram(fill='lightgrey', bins=20)+geom_vline(xintercept = deg.cis.single, color='red')+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps2.single = igraph::graph.adjacency(thetas.est.mono.single[[2]]!=0, mode='undirected')
  df.deg.lps2.single = data.frame(degree=igraph::degree(g.lps2.single))
  rownames(df.deg.lps2.single) = genes_id
  deg.cis.single = df.deg.lps2.single[cis.names[i],]
  g.deg.lps2.single =ggplot2::ggplot(df.deg.lps2.single, aes(x=degree))+geom_histogram(fill='lightgrey', bins=20)+geom_vline(xintercept = deg.cis.single, color='red')+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps24.single = igraph::graph.adjacency(thetas.est.mono.single[[3]]!=0, mode='undirected')
  df.deg.lps24.single = data.frame(degree=igraph::degree(g.lps24.single))
  rownames(df.deg.lps24.single) = genes_id
  deg.cis.single = df.deg.lps24.single[cis.names[i],]
  g.deg.lps24.single = ggplot2::ggplot(df.deg.lps24.single, aes(x=degree))+geom_histogram(fill='lightgrey', bins=20)+geom_vline(xintercept = deg.cis.single, color='red')+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.unstim.single = igraph::graph.adjacency(thetas.est.mono.single[[4]]!=0, mode='undirected')
  df.deg.unstim.single = data.frame(degree=igraph::degree(g.unstim.single))
  rownames(df.deg.unstim.single) = genes_id
  deg.cis.single = df.deg.unstim.single[cis.names[i],]
  g.deg.unstim.single =ggplot2::ggplot(df.deg.unstim.single , aes(x=degree))+geom_histogram(fill='lightgrey', bins=20)+geom_vline(xintercept = deg.cis.single, color='red')+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_", cis.names[i],"_single_combined_jointGHS.pdf"), 8, 8)
  gridExtra::grid.arrange(g.deg.ifn.single, g.deg.lps2.single, g.deg.lps24.single, g.deg.unstim.single)
  dev.off()
}

# Also make density plot
g.ifn.single = igraph::graph.adjacency(thetas.est.mono.single[[1]]!=0, mode='undirected')
g.lps2.single = igraph::graph.adjacency(thetas.est.mono.single[[2]]!=0, mode='undirected')
g.lps24.single = igraph::graph.adjacency(thetas.est.mono.single[[3]]!=0, mode='undirected')
g.unstim.single = igraph::graph.adjacency(thetas.est.mono.single[[4]]!=0, mode='undirected')

df.degree.single = data.frame(degree=c(igraph::degree(g.ifn.single), igraph::degree(g.lps2.single), igraph::degree(g.lps24.single), igraph::degree(g.unstim.single)), 
                       condition=factor(c(rep(condition.names[1], p), rep(condition.names[2], p),rep(condition.names[3], p),rep(condition.names[4], p))))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_single_jointGHS.pdf"), 8, 8)
ggplot2::ggplot(df.degree.single, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) 
dev.off()

# Only considering degree of genes controlled by top hotspot

# Get indices
ind.flagget.ifn = which(genes_id %in% genes.flagged[[1]])
ind.flagget.lps2 = which(genes_id %in% genes.flagged[[2]])
ind.flagget.lps24 = which(genes_id %in% genes.flagged[[3]])
ind.flagget.unstim = which(genes_id %in% genes.flagged[[4]])

df.degree.single = data.frame(degree=c(igraph::degree(g.ifn.single)[ind.flagget.ifn], igraph::degree(g.lps2.single)[ind.flagget.lps2], igraph::degree(g.lps24.single)[ind.flagget.lps24], igraph::degree(g.unstim.single)[ind.flagget.unstim]), 
                       condition=factor(c(rep(condition.names[1], length(ind.flagget.ifn)), rep(condition.names[2], length(ind.flagget.lps2)),rep(condition.names[3], length(ind.flagget.lps24)),rep(condition.names[4], length(ind.flagget.unstim)))))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_single_onlyflagged_jointGHS.pdf"), 8, 8)
ggplot2::ggplot(df.degree.single, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) 
dev.off()


# Degree of cis genes in single networks (from joint procedure) -----------------------------------------------------

g.ifn.single = igraph::graph.adjacency(abs(cov2cor(res.joint$theta_single[[1]]))>1e-5, mode='undirected')
df.deg.ifn.single = data.frame(degree=igraph::degree(g.ifn.single))
hist(df.deg.ifn.single$degree, breaks=20)
abline(v=df.deg.ifn.single$degree[which(genes_id=='LYZ')])
abline(v=df.deg.ifn.single$degree[which(genes_id=='YEATS4')], col='red')
abline(v=df.deg.ifn.single$degree[which(genes_id=='CREB1')], col='blue')

g.lps.single = igraph::graph.adjacency(abs(cov2cor(res.joint$theta_single[[2]]))>1e-5, mode='undirected')
df.deg.lps.single = data.frame(degree=igraph::degree(g.lps.single))
hist(df.deg.lps.single$degree, breaks=20)
abline(v=df.deg.lps.single$degree[which(genes_id=='LYZ')])
abline(v=df.deg.lps.single$degree[which(genes_id=='YEATS4')], col='red')
abline(v=df.deg.lps.single$degree[which(genes_id=='CREB1')], col='blue')

g.lps24.single = igraph::graph.adjacency(abs(cov2cor(res.joint$theta_single[[3]]))>1e-5, mode='undirected')
df.deg.lps24.single = data.frame(degree=igraph::degree(g.lps24.single))
hist(df.deg.lps24.single$degree, breaks=20)
abline(v=df.deg.lps24.single$degree[which(genes_id=='LYZ')])
abline(v=df.deg.lps24.single$degree[which(genes_id=='YEATS4')], col='red')
abline(v=df.deg.lps24.single$degree[which(genes_id=='CREB1')], col='blue')

g.unstim.single = igraph::graph.adjacency(abs(cov2cor(res.joint$theta_single[[4]]))>1e-5, mode='undirected')
df.deg.unstim.single = data.frame(degree=igraph::degree(g.unstim.single))
hist(df.deg.unstim.single$degree, breaks=20)
abline(v=df.deg.unstim.single$degree[which(genes_id=='LYZ')])
abline(v=df.deg.unstim.single$degree[which(genes_id=='YEATS4')], col='red')
abline(v=df.deg.unstim.single$degree[which(genes_id=='CREB1')], col='blue')


# What are the hubs? ------------------------------------------------------------

# Defining hubs as having top 10% highest degree

df.degree.ifn = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.ifn)))], rev(sort(igraph::degree(g.ifn)))))
names(df.degree.ifn)=c('Gene', 'Degree')
top.ifn = as.numeric(df.degree.ifn$Degree)>=quantile(as.numeric(df.degree.ifn$Degree), 0.9)
df.degree.ifn[top.ifn,]
length(df.degree.ifn[top.ifn,1]) # 56
df.degree.lps2 = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.lps2)))], rev(sort(igraph::degree(g.lps2)))))
names(df.degree.lps2)=c('Gene', 'Degree')
top.lps2 = as.numeric(df.degree.lps2$Degree)>=quantile(as.numeric(df.degree.lps2$Degree), 0.9)
df.degree.lps2[top.lps2,]
length(df.degree.lps2[top.lps2,1]) # 46
df.degree.lps24 = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.lps24)))], rev(sort(igraph::degree(g.lps24)))))
names(df.degree.lps24)=c('Gene', 'Degree')
top.lps24 = as.numeric(df.degree.lps24$Degree)>=quantile(as.numeric(df.degree.lps24$Degree), 0.9)
df.degree.lps24[top.lps24,]
length(df.degree.lps24[top.lps24,1]) # 49
df.degree.unstim = as.data.frame(cbind(genes_id[rev(order(igraph::degree(g.unstim)))], rev(sort(igraph::degree(g.unstim)))))
names(df.degree.unstim) =c('Gene', 'Degree')
top.unstim = as.numeric(df.degree.unstim$Degree)>=quantile(as.numeric(df.degree.unstim$Degree), 0.9)
df.degree.unstim[top.unstim,]
length(df.degree.unstim[top.unstim,1]) # 48

# How much do the hubs agree?
sum(df.degree.ifn[top.ifn,1]%in% df.degree.lps2[top.ifn,1])
# 33
sum(df.degree.ifn[top.lps2,1]%in% df.degree.lps24[top.lps24,1])
# 28
sum(df.degree.ifn[top.lps2,1]%in% df.degree.unstim[top.unstim,1])
# 27
sum(df.degree.lps2[top.ifn,1]%in% df.degree.lps24[top.lps24,1])
# 28
sum(df.degree.lps2[top.ifn,1]%in% df.degree.unstim[top.unstim,1])
# 30
sum(df.degree.lps24[top.lps24,1]%in% df.degree.unstim[top.unstim,1])
# 23

# Is YEATS4 is a hub in all? 
'YEATS4' %in% df.degree.ifn[top.ifn,1]
'YEATS4' %in% df.degree.lps2[top.lps2,1]
'YEATS4' %in% df.degree.lps24[top.lps24,1]
'YEATS4' %in% df.degree.unstim[top.unstim,1]
# All but LPS2.

# Is LYZ is a hub in all? 
'LYZ' %in% df.degree.ifn[top.ifn,1]
'LYZ' %in% df.degree.lps2[top.lps2,1]
'LYZ' %in% df.degree.lps24[top.lps24,1]
'LYZ' %in% df.degree.unstim[top.unstim,1]
# Only in unstim. 

# Show as heatmap (frac present in another.)
sim.list = list(df.degree.ifn[top.ifn,1], df.degree.lps2[top.lps2,1], df.degree.lps24[top.lps24,1], df.degree.unstim[top.unstim,1])
sim.heat = sapply(sim.list, function(x) sapply(sim.list, function(y) 100*sum(x %in% y)/(length(y)+length(x))))
diag(sim.heat) = NA
rownames(sim.heat)=colnames(sim.heat) = condition.names

png("Monocytes/plots/commonhubs_heatmap_jointGHS.png", width=1200, height =1200)
superheat(sim.heat, heat.col.scheme = "red", heat.lim = c(min(sim.heat,na.rm = T), max(sim.heat,na.rm = T)+0.2))
dev.off()

# Are any genes a hub in all conditions?
hub.in.all = genes_id[which(genes_id %in% df.degree.ifn[top.ifn,1] &  genes_id %in% df.degree.lps2[top.lps2,1] &  genes_id %in%  df.degree.lps24[top.lps24,1] &  genes_id %in%  df.degree.unstim[top.unstim,1])]
hub.in.all
# "ACBD5"    "AFMID"    "AFTPH"    "AIRE"     "COX6A1"   "GIMAP1"   "IMPDH1"   "KIAA0101" "LGALS3"   "RELB"     "SLC3A2"   "SORL1"    "STAG3L3" 


# Save as RData object 
top.hubs.ifn = df.degree.ifn[top.ifn,1]
top.hubs.lps2 = df.degree.lps2[top.lps2,1]
top.hubs.lps24 = df.degree.lps24[top.lps24,1]
top.hubs.unstim = df.degree.unstim[top.unstim,1]
save(top.hubs.ifn, top.hubs.lps2, top.hubs.lps24, top.hubs.unstim, file='Monocytes/data/top_hubs.RData')



# What are their neighbours? Looking at most influential hub

# IFN-GAMMA
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[1]),]!=0)],]

# lps2
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[1]),]!=0)],]
# YEATS4 a neighbour!

# lps24
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[1]),]!=0)],]

# Unstim
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[1]),]!=0)],]

# Only in lps2 is a cis gene a neighbour to the main hub.

# Second most influential

# IFN-GAMMA
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[2]),]!=0)],]
# YEATS4 a neighbour

# lps2
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[2]),]!=0)],]

# lps24
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[2]),]!=0)],]

# Unstim
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[2]),]!=0)],]
# LYZ a neighbour

# Third most influential

# IFN-GAMMA
df.degree.ifn$Gene[3] # AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[1]][which(genes_id==df.degree.ifn$Gene[3]),]!=0)],]
# LYZ a neighbour

# lps2
df.degree.lps2$Gene[3]  #AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[2]][which(genes_id==df.degree.lps2$Gene[3]),]!=0)],]
# LYZ a neighbour

# lps24
df.degree.lps24$Gene[3] # AFMID
annot_expr_full[genes_id[which(thetas.est.mono[[3]][which(genes_id==df.degree.lps24$Gene[3]),]!=0)],]
# LYZ a neighbour

# Unstim 
df.degree.unstim$Gene[3] # PHRF1
annot_expr_full[genes_id[which(thetas.est.mono[[4]][which(genes_id==df.degree.unstim$Gene[3]),]!=0)],]
# LYZ a neighbour

# So LYZ is a neighbour of one of the top three hubs in all conditions!

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

# LYZ has multiple top hubs as neighbours in each condition.
# YEATS4 has the most.

# Look at location of hubs
df.loc.ifn = data.frame(chromosome=factor(annot_expr[df.degree.ifn[top.ifn,1],]$Chr, levels=1:22))
g.loc.ifn = ggplot2::ggplot(df.loc.ifn,aes(x=chromosome))+geom_histogram(fill='lightblue', alpha=0.8, stat='count')+geom_vline(xintercept = factor(12), color='hotpink')+
  ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+scale_x_discrete(drop = FALSE)+ylim(0,8)
df.loc.lps2 = data.frame(chromosome=factor(annot_expr[df.degree.lps2[top.lps2,1],]$Chr, levels=1:22))
g.loc.lps2= ggplot2::ggplot(df.loc.lps2,aes(x=chromosome))+geom_histogram(fill='lightblue', alpha=0.8, stat='count')+geom_vline(xintercept = factor(12), color='hotpink')+
  ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+scale_x_discrete(drop = FALSE)+ylim(0,8)
df.loc.lps24 = data.frame(chromosome=factor(annot_expr[df.degree.lps24[top.lps24,1],]$Chr, levels=1:22))
g.loc.lps24= ggplot2::ggplot(df.loc.lps24,aes(x=chromosome))+geom_histogram(fill='lightblue', alpha=0.8, stat='count')+geom_vline(xintercept = factor(12), color='hotpink')+
  ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+scale_x_discrete(drop = FALSE)+ylim(0,8)
df.loc.unstim = data.frame(chromosome=factor(annot_expr[df.degree.unstim[top.unstim,1],]$Chr, levels=1:22))
g.loc.unstim= ggplot2::ggplot(df.loc.unstim,aes(x=chromosome))+geom_histogram(fill='lightblue', alpha=0.8, stat='count')+geom_vline(xintercept = factor(12), color='hotpink')+
  ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+scale_x_discrete(drop = FALSE)+ylim(0,8)

pdf("Monocytes/plots/tophubs_jointGHS.pdf", 8, 8)
gridExtra::grid.arrange(g.loc.ifn,g.loc.lps2, g.loc.lps24, g.loc.unstim)
dev.off()


# Order by chromosome -----------------------------------------------------------

# Order genes before plotting (Must sort this as order in annot_table does not match expr_all)
# Sorting according to (1) chromosome and (2) location
df.chrom.loc.tmp = data.frame(x=as.numeric(annot_expr_full$Chr), z=coords.all[,1])
ordered.ind = with(df.chrom.loc.tmp, order(x, z))
annot_expr_full_ordered = annot_expr_full[ordered.ind,]
annot_expr_full_ordered = annot_expr_full_ordered[rownames(annot_expr_full_ordered) %in% genes_id,] # Only the genes we have used for analysis
genes_id_ordered = rownames(annot_expr_full_ordered)
new.order.ind = match(genes_id_ordered, genes_id)
# Reorder estimated precision matrices
thetas.est.ordered = list()
for(i in 1:length(thetas.est.mono)){
  m.tmp = thetas.est.mono[[i]]
  thetas.est.ordered[[i]] = -cov2cor(m.tmp[new.order.ind, new.order.ind])
  diag(thetas.est.ordered[[i]]) = abs(diag(thetas.est.ordered[[i]]))
  thetas.est.ordered[[i]] = sign(thetas.est.ordered[[i]])
}

# Information about which networks an edge is present and about pos/neg correlation is shared though value/sign.


unique.list.ordered = list()
theta.joint.array.ordered = array(unlist(thetas.est.ordered), dim=c(ncol(thetas.est.ordered[[1]]),ncol(thetas.est.ordered[[1]]),length(thetas.est.ordered)))
theta.joint.array.ordered = theta.joint.array.ordered!=0
for(k in 1:length(thetas.est.ordered)){
  unique.list.ordered[[k]] = ifelse(theta.joint.array.ordered[,,k] !=0, 1, 0) # 1 means edge
  unique.list.ordered[[k]][which((theta.joint.array.ordered[,,k] !=0 ) & (apply(theta.joint.array.ordered[,,-k],c(1,2),FUN = function(s) sum(abs(s)) ) ==0))] = 2 # only present in network k
  unique.list.ordered[[k]][which(apply(theta.joint.array.ordered[,,],c(1,2),FUN = function(s) sum(abs(s)) )==length(thetas.est.ordered))] = 3 # present in all
  unique.list.ordered[[k]] = sign(thetas.est.ordered[[k]])*unique.list.ordered[[k]]
}


chrom.ordered = as.numeric(annot_expr_full_ordered$Chr)
chr.unique = unique(chrom.ordered)
cols.unique = rainbow(length(chr.unique))
palette.chrom = c(cols.unique)
names(palette.chrom) = chr.unique

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Plot neighbors of LYZ -----------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.lyz = diag(1,p)
  net.lyz[genes_id_ordered=='LYZ',] = unique.list.ordered[[i]][genes_id_ordered=='LYZ',]
  net.lyz[,genes_id_ordered=='LYZ'] = unique.list.ordered[[i]][,genes_id_ordered=='LYZ']
  net = network::network(net.lyz,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  names.net = rep('', p)
  names.net[genes_id_ordered=='LYZ'] = 'LYZ'
  network.vertex.names(net) = names.net
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.3, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
                             geom_text(aes(label =  names.net), size=3, hjust=1.2, vjust=1.2)+
                             ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/LYZ_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                                               legend.ordered, ncol=2, widths=c(8,1))
dev.off()


# Plot neighbors of YEATS4 -----------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.yeat = diag(1,p)
  net.yeat[genes_id_ordered=='YEATS4',] = unique.list.ordered[[i]][genes_id_ordered=='YEATS4',]
  net.yeat[,genes_id_ordered=='YEATS4'] = unique.list.ordered[[i]][,genes_id_ordered=='YEATS4']
  net = network::network(net.yeat,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  names.net = rep('', p)
  names.net[genes_id_ordered=='YEATS4'] = 'YEATS4'
  network.vertex.names(net) = names.net
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.3, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
    geom_text(aes(label =  names.net), size=3, hjust=1.2, vjust=1.2)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/YEATS4_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                                               legend.ordered, ncol=2, widths=c(8,1))
dev.off()

# Plot neighbors of CREB1 -----------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.creb = diag(1,p)
  net.creb[genes_id_ordered=='CREB1',] = unique.list.ordered[[i]][genes_id_ordered=='CREB1',]
  net.creb[,genes_id_ordered=='CREB1'] = unique.list.ordered[[i]][,genes_id_ordered=='CREB1']
  net = network::network(net.creb,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  names.net = rep('', p)
  names.net[genes_id_ordered=='CREB1'] = 'CREB1'
  network.vertex.names(net) = names.net
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.3, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
    geom_text(aes(label =  names.net), size=3, hjust=-0.5, vjust=-0.5)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/CREB1_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()

# Plot neighbors of all CIS genes -----------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.cis = diag(1,p)
  net.cis[annot_expr_full_ordered$Symbol %in% cis.genes.names,] = unique.list.ordered[[i]][annot_expr_full_ordered$Symbol %in% cis.genes.names,]
  net.cis[,annot_expr_full_ordered$Symbol %in% cis.genes.names] = unique.list.ordered[[i]][,annot_expr_full_ordered$Symbol %in% cis.genes.names]
  net = network::network(net.cis,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  names.net = rep('', p)
  names.net[annot_expr_full_ordered$Symbol %in% cis.genes.names] = genes_id_ordered[annot_expr_full_ordered$Symbol %in% cis.genes.names]
  network.vertex.names(net) = names.net
  h.just = rep(0,length(genes_id_ordered))
  v.just = rep(0,length(genes_id_ordered))
  h.just[annot_expr_full_ordered$Symbol %in% cis.genes.names]=  c(1.6,1.4) #c(1.3,1.2,3,1.4,1.4) # c(0.2,0.4,1.8,1.2,1.2) #rev((1:sum(annot_expr_full_ordered$cis_to_LYZ=='yes'))*0.4)
  v.just[annot_expr_full_ordered$Symbol %in% cis.genes.names]=  c(2,1)#rev(1:sum(annot_expr_full_ordered$Symbol %in% cis.genes.names))*1.6
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.3, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
    geom_text(aes(label =  names.net), size=2, hjust=h.just, vjust=v.just)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/CIS_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()


# Plot all edges ----------------------------------------------------------------- 


nets = list()
for(i in 1:length(unique.list.ordered)){
  net.all =  unique.list.ordered[[i]]
  net = network::network(net.all,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  #names.net = rep('', p)
  #names.net[genes_id_ordered=='YEATS4'] = 'YEATS4'
  #network.vertex.names(net) = names.net
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.1, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
    #geom_text(aes(label =  names.net), size=3, hjust=1.2, vjust=1.2)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/full_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()

# Plot shared edges ------------------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.all =  unique.list.ordered[[i]]
  net.all[which(abs(net.all)!=3)]=0
  net = network::network(net.all,directed=F, ignore.eval=F,names.eval='weights')
  set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
  # Assign color according to chromosome
  net %v% "chromosome" = chrom.ordered+0
  nets[[i]] = GGally::ggnet2(net,size = "degree", max_size = 6,edge.size = 0.1, alpha=0.7,mode = "circle",color = 'chromosome', palette=palette.chrom,
                             edge.color = 'color',edge.lty = "lty",legend.size = 10)+
    #geom_text(aes(label =  names.net), size=3, hjust=1.2, vjust=1.2)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/common_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]]+ theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()




# Investigating particular edges-------------------------------------------------


# Find the one edge that is negative in all but IFN-GAMMA

ind.creb1 = which(genes_id=='CREB1')
gene.neg = genes_id[which(sign(thetas.est.mono[[1]][ind.creb1,])==-1 &  sign(thetas.est.mono[[2]][ind.creb1,])==1 & sign(thetas.est.mono[[3]][ind.creb1,])==1 & sign(thetas.est.mono[[4]][ind.creb1,])==1)]
annot_expr[gene.neg,]


thetas.est.mono[[1]][which(genes_id==gene.neg),which(genes_id=='CREB1')]
thetas.est.mono[[2]][which(genes_id==gene.neg),which(genes_id=='CREB1')]
thetas.est.mono[[3]][which(genes_id==gene.neg),which(genes_id=='CREB1')]
thetas.est.mono[[4]][which(genes_id==gene.neg),which(genes_id=='CREB1')]



which.max(res.joint$E_NuInv[[1]],arr.ind=T)


# Checking for enrichment-------------------------------------------------------

top.hubs = list(df.degree.ifn[top.ifn,1], df.degree.lps2[top.lps2,1],df.degree.lps24[top.lps24,1],df.degree.unstim[top.unstim,1])
save(top.hubs, file='Monocytes/data/temp.RData')


library(enrichR)
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
head(dbs)
dbs$libraryName[grepl('Disease', dbs$libraryName, fixed = TRUE)]


dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
#dbs = c('Chromosome_Location')
enriched.ifn <- enrichr(top.hubs[[1]], dbs)
enriched.lps2 <- enrichr(top.hubs[[2]], dbs)
enriched.lps24 <- enrichr(top.hubs[[3]], dbs)
enriched.unstim <- enrichr(top.hubs[[4]], dbs)


enriched[["GO_Molecular_Function_2015"]]
plotEnrich(enriched.ifn[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps2[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps24[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.unstim[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

dbs = c("Disease_Signatures_from_GEO_up_2014", "OMIM_Disease", "Disease_Signatures_from_GEO_down_2014" )
#dbs = c('Chromosome_Location')
enriched.ifn <- enrichr(top.hubs[[1]], dbs)
enriched.lps2 <- enrichr(top.hubs[[2]], dbs)
enriched.lps24 <- enrichr(top.hubs[[3]], dbs)
enriched.unstim <- enrichr(top.hubs[[4]], dbs)
plotEnrich(enriched.ifn[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps2[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps24[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.unstim[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

plotEnrich(enriched.ifn[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps2[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.lps24[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched.unstim[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")


