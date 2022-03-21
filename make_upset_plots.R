rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_and_bcells_LYZ_region.RData")
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")
load('Monocytes/data/monocytes_jointGHS_smalleps_largeAIC.RData')
load('Monocytes/data/genes_flagged.RData')
load('Monocytes/data/top_hubs.RData')
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

# Make plot showing their agreement, with chromosome marked --------------------------------------------

# Get the right format

dat.edges = data.frame(IFN =thetas.est.mono[[1]][upper.tri(thetas.est.mono[[1]])]!=0, LPS2 =thetas.est.mono[[2]][upper.tri(thetas.est.mono[[2]])]!=0, 
                       LPS24 =thetas.est.mono[[3]][upper.tri(thetas.est.mono[[3]])]!=0, Unstim =thetas.est.mono[[4]][upper.tri(thetas.est.mono[[4]])]!=0)
names(dat.edges) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

plot.upset = upset(dat.edges, names(dat.edges)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=50000, sort_sets=F)


# Add chromosome information
dat.edges = rbind(dat.edges,dat.edges)
chr.info = annot_expr[genes_id, 'Chr']
chr.info.row = matrix(rep(chr.info,p), nrow=p, byrow=T)
chr.info.col = matrix(rep(chr.info,p), nrow=p, byrow=F)
chr.info.row.upper = chr.info.row[upper.tri(chr.info.row)]
chr.info.col.upper = chr.info.col[upper.tri(chr.info.col)]
dat.edges$`Gene` = c(factor(chr.info.row.upper, levels=1:22), factor(chr.info.col.upper, levels=1:22))

plot.upset.chr = upset(dat.edges, names(dat.edges)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=100000,
                       annotations = list(
                         'Chromosome'=(
                           ggplot(mapping=aes(fill=`Gene`))
                           + geom_bar(stat='count', position='fill')
                           + scale_y_continuous(labels=scales::percent_format())
                           + scale_fill_manual(values= rainbow(22))
                           + ylab('Chromosome')
                         )
                       ))$patches$plots[[2]]


layout <- "
##AAAAAAAAAAAA
BBBBBBBBBBBBBB
"

pdf("Monocytes/plots/intersection_jointGHS.pdf", 12, 10)
plot.upset.chr/plot.upset + plot_layout(design = layout, heights = c(1.3,2))
dev.off()

# Compare intersection of single and joint estimates (at same sparsities) ---------------------------


# Get intersection order
inter.order = list()
for(k in 1:length(plot.upset$data$intersection)){
  inter.order[[k]] = c("IFN-gamma", 'LPS-24h', 'LPS-2h','Unstim')[as.numeric(strsplit(toString(plot.upset$data$intersection[k]), '-')[[1]])]
}

dat.edges.single = data.frame(IFN =thetas.est.mono.single[[1]][upper.tri(thetas.est.mono.single[[1]])]!=0, LPS2 =thetas.est.mono.single[[2]][upper.tri(thetas.est.mono.single[[2]])]!=0, 
                       LPS24 =thetas.est.mono.single[[3]][upper.tri(thetas.est.mono.single[[3]])]!=0, Unstim =thetas.est.mono.single[[4]][upper.tri(thetas.est.mono.single[[4]])]!=0)
names(dat.edges.single) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

plot.upset.single = upset(dat.edges.single, names(dat.edges.single)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=50000, sort_sets=F, sort_intersections=F, 
                          intersections=rev(inter.order))

layout <- "
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB#
"

pdf("Monocytes/plots/intersection_withsingle_jointGHS.pdf", 12, 10)
plot.upset.single/plot.upset #+ plot_layout( heights = c(1,1.1))
dev.off()





# Same plot as first, but marked by whether it's controlled by a hotspot instead--------------------

# Add hotspot information (one row per condition)

# only TRUE if hotspot in all conditions in intersection

for(i in 1:length(condition.names)){
  dat.edges = rbind(dat.edges,dat.edges) # Doubling up
  hotspot.info = genes_id %in% genes.flagged[[i]]
  hotspot.info.row = matrix(rep(hotspot.info,p), nrow=p, byrow=T)
  hotspot.info.col = matrix(rep(hotspot.info,p), nrow=p, byrow=F)
  hotspot.info.row.upper = hotspot.info.row[upper.tri(hotspot.info.row)]
  hotspot.info.col.upper = hotspot.info.col[upper.tri(hotspot.info.col)]
  titl = paste0('Hotspot in ', condition.names[i])
  dat.edges$tmp = c(hotspot.info.row.upper, hotspot.info.col.upper)
  names(dat.edges)[5+i] = titl
}


# Add Hotspot information

# Loop over edges 
dat.edges.full = dat.edges
dat.edges.full$hotspot = rep('no', nrow(dat.edges.full))
for (i in 1:nrow(dat.edges)){
  which.true = which(dat.edges[i,1:4]==1)
  if(length(which.true)!=0){
    dat.edges.full$hotspot[i] = ifelse(sum(dat.edges[i,which.true+5])==length(which.true), 'yes', 'no')
  }
}
names(dat.edges.full)[10] = 'Top hotspot controlled'


# Loop over rows in dat.edges
#dat.edges.full = dat.edges
#for (i in 1:nrow(dat.edges)){
#  which.true = which(dat.edges[i,1:4]==1)
#  if(length(which.true)!=0){
#    # For each condition with the hub, add info
#    for(j in 1:length(which.true)){
#      if(j==1){
#        # Each row corresponds to one of the two genes in an edge
#        # Check if genes corresponding to row i, condition j is flagged using existing variables
#        dat.edges.full$hotspot[i] = ifelse(dat.edges[i,j+5]==T, 'yes', 'no')
#      }
#      else{
#        # Add row
#        dat.edges.full=rbind(dat.edges.full, dat.edges.full[i,])
#        dat.edges.full$hotspot[nrow(dat.edges.full)] = ifelse(dat.edges[i,j+5]==T, 'yes', 'no')
#      }
#    }
#  }
#}
#names(dat.edges.full)[10] = 'Top hotspot controlled'

plot.upset.hotspot = upset(dat.edges.full, names(dat.edges.full)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=100000,
                       annotations = list(
                         'Controlled by top hotspot'=(
                           ggplot(mapping=aes(fill=`Top hotspot controlled`))
                           + geom_bar(stat='count', position='fill')
                           + scale_y_continuous(labels=scales::percent_format())
                           + scale_fill_manual(values=c('ivory4','red2'))
                           + ylab('Controlled buy top hotspot')
                         )
                       ))$patches$plots[[2]]

layout <- "
##AAAAAAAAAAAA
BBBBBBBBBBBBBB
"

pdf("Monocytes/plots/intersection_hotspotmarked_jointGHS.pdf", 12, 10)
plot.upset.hotspot/plot.upset + plot_layout(design = layout, heights = c(1.3,2))
dev.off()


# Upset plot for common hubs ------------------------------------------------------


# Get the right format

dat.hubs = data.frame(IFN = genes_id %in% top.hubs.ifn, LPS2 = genes_id %in% top.hubs.lps2, LPS24 = genes_id %in% top.hubs.lps24, UNSTIM = genes_id %in% top.hubs.unstim)
names(dat.hubs) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

plot.upset = upset(dat.hubs, names(dat.hubs)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=200, sort_sets=F)


# Add Hotspot information
dat.hubs$hotspot = rep('no', nrow(dat.hubs))
  
# Loop over rows in dat.hubs (Loop over genes)
dat.hubs.full = dat.hubs
for (i in 1:nrow(dat.hubs)){
  # Iteration i = gene_id[i]
  which.true = which(dat.hubs[i,1:4]==1)
  if(length(which.true)!=0){
    dat.hubs.full$hotspot[i] = ifelse(sum(unlist(lapply(1:length(which.true), FUN=  function(j) genes_id[i] %in% genes.flagged[[which.true[j]]])))==length(which.true), 'yes', 'no')
  }
}

names(dat.hubs.full)[5] = 'Top hotspot controlled'

# Loop over rows in dat.hubs
#dat.hubs.full = dat.hubs
#for (i in 1:nrow(dat.hubs)){
#  which.true = which(dat.hubs[i,1:4]==1)
#  if(length(which.true)!=0){
#    # For each condition with the hub, add info
#    for(j in 1:length(which.true)){
#      if(j==1){
#        # Check if gene corresponding to row i is flagged
#        dat.hubs.full$hotspot[i] = ifelse(genes_id[i] %in% genes.flagged[[which.true[j]]], 'yes', 'no')
#      }
#      else{
#        # Add row
#        dat.hubs.full=rbind(dat.hubs.full, dat.hubs.full[i,])
#        dat.hubs.full$hotspot[nrow(dat.hubs.full)] = ifelse(genes_id[i] %in% genes.flagged[[which.true[j]]], 'yes', 'no')
#      }
#    }
#  }
#}

#names(dat.hubs.full)[5] = 'Top hotspot controlled'



plot.upset.hubs = upset(dat.hubs.full, names(dat.hubs.full)[1:4], name='condition', width_ratio=0.1,keep_empty_groups=TRUE, max_size=200,
                       annotations = list(
                         'Chromosome'=(
                           ggplot(mapping=aes(fill=`Top hotspot controlled`))
                           + geom_bar(stat='count', position='fill')
                           + scale_y_continuous(labels=scales::percent_format())
                           + scale_fill_manual(values=c('ivory4','red2'))
                           + ylab('')
                         )
                       ))$patches$plots[[2]]


layout <- "
##AAAAAAAAAAAA
BBBBBBBBBBBBBB
"

pdf("Monocytes/plots/intersection_hubs_jointGHS.pdf", 12, 10)
plot.upset.hubs/plot.upset + plot_layout(design = layout, heights = c(1.3,2))
dev.off()


# Histrogram of hubs and all nodes, differentiating between 1st, 2nd and 3rd order neighbours of cis genes (LYZ and YEATS4).----------------------

library(igraph)
library(RColorBrewer)


# Histogram of hubs (Not upset plot)

dat.hubs.neigh = data.frame(IFN = genes_id %in% top.hubs.ifn, LPS2 = genes_id %in% top.hubs.lps2, LPS24 = genes_id %in% top.hubs.lps24, UNSTIM = genes_id %in% top.hubs.unstim)
names(dat.hubs) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

lyz.ind = which(genes_id=='LYZ')
yeats4.ind = which(genes_id=='YEATS4')
#creb1.ind = which(genes_id=='CREB1')
gg.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected')
dists.lyz.ifn = igraph::distances(gg.ifn, to=lyz.ind) # Get distances
dists.y4.ifn = igraph::distances(gg.ifn, to=yeats4.ind) # Get distances
gg.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected')
dists.lyz.lps2 = igraph::distances(gg.lps2, to=lyz.ind) # Get distances
dists.y4.lps2 = igraph::distances(gg.lps2, to=yeats4.ind) # Get distances
gg.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected')
dists.lyz.lps24 = igraph::distances(gg.lps24, to=lyz.ind) # Get distances
dists.y4.lps24 = igraph::distances(gg.lps24, to=yeats4.ind) # Get distances
gg.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected')
dists.lyz.unstim = igraph::distances(gg.unstim, to=lyz.ind) # Get distances
dists.y4.unstim = igraph::distances(gg.unstim, to=yeats4.ind) # Get distances
#dists.creb1 = igraph::distances(gg.ifn, to=creb1.ind) # Get distances
dists.to.cis.ifn = apply(cbind(dists.lyz.ifn, dists.y4.ifn), 1, min)
dists.to.cis.lps2 = apply(cbind(dists.lyz.lps2, dists.y4.lps2), 1, min)
dists.to.cis.lps24 = apply(cbind(dists.lyz.lps24, dists.y4.lps24), 1, min)
dists.to.cis.unstim = apply(cbind(dists.lyz.unstim, dists.y4.unstim), 1, min)
dists.cis.ifn = dists.to.cis.ifn[which(genes_id %in% top.hubs.ifn)]
dists.cis.lps2 = dists.to.cis.lps2[which(genes_id %in% top.hubs.lps2)]
dists.cis.lps24 = dists.to.cis.lps24[which(genes_id %in% top.hubs.lps24)]
dists.cis.unstim = dists.to.cis.unstim[which(genes_id %in% top.hubs.unstim)]
dists.cis.all = c(dists.cis.ifn, dists.cis.lps2, dists.cis.lps24, dists.cis.unstim)
dists.cis.all = dists.cis.all[-which(dists.cis.all==0)]
dists.cis.all = c('1st','2nd','3rd','4th')[dists.cis.all]
dists.factor = factor(dists.cis.all,levels=c('1st','2nd','3rd','4th'))

dat.hubs.neigh = data.frame(Neighbour = dists.factor,
                            Condition= c(rep("IFN-gamma", length(dists.cis.ifn)-1),rep("LPS-2h", length(dists.cis.lps2)-1),
                                         rep("LPS-24h", length(dists.cis.lps24)-1),rep("Unstim", length(dists.cis.unstim)-1)))


pdf("Monocytes/plots/hist_cisorder_hubs_jointGHS.pdf", 10, 12)
ggplot(dat.hubs.neigh, aes(x=Condition, fill=Neighbour))+geom_bar(position='fill')+ guides(fill=guide_legend(title="Cis neighbour order"))+
      scale_fill_manual(values=RColorBrewer::brewer.pal(3,'YlOrRd')[3:1])
dev.off()

# And then for ALL nodes:

dists.to.cis.all = c(dists.to.cis.ifn, dists.to.cis.lps2, dists.to.cis.lps24, dists.to.cis.unstim)
dists.to.cis.all = dists.to.cis.all[-which(dists.to.cis.all==0)]
dists.to.cis.all[which(dists.to.cis.all==Inf)] = 6
dists.to.cis.all = c('1st','2nd','3rd','4th', '5th', 'Not a neighbour')[dists.to.cis.all]
dists.factor = factor(dists.to.cis.all,levels=c('1st','2nd','3rd','4th', '5th', 'Not a neighbour'))

dat.hubs.neigh = data.frame(Neighbour = dists.factor,
                            Condition= c(rep("IFN-gamma", p-2),rep("LPS-2h", p-2),
                                         rep("LPS-24h", p-2),rep("Unstim", p-2)))


pdf("Monocytes/plots/hist_cisorder_all_jointGHS.pdf", 10, 12)
ggplot(dat.hubs.neigh, aes(x=Condition, fill=Neighbour))+geom_bar(position='fill')+ guides(fill=guide_legend(title="Cis neighbour order"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(6,'YlOrRd')[6:1])
dev.off()


# NOT DONE since we can't average over conditions.
# Make upset plot, marked by the smallest distance to a cis gene among the two genes.
#dat.edges = data.frame(IFN =thetas.est.mono[[1]][upper.tri(thetas.est.mono[[1]])]!=0, LPS2 =thetas.est.mono[[2]][upper.tri(thetas.est.mono[[2]])]!=0, 
#                      LPS24 =thetas.est.mono[[3]][upper.tri(thetas.est.mono[[3]])]!=0, Unstim =thetas.est.mono[[4]][upper.tri(thetas.est.mono[[4]])]!=0)
#names(dat.edges) = c("IFN-gamma", 'LPS-2h', 'LPS-24h', 'Unstim')

#lyz.ind = which(genes_id=='LYZ')
#yeats4.ind = which(genes_id=='YEATS4')
#gg.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected')
#dists.lyz = igraph::distances(gg.ifn, to=lyz.ind) # Get distances
#dists.y4 = igraph::distances(gg.ifn, to=yeats4.ind) # Get distances
#dists.to.cis = apply(cbind(dists.lyz, dists.y4), 1, min)
#dat.edges.neigh = dat.edges[,1:4]

# Add info on neighbours
#dat.edges.neigh = rbind(dat.edges.neigh,dat.edges.neigh) # one row per gene
#neigh.info.row = matrix(rep(dists.to.cis,p), nrow=p, byrow=T)
#neigh.info.col = matrix(rep(dists.to.cis,p), nrow=p, byrow=F)
#neigh.info.row.upper = neigh.info.row[upper.tri(neigh.info.row)]
#neigh.info.col.upper = neigh.info.col[upper.tri(neigh.info.col)]
#dat.edges.neigh$Neighbour = c(factor(chr.info.row.upper, levels=c(1,2,3,4,Inf)), factor(chr.info.col.upper, levels=c(1,2,3,4,Inf)))







# Evaluate significance of density of top hotspot mediated subnetworks ------------------------

# IFN
spars.ifn = tailoredGlasso::sparsity(thetas.est.mono[[1]][which(genes_id %in% genes.flagged[[1]]),which(genes_id %in% genes.flagged[[1]])])
spars.ifn 
#0.01852755
# LPS2
spars.lps2 = tailoredGlasso::sparsity(thetas.est.mono[[2]][which(genes_id %in% genes.flagged[[2]]),which(genes_id %in% genes.flagged[[2]])])
spars.lps2
# 0.04205852
# LPS24
spars.lps24 = tailoredGlasso::sparsity(thetas.est.mono[[3]][which(genes_id %in% genes.flagged[[3]]),which(genes_id %in% genes.flagged[[3]])])
spars.lps24
# 0.1833333
# Unstim
spars.unstim = tailoredGlasso::sparsity(thetas.est.mono[[4]][which(genes_id %in% genes.flagged[[4]]),which(genes_id %in% genes.flagged[[4]])])
spars.unstim
# 0.02242991

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
# 0.08279172
(sum(spars.vals.lps2>spars.lps2)+1)/(n.sample+1)
# 9.999e-05
(sum(spars.vals.lps24>spars.lps24)+1)/(n.sample+1)
# 9.999e-05
(sum(spars.vals.unstim>spars.unstim)+1)/(n.sample+1)
# 9.999e-05


# Evaluate siginificance of neighbours of LYZ and YEATS4 being top hotspot mediated in different conditions ------------------------------

# IFN
neigh.lyz.ifn = which(thetas.est.mono[[1]][which(genes_id == 'LYZ'),]!=0)
neigh.lyz.ifn = neigh.lyz.ifn[-which(neigh.lyz.ifn== which(genes_id == 'LYZ'))]
mediated.lyz.ifn = neigh.lyz.ifn[which(genes_id[neigh.lyz.ifn] %in% genes.flagged[[1]])]
neigh.y4.ifn = which(thetas.est.mono[[1]][which(genes_id == 'YEATS4'),]!=0)
neigh.y4.ifn = neigh.y4.ifn[-which(neigh.y4.ifn== which(genes_id == 'YEATS4'))]
mediated.y4.ifn = neigh.y4.ifn[which(genes_id[neigh.y4.ifn] %in% genes.flagged[[1]])]
neigh.creb1.ifn = which(thetas.est.mono[[1]][which(genes_id == 'CREB1'),]!=0)
neigh.creb1.ifn = neigh.creb1.ifn[-which(neigh.creb1.ifn== which(genes_id == 'CREB1'))]
mediated.creb1.ifn = neigh.creb1.ifn[which(genes_id[neigh.creb1.ifn] %in% genes.flagged[[1]])]
# LPS2
neigh.lyz.lps2 = which(thetas.est.mono[[2]][which(genes_id == 'LYZ'),]!=0)
neigh.lyz.lps2 = neigh.lyz.lps2[-which(neigh.lyz.lps2== which(genes_id == 'LYZ'))]
mediated.lyz.lps2 = neigh.lyz.lps2[which(genes_id[neigh.lyz.lps2] %in% genes.flagged[[2]])]
neigh.y4.lps2  = which(thetas.est.mono[[2]][which(genes_id == 'YEATS4'),]!=0)
neigh.y4.lps2 = neigh.y4.lps2[-which(neigh.y4.lps2== which(genes_id == 'YEATS4'))]
mediated.y4.lps2  = neigh.y4.lps2[which(genes_id[neigh.y4.lps2] %in% genes.flagged[[2]])]
neigh.creb1.lps2 = which(thetas.est.mono[[2]][which(genes_id == 'CREB1'),]!=0)
neigh.creb1.lps2 = neigh.creb1.lps2[-which(neigh.creb1.lps2== which(genes_id == 'CREB1'))]
mediated.creb1.lps2 = neigh.creb1.lps2[which(genes_id[neigh.creb1.lps2] %in% genes.flagged[[2]])]
# LPS24
neigh.lyz.lps24 = which(thetas.est.mono[[3]][which(genes_id == 'LYZ'),]!=0)
neigh.lyz.lps24 = neigh.lyz.lps24[-which(neigh.lyz.lps24== which(genes_id == 'LYZ'))]
mediated.lyz.lps24 = neigh.lyz.lps24[which(genes_id[neigh.lyz.lps24] %in% genes.flagged[[3]])]
neigh.y4.lps24  = which(thetas.est.mono[[3]][which(genes_id == 'YEATS4'),]!=0)
neigh.y4.lps24 = neigh.y4.lps24[-which(neigh.y4.lps24== which(genes_id == 'YEATS4'))]
mediated.y4.lps24  = neigh.y4.lps24[which(genes_id[neigh.y4.lps24] %in% genes.flagged[[3]])]
neigh.creb1.lps24 = which(thetas.est.mono[[3]][which(genes_id == 'CREB1'),]!=0)
neigh.creb1.lps24 = neigh.creb1.lps24[-which(neigh.creb1.lps24== which(genes_id == 'CREB1'))]
mediated.creb1.lps24 = neigh.creb1.lps24[which(genes_id[neigh.creb1.lps24] %in% genes.flagged[[3]])]
# Unstim
neigh.lyz.unstim = which(thetas.est.mono[[4]][which(genes_id == 'LYZ'),]!=0)
neigh.lyz.unstim = neigh.lyz.unstim[-which(neigh.lyz.unstim== which(genes_id == 'LYZ'))]
mediated.lyz.unstim  = neigh.lyz.unstim[which(genes_id[neigh.lyz.unstim] %in% genes.flagged[[4]])]
neigh.y4.unstim   = which(thetas.est.mono[[4]][which(genes_id == 'YEATS4'),]!=0)
neigh.y4.unstim = neigh.y4.unstim[-which(neigh.y4.unstim== which(genes_id == 'YEATS4'))]
mediated.y4.unstim   = neigh.y4.unstim[which(genes_id[neigh.y4.unstim] %in% genes.flagged[[4]])]
neigh.creb1.unstim = which(thetas.est.mono[[4]][which(genes_id == 'CREB1'),]!=0)
neigh.creb1.unstim = neigh.creb1.unstim[-which(neigh.creb1.unstim== which(genes_id == 'CREB1'))]
mediated.creb1.unstim = neigh.creb1.unstim[which(genes_id[neigh.creb1.unstim] %in% genes.flagged[[4]])]


n.sample = 10000
n.vals.lyz.ifn = c()
n.vals.y4.ifn = c()
n.vals.creb1.ifn = c()
n.vals.lyz.lps2 = c()
n.vals.y4.lps2 = c()
n.vals.creb1.lps2 = c()
n.vals.lyz.lps24 = c()
n.vals.y4.lps24 = c()
n.vals.creb1.lps24 = c()
n.vals.lyz.unstim = c()
n.vals.y4.unstim = c()
n.vals.creb1.unstim = c()
set.seed(123)
for(i in 1:n.sample){
  # Sample a set of genes, and count how many are in flagged
  n.vals.lyz.ifn[i] = length(which(genes_id[sample(1:p, length(neigh.lyz.ifn))]%in% genes.flagged[[1]]))
  n.vals.y4.ifn[i] = length(which(genes_id[sample(1:p, length(neigh.y4.ifn))]%in% genes.flagged[[1]]))
  n.vals.creb1.ifn[i] = length(which(genes_id[sample(1:p, length(neigh.creb1.ifn))]%in% genes.flagged[[1]]))
  n.vals.lyz.lps2[i] = length(which(genes_id[sample(1:p, length(neigh.lyz.lps2))]%in% genes.flagged[[2]]))
  n.vals.y4.lps2[i] = length(which(genes_id[sample(1:p, length(neigh.y4.lps2))]%in% genes.flagged[[2]]))
  n.vals.creb1.lps2[i] = length(which(genes_id[sample(1:p, length(neigh.creb1.lps2))]%in% genes.flagged[[2]]))
  n.vals.lyz.lps24[i] = length(which(genes_id[sample(1:p, length(neigh.lyz.lps24))]%in% genes.flagged[[3]]))
  n.vals.y4.lps24[i] = length(which(genes_id[sample(1:p, length(neigh.y4.lps24))]%in% genes.flagged[[3]]))
  n.vals.creb1.lps24[i] = length(which(genes_id[sample(1:p, length(neigh.creb1.lps24))]%in% genes.flagged[[3]]))
  n.vals.lyz.unstim[i] = length(which(genes_id[sample(1:p, length(neigh.lyz.unstim))]%in% genes.flagged[[4]]))
  n.vals.y4.unstim[i] = length(which(genes_id[sample(1:p, length(neigh.y4.unstim))]%in% genes.flagged[[4]]))
  n.vals.creb1.unstim[i] = length(which(genes_id[sample(1:p, length(neigh.creb1.unstim))]%in% genes.flagged[[4]]))
}

(sum(n.vals.lyz.ifn>length(mediated.lyz.ifn))+1)/(n.sample+1)
# 0.3520648
(sum(n.vals.y4.ifn>length(mediated.y4.ifn))+1)/(n.sample+1)
# 0.870313
(sum(n.vals.creb1.ifn>length(mediated.creb1.ifn))+1)/(n.sample+1)
#  9.999e-05

(sum(n.vals.lyz.lps2>length(mediated.lyz.lps2))+1)/(n.sample+1)
# 0.03389661
(sum(n.vals.y4.lps2>length(mediated.y4.lps2))+1)/(n.sample+1)
# 0.5551445
(sum(n.vals.creb1.lps2>length(mediated.creb1.lps2))+1)/(n.sample+1)
# 0.00229977

(sum(n.vals.lyz.lps24>length(mediated.lyz.lps24))+1)/(n.sample+1)
# 0.00039996
(sum(n.vals.y4.lps24>length(mediated.y4.lps24))+1)/(n.sample+1)
# 0.00049995
(sum(n.vals.creb1.lps24>length(mediated.creb1.lps24))+1)/(n.sample+1)
# 9.999e-05

(sum(n.vals.lyz.unstim>length(mediated.lyz.unstim))+1)/(n.sample+1)
# 0.03269673
(sum(n.vals.y4.unstim>length(mediated.y4.unstim))+1)/(n.sample+1)
# 0.7830217
(sum(n.vals.creb1.unstim>length(mediated.creb1.unstim))+1)/(n.sample+1)
# 9.999e-05

# We have too many top hotspot mediated genes in IFNg and Unstim to be able to interpret.
# Notably the neighbours of CREB1 are significantly mediated in all.
# Significant in LPS-24h, but very few hotspot mediated here to begin with.



# Just looking at the numbers instead: fraction of neighbour that are mediated: 

frac.lyz.ifn = length(mediated.lyz.ifn)/length(neigh.lyz.ifn)
frac.y4.ifn = length(mediated.y4.ifn)/length(neigh.y4.ifn)
frac.creb1.ifn = length(mediated.creb1.ifn)/length(neigh.creb1.ifn)

frac.lyz.lps2= length(mediated.lyz.lps2)/length(neigh.lyz.lps2)
frac.y4.lps2 = length(mediated.y4.lps2)/length(neigh.y4.lps2)
frac.creb1.lps2 =length(mediated.creb1.lps2)/length(neigh.creb1.lps2)

frac.lyz.lps24 = length(mediated.lyz.lps24)/length(neigh.lyz.lps24)
frac.y4.lps24 =length(mediated.y4.lps24)/length(neigh.y4.lps24)
frac.creb1.lps24 = length(mediated.creb1.lps24)/length(neigh.creb1.lps24)

frac.lyz.unstim = length(mediated.lyz.unstim)/length(neigh.lyz.unstim)
frac.y4.unstim = length(mediated.y4.unstim)/length(neigh.y4.unstim)
frac.creb1.unstim = length(mediated.creb1.unstim)/length(neigh.creb1.unstim)

# High proportion, particularly in IFN

# Compare to overall fraction: 

frac.mediated.ifn = length(genes.flagged[[1]])/p
frac.mediated.lps2 = length(genes.flagged[[2]])/p
frac.mediated.lps24 = length(genes.flagged[[3]])/p
frac.mediated.unstim = length(genes.flagged[[4]])/p

df.mediated = data.frame(IFN=c(frac.mediated.ifn, frac.lyz.ifn, frac.y4.ifn, frac.creb1.ifn), LPS2 = c(frac.mediated.lps2, frac.lyz.lps2, frac.y4.lps2, frac.creb1.lps2), 
                         LPS24 = c(frac.mediated.lps24, frac.lyz.lps24, frac.y4.lps24, frac.creb1.lps24),
                         Unstim = c(frac.mediated.unstim, frac.lyz.unstim, frac.y4.unstim, frac.creb1.unstim)) 
names(df.mediated) = condition.names
rownames(df.mediated) = c('Overall', 'LYZ neighbourhood', 'YEATS4 neighbourhood', 'CREB1 neighbourhood')
df.mediated= round(df.mediated,3)
df.mediated
#                     IFN-gamma LPS-2h LPS-24h Unstim
#Overall                  0.772  0.231   0.042  0.564
#LYZ neighbourhood        0.778  0.417   0.333  0.769
#YEATS4 neighbourhood     0.625  0.167   0.222  0.429
#CREB1 neighbourhood      1.000  0.600   0.333  1.000

# YEATS4 not really more mediated neighbours except in LPS-24h. 
# CREB1 more mediated neighbours in all! And all its neighbours are mediated in IFNg and Unstim!
# LYZ more in all, though not very much more in IFNg.

# In LPS-24h, all have more mediated neighbours. 

