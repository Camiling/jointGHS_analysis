rm(list=ls())
source('utils.R')
load("Monocytes/data/expression_data_monocytes_and_bcells_LYZ_region_notOnGit.RData")
load("Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5_notOnGit.RData")
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

# Check edge agreement
100*round(tailoredGlasso::precision(thetas.est.mono.single[[1]]!=0, thetas.est.mono[[1]]!=0), 3)
100*round(tailoredGlasso::precision(thetas.est.mono.single[[2]]!=0, thetas.est.mono[[2]]!=0), 3)
100*round(tailoredGlasso::precision(thetas.est.mono.single[[3]]!=0, thetas.est.mono[[3]]!=0), 3)
100*round(tailoredGlasso::precision(thetas.est.mono.single[[4]]!=0, thetas.est.mono[[4]]!=0), 3)

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
  g.deg.ifn.noline = ggplot2::ggplot(df.deg.ifn, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.ifn = ggplot2::ggplot(df.deg.ifn, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[1]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected')
  df.deg.lps2 = data.frame(degree=igraph::degree(g.lps2))
  rownames(df.deg.lps2) = genes_id
  deg.cis = df.deg.lps2[cis.names[i],]
  g.deg.lps2.noline =ggplot2::ggplot(df.deg.lps2, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.lps2 =ggplot2::ggplot(df.deg.lps2, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[2]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected')
  df.deg.lps24 = data.frame(degree=igraph::degree(g.lps24))
  rownames(df.deg.lps24) = genes_id
  deg.cis = df.deg.lps24[cis.names[i],]
  g.deg.lps24.noline = ggplot2::ggplot(df.deg.lps24, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.lps24 = ggplot2::ggplot(df.deg.lps24, aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[3]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected')
  df.deg.unstim = data.frame(degree=igraph::degree(g.unstim ))
  rownames(df.deg.unstim) = genes_id
  deg.cis = df.deg.unstim[cis.names[i],]
  g.deg.unstim.noline =ggplot2::ggplot(df.deg.unstim , aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  g.deg.unstim =ggplot2::ggplot(df.deg.unstim , aes(x=degree))+geom_histogram(fill='lightgrey', bins=30)+geom_vline(xintercept = deg.cis, color='red')+
    ggplot2::ggtitle(condition.names[4]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_", cis.names[i],"_combined_jointGHS.pdf"), 8, 8)
  gridExtra::grid.arrange(g.deg.ifn, g.deg.lps2, g.deg.lps24, g.deg.unstim)
  dev.off()
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

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_hist_allcis_jointGHS.pdf"), 8, 8)
gridExtra::grid.arrange(g.deg.ifn.noline+geom_vline(xintercept = degs.cis.ifn[1]+0.2, color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.ifn[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.ifn[3], color='red', alpha=0.7), 
                        g.deg.lps2.noline+geom_vline(xintercept = degs.cis.lps2[1]+0.2, color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.lps2[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.lps2[3], color='red', alpha=0.7),
                        legend,
                        g.deg.lps24.noline+geom_vline(xintercept = degs.cis.lps24[1], color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.lps24[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.lps24[3], color='red', alpha=0.7),
                        g.deg.unstim.noline+geom_vline(xintercept = degs.cis.unstim[1], color='green', alpha=0.7) + geom_vline(xintercept = degs.cis.unstim[2], color='blue', alpha=0.7) +geom_vline(xintercept = degs.cis.unstim[3], color='red', alpha=0.7),
                        ncol=3,widths=c(1,1,0.3))
dev.off()


# Also make density plot
g.ifn = igraph::graph.adjacency(thetas.est.mono[[1]]!=0, mode='undirected')
g.lps2 = igraph::graph.adjacency(thetas.est.mono[[2]]!=0, mode='undirected')
g.lps24 = igraph::graph.adjacency(thetas.est.mono[[3]]!=0, mode='undirected')
g.unstim = igraph::graph.adjacency(thetas.est.mono[[4]]!=0, mode='undirected')

df.degree = data.frame(degree=c(igraph::degree(g.ifn), igraph::degree(g.lps2), igraph::degree(g.lps24), igraph::degree(g.unstim)), 
                       condition=factor(c(rep(condition.names[1], p), rep(condition.names[2], p),rep(condition.names[3], p),rep(condition.names[4], p))))

gg.dens = ggplot2::ggplot(df.degree, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))  +scale_x_continuous(limits = c(0,40))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_jointGHS.pdf"), 8, 8)
gg.dens+theme_bw()
dev.off()

# Save for later
df.degree.jointGHS = df.degree

# Divided up

gg.dens.split = ggplot2::ggplot(df.degree, aes(degree))+geom_density(color='darkgray', fill='darkgray', alpha=0.2)+scale_fill_manual(values='darkgray')+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))  + scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition)+
  geom_vline(data=dplyr::filter(df.degree, condition=="IFN-gamma"), aes(xintercept=degs.cis.ifn[1]+0.2), colour="green",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="IFN-gamma"), aes(xintercept=degs.cis.ifn[2]), colour="blue",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="IFN-gamma"), aes(xintercept=degs.cis.ifn[3]), colour="red",alpha=0.7)+
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-2h"), aes(xintercept=degs.cis.lps2[1]+0.2), colour="green",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-2h"), aes(xintercept=degs.cis.lps2[2]), colour="blue",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-2h"), aes(xintercept=degs.cis.lps2[3]), colour="red",alpha=0.7)  +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-24h"), aes(xintercept=degs.cis.lps24[1]), colour="green",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-24h"), aes(xintercept=degs.cis.lps24[2]), colour="blue",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="LPS-24h"), aes(xintercept=degs.cis.lps24[3]), colour="red",alpha=0.7)+
  geom_vline(data=dplyr::filter(df.degree, condition=="Unstim"), aes(xintercept=degs.cis.unstim[1]), colour="green",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="Unstim"), aes(xintercept=degs.cis.unstim[2]), colour="blue",alpha=0.7) +
  geom_vline(data=dplyr::filter(df.degree, condition=="Unstim"), aes(xintercept=degs.cis.unstim[3]), colour="red",alpha=0.7)    

pdf("Monocytes/plots/degreedistributions/degreedistr_density_cisgenes_CREB1_jointGHS.pdf", 10, 8)
gridExtra::grid.arrange(gg.dens.split+theme(legend.position = 'none')+theme_bw(), legend,ncol=2,widths=c(1,0.1))
dev.off()

# Not with CREB1, final plot

gg.dens.split = ggplot2::ggplot(df.degree, aes(degree))+geom_density(color='darkgray', fill='darkgray', alpha=0.2)+scale_fill_manual(values='darkgray')+
  ggplot2::theme(legend.position = 'none')  + scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition)+
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


pdf("Monocytes/plots/degreedistributions/degreedistr_density_cisgenes_jointGHS.pdf", 10, 8)
gridExtra::grid.arrange(gg.dens.split+theme_bw()+ggplot2::theme(text = element_text(size = 15)), legend ,ncol=2,widths=c(1,0.15))
dev.off()



# Look at clustering coefficient
clustering.coef.all = data.frame(IFN=igraph::transitivity(g.ifn), LPS2=igraph::transitivity(g.lps2), LPS24=igraph::transitivity(g.lps24), Unstim=igraph::transitivity(g.unstim))
names(clustering.coef.all) = condition.names
igraph::transitivity(g.ifn)

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
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))  +scale_x_continuous(limits = c(0,40))
dev.off()

# Flagged vs all 

df.degree.flagvsall = rbind(df.degree.jointGHS, df.degree)
df.degree.flagvsall$genes = c(rep('all', nrow(df.degree.jointGHS)),rep('top hotspot controlled', nrow(df.degree)))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_flaggedVSall_jointGHS.pdf"), 12, 8)
ggplot2::ggplot(df.degree.flagvsall, aes(degree, group=genes, colour=genes,fill=genes))+geom_density(alpha=0.4)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition)+theme_bw()
dev.off()


# Add info to clustering coef df
g.all.hotspot = list(igraph::subgraph(g.ifn,ind.flagget.ifn), igraph::subgraph(g.lps2,ind.flagget.lps2),igraph::subgraph(g.lps24,ind.flagget.lps24),igraph::subgraph(g.unstim,ind.flagget.unstim))
clustering.coef.all[2,] = unlist(lapply(g.all.hotspot, igraph::transitivity))

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

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_single.pdf"), 8, 8)
ggplot2::ggplot(df.degree.single, aes(degree, group=condition, colour=condition))+geom_density()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))
dev.off()


# Compare to joint density 
df.dens.total = rbind(df.degree.jointGHS,df.degree.single)
df.dens.total$method = c(rep('jointGHS', nrow(df.degree.jointGHS)),rep('fastGHS', nrow(df.degree.single)))

pdf("Monocytes/plots/degreedistributions/degreedistr_density_single_vs_joint.pdf", 8, 8)
ggplot2::ggplot(df.dens.total, aes(degree, group=method, colour=method, fill=method))+geom_density(alpha=0.4)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition)
dev.off()

# Only tail
#pdf("Monocytes/plots/degreedistributions/degreedistr_density_single_vs_joint.pdf", 8, 8)
#ggplot2::ggplot(df.dens.total, aes(degree, group=method, colour=method))+geom_histogram()+
 # ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(20,40))+facet_wrap(~ condition)
#dev.off()
# Not really looking great...skip


# Add clustering coef info
clustering.coef.all[3,] = c(igraph::transitivity(g.ifn.single), igraph::transitivity(g.lps2.single), igraph::transitivity(g.lps24.single),igraph::transitivity(g.unstim.single))


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
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))  +scale_x_continuous(limits = c(0,40))
dev.off()


# Single vs joint, only flagged

df.degree.flag.jointVSsingle = rbind(df.degree, df.degree.single)
df.degree.flag.jointVSsingle$method = c(rep('jointGHS', nrow(df.degree)),rep('fastGHS', nrow(df.degree.single)))

pdf(paste0("Monocytes/plots/degreedistributions/degreedistr_density_flagged_single_vs_jointGHS.pdf"), 8, 8)
ggplot2::ggplot(df.degree.flag.jointVSsingle, aes(degree, group=method, colour=method,fill=method))+geom_density(alpha=0.4)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition)+theme_bw()
dev.off()



# Add info to clustering coef df
g.all.hotspot = list(igraph::subgraph(g.ifn.single,ind.flagget.ifn), igraph::subgraph(g.lps2.single,ind.flagget.lps2),igraph::subgraph(g.lps24.single,ind.flagget.lps24),igraph::subgraph(g.unstim.single,ind.flagget.unstim))
clustering.coef.all[4,] = unlist(lapply(g.all.hotspot, igraph::transitivity))
rownames(clustering.coef.all) = c('JointGHS', 'jointGHS top hotspot mediated', 'fastGHS', 'fastGHS top hotspot mediated')

clustering.coef.all

# Make one large combined plot ----------------------------------------------------
load('Monocytes/data/neigh_common_unique_singleVSjoint.RData') # dat.edge.combined

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg.cols = gg_color_hue(2)

median.deg.joint = dplyr::filter(df.degree.flag.jointVSsingle,method=='jointGHS') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Median= median(degree))
median.deg.single = dplyr::filter(df.degree.flag.jointVSsingle,method=='fastGHS') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Median= median(degree))
mean.deg.joint = dplyr::filter(df.degree.flag.jointVSsingle,method=='jointGHS') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Mean= mean(degree))
mean.deg.single = dplyr::filter(df.degree.flag.jointVSsingle,method=='fastGHS') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Mean= mean(degree))

names(df.degree.flag.jointVSsingle) = c(names(df.degree.flag.jointVSsingle)[1:2], 'Method')

gg.comb.flag = ggplot2::ggplot(df.degree.flag.jointVSsingle, aes(degree, group=Method, colour=Method,fill=Method))+geom_density(alpha=0.4)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition, ncol=1)+theme_bw()+
  geom_vline(data=median.deg.joint,aes(xintercept=Median,linetype="median"),color=gg.cols[2],size=0.5)+
  geom_vline(data=median.deg.single,aes(xintercept=Median,linetype="median"),color=gg.cols[1],size=0.5)+
  geom_vline(data=mean.deg.joint,aes(xintercept=Mean,linetype="mean"),color=gg.cols[2],size=0.5)+
  geom_vline(data=mean.deg.single,aes(xintercept=Mean,linetype="mean"),color=gg.cols[1],size=0.5)+
  scale_linetype_manual(name='Statistic', values=c(mean='solid',median='dashed'))+theme(legend.position = 'bottom')
  
gg.comb.order = ggplot(dat.edge.combined, aes(x=method, fill=Neighbour))+geom_bar(position='fill')+ guides(fill=guide_legend(title="Edge order"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(5,'YlOrRd')[5:1])+facet_wrap(~ Condition, ncol=1)+ylab('Fraction of edges')+theme_bw()+theme(legend.position = 'bottom')

pdf("Monocytes/plots/degreedistributions/combinedplot_single_vs_jointGHS.pdf", 13, 15)
(gg.comb.flag+ggplot2::theme(text = element_text(size = 15)))+(gg.comb.order+ggplot2::theme(text = element_text(size = 15)))
dev.off()

# Including the non-flagged genes in the density plot as well

ind.other.ifn = which(! genes_id %in% genes.flagged[[1]])
ind.other.lps2 = which(! genes_id %in% genes.flagged[[2]])
ind.other.lps24 = which(! genes_id %in% genes.flagged[[3]])
ind.other.unstim = which(! genes_id %in% genes.flagged[[4]])

df.degree.other = data.frame(degree=c(igraph::degree(g.ifn)[ind.other.ifn], igraph::degree(g.lps2)[ind.other.lps2], igraph::degree(g.lps24)[ind.other.lps24], igraph::degree(g.unstim)[ind.other.unstim]), 
                       condition=factor(c(rep(condition.names[1], length(ind.other.ifn)), rep(condition.names[2], length(ind.other.lps2)),rep(condition.names[3], length(ind.other.lps24)),rep(condition.names[4], length(ind.other.unstim)))))
df.degree.other$Method=rep('Not controlled', nrow(df.degree.other))

df.degree.flag.jointVSsingle.other = rbind(df.degree.flag.jointVSsingle, df.degree.other)

median.deg.other = dplyr::filter(df.degree.flag.jointVSsingle.other,Method=='Not controlled') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Median= median(degree))
mean.deg.other = dplyr::filter(df.degree.flag.jointVSsingle.other,Method=='Not controlled') %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(Mean= mean(degree))


gg.comb.flag.other = ggplot2::ggplot(df.degree.flag.jointVSsingle.other, aes(degree, group=Method, colour=Method,fill=Method))+geom_density(alpha=0.4)+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition, ncol=1)+theme_bw()+
  geom_vline(data=median.deg.joint,aes(xintercept=Median,linetype="median"),color=gg.cols[2],size=0.5)+
  geom_vline(data=median.deg.single,aes(xintercept=Median,linetype="median"),color=gg.cols[1],size=0.5)+
  geom_vline(data=mean.deg.joint,aes(xintercept=Mean,linetype="mean"),color=gg.cols[2],size=0.5)+
  geom_vline(data=mean.deg.single,aes(xintercept=Mean,linetype="mean"),color=gg.cols[1],size=0.5)+
  scale_linetype_manual(name='Statistic', values=c(mean='solid',median='dashed'))+theme(legend.position = 'bottom')+
  scale_color_manual(values=c(gg.cols,'grey'))+scale_fill_manual(values=c(gg.cols,'grey'))+ theme(legend.title = element_blank())
  

pdf("Monocytes/plots/degreedistributions/combinedplot_with_unflagged_single_vs_jointGHS.pdf", 13, 15)
(gg.comb.flag.other+ggplot2::theme(text = element_text(size = 15)))+(gg.comb.order+ggplot2::theme(text = element_text(size = 15)))
dev.off()

# Add mean and median for unflagged too
pdf("Monocytes/plots/degreedistributions/combinedplot_with_unflagged_allLines_single_vs_jointGHS.pdf", 13, 15)
(gg.comb.flag.other+
  geom_vline(data=mean.deg.other,aes(xintercept=Mean,linetype="mean"),color='grey',size=0.5)+
    geom_vline(data=median.deg.other,aes(xintercept=Median,linetype="median"),color='grey',size=0.5)+
    ggplot2::theme(text = element_text(size = 15)))+(gg.comb.order+ggplot2::theme(text = element_text(size = 15)))
dev.off()

# Also with LYZ and YEATS4 degree marked
df.degree.jointGHS.cis = df.degree.jointGHS
df.degree.jointGHS.cis$Method = rep('jointGHS', nrow(df.degree.jointGHS.cis))

gg.comb.flag.other.withcis = gg.comb.flag.other+
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="IFN-gamma"), aes(x=degs.cis.ifn[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="IFN-gamma"), aes(x=degs.cis.ifn[2],y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-2h"), aes(x=degs.cis.lps2[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-2h"), aes(x=degs.cis.lps2[2]+0.2,y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-24h"), aes(x=degs.cis.lps24[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-24h"), aes(x=degs.cis.lps24[2],y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="Unstim"), aes(x=degs.cis.unstim[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="Unstim"), aes(x=degs.cis.unstim[2],y=0), colour="blue",alpha=0.7, show.legend = F)

pdf("Monocytes/plots/degreedistributions/combinedplot_with_unflagged_cismarked_single_vs_jointGHS.pdf", 13, 15)
((gg.comb.flag.other.withcis+ggplot2::theme(text = element_text(size = 15))))+(gg.comb.order+ggplot2::theme(text = element_text(size = 15)))
dev.off()


# And finally, without the unflagged density but with the cis genes degree marked

gg.comb.flag.withcis = gg.comb.flag+
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="IFN-gamma"), aes(x=degs.cis.ifn[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="IFN-gamma"), aes(x=degs.cis.ifn[2],y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-2h"), aes(x=degs.cis.lps2[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-2h"), aes(x=degs.cis.lps2[2]+0.2,y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-24h"), aes(x=degs.cis.lps24[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="LPS-24h"), aes(x=degs.cis.lps24[2],y=0), colour="blue",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="Unstim"), aes(x=degs.cis.unstim[1],y=0), colour="orange",alpha=0.7, show.legend = F) +
  geom_point(data=dplyr::filter(df.degree.jointGHS.cis, condition=="Unstim"), aes(x=degs.cis.unstim[2],y=0), colour="blue",alpha=0.7, show.legend = F)

pdf("Monocytes/plots/degreedistributions/combinedplot_cismarked_single_vs_jointGHS.pdf", 13, 15)
((gg.comb.flag.withcis+ggplot2::theme(text = element_text(size = 15))))+(gg.comb.order+ggplot2::theme(text = element_text(size = 15)))
dev.off()



# Same with ALL edges common, not compliment
load('Monocytes/data/neigh_common_all_singleVSjoint.RData')

gg.comb.order.all = ggplot(dat.edge.combined.all, aes(x=method, fill=Neighbour))+geom_bar(position='fill')+ guides(fill=guide_legend(title="Edge order"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(7,'YlOrRd')[7:1])+facet_wrap(~ Condition, ncol=1)+ylab('Fraction of edges')+theme_bw()+theme(legend.position = 'bottom')

pdf("Monocytes/plots/degreedistributions/combinedplot_all_single_vs_jointGHS.pdf", 13, 15)
(gg.comb.flag+ggplot2::theme(text = element_text(size = 15)))+(gg.comb.order.all+ggplot2::theme(text = element_text(size = 15)))
dev.off()

# Also using count

gg.comb.order.all.count = ggplot(dat.edge.combined.all, aes(x=method, fill=Neighbour))+geom_histogram(stat='count')+ guides(fill=guide_legend(title="Edge order"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(7,'YlOrRd')[7:1])+facet_wrap(~ Condition, ncol=1)+ylab('Number of edges')+theme_bw()+theme(legend.position = 'bottom')


pdf("Monocytes/plots/degreedistributions/combinedplot_count_all_single_vs_jointGHS.pdf", 13, 15)
(gg.comb.flag+ggplot2::theme(text = element_text(size = 15)))+(gg.comb.order.all.count+ggplot2::theme(text = element_text(size = 15)))
dev.off()





# only joint, adding cis degrees (already did this in another fig, so not neccesary here...)

#median.deg.joint = dplyr::filter(df.dens.total,method=='jointGHS') %>%
#  dplyr::group_by(condition) %>%
#  dplyr::summarize(Median= median(degree))
#mean.deg.joint = dplyr::filter(df.dens.total,method=='jointGHS') %>%
#  dplyr::group_by(condition) %>%
#  dplyr::summarize(Mean= mean(degree))#

#USE df.degree.jointGHS
#names(df.degree.jointGHS) = c(names(df.degree.jointGHS)[1:2],'Method')
#ggplot2::ggplot(df.dens.total, aes(degree, group=Method, colour=Method,fill=Method))+geom_density(alpha=0.4)+
#  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +scale_x_continuous(limits = c(0,40))+facet_wrap(~ condition, ncol=1)+theme_bw()+
#  geom_vline(data=median.deg.joint,aes(xintercept=Median,linetype="median"),color=gg.cols[2],size=0.5)+
#  geom_vline(data=median.deg.single,aes(xintercept=Median,linetype="median"),color=gg.cols[1],size=0.5)+
#  geom_vline(data=mean.deg.joint,aes(xintercept=Mean,linetype="mean"),color=gg.cols[2],size=0.5)+
#  geom_vline(data=mean.deg.single,aes(xintercept=Mean,linetype="mean"),color=gg.cols[1],size=0.5)+
#  scale_linetype_manual(name='Statistic', values=c(mean='solid',median='dashed'))+theme(legend.position = 'bottom')





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
#IFN:  LYZ :  ADAMTS4 C1orf210 C21orf24 C21orf55 C5orf28 
#LPS2h:  LYZ :  ADAMTS4 C19orf12 C1orf151 C21orf24 
#LPS24h:  LYZ :  ACBD5 C19orf31 C1QA 
#LPS2h:  LYZ :  ADAMTS4 C15orf63 C19orf31 C1orf151 
#IFN:  YEATS4 :  AFMID AGTRAP APOBEC3D C19orf12 C1orf210 C9orf116 CATSPER2 CCR6 CDK5RAP2 
#LPS2h:  YEATS4 :  AFMID AGTRAP C12orf43 C19orf12 C4orf34 CAPS2 CATSPER2 
#LPS24h:  YEATS4 :  ADAMTS4 AFMID ANKRD30B C15orf63 C19orf31 C4orf34 CAPS2 CATSPER2 CCNB1IP1 
#LPS2h:  YEATS4 :  AFMID AGTRAP APOBEC3D BMS1P5 C15orf63 C5orf28 C9orf80 CAPS2 CCDC125 
#IFN:  CREB1 :  ACBD5 
#LPS2h:  CREB1 :  ACBD5 BATF3 
#LPS24h:  CREB1 :   
#  LPS2h:  CREB1 :  ACBD5 

# LYZ has multiple top hubs as neighbours in each condition.
# YEATS4 has the most.

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


# LYZ 
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


# Which are common to all conditions? LYZ
genes_id[which(genes_id %in% degree.list.ifn[[1]]$gene &  genes_id %in% degree.list.lps2[[1]]$gene & genes_id %in% degree.list.lps24[[1]]$gene & genes_id %in% degree.list.unstim[[1]]$gene)]
# only ifn
genes_id[which(genes_id %in% degree.list.ifn[[1]]$gene &  ! genes_id %in% degree.list.lps2[[1]]$gene & ! genes_id %in% degree.list.lps24[[1]]$gene & ! genes_id %in% degree.list.unstim[[1]]$gene)]
# only lps2
genes_id[which(! genes_id %in% degree.list.ifn[[1]]$gene &  genes_id %in% degree.list.lps2[[1]]$gene & ! genes_id %in% degree.list.lps24[[1]]$gene & ! genes_id %in% degree.list.unstim[[1]]$gene)]
# only lps24
genes_id[which(! genes_id %in% degree.list.ifn[[1]]$gene &  !genes_id %in% degree.list.lps2[[1]]$gene &  genes_id %in% degree.list.lps24[[1]]$gene & ! genes_id %in% degree.list.unstim[[1]]$gene)]
# only unstim
genes_id[which(! genes_id %in% degree.list.ifn[[1]]$gene & ! genes_id %in% degree.list.lps2[[1]]$gene & ! genes_id %in% degree.list.lps24[[1]]$gene &  genes_id %in% degree.list.unstim[[1]]$gene)]

# Which are common to all conditions? YEATS4
genes_id[which(genes_id %in% degree.list.ifn[[2]]$gene &  genes_id %in% degree.list.lps2[[2]]$gene & genes_id %in% degree.list.lps24[[2]]$gene & genes_id %in% degree.list.unstim[[2]]$gene)]
# only ifn
genes_id[which(genes_id %in% degree.list.ifn[[2]]$gene &  ! genes_id %in% degree.list.lps2[[2]]$gene & ! genes_id %in% degree.list.lps24[[2]]$gene & ! genes_id %in% degree.list.unstim[[2]]$gene)]
# only lps2
genes_id[which(! genes_id %in% degree.list.ifn[[2]]$gene &  genes_id %in% degree.list.lps2[[2]]$gene & ! genes_id %in% degree.list.lps24[[2]]$gene & ! genes_id %in% degree.list.unstim[[2]]$gene)]
# only lps24
genes_id[which(! genes_id %in% degree.list.ifn[[2]]$gene &  !genes_id %in% degree.list.lps2[[2]]$gene &  genes_id %in% degree.list.lps24[[2]]$gene & ! genes_id %in% degree.list.unstim[[2]]$gene)]
# only unstim
genes_id[which(! genes_id %in% degree.list.ifn[[2]]$gene & ! genes_id %in% degree.list.lps2[[2]]$gene & ! genes_id %in% degree.list.lps24[[2]]$gene &  genes_id %in% degree.list.unstim[[2]]$gene)]


# Which are controlled by the top hotspot? 

# LYZ
degree.list.ifn[[1]]$gene[which(degree.list.ifn[[1]]$gene %in% genes.flagged[[1]])]
degree.list.lps2[[1]]$gene[which(degree.list.lps2[[1]]$gene %in% genes.flagged[[2]])]
degree.list.lps24[[1]]$gene[which(degree.list.lps24[[1]]$gene %in% genes.flagged[[3]])]
degree.list.unstim[[1]]$gene[which(degree.list.unstim[[1]]$gene %in% genes.flagged[[4]])]

# YEATS4
degree.list.ifn[[2]]$gene[which(degree.list.ifn[[2]]$gene %in% genes.flagged[[1]])]
degree.list.lps2[[2]]$gene[which(degree.list.lps2[[2]]$gene %in% genes.flagged[[2]])]
degree.list.lps24[[2]]$gene[which(degree.list.lps24[[2]]$gene %in% genes.flagged[[3]])]
degree.list.unstim[[2]]$gene[which(degree.list.unstim[[2]]$gene %in% genes.flagged[[4]])]



# Quick look at top degree neighbours and their partial correlations
ind.a = which(genes_id == 'AFMID')
hist(cov2cor(thetas.est.mono[[1]])[ind.a,][-which(thetas.est.mono[[1]][ind.a,]%in% c(0,1))], breaks=100)

ind.a = which(genes_id == 'AFTPH')
hist(cov2cor(thetas.est.mono[[1]])[ind.a,][-which(thetas.est.mono[[1]][ind.a,]%in% c(0,1))], breaks=100)


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


# Print for table
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


# What are their hotspot statuses?
df.deg.ifn.ordered.top$genes[which(df.deg.ifn.ordered.top$genes %in% genes.flagged[[1]])]
df.deg.lps2.ordered.top$genes[which(df.deg.lps2.ordered.top$genes %in% genes.flagged[[2]])]
df.deg.lps24.ordered.top$genes[which(df.deg.lps24.ordered.top$genes %in% genes.flagged[[3]])]
df.deg.unstim.ordered.top$genes[which(df.deg.unstim.ordered.top$genes %in% genes.flagged[[4]])]


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
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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
                             edge.color = 'color',edge.lty = "lty",legend.size = 20)+
    geom_text(aes(label =  names.net), size=4, hjust=h.just, vjust=v.just)+
    ggplot2::ggtitle(condition.names[i]) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

legend.ordered = get_legend(nets[[1]]+ theme(legend.position="right")+ guides(size ='none'))

pdf("Monocytes/plots/CIS_net_circle_jointGHS.pdf", 15, 12)
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]]+ theme(legend.position="none",text = element_text(size = 20)),nets[[2]]+ theme(legend.position="none",text = element_text(size = 20)), 
                                               nets[[3]]+ theme(legend.position="none",text = element_text(size = 20)), nets[[4]]+ theme(legend.position="none",text = element_text(size =20)), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()



# Plot all edges ----------------------------------------------------------------- 


nets = list()
for(i in 1:length(unique.list.ordered)){
  net.all =  unique.list.ordered[[i]]
  net = network::network(net.all,directed=F, ignore.eval=F,names.eval='weights')
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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
gridExtra::grid.arrange(gridExtra::arrangeGrob(nets[[1]] + theme(legend.position="none"),nets[[2]]+ theme(legend.position="none"), 
                                               nets[[3]] + theme(legend.position="none"), nets[[4]]+ theme(legend.position="none"), ncol=2),
                        legend.ordered, ncol=2, widths=c(8,1))
dev.off()

# Plot shared edges ------------------------------------------------------------

nets = list()
for(i in 1:length(unique.list.ordered)){
  net.all =  unique.list.ordered[[i]]
  net.all[which(abs(net.all)!=3)]=0
  net = network::network(net.all,directed=F, ignore.eval=F,names.eval='weights')
  network::set.edge.attribute(net, "lty", ifelse(sign(net %e% "weights") == 1, 1, 2))
  network::set.edge.attribute(net, "color", c("black", "grey75","red", "blue")[abs(net %e% "weights")+1])
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


# Write list of edges to file: -------------------------------------------

# We now write a list of all edges indentified by jointGHS to a .csv file.

get_and_print_edges <- function(a.mat, col.names, theta.mat=NULL) {
  # Function for printing all the edges in a graph with a layout that can be inserted into a latex table.
  # Also returns a data frame containing the edges
  # a.mat:          the adjacency matrix
  # col.names:      the names of the nodes in the graph
  # theta.mat:      the precision matrix. If included, the size the partial correlations are included in the table as well
  a.mat[which(diag(rep(1, ncol(a.mat))) == 1, arr.ind = T)] <- 0 # Make diagonal zero
  pairs <- which(a.mat[, ] == 1, arr.ind = T)
  df <- data.frame(t(apply(pairs, 1, sort))) # Sort so that the node in the pair whose name is first in the alphabet is first.
  df <- unique(df)
  names <- cbind(col.names[df[, 1]], col.names[df[, 2]])
  if (!is.null(theta.mat)){
    effect = round(-cov2cor(theta.mat)[cbind(df[,1], df[,2])], 5)
    names = cbind(names, effect)
    return(names)
  }
  for (i in 1:nrow(names)) {
    cat(names[i, 1], " & ", names[i, 2], " \\\\ \n")
  }
  return(names)
}

edges.ifn = get_and_print_edges(thetas.est.mono[[1]]!=0,genes_id, thetas.est.mono[[1]])
edges.lps2 = get_and_print_edges(thetas.est.mono[[2]]!=0,genes_id, thetas.est.mono[[2]])
edges.lps24 = get_and_print_edges(thetas.est.mono[[3]]!=0,genes_id, thetas.est.mono[[3]])
edges.unstim = get_and_print_edges(thetas.est.mono[[4]]!=0,genes_id, thetas.est.mono[[4]])
colnames(edges.ifn) <- c("Gene1", "Gene2", "PartialCor")
colnames(edges.lps2) <- c("Gene1", "Gene2", "PartialCor")
colnames(edges.lps24) <- c("Gene1", "Gene2", "PartialCor")
colnames(edges.unstim) <- c("Gene1", "Gene2", "PartialCor")

# Write to file
write.csv(edges.ifn, file = "Monocytes/edge_lists/edgesIFNg.csv", row.names = F,quote=F)
write.csv(edges.lps2, file = "Monocytes/edge_lists/edgesLps2h.csv", row.names = F,quote=F)
write.csv(edges.lps24, file = "Monocytes/edge_lists/edgesLps24h.csv", row.names = F,quote=F)
write.csv(edges.unstim, file = "Monocytes/edge_lists/edgesUnstim.csv", row.names = F,quote=F)

# Write to file, pooled
edges.comb = rbind(edges.ifn, edges.lps2, edges.lps24, edges.unstim)
edges.comb = cbind(edges.comb,c(rep(condition.names[1], nrow(edges.ifn)), rep(condition.names[2], nrow(edges.lps2)), rep(condition.names[3], nrow(edges.lps24)),
                   rep(condition.names[4], nrow(edges.unstim))))
colnames(edges.comb) <- c("Gene1", "Gene2", "PartialCor", "Condition")
write.csv(edges.comb, file = "Monocytes/edge_lists/edgesComb.csv", row.names = F,quote=F)


# Write list of node degree for all genes -------------------------------------------

write.csv(df.degree.ifn, file = "Monocytes/Edge_lists/degreesIFNg.csv", row.names = T,quote=F)
write.csv(df.degree.lps2, file = "Monocytes/Edge_lists/degreesLps2h.csv", row.names = T,quote=F)
write.csv(df.degree.lps24, file = "Monocytes/Edge_lists/degreesLps24h.csv", row.names = T,quote=F)
write.csv(df.degree.unstim, file = "Monocytes/Edge_lists/degreesUnstim.csv", row.names = T,quote=F)

# Write to file for all conditions, pooled
df.degree.comb = rbind(df.degree.ifn, df.degree.lps2, df.degree.lps24, df.degree.unstim)
df.degree.comb = cbind(df.degree.comb,c(rep(condition.names[1], p), rep(condition.names[2], p), rep(condition.names[3], p),
                                        rep(condition.names[4], p)))
colnames(df.degree.comb) <- c("Gene", "Degree", "Condition")
write.csv(df.degree.comb, file = "Monocytes/Edge_lists/degreesComb.csv", row.names = T,quote=F)

# Not in long format

df.degree.wide = cbind(sort(df.degree.ifn), sort(df.degree.lps2)$Degree,sort(df.degree.lps24)$Degree,sort(df.degree.unstim)$Degree)
colnames(df.degree.wide) = c('Gene', condition.names)
write.csv(df.degree.wide, file = "Monocytes/Edge_lists/degreesCombWide.csv", row.names = F,quote=F)


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


