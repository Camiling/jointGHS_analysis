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

# Write list of node degree for all genes -------------------------------------------

write.csv(df.degree.ifn, file = "Monocytes/Edge_lists/degreesIFNg.csv", row.names = T,quote=F)
write.csv(df.degree.lps2, file = "Monocytes/Edge_lists/degreesLps2h.csv", row.names = T,quote=F)
write.csv(df.degree.lps24, file = "Monocytes/Edge_lists/degreesLps24h.csv", row.names = T,quote=F)
write.csv(df.degree.unstim, file = "Monocytes/Edge_lists/degreesUnstim.csv", row.names = T,quote=F)
