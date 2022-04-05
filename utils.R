library(jointGHS)
library(fastGHS)
library(foreach)
library(tailoredGlasso)
library(superheat)
library(GGally)
library(sna)
library(gridExtra)
library(ggplot2)
library(ComplexUpset)
library(patchwork)
library(igraph)
library(network)
library(ggalluvial)

# Find Matthews correlation coefficient for estimated graphs
MCC = function(g,g.hat){
  p = nrow(g[,])
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat) = rep(0,p)
  tp = sum(g.hat ==1 & g ==1)/10 # True positives. Divide by 10 to avoid integer overflow.
  fp = sum(g.hat ==1 & g ==0)/10 # False positives
  tn = (sum(g.hat == 0 & g == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat == 0 & g == 1)/10 # False negatives
  return((tp*tn - fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
}

gaussianloglik_util = function (sample.cov, theta, n) { # Avoid overflow
  if (mean(dim(sample.cov) == dim(theta)) != 1) 
    stop("matrices must have the same dimension")
  if (det(theta) <= 0) 
    stop("precision matrix must be positive definite.")
  if (!isSymmetric(sample.cov)) 
    stop("sample covariance matrix must be symmetric")
  if (n <= 0) 
    stop("number of observations n must be positive")
  p <- nrow(theta)
  return(-p * n * log(2 * pi)/2 + n * log(det(theta/5))/2 + n*log(5^p)/2 - n * 
           sum(diag(sample.cov %*% theta))/2)
}