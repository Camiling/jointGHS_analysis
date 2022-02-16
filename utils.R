library(jointGHS)
library(fastGHS)
library(foreach)
library(tailoredGlasso)
library(superheat)
library(GGally)
library(sna)
library(gridExtra)
library(ggplot2)

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
