library("clime")
library("e1071")
library("glasso")
setwd("~/SLFA")
X = read.table("./data/bc/X.txt")
Y = read.table("./data/bc/Y.txt")
X = as.matrix(X)
Y = as.matrix(Y)
X = apply(X,2,function(x) x - mean(x))
q = 100
rho = 0.1
lambda = 0.01
theta = sqrt(0.1)