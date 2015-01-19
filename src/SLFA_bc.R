#----------------------------------------#
# Name: Beilun Wang, Date: 06/18/2014
# 
# For paper: 

# @INPROCEEDINGS{yhe12NIPS,
#   title={Learning the Dependency Structure of Latent Factors},
#   author={Y. He and Y. Qi and K. Kavukcuoglu and H. Park},
#   booktitle={Proceedings of Advances in Neural Information Processing Systems (NIPS)},
#   year={2012},
# }

#----------------------------------------#

setwd("~/bioNet")
X = read.table("./test-data/bc/X.txt")
Y = read.table("./test-data/bc/Y.txt")
X = as.matrix(X)
Y = as.matrix(Y)
X = apply(X,2,function(x) x - mean(x))
q = 100
rho = 0.1
lambda = 0.01
theta = sqrt(0.1)
#LSFA <- function(X, q, lambda, theta){
    n = nrow(X)
    p = ncol(X)
    X = t(X)
    B = diag(1,p,q)
    S = matrix(0,q,n)
    Z = matrix(0.1,q,n)
    U = matrix(0.1,q,n)
    Phi = diag(1,q,q)
    Z_last = Z
    while((norm((Z-S),'F') > 10^(-4) * sqrt(q*n) + 10^(-4) * max(norm(Z,'F'),norm(S,'F'))) || (norm((Z-Z_last),'F') > 10^(-4) * sqrt(q*n) + 10^(-4) * norm(Z,'F'))){
        Z_last = Z
        S = solve(((1/n)*t(B)%*%B + (rho/2)*diag(1,q,q)),(1/n) * t(B) %*% X + (rho/2)*(Z - U))
        Z = solve(((theta^2/n)*Phi + rho/2*diag(1,q,q)),((rho/2)*(S + U)))
        #re.clime <- clime(Z %*% t(Z), standardize=FALSE, linsolver="simplex")
        #re.cv <- cv.clime(re.clime)
        #Phi <-  clime(Z %*% t(Z), standardize=FALSE, re.cv$lambdaopt)$Omega[[1]]
        glasso(Z %*% t(Z), lambda)
        B = t(solve((S %*% t(S)), (S%*%t(X))))
        U = U + S - Z
        print(norm(Z - Z_last,'F'))
    }
#}
error_rate = rep(0,50)
for(i in 1:50){
    index = 1 : n
    testindex = sample(index, trunc(length(index)/3))
    testset = t(S)[testindex,] 
    trainset = t(S)[-testindex,]
    svm.model = svm( Y[-testindex,2],x = trainset, cost = 100, kernel = 'linear',cross = 10)
    svm.pre = predict(svm.model,testset)
    error = (Y[testindex,2] > 0) - (as.vector(svm.pre) > 0)
    error_rate[i] = sum(abs(error))/98
}

error_rate_PCA = rep(0,50)
for(i in 1:50){
    index = 1 : n
    testindex = sample(index, trunc(length(index)/3))
    testset = t(S2)[testindex,] 
    trainset = t(S2)[-testindex,]
    svm.model = svm( Y[-testindex,2],x = trainset, cost = 100, kernel = 'linear',cross = 10)
    svm.pre = predict(svm.model,testset)
    error = (Y[testindex,2] > 0) - (as.vector(svm.pre) > 0)
    error_rate_PCA[i] = sum(abs(error))/98
}
prin = svd(X%*%t(X))
S_PCA = solve(prin$u,X)

error_rate_o = rep(0,50)
for(i in 1:50){
    index = 1 : n
    testindex = sample(index, trunc(length(index)/3))
    testset = t(X)[testindex, ] 
    trainset = t(X)[-testindex, ]
    svm.model = svm( Y[-testindex, 2], x = trainset, cost = 100, kernel = 'linear',cross = 10)
    svm.pre = predict(svm.model, testset)
    error = (Y[testindex, 2] > 0) - (as.vector(svm.pre) > 0)
    error_rate_o[i] = sum(abs(error))/98
}
