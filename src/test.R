setwd("~/SLFA")
source("./src/Input.R")
source("./src/LSFA.R")
S = LSFA(X, q , lambda , theta )[[1]]


#####################test case for breast cancer data##########################
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

prin = svd(X%*%t(X))
S2 = solve(prin$u,X)
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

