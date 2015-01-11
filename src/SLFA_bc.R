library("QUIC")
setwd("~/bioNet")
X = read.table("./test-data/X.txt")
Y = read.table("./test-data/Y.txt")
X = as.matrix(X)
Y = as.matrix(Y)
k = 100

sel_bc = which(Y[, 2] == -1)
sel_normal = which(Y[, 2] == 1)
X1 = X[sel_bc, ]
means <- colMeans(X1)
for (i in 1:length(means)){
    X1[,i] <- X1[,i]-means[i] 
}
X2 = X[sel_normal, ]
means <- colMeans(X2)
for (i in 1:length(means)){
    X2[,i] <- X2[,i]-means[i] 
}

    #initialize B, S, Phi and theta
    B = (eigen(cov(X)))$vectors[,1:k]
    B1= B
    B = B1
    X1 = t(X1)
    X2 = t(X2)
    S1 = solve((t(B) %*% B), t(B) %*% X1)
    S2 = solve((t(B) %*% B), t(B) %*% X2)
    theta = 10^(-2)
    n1 = ncol(X1)
    n2 = ncol(X2)
    Phi1 = QUIC(S1 %*% t(S1) / n1, rho)$W
    Phi2 = QUIC(S2 %*% t(S2) / n2, rho)$W
    Phi_last1 = diag(0, k)
    Phi_last2 = diag(0, k)
    sigma = 0.1
    rho = 1
    
    #interate until stopping
    while((norm((Phi1-Phi_last1), 'F') > theta) && (norm((Phi2-Phi_last2), 'F') > theta)){
        
        Phi_last1 = Phi1
        Phi_last2 = Phi2
        
        S1 = solve((t(B) %*% B + sigma * sigma * Phi1), t(B) %*% X1)
        S2 = solve((t(B) %*% B + sigma * sigma * Phi2), t(B) %*% X2)
        
        B = t(solve((S1 %*% t(S1)) / n1 + (S2 %*% t(S2)) / n2, (S1 %*% t(X1)) / n1 + (S2 %*% t(X2)) / n2))
        Phi1 = QUIC(S1 %*% t(S1) / n1, rho)$W
        Phi2 = QUIC(S2 %*% t(S2) / n2, rho)$W
            
        for(i in 1:k){
            B[,i] = B[,i] / norm(B[,i], '2')
        }
    }