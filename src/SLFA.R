SLFA = function(X, k = 10, sigma = 0.1){
    
    #initialize B, S, Phi and theta
    B = (princomp(X)$loadings)[1:k]
    S = solve((t(B) %*% B), t(B) %*% X)
    theta = 10^(-4)
    Phi = diag(1, k)
    Phi_last = diag(0, k)
    
    #interate until stopping
    while(norm((Phi-Phi_last), 'F') < theta){
        Phi_last = Phi
        S = solve((t(B) %*% B + sigma * sigma * Phi), t(B) %*% X)
        B = t(solve(S %*% t(S), S %*% t(X)))
        Phi = QUIC(S %*% t(S))
        for(i in 1:k){
            B[,i] = B[,i] / norm(B[,i], '2')
        }
    }
    
}