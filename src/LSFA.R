#----------------------------------------#
# Name: Beilun Wang, Date: 06/18/2014
# This function is used for implementing SLFA 
# 
# Please refer 

# # @INPROCEEDINGS{yhe12NIPS,
  # title={Learning the Dependency Structure of Latent Factors},
  # author={Y. He and Y. Qi and K. Kavukcuoglu and H. Park},
  # booktitle={Proceedings of Advances in Neural Information Processing Systems (NIPS)},
  # year={2012},
# }

# 
#----------------------------------------#
LSFA <- function(X, q = 100, lambda = 0.01, theta = sqrt(0.1)){
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
        Phi = glasso(Z %*% t(Z), lambda)$wi
        B = t(solve((S %*% t(S)), (S%*%t(X))))
        U = U + S - Z
        #print(norm(Z - Z_last,'F'))
    }
    result <- list()
    result[[1]] <- S
    result[[2]] <- Phi
    result[[3]] <- B
    return(result)
}
