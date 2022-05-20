gmm2s_iter <- function(Ti, N, y, X, Z, yy, k_t,tol, guess = -99, maxiter, second) {  #tol, guess, maxiter optional
  
  
  K <- NCOL(Z)
  num_var <- NCOL(X)
  #print(dim(z))
  dof <- N*NCOL(Z) - NCOL(X) # ncol gives number of columns of a dataframe or matrix
  #NOT SURE IF THIS WORKS CORRECTLY FOR K>1
  
  
  bigZ  <- kronecker(matrix(1,1,N),Z)
  zero_out <- kronecker(diag(N),matrix(1,Ti,K)) 
  bigZ <- hadamard.prod(bigZ,zero_out) #element by element multiplication of matrices
  aux_yz2 <- t(bigZ)%*%y
  aux_xz2 <- t(bigZ)%*%X
  
  Omega_hat <- matrix(0,N*K,N*K)
  Omega_u_hat <- matrix(0,N*K,N*K)
    
  W_0 <- diag(N*K)
  phi_1st <- base::solve(t(aux_xz2)%*%W_0%*%aux_xz2) %*% t(aux_xz2)%*%W_0%*%aux_yz2 
  
  phi_hat <- phi_1st
  phi_guess = phi_1st
  
  iter <- 0
  
  while (iter < maxiter) {
  
    iter <- iter + 1 
    u_hat <- y - X%*%phi_guess 
  
    Omega_hat_temp <- matrix(0,N*K,N*K)
    Omega_u_hat_temp <- matrix(0,N*K,N*K)
  
    u2 <- Reshape(u_hat,Ti,N)
    bigU <- kronecker(u2,matrix(1,1,K))
  
    z_re <- matrix(0,Ti,N*K)
    if (N>1) {
      for (jj in 1:N) {
        z_re[1:Ti,K*(jj-1)+1:K] <- Z[Ti*(jj-1)+(1:Ti),1:K]
      } 
    } else {
      z_re[1:Ti,1:K] <- Z
    }
    mom_re  <- hadamard.prod(z_re,bigU)
    Omega_hat <- matrix(0,N*K,N*K)
    #Omega_hat <- (t(mom_re)%*%mom_re)/Ti
  

    for (k in -yy:yy) { # yy to 0
      start <- max(1,1-k)
      finish <- min(Ti,Ti-k)
    
      w <- 1 - abs(k)/(yy+1)
    
      if (second == 0) {
        Omega_hat_temp <- (1/Ti)*w*(t(mom_re[start:finish,])%*%mom_re[(start+k):(finish+k),]) + Omega_hat_temp 
    } else { 
        Omega_U <-(1/Ti)*t(bigU[start:finish,])%*%bigU[(start+k):(finish+k),]
        Omega_Z <- (1/Ti)*t(z_re[start:finish,])%*%z_re[(start+k):(finish+k),]
        Omega_hat_temp <- w*(hadamard.prod(Omega_U,Omega_Z)) + Omega_hat_temp
      }
    }
  
    W <- base::solve(Omega_hat_temp)
    #W <- diag(N*K)*W
    phi_hat <- base::solve(t(aux_xz2)%*%W%*%aux_xz2) %*% t(aux_xz2)%*%W%*%aux_yz2
    G <- (1/Ti)*aux_xz2
    temp <- base::solve(t(G)%*%W%*%G)
    #var <- (1/Ti)*temp%*%t(G)%*%W%*%Omega_hat%*%t(W)%*%G%*%temp
    var <- (1/Ti)*base::solve((1/Ti)*t(aux_xz2)%*%W%*%aux_xz2*(1/Ti))
    phi_guess <- phi_hat
    Omega_hat <- Omega_hat_temp
    Omega_u_hat <- Omega_u_hat_temp
  
  }

  m <- (1/sqrt(Ti))*t(mom_re)%*%matrix(1,Ti,1)
  #J <- t(m)%*%solve(Omega_hat)%*%m
  J <- t(m)%*%W%*%m
  
  pval = pchisq(J, df=dof, lower.tail=FALSE)
  
  t = phi_hat/sqrt(diag(var))
  p = 2*pt(-abs(t),1000000)
  
  return_list <- list("phi_hat"=phi_hat,"phi_1st"=phi_1st,"Var"=var,"J"=J,"pval"=pval,
                      "nobs"=(Ti*N),"p"=p,"Omega_hat"=Omega_hat,"W"=W)
  return(return_list)
  
}