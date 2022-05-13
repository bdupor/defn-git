anteAgg <- function(fips, Ti, N, y, shock, yy, k_t,tol, guess = -99, maxiter, second) {  #tol, guess, maxiter optional
  
  #Hd is the weighted vector for the group-level coefficients
  #Common parameter sectors must be ordered first of all sectors in design matrix
  Hd <- matrix(1,1,N)
  
  #create matrix of possible combinations of subset differences
  some_diff <- matrix(,1,6)
  colnames(some_diff) <- c("est","se","ci_lo","ci_hi","pval_d","J_d")
  
  #run through all possible combinations of different coefficient subsets
  for (ii in 1:NROW(comb)) {
    Fips = as.numeric(unique(fips_seq))
    Same = subset(Fips,Fips  %notin% comb[ii,])
    fips2 = matrix(fips_seq,Ti,N)
    fips2[,Same] = 0
    X <- model.matrix(~factor(fips) + factor(fips2)*shock - shock - factor(fips2))
    Z <- model.matrix(~shock)
    est_dL <- gmm2S(Ti,N,y,X,Z, yy, k_t,tol, guess, maxiter,second)
    phi_d <- est_dL$phi_hat
    V_d <- est_dL$Var
    nbeta= NROW(phi_d)
    some_diff[ii,c("est","se")] =  c(Hd%*%(phi_d[(nbeta-Nd+1):nbeta]),
                                     sqrt(Hd%*%V_d[(nbeta-Nd+1):nbeta,(nbeta-Nd+1):nbeta]%*%t(Hd)))
    some_diff[ii,c("ci_lo","ci_hi")] <- some_diff[ii,"est"] + some_diff[ii,"se"]*c(-1.7,1.7)
    some_diff[ii,c("pval_d","J_d")] <- c(est_dL$pval,est_dL$J)
  }
  
  not_reject_diff <- some_diff[some_diff[,"pval_d"]>0.1,]
  
  ci.results <- matrix(,3,3)
  rownames(ci.results) <- c("aggregate","common","some diff")
  colnames(ci.results) <- c("ci_lo","ci_hi","J-pval")
  ci.results[c("some diff"),c("ci_lo","ci_hi")] <- 
            c(min(not_reject_diff[,"ci_lo"]),max(not_reject_diff[,"ci_hi"]))
  
  return_list <- list("some_diff"=some_diff,"ci.results"=ci.results)
  return(return_list)
  
}