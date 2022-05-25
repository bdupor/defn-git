
/*******************************************************************************
**Efficient Iterative GMM:
********************************************************************************

Requires
	t     		- length of panel
	i     		- size of cross-section
	dep   		- dependent variable (as a vector)
	indep       - independent varibles (as a matrix)
	instrument  - instrumental variables (as a matrix)
	yy  		- serial correlation bandwidth (t_0 in the paper)
	partition   - N x 1 vector id-ing partition assignments when combining moments
 
Optional
	tol     - tolerance
	guess   - replace initial guess with user entered values
	maxiter - maximum number of iterations
	geniv   - use generalize IV as the initial guess weighting matrix
	second  - use second moment independence 

*********************************************************************************/


program define PartialGMM2S, eclass
	version 10
	syntax, t(string) i(string) dep(string) indep(string) instrument(string) ///
		    yy(int) partition(string)  ///
			 [ tol(real 1e-6) guess(string) maxiter(int 1000) geniv(int 0) second(int 1) ]
			 
	matrix indep = `indep'
	matrix dep 	 = `dep'
	matrix instrument = `instrument'
	matrix part = `partition'
	
	if "`guess'" == ""{
		matrix guess = (-99)
	}
	else{
		matrix guess = `guess'
	}
	
mata{
		****  Setup	
		T = strtoreal(st_local("t"))
		N = strtoreal(st_local("i"))
		yy = strtoreal(st_local("yy"))
		tol = strtoreal(st_local("tol"))
		maxiter = strtoreal(st_local("maxiter"))
		geniv = strtoreal(st_local("geniv"))
		second = strtoreal(st_local("second"))
		
		y = st_matrix("dep")
		x = st_matrix("indep")
		z = st_matrix("instrument")
		part = st_matrix("part")
		guess = st_matrix("guess")
		K = cols(z)
		unique_part = uniqrows(part)
		pN = rows(unique_part)
		num_var = cols(x)

		
		/*** Setup 
		Following allows for moment conditions and derivatives of moment 
		conditions to be combined along the partition used when multiplying 
		x or u from the right.
		*/
		z2 = J(N*T,pN*K,0)
		for (p=1; p<=pN; p++){
			p_r = unique_part[p,1]
			for (i = 1; i<=N; i++){
				ii = T*(i-1) +1
				jj = K*(p-1) +1
				if (p_r == part[i,1]){
					aux = sum(part:==part[i,1])
					z2[(ii::ii+T-1),(jj::jj+K-1)] = z[(ii::ii+T-1),1...]  / aux
				}
			}
		}
		

		**** Initial Guess: W=I or Gen IV or Entry
		if (guess == (-99)){
			W_0 = I(pN*K)
			if (geniv==1){
				W_0 = invsym(z2'*z2)
			}
			phi_1st = invsym(x'*z2*W_0*z2'*x)*x'*z2*W_0*z2'*y
		}
		else {
			phi_1st = guess
		}
		
		
		**** Covariance estimation
		Omega_hat   = J(pN*K,pN*K,0)
		W           = J(pN*K,pN*K,0)
		phi_hat     = phi_1st
		phi_guess   = phi_1st
		var         = J(cols(x),cols(x),0)
		Omega_u_hat = J(pN,pN,0)
		iter = 0
		diff = 1
		
		
		
		while (diff > tol & iter < maxiter){
			**** First stage Residuals
			u_hat = y - x*phi_guess
			
			**** Covariance estimation
			Omega_hat_temp   = J(pN*K,pN*K,0)
			Omega_u_hat_temp = J(pN,pN,0)
			
			for (k = -yy; k<= yy; k ++){
				start = max((1,1-k))
				finish = min((T, T-k))
				
				w  = 1 - abs(k)/(yy+1) 
				
				if (second == 1){
					Omega_u_aux = J(pN,pN,0)
					Omega_z_aux = J(pN*K,pN*K,0)
				}
				
				for (t = start; t<= finish; t++){
					if (second == 1) {
						u_t  = m2_t(t,u_hat,N,T,part)
						u_tk = m2_t(t+k,u_hat,N,T,part)
						z_t  = m2_t(t,z,N,T,part)
						z_tk = m2_t(t+k,z,N,T,part)
						Omega_u_aux = Omega_u_aux + u_t*u_tk'
						Omega_z_aux = Omega_z_aux + z_t*z_tk'
					}
					else{
						mom_t   = m_t(t,u_hat,z,N,T,part)
						mom_tk  = m_t(t+k,u_hat,z,N,T,part)
						Omega_hat_temp = (1/T)*w*mom_t*mom_tk' + Omega_hat_temp
					}
				}
				
				if (second == 1){
					Omega_u_hat_temp = (1/T)*w*Omega_u_aux + Omega_u_hat_temp
					
					Omega_u = (1/T)*Omega_u_aux#J(K,K,1)
					Omega_z = (1/T)*Omega_z_aux
					Omega_hat_temp = w*Omega_u :* Omega_z + Omega_hat_temp
				}
				
				
			}

			**** Put Together
			W = invsym(Omega_hat_temp)
			phi_hat = invsym(x'*z2*W*z2'*x)*x'*z2*W*z2'*y
			var = (1/T)*invsym( (1/T)*x'*z2 * W * z2'*x*(1/T))
			
			diff = max(abs(phi_hat - phi_guess))
			iter = iter + 1
			
			phi_guess   = phi_hat
			Omega_hat   = Omega_hat_temp
			Omega_u_hat = Omega_u_hat_temp
		}
		
		**** J-stat
		u_hat = y - x*phi_hat
		m_bar = J(pN*K,1,0)
		for (t = 1; t<= T; t++){
			m = m_t(t,u_hat,z,N,T,part)
			m_bar = m_bar + m
		}
		m_bar = (1/sqrt(T))*m_bar
		J = m_bar'*W*m_bar
		
		
		**** Wald Statistic
		Wald = phi_hat'*invsym(var)*phi_hat
	
		**** Collect Results
		st_matrix("Omega", Omega_hat)
		st_matrix("Omega_z", Omega_z)
		st_matrix("phi",phi_hat)
		st_matrix("phi_1st",phi_1st)
		st_matrix("Sigma", var)
		st_matrix("W", W)
		st_matrix("J",J)
		st_matrix("Wald", Wald)
		st_matrix("Omega_u", Omega_u_hat)
		st_matrix("iter", iter)
		st_matrix("diff", diff)

		printf("Exited with:")
		printf(" %g iterations \n",iter)
		printf(" and %g tolerance\n",diff)
	
}
		**** Save Results
		ereturn clear
		ereturn matrix phi       = phi
		ereturn matrix phi_1st   = phi_1st
		ereturn matrix Omega     = Omega
		ereturn matrix Sigma     = Sigma
		ereturn matrix W         = W
		ereturn matrix J         = J
		ereturn matrix Wald		 = Wald
		ereturn matrix Omega_u   = Omega_u
		ereturn matrix iter      = iter
		ereturn matrix diff      = diff
		ereturn matrix Omega_z   = Omega_z

end

*** Functions
mata:

	// See notes for details on creation of P matrix below

	// uz moment at t
	function m_t(t,u,z,N,T,part){
		K = cols(z)
		unique_part = uniqrows(part)
		pN = rows(unique_part)
		
		E_k = I(K) 
		P_aux = J(pN,N,0)
		for (i=1; i<=N; i++){
			p_r = part[i,1]
			for (r=1; r<=pN; r++){
				if (r == p_r){
					aux = sum(p_r:==part)
					P_aux[r,i] = 1  / aux
				}
			}
		}
		
		P = P_aux#E_k
		
		m = J(N*K,1,0) 
		for (i = 1; i<=N; i++) {
			for (r = 1; r<=K; r++){
				idx = K*(i - 1) + r 
				j = t+T*(i-1)
				m[idx,1] =  u[j,1]*z[j,r]  
			}
		}
		
		
		return(P*m)
	}
	
	// u or z moment at t (for second mom. ind.)
	function m2_t(t,u,N,T,part){
		K = cols(u)
		unique_part = uniqrows(part)
		pN = rows(unique_part)
		
		E_k = I(K) 
		P_aux = J(pN,N,0)
		for (i=1; i<=N; i++){
			p_r = part[i,1]
			for (r=1; r<=pN; r++){
				if (r == p_r){
					aux = sum(p_r:==part)
					P_aux[r,i] = 1   / aux
				}
			}
		}
		
		P = P_aux#E_k
		
		m = J(N*K,1,0) 
		for (i = 1; i<=N; i++) {
			for (r = 1; r<=K; r++){
				idx = K*(i - 1) + r 
				j = t+T*(i-1)
				m[idx,1] = u[j,r]  
			}
		}
		
		
		return(P*m)
	}
	
	
	
	
end
