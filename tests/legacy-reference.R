# Original supplied simulation implementation, with only the function name
# and candidate grid changed for a compact regression test.
legacy_reference = function(Y, X, C, eigen_values_positive, eigen_vector_sigma, BBB = 10000){
        
        p = nrow(Y)
        q = ncol(C)
        k = nrow(X)
        n = ncol(Y) -k
        gamma = p/(n)
        
        XXTinverse = solve(X%*%t(X))
        
        PX = t(X) %*% XXTinverse %*% X
        
        eigen_CX = eigen(t(C) %*% XXTinverse %*% C)
        
        CXinverse = eigen_CX$vectors %*% diag( (eigen_CX$values)^(-1/2),  nrow = length(eigen_CX$values)) %*% t(eigen_CX$vectors)
        
        Q_n =  t(X) %*% XXTinverse %*% C %*% CXinverse 
        
        #Sigma_hat = 1/(n) * Y %*% (diag(ncol(Y)) - PX) %*% t(Y)

        #half_Sigma =  1/sqrt(n) * Y %*% (diag(ncol(Y)) - PX)
        #svd_half_Sigma_hat = svd( half_Sigma, nu = nrow(half_Sigma), nv = 0 )
        #svdofSigma_hat =  svd(Sigma_hat)
        
        emp_eigen_positive = eigen_values_positive # svd_half_Sigma_hat$d^2
        
        no_negative_eigen = p - length(emp_eigen_positive)

        emp_eigen = c(emp_eigen_positive, rep(0, times =no_negative_eigen ))   ## shrink small eigenvalues to exactly 0. 
        
        #Sigma_hat = svd_half_Sigma_hat$u %*% diag(emp_eigen) %*% t(svd_half_Sigma_hat$u)
        projected_Y =  1/sqrt(n) * t(eigen_vector_sigma) %*% Y %*% Q_n ## M(f) = t(projected_Y) %*% diag(shrunk eigenvalues) %*% projected_Y
        
        ### construct optimization region 
        
        frakM = mean(emp_eigen)
        
        lambda <- c(-0.1, -1, -10) # Small shared grid for regression testing
      

        m_n = 1/p* ( rowSums( -1/ outer(lambda, emp_eigen_positive, FUN = "-") ) - no_negative_eigen / lambda )
        
        m_n_prime = 1/p* ( rowSums( 1 / (outer(lambda, emp_eigen_positive, FUN = "-"))^2 ) + no_negative_eigen / lambda^2 )
                
        Theta =  1/( 1 - gamma - gamma * lambda * m_n )
        
        tau =   (Theta -1)
        
        var11 = gamma * (1 + lambda * m_n) * Theta^3 + gamma * lambda * ( m_n + lambda * m_n_prime) * Theta^4
        Delta = 2 * var11
        
        
        rho0 = m_n
        rho1 = Theta * (1 + lambda * m_n)
        rho2 = Theta * (frakM + lambda * rho1)
        
        opt0_index = which.max(rho0 / sqrt(Delta))
        
        opt1_index = which.max(rho1 / sqrt(Delta))
        
        opt2_index = which.max(rho2 / sqrt(Delta))
                
        opt0_lambda = lambda[opt0_index]
        opt1_lambda = lambda[opt1_index]
        opt2_lambda = lambda[opt2_index]
        
        opt0_mean = tau[opt0_index]
        opt1_mean = tau[opt1_index]
        opt2_mean = tau[opt2_index]
        
        opt0_var = Delta[opt0_index]
        opt1_var = Delta[opt1_index]
        opt2_var = Delta[opt2_index]

        opt0_Theta = Theta[opt0_index]
        opt1_Theta = Theta[opt1_index]
        opt2_Theta = Theta[opt2_index]

        if(opt0_lambda == opt1_lambda){
            cor01 = 1
        }else{
            cor01 = (opt0_var * opt1_var)^(-1/2) * 2 * opt0_Theta * opt1_Theta * 
                    ((opt0_lambda * opt0_Theta - opt1_lambda * opt1_Theta)/(opt0_lambda - opt1_lambda)  - 1) 
        }
        if(opt0_lambda == opt2_lambda){
            cor02 = 1
        }else{
            cor02 = (opt0_var * opt2_var)^(-1/2) * 2 * opt0_Theta * opt2_Theta * 
                    ((opt0_lambda * opt0_Theta - opt2_lambda * opt2_Theta)/(opt0_lambda - opt2_lambda)  - 1) 
        }
        if(opt1_lambda == opt2_lambda){
            cor12 = 1
        }else{
            cor12 = (opt1_var * opt2_var)^(-1/2) * 2 * opt1_Theta * opt2_Theta * 
                    ((opt1_lambda * opt1_Theta - opt2_lambda * opt2_Theta)/(opt1_lambda - opt2_lambda)  - 1) 
        }
        Cov_Kernel = matrix( c(1, cor01, cor02, cor01, 1, cor12, cor02, cor12, 1), 3,3)
        eigen_Cov_Kernel = eigen(Cov_Kernel)
        eigen_Cov_Kernel_plus =  eigen_Cov_Kernel$values * ( eigen_Cov_Kernel$values >= 0) 
        Cov_Kernel_sqrt = eigen_Cov_Kernel$vectors %*% diag(sqrt(eigen_Cov_Kernel_plus))
        

        #cor01 = ifelse(opt0_lambda != opt1_lambda,  
        #                ( opt0_var * opt1_var)^(-1/2) * 2 * Theta[opt0_index] * Theta[opt1_index] * 
        #                ((opt0_lambda * Theta[opt0_index] - opt1_lambda * Theta[opt1_index])/( opt0_lambda - opt1_lambda)  -1),
        #                1)
        #cor02 = ifelse(opt0_lambda != opt2_lambda,  ( opt0_var * opt2_var)^(-1/2) * 2 * Theta[opt0_index] * Theta[opt2_index] * 
        #        ( (opt0_lambda * Theta[opt0_index] - opt2_lambda * Theta[opt2_index])/( opt0_lambda - opt2_lambda)  -1  ), 1) 
        #cor12 = ifelse( opt1_lambda != opt2_lambda, ( opt1_var * opt2_var)^(-1/2) * 2 * Theta[opt1_index] * Theta[opt2_index] * 
        #        ( (opt1_lambda * Theta[opt1_index] - opt2_lambda * Theta[opt2_index])/( opt1_lambda - opt2_lambda)  -1  ), 1)


        
        opt0_shrunk_precision = 1 /(emp_eigen - opt0_lambda)
        opt1_shrunk_precision = 1 /(emp_eigen - opt1_lambda)
        opt2_shrunk_precision = 1 /(emp_eigen - opt2_lambda)
        
        M_f0 = eigen( t(projected_Y) %*% diag(opt0_shrunk_precision) %*% projected_Y )$values
        M_f1 = eigen( t(projected_Y) %*% diag(opt1_shrunk_precision) %*% projected_Y )$values
        M_f2 = eigen( t(projected_Y) %*% diag(opt2_shrunk_precision) %*% projected_Y )$values
        
        T_LR0 = sum(log(1+M_f0))
        T_LR1 = sum(log(1+M_f1))
        T_LR2 = sum(log(1+M_f2))
        
        T_LH0 = sum(M_f0)
        T_LH1 = sum(M_f1)
        T_LH2 = sum(M_f2)
        
        T_BNP0 = sum(M_f0/(1+M_f0))
        T_BNP1 = sum(M_f1/(1+M_f1))
        T_BNP2 = sum(M_f2/(1+M_f2))
        
        stand_T0_LR = sqrt(n) * (1+ opt0_mean) / q^(1/2) / opt0_var^(1/2) * ( T_LR0 - q* log( 1+ opt0_mean) ) 
        stand_T1_LR = sqrt(n) * (1+ opt1_mean) / q^(1/2) / opt1_var^(1/2) * ( T_LR1 - q* log( 1+ opt1_mean) )
        stand_T2_LR = sqrt(n) * (1+ opt2_mean) / q^(1/2) / opt2_var^(1/2) * ( T_LR2 - q* log( 1+ opt2_mean) )
        
        comp_LR = max(stand_T0_LR, stand_T1_LR, stand_T2_LR)

        stand_T0_LH = sqrt(n) / q^(1/2) / opt0_var^(1/2) *(T_LH0 - q* opt0_mean)
        stand_T1_LH = sqrt(n) / q^(1/2) / opt1_var^(1/2) *(T_LH1 - q* opt1_mean)
        stand_T2_LH = sqrt(n) / q^(1/2) / opt2_var^(1/2) *(T_LH2 - q* opt2_mean)
        
        comp_LH = max(stand_T0_LH, stand_T1_LH, stand_T2_LH)

        stand_T0_BNP = sqrt(n) * (1 + opt0_mean)^2 / q^(1/2) / opt0_var^(1/2) * ( T_BNP0 - q * opt0_mean / (1+opt0_mean) )
        stand_T1_BNP = sqrt(n) * (1 + opt1_mean)^2 / q^(1/2) / opt1_var^(1/2) * ( T_BNP1 - q * opt1_mean / (1+opt1_mean) )
        stand_T2_BNP = sqrt(n) * (1 + opt2_mean)^2 / q^(1/2) / opt2_var^(1/2) * ( T_BNP2 - q * opt2_mean / (1+opt2_mean) )
        
        comp_BNP = max(stand_T0_BNP, stand_T1_BNP, stand_T2_BNP)

        stat = cbind(comp_LR = comp_LR, comp_LH = comp_LH, comp_BNP = comp_BNP)

        Tmax =  apply(Cov_Kernel_sqrt %*% matrix(rnorm(3*BBB), 3, BBB), 2, max)

        pvalue = cbind(comp_LR = mean(comp_LR <= Tmax),
                       comp_LH = mean(comp_LH <= Tmax),
                       comp_BNP = mean(comp_BNP <= Tmax)
            )
        reject_Ind = pvalue < 0.05  
        
        return(list(stat = stat, reject_Ind = reject_Ind, pvalue =pvalue))    
} 


