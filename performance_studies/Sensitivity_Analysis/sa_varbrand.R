library(rrcov)
library(ramify)

run_sa <- function(generated_data) {
  Y_tot <- generated_data[generated_data$train==0, 1:2]
  Y_training <- generated_data[generated_data$train==1, 1:2]
  
  labels_tot <- generated_data[generated_data$train==0, 3]
  labels_training <- generated_data[generated_data$train==1, 3]
   
  #SPECIFY USER INPUT
  #Number of components
  p = dim(Y_tot)[2]
  
  #Number of classes in the training set
  J = length(unique(labels_training))[1]
  
  #Max number of classes in the Dirichlet Process
  Tr = 10
  
  #Number of iterations and tolerance -  CAVI 
  n_iter = 1000
  tol = 1e-5
  
  #HYPERPARAMETERS
  #a_dir_k -> Dirichlet components -> vector of (J+1) components
  a_dir_k = rep(1, J + 1)
  
  #nIW_DP_0
  # 	mu_0_DP -> vector of (p) components
  # 	nu_0_DP -> scalar
  # 	lambda_0_DP -> scalar
  # 	PSI_0_DP -> matrix (pxp)
  
  mu_0_DP = rep(0, p)
  
  PSI_0_DP = matrix(diag(p), nrow = p, ncol = p)
  
  #NIW_MIX_0
  # mu_0_MIX -> matrix (pxJ)
  #           -> in every column, the mean components of the mixture
  #               (mu_0_MIX[:,i] -> mean component of the (i+1)-th NIW of the mixture)
  #
  # nu_0_MIX -> vector of (J) components
  #               (nu_0_MIX[i] = nu_0 of the (i+1)-th NIW of the mixture)
  #
  # lambda_0_MIX -> vector (J) components
  #               (lambda_0_MIX[i] = lambda_0 of the (i+1)-th NIW of the mixture)
  #
  # PSI_0_MIX -> vector of matrices (Jxpxp), essentially an ndarray
  #           -> PSI_0_MIX[i,:,:] =  PSI Matrix of the (i+1)-th NIW of the mixture
  
  mu_0_MIX = NULL 
  # mu_0_MIX set automatically post robust parameters calculation:
  # if mu_0_MIX == NULL:
  #   mu_0_MIX = robust_mean 
  
  nu_0_MIX = rep(p+2,J)
  
  lambda_0_MIX = rep(1,J)
  
  PSI_0_MIX = NULL
  PSI_0_MIX_FACTOR = 1
  
  # PSI_0_MIX set automatically post robust parameters calculation:
  # if PSI_0_MIX == NULL =>
  #   PSI_0_MIX = robust_inv_cov_mat 
  #   for(i in 1:length(PSI_0_MIX)){
  #     PSI_0_MIX[[i]] = PSI_0_MIX_FACTOR*PSI_0_MIX[[i]]
  #   } 
  
  
  #VARIATIONAL PARAMETERS
  M = dim(Y_tot)[1]
  Phi_m_k = matrix(1, M, J + Tr) * (1 / (J + T))
  
  eta_k = a_dir_k
  
  a_k_beta = rep(1, Tr - 1)
  
  # NIW_DP_VAR
  # mu_var_DP -> matrix (Txp)
  # -> in every column,  the mean components of the mixtureDP
  #    (mu_var_DP[:,i] -> mean component of the (i+1)-th NIW of the mixtureDP)
  
  # nu_var_DP -> vector of (T) components
  # (nu_var_MIX[i] = nu_var of the (i+1)-th NIW of the mixtureDP)
  
  # lambda_var_DP -> vector of (T) components
  # (lambda_var_DP[i] = lambda_var of the (i+1)-th NIW of the mixtureDP)
  
  # PSI_var_DP -> vector of matrices (Txpxp),  essentially an ndarray
  # -> PSI_var_DP[i,:,:] = PSI Matrix of the (i+1)-th NIW of the mixtureDP
  
  mu_var_DP = t(kmeans(Y_tot, Tr)$center) 
  # Alternative
  # mu_var_DP <- matrix(1, Tr, p) 
  # 
  # for(i in 1:Tr){
  #   for(j in 1:p){
  #     mu_var_DP[i,j] = rnorm(n=1, mean = 0, sd = 5)
  #   }
  # }
  
  PSI_var_DP = list()
  for (i in 1:Tr) PSI_var_DP[[i]] = PSI_0_DP
  
  #NIW_MIX_VAR:
  # mu_VAR_MIX -> matrix (Jxp)
  #           -> in every column,  the mean components of the mixture
  #               (mu_VAR_MIX[:,i] -> mean of the (i+1)-th NIW of the mixture)
  #
  # nu_VAR_MIX -> vector of (J) components
  #               (nu_VAR_MIX[i] = nu_0 of the (i+1)-th NIW of the mixture)
  #
  # lambda_VAR_MIX ->  vector of  (J) components
  #               (lambda_VAR_MIX[i] = lambda_0 of the (i+1)-th NIW of the mixture)
  #
  # PSI_VAR_MIX ->  vector of matrices (Jxpxp), essentially an ndarray
  #           -> PSI_VAR_MIX[i,:,:] = PSI Matrix of the (i+1)-th NIW of the mixture
  
  mu_VAR_MIX = NULL
  # mu_VAR_MIX set automatically post robust parameters calculation:
  # if mu_VAR_MIX == NULL
  #   mu_VAR_MIX = mu_0_MIX 
  
  nu_VAR_MIX = nu_0_MIX
  
  lambda_VAR_MIX = lambda_0_MIX
  
  PSI_VAR_MIX = NULL 
  # PSI_VAR_MIX set automatically post robust parameters calculation:
  # if PSI_VAR_MIX == NULL, 
  #   PSI_VAR_MIX = PSI_0_MIX 
  
  ####  FIT - STEP 1  ####
  # Generate Robust Parameters
  robust_mean <- c()
  robust_inv_cov_mat <- list()
  
  i = 1;
  while(i<=J){
    index = which(labels_training==i)
    robust_mean <- rbind(robust_mean, CovMcd(Y_training[index,])$center)
    robust_inv_cov_mat[[i]] <- CovMcd(Y_training[index,])$cov
    i = i+1
  }
  robust_mean
  robust_inv_cov_mat
  
  # Update model
  # mu_0_MIX
  if (is.null(mu_0_MIX)){
    mu_0_MIX = t(robust_mean)
  }#EVENTUALMENTE TRASPOSTA
  
  # PSI_0_MIX
  if (is.null(PSI_0_MIX)){
    PSI_0_MIX = robust_inv_cov_mat
    for(i in 1:length(PSI_0_MIX)){
      PSI_0_MIX[[i]] = PSI_0_MIX_FACTOR *PSI_0_MIX[[i]]
    }
  }
  
  # mu_VAR_MIX
  if (is.null(mu_VAR_MIX)){
    mu_VAR_MIX = mu_0_MIX
  }
  
  # PSI_VAR_MIX
  if (is.null(PSI_VAR_MIX)){
    PSI_VAR_MIX = PSI_0_MIX
  }
  
  #Number of points in the interval of the single variable
  n_points <- 3
  
  result <- explore_space(Y_tot, Y_training, labels_tot, labels_training, p, J, Tr, n_iter, tol, gamma, 
                                a_dir_k, mu_0_DP, PSI_0_DP, mu_0_MIX, nu_0_MIX, lambda_0_MIX, PSI_0_MIX, 
                                PSI_0_MIX_FACTOR, M, Phi_m_k, eta_k, a_k_beta, b_k_beta, mu_var_DP, 
                                PSI_var_DP, mu_VAR_MIX, nu_VAR_MIX, lambda_VAR_MIX, PSI_VAR_MIX, n_points)
  
  return(result)
}
