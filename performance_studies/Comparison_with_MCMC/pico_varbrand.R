## Install Packages ##
#install.packages("Rcpp")
#install.packages("RcppEigen")
#install.packages("foobar2_1.0.tar.gz",repos = NULL, type = "source")
#install.packages("ramify")

## Load Libraries ##
library(Rcpp)
library(RcppEigen)
library(foobar2)


#PICO_VARBRAND
pico_varbrand <- function(X,Y,p,cl_train_label_noise, G, cl_tot){ 
  start_time = Sys.time()
  
  # LIBRARIES
  #install.packages("Rcpp")
  #install.packages("rrcov")
  # install.packages("foobar2")
  library(foobar2)
  library(Rcpp)
  
  library(dendextend)
  library(aricode)
  
  
  # LOAD DATASETS
  labels_tot = cl_tot
  labels_training = cl_train_label_noise
  Y_training = X
  Y_tot = Y
  
  #SPECIFY USER INPUT
  #Number of components
  p = dim(Y_tot)[2]
  
  #Number of classes in the training set
  J = G #dim(unique(labels_training))[1]
  
  #Max number of classes in the Dirichlet Process
  Tr = 20
  
  #Number of iterations and tolerance -  CAVI 
  n_iter = 1000
  tol = 1e-5
  
  #HYPERPARAMETERS
  #gamma -> Stick Breaking parameter -> scalar
  gamma = 5
  
  #a_dir_k -> Dirichlet components -> vector of (J+1) components
  a_dir_k = rep(0.1, J + 1)
  
  #nIW_DP_0
  # 	mu_0_DP -> vector of (p) components
  # 	nu_0_DP -> scalar
  # 	lambda_0_DP -> scalar
  # 	PSI_0_DP -> matrix (pxp)
  
  mu_0_DP = rep(0, p)
  
  nu_0_DP = c(p+1+1)
  
  lambda_0_DP = c(0.1)
  
  PSI_0_DP = 3*matrix(diag(p), nrow = p, ncol = p)
  
  
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
  
  nu_0_MIX = rep(10 + p + 1,J)
  
  lambda_0_MIX = rep(20,J)
  
  PSI_0_MIX = NULL
  PSI_0_MIX_FACTOR = 1
  
  # PSI_0_MIX set automatically post robust parameters calculation:
  # if (PSI_0_MIX == NULL) =>
  #   PSI_0_MIX = robust_inv_cov_mat 
  #   for(i in 1:length(PSI_0_MIX)){
  #     PSI_0_MIX[[i]] = PSI_0_MIX_FACTOR*PSI_0_MIX[[i]]
  #   } 
  
  
  #VARIATIONAL PARAMETERS
  M = dim(Y_tot)[1]
  Phi_m_k = matrix(1, M, J + Tr) * (1 / (J + T))
  
  eta_k = a_dir_k
  
  a_k_beta = rep(1, Tr - 1)
  
  b_k_beta = rep(gamma, Tr - 1)
  
  
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
  
  nu_var_DP = rep(nu_0_DP, Tr)
  
  lambda_var_DP = rep(lambda_0_DP, Tr)
  
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
  
  result = foobar2::var_brand(Y_tot, Y_training, labels_tot, labels_training, p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                              lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                              PSI_0_MIX,PSI_0_MIX_FACTOR,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                              lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX)
  
  # print(result$elbo_values)
  # #REPORT 
  labels_pred = result$labels_pred
  n = length(labels_pred)
  
  tot_clusters = length(unique(cl_tot))
  cluster_found = length(unique(labels_pred))
  ARI = result$ari
  
  FMI = FM_index(labels_tot, labels_pred)[1]
  ami = AMI(labels_tot, labels_pred)
  
  
  end_time = Sys.time()
  elapsed_time_mins = end_time - start_time
  
  output = list("elapsed_time_mins" = elapsed_time_mins, "tot_clusters" = tot_clusters, "cluster_found" = cluster_found, "ARI" = ARI, "FMI" = FMI, "ami" = ami)
  return(output)
}

## Test ##
# N_FACTOR = 0.1
# p = 2
# easy = 0
# generated_data = generate_data(N_FACTOR, p, easy)
# output_varbrand = pico_varbrand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
# print(output_varbrand)