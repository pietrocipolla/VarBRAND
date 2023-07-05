wd <- "***//VarBRAND"

#Packages

install.packages("packages//variational//foobar2_1.0_nostdout.tar.gz",repos = NULL, type = "source")
#package installation procedure may be machine-dependent due safety options
#many issues can be solved by following the procedure reported in https://stackoverflow.com/questions/42807247/installing-package-cannot-open-file-permission-denied
library(foobar2)
library(FME)
library(rrcov)
library(Rcpp)
library(dendextend)
library(aricode)
library("ramify")
library(data.table)

setwd(paste(wd, "//application", sep = ""))

#PICO_VARBRAND
pico_varbrand_soil <- function(training_set, test_set, rand_par) {
  # SET SEED FOR REPRODUCIBILITY
  set.seed(777)
  
  start_time <- Sys.time()
  
  # LOAD DATASETS
  labels_tot = test_set[,1]
  labels_training = training_set[,1]
  Y_training = training_set[,-1]
  Y_tot = test_set[,-1]
  
  
  #Number of components
  p = dim(Y_tot)[2]
  
  #Number of classes in the training set
  J = length(unique(labels_training))
  
  #Max number of classes in the Dirichlet Process
  Tr = 10
  
  #Number of iterations and tolerance -  CAVI 
  n_iter = 1000
  tol = 1e-5
  
  # Prior specification
  gamma = 10
  
  a_dir_k = rep(0.1, J + 1)
  
  mu_0_DP = colMeans(Y_training)
  
  nu_0_DP = c(p+2)
  
  lambda_0_DP = c(0.1)
  
  PSI_0_DP = (p+1)*var(Y_training)
  
  
  mu_0_MIX = NULL 
  nu_0_MIX = rep(200+p+1,J)
  
  lambda_0_MIX = rep(200,J)
  
  PSI_0_MIX = NULL
  PSI_0_MIX_FACTOR = 200
  
  
  
  # variational parameters
  M = dim(Y_tot)[1]
  Phi_m_k = matrix(1, M, J + Tr) * (1 / (J + T))
  
  eta_k = a_dir_k*rand_par[1]
  
  a_k_beta = rep(1, Tr - 1)
  
  b_k_beta = rep(gamma, Tr - 1)
  
  
  mu_var_DP = t(kmeans(Y_tot, Tr, iter.max = 400)$center) 
  
  nu_var_DP = rep(nu_0_DP, Tr)*rand_par[2]
  
  lambda_var_DP = rep(lambda_0_DP, Tr)*rand_par[3]
  
  PSI_var_DP = list()
  for (i in 1:Tr) PSI_var_DP[[i]] = PSI_0_DP
  
  mu_VAR_MIX = NULL # later initialized as robust parameters
  
  nu_VAR_MIX = nu_0_MIX
  
  lambda_VAR_MIX = lambda_0_MIX
  
  PSI_VAR_MIX = NULL # later initialized as robust parameters
  
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
  
  # Update model
  # mu_0_MIX
  if (is.null(mu_0_MIX)){
    mu_0_MIX = t(robust_mean)
  }
  
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
  
  ####  FIT - STEP 2  ####
  result = foobar2::VarBrand(as.matrix(Y_tot), p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                             lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                             PSI_0_MIX,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                             lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX)
  
  
  Phi_m_k_l = result$Phi_m_k
  mu_VAR_MIX_updated  = result$mu_VAR_MIX
  mu_var_DP_updated  = result$mu_var_DP
  elbo_values = result$elbo_values
  
  
  ####  RESULTS  ####
  #Generate labels_pred
  labels_pred = argmax(Phi_m_k_l, rows = TRUE)
  
  ## generate report (if CAVI converged)
  if(sum(Phi_m_k_l) != 0) {
    # ARI
    n = length(labels_pred)
    
    
    tot_clusters = length(unique(labels_tot))
    cluster_found = length(unique(labels_pred))
    
    ARI = ARI(labels_pred,labels_tot)
    print("ARI:")
    print(ARI)
    FMI = FM_index(labels_tot, labels_pred)[1]
    max_elbo <- elbo_values[length(elbo_values)]
  } else {
    ARI = NaN
    FMI = NaN
    max_elbo = NaN
  }
  
  end_time = Sys.time()
  elapsed_time_mins = end_time - start_time
  
  ret <- list()
  ret[["output"]] <- c(rand_par, ARI, FMI, max_elbo)
  ret[["elbo"]] <- data.frame(elbo_values)
  
  return(ret)
}

run_mc_varbrand <- function(training_set, test_set) {
    set.seed(42)

    N <- 4

    output <- data.frame(matrix(ncol = 6, nrow = N))
    colnames(output) <- c("a_dir_k", "nu_0_DP", "lambda_0_DP", "ARI", "FMI", "max_elbo")
    elbos <- list()

    #define grid for random parameters
    a_dir_k <- c(0.1, 1)
    nu_0_DP <- c(1, 10)
    lambda_0_DP <- c(1, 10)

    parRange <- rbind(a_dir_k, nu_0_DP, lambda_0_DP)
    input <- Latinhyper(parRange = parRange, N)

    #run
    for (i in 1:N) {
        temp.out <- pico_varbrand_soil(training_set, test_set, input[i, ])
        output[i, ] <- temp.out$output
        elbos[[i]] <- temp.out$elbo
    }

    elbos <- rbindlist(elbos, idcol = TRUE)

    ret <- list()
    ret[["indices"]] <- output
    ret[["elbo"]] <- elbos

    return(ret)
}