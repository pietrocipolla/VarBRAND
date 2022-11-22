### GENERATE_DATA.R ###
## Install Packages ##
#install.packages("mvnfast")
#install.packages("purrr")
#install.packages("dplyr")

## Load Libraries ##
library(mvnfast)
library(purrr)
library(dplyr)

#GENERATE DATA
generate_data <- function(N_FACTOR, p, easy) {
  # training set realizations
  N <- 1000 * N_FACTOR
  
  # full dataset realizations
  M <- 1000 * N_FACTOR
  
  # training set cardinalities
  N_vec <- list(500, 500)
  N_vec = lapply(N_vec,"*",N_FACTOR)
  
  # full dataset cardinalities
  #M_vec <- list(200, 200, 250, 90, 100, 100, 10) <= paper_denti_cappozzo
  M_vec <- list(200, 200, 200, 200, 200)
  M_vec = lapply(M_vec,"*",N_FACTOR)
  
  #alternative cardinalities
  # "small" = list(350, 250, 250, 49, 50, 50, 1)
  # "not_small" = list(200, 200, 250, 90, 100, 100, 10)
  
  # clusters in the training set
  G <- 2
  
  # training sets labels
  cl_train <- rep(1:G, N_vec)
  
  #full labels
  cl_tot <- rep(1:length(M_vec), M_vec)
  
  #mu
  MU <-
    list(c(2, 2),
         c(-2, -2),
         c(-2, 2),
         c(2, -2),
         c(0, 0))
  
  #sigma
  if(easy){
    SIGMA <-
      list(
        diag(0.2, 2),
        diag(0.2, 2),
        diag(0.2, 2),
        diag(0.2, 2),
        diag(0.2, 2)
      )
  }
  else{
    SIGMA <-
      list(
        diag(0.75, 2),
        diag(0.75, 2),
        diag(0.75, 2),
        diag(0.75, 2),
        0.5*diag(0.75, 2)
      )
  }
  
  
  #multidimensional MU
  desired_length <- length(MU) 
  updated_mu <- vector(mode = "list", length = desired_length)
  
  i = 1
  while(i<length(MU)+1){
    mu = matrix(0,1,p)
    mu[1:2]=MU[[i]]
    updated_mu[[i]]=mu
    i=i+1
  }
  
  MU = updated_mu
  
  #multidimensional SIGMA
  desired_length <- length(SIGMA) 
  updated_SIGMA <- vector(mode = "list", length = desired_length)
  
  i = 1
  while(i<length(SIGMA)+1){
    sigma = 0.2*diag(p)
    sigma[1:2,1:2] = SIGMA[[i]]
    updated_SIGMA[[i]]=sigma
    i=i+1
  }
  
  SIGMA = updated_SIGMA
  
  #training set labels
  cl_test <- rep(1:length(M_vec), M_vec)
  
  #optionally add label noise
  label_noise <- FALSE
  cl_train_label_noise <- cl_train
  if (label_noise == TRUE) {
    obs_LB <-
      sample(which(cl_train_label_noise %in% c(2, 3)), size = 120)
    cl_train_label_noise[obs_LB] <-
      ifelse(test = cl_train_label_noise[obs_LB] == 2, 3, 2) # add label noise
  }
  
  #Y_TRAINING
  X <- 
    purrr::map_dfr(1:G,  ~ as.data.frame(mvnfast::rmvn(
      n = N_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]],
    )))
  
  #Y_ALL
  Y <- 
    purrr::map_dfr(1:length(M_vec),  ~ as.data.frame(mvnfast::rmvn(
      n = M_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
    ))) #NEW: 1:length(M_vec), OLD 1:H, BUG/svista?
  
  #updated labels_training
  cl_train_label_noise
  
  #PLOTS
  #plots training and all dataset
  # x11()
  # plot(X, col=cl_train_label_noise)
  # x11()
  # plot(Y, col=cl_tot)
  
  generated_data <- list("X" = X,"Y" = Y, "p" = p, "cl_train_label_noise" = cl_train_label_noise, "G" = G , "cl_tot" = cl_tot)
  
  return(generated_data)
}

## Test ##
# N_FACTOR = 1
# p = 2
# easy = 0
# generated_data = generate_data(N_FACTOR, p, easy)
# print(generated_data)




### SAVE.R ###

## Save Generated Data ##
save_generated_data <- function(SEED, N_FACTOR, p, easy, generated_data) {
  base_filename = paste('_',SEED,'_',N_FACTOR,'_',p,'_',easy)
  
  write.table(generated_data$Y,paste("Y",base_filename,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(generated_data$X,paste("Y_training",base_filename,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(generated_data$cl_train_label_noise, paste("labels_training",base_filename,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(generated_data$cl_tot,paste("labels_tot",base_filename,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
}

## Test ##
# setwd("/home/eb/Desktop/comparison/results/test")
# SEED = 0
# N_FACTOR = 1
# p = 2
# easy = 0
# generated_data = generate_data(N_FACTOR, p, easy)
# save_generated_data(SEED, N_FACTOR, p, easy, generated_data)


## Save output ##
save_output <- function(algo, SEED, n_factor, p, easy, output) {
  title = paste(algo,':','SEED=',SEED,'N_FACTOR=',n_factor,'p=',p,'easy=', easy,sep=' ')
  subtitle = paste("SEED", "algo","N_FACTOR","p"," easy", "elapsed_time_mins", "tot_clusters", "cluster_found","ARI","FMI","AMI",sep='\t' )
  txt = paste(title,subtitle,sep='\n')
  
  base_filename = paste('_',algo,'_',SEED,'_',n_factor,'_',p,'_',easy)
  
  writeLines(txt, paste("outfile",base_filename,'.txt',  sep = ""))
  
  intro = list("SEED"= SEED, "algo"= algo, "N_FACTOR"= n_factor, "p"= p, "easy" = easy)
  final_output = append(intro, output)
  
  write.table(final_output, paste("outfile",base_filename,'.txt',  sep = ""),sep="\t",
              row.names=FALSE,col.names=FALSE, append = TRUE)
  
  
  write.table(final_output, paste("output_csv",base_filename,'.csv',  sep = ""),sep="\n",
              row.names=FALSE,col.names=FALSE, append = FALSE)
}

## Test ##
# setwd("/home/eb/Desktop/comparison/results/test")
# algo = "BRAND"
# SEED = 0
# N_FACTOR = 0.1
# p = 2
# easy = 0
# 
# generated_data = generate_data(N_FACTOR, p, easy)
# output_brand = pico_brand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
# save_output(algo, SEED, N_FACTOR, p, easy, output_brand)



### PICO_BRAND.r ###

## Install Packages ##
#install.packages("adaptDA_1.0.tar.gz",repos = NULL, type = "source")
#install.packages("remotes")
#remotes::install_github("AndreaCappozzo/raedda")
#remotes::install_github("Fradenti/Brand")
#install.packages("mcclust")
#install.packages("mclust")
#remotes::install_github("sarawade/mcclust.ext")
#install.packages("dendextend")
#install.packages("aricode")

## Load Libraries ##
library(adaptDA) 
library(raedda) 
library(brand) 
library(mcclust) 
library(mclust) 
library(mcclust.ext) 
library(dendextend)
library(aricode)


#PICO_BRAND
pico_brand <- function(X,Y,p,cl_train_label_noise, G, cl_tot, SEED, easy) {
  #MCMC
  # truncation
  H <- 20 #97% of SB in E value, conti da controllare
  
  BURN_IN <- 20000 
  length_chain <- 60000 
  
  start_time = Sys.time()
  
  #num clusters after noise
  J <- length(unique(cl_train_label_noise))
  
  #hyperparameters
  alpha_MCD = .75 #c(.75, 1)
  k_tilde_train = 10+p+1 #c(10, 1000)
  
  prior <- list(
    aDir = c(rep(1, J), .1),
    aDP = .1,
    #(.5?)
    m_H = rep(0, p),
    k_H = .1,
    v_H = p+1+1, #smallest feasible
    S_H = diag(p) * 3,
    k_g = k_tilde_train,
    v_g = 20,
    a_alphaDP = 1,
    b_alphaDP = 5
  )
  
  #brand
  fit_adapt_RDA_bnp <-
    Brand_mlvt(
      Y = Y,
      X = X,
      categ = cl_train_label_noise,
      prior = prior,
      L = 20,
      burn_in = BURN_IN,
      thinning = 1,
      nsim = length_chain,
      fixed_alphaDP = FALSE,
      h_MCD = alpha_MCD,
      raw_MCD = FALSE,
      kappa = .25,
      learning_type = "transductive",
      light = TRUE
    )
  
  cl_adapt_BNP <-
    apply(fit_adapt_RDA_bnp$AB[, 1, ], 1, major.vote)
  
  novelty_adapt_BNP <- which(cl_adapt_BNP == 0)
  
  BET    <- fit_adapt_RDA_bnp$AB[novelty_adapt_BNP, 2,]
  psmBET <- comp.psm(t(BET) + 1)
  
  cl_beta_VI <- minVI(psmBET)$cl
  cl_beta_adapt_BNP <- cl_adapt_BNP
  cl_beta_adapt_BNP[cl_adapt_BNP == 0] <-
    cl_beta_VI + G # I add G to separate from the original tr labels
  
  a_posteriori_prob_novelty <-
    apply(fit_adapt_RDA_bnp$AB[, 1, ] == 0, 1, mean)
  
  result <-
    list(
      cl_alpha = cl_adapt_BNP,
      cl_beta = cl_beta_adapt_BNP,
      a_posteriori_prob_novelty = a_posteriori_prob_novelty,
      ari = mclust::adjustedRandIndex(cl_tot, cl_beta_adapt_BNP),
      cluster_found = length(unique(cl_beta_adapt_BNP))
    )
  
  
  end_time = Sys.time()
  elapsed_time_mins = end_time - start_time
  
  
  #REPORT 
  n = length(result$cl_beta)
  
  tot_clusters = length(unique(cl_tot))
  cluster_found = length(unique(result$cl_beta))
  ARI = result$ari
  FMI = FM_index(cl_tot, cl_beta_adapt_BNP)[1]
  ami = AMI(cl_tot, cl_beta_adapt_BNP)
  
  output = list("elapsed_time_mins" = elapsed_time_mins, "tot_clusters" = tot_clusters, "cluster_found" = cluster_found, "ARI" = ARI, "FMI" = FMI, "ami" = ami)
  return(output)
}

## Test ##
# N_FACTOR = 0.1
# p = 2
# easy = 0
# generated_data = generate_data(N_FACTOR, p, easy)
# output_brand = pico_brand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
# print(output_brand)




### PICO_VARBRAND.R ###

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
  
  print(result$elbo_values)
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



#install.packages("matrixStats")
library(matrixStats)

# REPORT 
avg_results = function(timestamp){
  filenames = list.files(pattern="*output_csv")
  
  data <- data.frame(matrix(ncol = 11, nrow = length(filenames)))
  colnames(data) <- c("SEED","algo","N_FACTOR", "p", "easy", "elapsed_time_secs", "tot_clusters", "cluster_found","ARI" ,"FMI","AMI")
  
  for (i in 1:length(filenames)){
    temp = read.table(filenames[i], header = FALSE)
    print(as.vector(temp$V1))
    data[i,] = as.vector(temp$V1)
  }
  print("filenames")
  print(filenames)
  
  print("data")
  print(data)
  
  combinations = unique(data[,c('algo','N_FACTOR', 'p','easy')])
  print("combinations")
  print(combinations)
  
  report <- data.frame(matrix(ncol = 13, nrow = nrow(combinations)))
  
  colnames(report) <- c("n_iter",'algo',"N_FACTOR", "p", "easy", "elapsed_time_secs_mean","elapsed_time_secs_sd", "mean_ARI","sd_ARI","mean_FMI", "sd_FMI","mean_AMI", "sd_AMI")
  
  i= 1 
  while(i <= nrow(combinations)){
    indexes = which((data$algo == combinations$algo[i]) & (data$N_FACTOR == combinations$N_FACTOR[i]) & (data$p == combinations$p[i]) & (data$easy == combinations$easy[i]))
    print("indexes")
    print(indexes)
    
    report[i,1] = nrow(data[indexes,])
    report[i,2:5] = combinations[i,1:4]
    report[i,6] = mean(as.numeric(data[indexes,6]))
    report[i,7] = sd(as.numeric(data[indexes,6]))
    report[i,8] = mean(as.numeric(data[indexes,9]))
    report[i,9] = sd(as.numeric(data[indexes, 9]))
    report[i,10] = mean(as.numeric(data[indexes,10]))
    report[i,11] = sd(as.numeric(data[indexes, 10]))
    report[i,12] = mean(as.numeric(data[indexes,11]))
    report[i,13] = sd(as.numeric(data[indexes, 11]))
    
    
    i=i+1
  }
  
  print("report")
  print(report)
  
  #setwd('C:/Users/user/Desktop/VariationalBRAND-cpp/tests/varbrand')
  report_filename = paste('report','_',timestamp,'.csv',sep = "")
  write.csv(report,report_filename, row.names = FALSE)
  
}


## Test ##
# setwd("/home/eb/Desktop/comparison/results/test_avg")
# 
# i = 0
# while(i<2){
#   seed = i
#   set.seed(seed)
#   n_factor_list = c(0.5) # c(0.5, 1, 2.5, 5, 10)
#   p_list= c(2) #c(2, 3, 5, 7, 10)
#   easy_list = c(0)# c(0, 1)
# 
#   for(n_factor in n_factor_list){
#     for(p in p_list){
#       for(easy in easy_list){
#         generated_data = generate_data(n_factor, p, easy)
#         save_generated_data(seed, n_factor, p, easy, generated_data)
# 
#         output_brand = pico_brand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
#         save_output("BRAND", seed, n_factor, p, easy, output_brand)
# 
#         output_varbrand = pico_varbrand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
#         save_output("VARBRAND", seed, n_factor, p, easy, output_varbrand)
# 
#       }
#     }
#   }
#   i=i+1
# }

#avg_results("00:00-test")









### PARALLEL_SIMULATIONS.R
## Install Packages ##
#-> parallel_simulations.R
#install.packages("foreach")
#install.packages("parallel")

#=> generate_data.R <=
#install.packages("mvnfast")
#install.packages("purrr")
#install.packages("dplyr")

#=> pico_brand.R <=
#install.packages("adaptDA_1.0.tar.gz",repos = NULL, type = "source")
#install.packages("remotes")
#remotes::install_github("AndreaCappozzo/raedda")
#remotes::install_github("Fradenti/Brand")
#install.packages("mcclust")
#install.packages("mclust")
#remotes::install_github("sarawade/mcclust.ext")
#install.packages("dendextend")
#install.packages("aricode")

#=> pico_varbrand.R <=
#install.packages("Rcpp")
#install.packages("RcppEigen")
#install.packages("foobar2_1.0.tar.gz",repos = NULL, type = "source")
#install.packages("ramify")

##=> avg_results.R <=
#install.packages("matrixStats")


## Load parallel libraries ##
library(foreach)
library(parallel)

## Init ##
base = "/home/eb/Desktop/comparison" # <= ! change with your folder ! <==
setwd(base) 

#- Load functions
# source("generate_data.R")
# source("save.R")
# source("pico_brand.R")
# source("pico_varbrand.R")

#- Create results' folder for current test
timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
dir.create("results")
setwd("results")
dir.create(timestamp)
setwd(timestamp)

## Parallel Simulations ##  <= ! change with your parallel parameters ! <==
NSIM <- 2 #50
i <- NULL # dummy to trick R CMD check
cl <- makeCluster(2, outfile="") #makeCluster(32, outfile="") # # of cores
doParallel::registerDoParallel(cl)

placeholder <- foreach(
  seed = 1:NSIM,
  # packages that need to be imported
  .packages = c(
    "mvnfast",
    "purrr",
    "dplyr",
    "adaptDA",
    "raedda",
    "brand",
    "mcclust",
    "mclust",
    "mcclust.ext",
    "dendextend",
    "aricode",
    "Rcpp",
    "RcppEigen",
    "foobar2",
    "matrixStats"
  )
) %dopar% {
  
  ## Main Simulation Function ##  <= ! change with your simulation parameters ! <==
  set.seed(seed)
  n_factor_list = c(0.5,1) # c(0.5, 1, 2.5, 5, 10)
  p_list= c(2,3) #c(2, 3, 5, 7, 10)
  easy_list = c(0)# c(0, 1)
  
  for(n_factor in n_factor_list){
    for(p in p_list){
      for(easy in easy_list){
        generated_data = generate_data(n_factor, p, easy)
        save_generated_data(seed, n_factor, p, easy, generated_data)
        
        output_brand = pico_brand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
        save_output("BRAND", seed, n_factor, p, easy, output_brand)
        
        output_varbrand = pico_varbrand(generated_data$X, generated_data$Y, generated_data$p ,generated_data$cl_train_label_noise, generated_data$G, generated_data$cl_tot)
        save_output("VARBRAND", seed, n_factor, p, easy, output_varbrand)
        
      }
    }
  }
}
# Generate avg results
avg_results(timestamp)

# reset working directory to base
setwd(base) 
