#####install required packages

#=> generate_data.R <=
#install.packages("mvnfast")
#install.packages("purrr")
#install.packages("dplyr")

#=> pico_varbrand.R <=
#install.packages("Rcpp")
#install.packages("RcppEigen")
#install.packages("ramify")
#install.packages("rrcov")

#install.packages("packages//variational//foobar2_1.0.tar.gz",repos = NULL, type = "source")
#package installation procedure may be machine-dependent due safety options
#many issues can be solved by following the procedure reported in https://stackoverflow.com/questions/42807247/installing-package-cannot-open-file-permission-denied

##=> avg_results.R <=
#install.packages("matrixStats")
#####


library(foobar2)

#GENERATE DATA
#p = number of dimensions, N_FACTOR = inflation factor for dataset cardinalities
generate_data <- function(N_FACTOR, p) {
  # number of training set realizations
  N <- 1000 * N_FACTOR
  
  # number of test set realizations
  M <- 1000 * N_FACTOR
  
  # training set class cardinalities
  N_vec <- list(300, 300, 400)
  N_vec = lapply(N_vec,"*",N_FACTOR)
  
  # test set class cardinalities
  M_vec <- list(200, 200, 250, 90, 100, 100, 10+50)
  M_vec = lapply(M_vec,"*",N_FACTOR)
  
  # clusters in the training set
  G <- 3
  
  # training set labels
  cl_train <- rep(1:G, N_vec)
  
  # test set labels
  cl_tot <- rep(1:length(M_vec), M_vec)
  
  # class centers and variances
  MU <-
    list(c(-5, 5),
         c(-4, -4),
         c(4, 4),
         c(-0, 0),
         c(5, -10),
         c(5, -10),
         c(-10, -10))
  SIGMA <-
    list(
      matrix(c(1, .9, .9, 1), 2, 2),
      diag(2),
      diag(2),
      matrix(c(1, -.75, -.75, 1), 2, 2),
      matrix(c(1, .9, .9, 1), 2, 2),
      matrix(c(1, -.9, -.9, 1), 2, 2),
      diag(.01, 2)
    )
  
  
  # dimensionality augmentation for p > 2
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
  
  desired_length <- length(SIGMA) 
  updated_SIGMA <- vector(mode = "list", length = desired_length)
  
  i = 1
  while(i<length(SIGMA)+1){
    sigma = diag(p)
    sigma[1:2,1:2] = SIGMA[[i]]
    updated_SIGMA[[i]]=sigma
    i=i+1
  }
  
  SIGMA = updated_SIGMA
  
  # training set labels
  cl_test <- rep(1:length(M_vec), M_vec)
  
  # test set label noise
  label_noise <- FALSE
  cl_train_label_noise <- cl_train
  if (label_noise == TRUE) {
    obs_LB <-
      sample(which(cl_train_label_noise %in% c(2, 3)), size = 120)
    cl_train_label_noise[obs_LB] <-
      ifelse(test = cl_train_label_noise[obs_LB] == 2, 3, 2) # add label noise
  }
  
  # generate training set
  X <- 
    purrr::map_dfr(1:G,  ~ as.data.frame(mvnfast::rmvn(
      n = N_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]],
    )))
  
  # generate test set
  Y <- 
    purrr::map_dfr(1:length(M_vec),  ~ as.data.frame(mvnfast::rmvn(
      n = M_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
    )))
  
  output_df <- list("X" = X,"Y" = Y, "p" = p, "cl_train_label_noise" = cl_train_label_noise, "G" = G , "cl_tot" = cl_tot)
  
  return(output_df)
}



#SAVE DATA
save_data <- function(X,Y,p,cl_train_label_noise, G, cl_tot, SEED, N_FACTOR) {
  working_directory <- getwd()
  setwd(paste(working_directory, "//bin", sep = ""))
  #save data
  write.table(Y,paste("Y_",SEED,'_',N_FACTOR,'_',p,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(X,paste("Y_training_",SEED,'_',N_FACTOR,'_',p,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(cl_train_label_noise, paste("labels_training_",SEED,'_',N_FACTOR,'_',p,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  write.table(cl_tot,paste("labels_tot_",SEED,'_',N_FACTOR,'_',p,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  setwd(working_directory)
}


#PICO_VARBRAND
pico_varbrand <- function(X,Y,p,cl_train_label_noise, G, cl_tot, SEED) {

  # libraries
  library(foobar2)
  library(Rcpp)
  
  library(dendextend)
  library(aricode)
  
  #start timer
  start_time = Sys.time()
  
  # load datasets
  labels_tot = cl_tot
  labels_training = cl_train_label_noise
  Y_training = X
  Y_tot = Y
  
  # simulation parameter
  p = p #dim(Y_tot)[2]
  
  #Number of classes in the training set
  J = G #dim(unique(labels_training))[1]
  
  # truncation parameter in the Dirichlet Process
  Tr = 20
  
  #Number of iterations and tolerance for CAVI 
  n_iter = 1000
  tol = 1e-5
  
  #Prior hyperparameters
  #gamma -> Stick Breaking parameter
  gamma = 5
  
  #a_dir_k -> Dirichlet hyperparameters
  a_dir_k = rep(0.1, J + 1)
  
  #nIW_DP_0 -> novelty prior distribution
  # 	mu_0_DP -> location
  # 	nu_0_DP -> degrees of freedom for IW
  # 	lambda_0_DP -> scaling coefficient for normal variance
  # 	PSI_0_DP -> inverse scale matrix
  
  mu_0_DP = rep(0, p)
  
  nu_0_DP = c(p+1+1)
  
  lambda_0_DP = c(0.1)
  
  PSI_0_DP = 10*matrix(diag(p), nrow = p, ncol = p)
  
  
  #NIW_MIX_0 -> Lists of known mixture components' priors hyperparameter
  mu_0_MIX = NULL #set to robust centers after initialization

  nu_0_MIX = rep(50+p+1,J)
  
  lambda_0_MIX = rep(50,J)
  
  PSI_0_MIX = NULL #set to robust variances after initialization
  PSI_0_MIX_FACTOR = 50 #inflation factor on robust estimates to keep IW expected value truthful
  

  #Variational parameters initialization
  M = dim(Y_tot)[1]
  Phi_m_k = matrix(1, M, J + Tr) * (1 / (J + T))
  
  eta_k = a_dir_k
  
  a_k_beta = rep(1, Tr - 1)
  
  b_k_beta = rep(gamma, Tr - 1)
  
  mu_var_DP = t(kmeans(Y_tot, Tr)$center) 
  nu_var_DP = rep(nu_0_DP, Tr)
  lambda_var_DP = rep(lambda_0_DP, Tr)
  PSI_var_DP = list()
  for (i in 1:Tr) PSI_var_DP[[i]] = PSI_0_DP

  mu_VAR_MIX = NULL
  nu_VAR_MIX = nu_0_MIX
  lambda_VAR_MIX = lambda_0_MIX
  PSI_VAR_MIX = NULL 
  
  
  #call VarBRAND, from initialization returns values
  result = foobar2::var_brand(Y_tot, Y_training, labels_tot, labels_training, p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                              lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                              PSI_0_MIX,PSI_0_MIX_FACTOR,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                              lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX)
  
  # Collect results 
  #predicted labels
  labels_pred = result$labels_pred
  n = length(labels_pred)

  # clusters found on test set
  tot_clusters = length(unique(cl_tot))
  cluster_found = length(unique(labels_pred))
  
  # performance metrics
  ARI = result$ari
  FMI = FM_index(labels_tot, labels_pred)[1]
  ami = AMI(labels_tot, labels_pred)
  
  
  #stop timer
  end_time = Sys.time()
  elapsed_time_mins = end_time - start_time

  #output
  out = cbind(SEED, n, p, elapsed_time_mins, tot_clusters, cluster_found,ARI, FMI, ami )

  # save output
  working_directory <- getwd()
  setwd(paste(working_directory, "//bin", sep = ""))
  title = paste('VARBRAND',':','n=',n,'p=', p,sep=' ')
  subtitle = paste("SEED", "n"," p", "elapsed_time_mins", "tot_clusters", "cluster_found","ARI","FMI","AMI",sep='\t' )
  txt = paste(title,subtitle,sep='\n')

  writeLines(txt, paste("outfile_",SEED,'_',n,'_',p,'.txt',  sep = ""))
  write.table(out, paste("outfile_",SEED,'_',n,'_',p,'.txt',  sep = ""),sep="\t",
              row.names=FALSE,col.names=FALSE, append = TRUE)

  write.table(out, paste("output_csv",SEED,'_',n,'_',p,'.csv',  sep = ""),sep="\n",
              row.names=FALSE,col.names=FALSE, append = FALSE)

  #jpeg(file="training_dataset.jpeg")
  jpeg(file=paste("training_dataset_",SEED,'_',n,'_',p,'.jpeg',  sep = ""))
  plot(X[,c(1,2)], col=cl_train_label_noise)
  dev.off()

  #jpeg(file="test_dataset.jpeg")
  jpeg(file=paste("test_dataset_",SEED,'_',n,'_',p,'.jpeg',  sep = ""))
  plot(Y[,c(1,2)], col=cl_tot)
  dev.off()

  #jpeg(file="output_brand_.jpeg")
  jpeg(file=paste("output_brand_",SEED,'_',n,'_',p,'.jpeg',  sep = ""))
  plot(Y[,c(1,2)], col=labels_pred)
  dev.off()

  write.table(labels_pred, paste("labels_pred_",SEED,'_',n,'_',p,'.csv',  sep = ""), row.names = FALSE, col.names=FALSE, sep=',')
  
  setwd(working_directory)
}

# average simulations results
avg_results = function(){
  working_directory <- getwd()
  setwd(paste(working_directory, "//bin", sep = ""))
  
  filenames = list.files(pattern="*output_csv")
  
  data <- data.frame(matrix(ncol = 9, nrow = length(filenames)))
  colnames(data) <- c("SEED", "n", "p", "elapsed_time_secs", "tot_clusters", "cluster_found","ARI" ,"FMI","AMI")
  
  for (i in 1:length(filenames)){
      temp = read.table(filenames[i], header = FALSE)
      print(as.vector(temp$V1))
      data[i,] = as.vector(temp$V1)
  }
  print(filenames)
  
  print(data)
  
  combinations = unique(data[,c('n','p')])
  print(combinations)
  
  library(matrixStats)
  report <- data.frame(matrix(ncol = 11, nrow = nrow(combinations)))

  colnames(report) <- c("n_trials","n", "p", "elapsed_time_secs_mean","elapsed_time_secs_sd", "mean_ARI","sd_ARI","mean_FMI", "sd_FMI","mean_AMI", "sd_AMI")

  i= 1 
  while(i <= nrow(combinations)){
    indexes = which((data$n == combinations$n[i]) & (data$p == combinations$p[i]))
    
    report[i,1] = nrow(data[indexes,])
    report[i,2:3] = combinations[i,1:2]
    report[i,4] = mean(data[indexes,4])
    report[i,5] = sd(data[indexes,4])
    report[i,6] = mean(data[indexes,7])
    report[i,7] = sd(data[indexes, 7])
    report[i,8] = mean(data[indexes,8])
    report[i,9] = sd(data[indexes, 8])
    report[i,10] = mean(data[indexes,9])
    report[i,11] = sd(data[indexes, 9])
    
    
    i=i+1
  }

  

  setwd(paste(working_directory, "//bin//reports", sep = ""))
  write.csv(report,"full_VB_report.csv", row.names = FALSE)
  setwd(working_directory)
}
