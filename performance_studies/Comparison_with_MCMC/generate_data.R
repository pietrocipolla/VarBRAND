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



