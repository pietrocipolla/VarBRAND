setwd("~/Documents/Polimi/VariationalBRAND-cpp/results/soil")
wd <- "***//VarBRAND"

#Packages

#install.packages("packages//variational//foobar2_1.0.tar.gz",repos = NULL, type = "source")
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
library(tidyverse)
library(viridis)
library(caret)
library(Rtsne)
library(tidyverse)
library(ggplot2)

setwd(paste(wd, "//application", sep =""))
working_directory <- getwd()
set.seed(27011999)

# Load stage I 'n' stage II data
test_set <- read.table("sat_test.tst", quote="\"", comment.char="")
training_set <- read.table("sat_training.trn", quote="\"", comment.char="")
test_set<-test_set[,c(37, c(1:36))]
training_set<-training_set[,c(37, c(1:36))]

# Apply correction
test_set[,-1] <- test_set[,-1]/4.5
training_set[,-1] <- training_set[,-1]/4.5

# Define the classes for the clustering and map them in numbers
colnames(test_set)[1] <- "soil_type"
colnames(training_set)[1] <- "soil_type"

mapper <- cbind(c("red_soil", "cotton_crop", "grey_soil", "dump_gray_soil", "soil_with_veg", "very_dump_g_s"), c(1:5,7))
mapper <- data.frame(mapper)
colnames(mapper) <- c("soil", "code")

for(s in mapper$soil){
  test_set$soil_type[which(test_set$soil_type == mapper$code[which(mapper$soil==s)])] <- s
  training_set$soil_type[which(training_set$soil_type == mapper$code[which(mapper$soil==s)])] <- s
}

recoder <- data.frame(cbind(c("red_soil", "grey_soil", "dump_gray_soil", "very_dump_g_s", "cotton_crop", "soil_with_veg"), 1:6))
colnames(recoder) <- c("soil", "code")

for(s in recoder$soil){
  test_set$soil_type[which(test_set$soil_type == s)] <- recoder$code[which(recoder$soil==s)]
  training_set$soil_type[which(training_set$soil_type ==s)] <-recoder$code[which(recoder$soil==s)]
}

training_set <- training_set[-which(training_set$soil_type>=5),]

# Read the parameters from the hyperparameter search and choose the ones with the best ELBO
tot_par <- read.table("Hyperpar_search/output_parsearch_Leo.csv")
best_par <- best_par[tot_par$Best_elbo==max(tot_par$Best_elbo),1:4]

#PICO_VARBRAND
pico_varbrand_soil <- function(training_set, test_set, ID, best_par, rand_par, Tr) {
  set.seed(777)
  
  # LIBRARIES
  library(foobar2)
  library(Rcpp)
  
  library(dendextend)
  library(aricode)
  
  start_time = Sys.time()
  # LOAD DATASETS
  labels_tot = test_set[,1]
  labels_training = training_set[,1]
  Y_training = training_set[,-1]
  Y_tot = test_set[,-1]
  
  
  #SPECIFY USER INPUT
  #Number of components
  p = dim(Y_tot)[2]
  
  #Number of classes in the training set
  J = length(unique(labels_training))
  
  #Number of iterations and tolerance -  CAVI 
  n_iter = 1000
  tol = 1e-5
  
  gamma = 10

  
  a_dir_k = rep(0.1, J + 1)
  
  
  mu_0_DP = colMeans(Y_training)
  
  nu_0_DP = c(p+2)
  
  lambda_0_DP = c(0.1)
  
  #PSI_0_DP = 15*matrix(diag(p), nrow = p, ncol = p)
  PSI_0_DP = (p+1)*var(Y_training)
  
  
  mu_0_MIX = NULL 
  nu_0_MIX = rep(200+p+1,J)
  
  lambda_0_MIX = rep(200,J)
  
  PSI_0_MIX = NULL
  PSI_0_MIX_FACTOR = 200
  
  
  
  #VARIATIONAL PARAMETERS
  M = dim(Y_tot)[1]
  Phi_m_k = matrix(1, M, J + Tr) * (1 / (J + T))
  
  eta_k = a_dir_k*best_par$a_dir_k
  
  a_k_beta = rep(1, Tr - 1)
  
  b_k_beta = rep(gamma, Tr - 1)
  
  
  mu_var_DP = t(rand_par) 
  
  nu_var_DP = rep(nu_0_DP, Tr)*best_par$nu_0_DP
  
  lambda_var_DP = rep(lambda_0_DP, Tr)*best_par$lambda_0_DP
  
  PSI_var_DP = list()
  for (i in 1:Tr) PSI_var_DP[[i]] = PSI_0_DP
  
  mu_VAR_MIX = NULL
  
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
  
  
  tot_clusters = length(unique(labels_tot))
  cluster_found = length(unique(labels_pred))
  ARI = result$ari
  
  FMI = FM_index(labels_tot, labels_pred)[1]
  ami = AMI(labels_tot, labels_pred)
  
  
  end_time = Sys.time()
  elapsed_time_mins = end_time - start_time
  
  out = cbind(elapsed_time_mins, tot_clusters, cluster_found,ARI, FMI, ami )
  
  labels_tot <- as.factor(labels_tot)
  table(labels_tot, labels_pred)
  
  rand_par <- array(rand_par)
  return(c(ARI,FMI,rand_par))
}

#extract partition obtained by best run
N <- 1
output <- data.frame(matrix(ncol = 6, nrow = N))

Tr <- 10
kmeans_centroids <- kmeans(test_set[,-1], Tr, iter.max = 400)$center

for(i in 1:N) {
  rand_par <- kmeans_centroids
  output[i,] <- pico_varbrand_soil(training_set, test_set, "Mean_search", best_par, kmeans_centroids, Tr)
}
  

# tSNE
tSNE_fit <- Rtsne(Y_tot)

tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number()) %>%
  mutate(label = labels_tot)

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = labels_tot))+
  geom_point()+
  theme_gray() +
  theme(legend.position="bottom")

write.table(tSNE_df, file = "tSNE_fit_dataset.csv")

tot_par <- tot_par %>% 
  mutate(Run=row_number())

tot_par %>%
  ggplot(aes(x = Run, y = Best_elbo)) +
  geom_point() +
  theme_bw() + 
  geom_abline(slope = 0, intercept = max(tot_par$Best_elbo), color = "Red")

labelz <- data.frame(labels_tot, labels_pred)

write.table(output, file = "output_labelsbestpar_Leo.csv")
