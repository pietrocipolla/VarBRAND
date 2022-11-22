## Install Packages ##
#install.packages("packages//MCMC//adaptDA_1.0.tar.gz",repos = NULL, type = "source")
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
  
  BURN_IN <- 10000 
  length_chain <- 10000 
  
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
      light = TRUE,
      verbose=FALSE
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


