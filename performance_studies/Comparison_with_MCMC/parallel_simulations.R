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
library(here)

## Init ##
base = here("")
setwd(base) 

#- Load functions
source("generate_data.R")
source("save.R")
source("pico_brand.R")
source("pico_varbrand.R")

#- Create results' folder for current test
timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
dir.create("results")
setwd("results")
dir.create(timestamp)
setwd(timestamp)

## Parallel Simulations ##  <= ! change with your parallel parameters ! <==
# NSIM <- 2 #50
NSIM <- 50
i <- NULL # dummy to trick R CMD check
# cl <- makeCluster(2, outfile="") #makeCluster(32, outfile="") # # of cores
cl <- makeCluster(32, outfile="") # # of cores
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
  # n_factor_list = c(0.5,1) 
  n_factor_list =  c(0.5, 1, 2.5, 5, 10)
  # p_list= c(2,3) #c(2, 3, 5, 7, 10)
  p_list= c(2, 3, 5, 7, 10)
  # easy_list = c(0)# c(0, 1)
  easy_list = c(0, 1)
  
  for(n_factor in n_factor_list){
    for(p in p_list){
      for(easy in easy_list){
        generated_data = generate_data(n_factor, p, easy)
        # save_generated_data(seed, n_factor, p, easy, generated_data)
        
        output_varbrand = pico_varbrand(
          generated_data$X,
          generated_data$Y,
          generated_data$p ,
          generated_data$cl_train_label_noise,
          generated_data$G,
          generated_data$cl_tot
        )
        save_output("VARBRAND", seed, n_factor, p, easy, output_varbrand)
        
        output_brand = tryCatch(pico_brand(
          generated_data$X,
          generated_data$Y,
          generated_data$p ,
          generated_data$cl_train_label_noise,
          generated_data$G,
          generated_data$cl_tot
        ),error=function(e) "MCMC did not converge")
        save_output("BRAND", seed, n_factor, p, easy, output_brand)
        
      }
    }
  }
}
# Generate avg results
avg_results(timestamp)

# reset working directory to base
setwd(base) 
