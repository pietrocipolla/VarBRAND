## Install Packages ##
#-> parallel_simulations.R
#install.packages("foreach")
#install.packages("parallel")

#=> generate_data.R <=
#install.library("MASS")

#=> sa_varbrand.R <=
#install.packages("Rcpp")
#install.packages("RcppEigen")
#install.packages("foobar2_1.0.tar.gz",repos = NULL, type = "source")
#install.packages("ramify")
#install.packages("rrcov")
#install.library("aricode")
#install.library("FME")

##=> avg_results.R <=
#install.packages("matrixStats")


## Load parallel libraries ##
library(foreach)
library(parallel)
library(aricode)
library(here)

## Init ##
base = here("SA") # <= ! change with your folder ! <==
setwd(base) 

#- Load functions
source("generate_data.R")
source("save.R")
source("sa_varbrand.R")
source("space_exploration.R")

#- Create results' folder for current test
timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
dir.create("results")
setwd("results")
dir.create(timestamp)
setwd(timestamp)

## Parallel Simulations ##  <= ! change with your parallel parameters ! <==
# NSIM <- 4 #50
i <- NULL # dummy to trick R CMD check
# cl <- makeCluster(4, outfile="") #makeCluster(32, outfile="") # # of cores
# doParallel::registerDoParallel(cl)

placeholder <- foreach(
  d = c(1:4),
  # packages that need to be imported
  .packages = c(
    "MASS",
    "rrcov",
    "ramify",
    "aricode",
    "Rcpp",
    "RcppEigen",
    "foobar2",
    "matrixStats"
  )
) %dopar% {
  ## Main Simulation Function ##  <= ! change with your simulation parameters ! <==
  set.seed(d+42)
  var_dummy = d %% 2
  numerosity_dummy = as.integer(d > 2)
  N_rep <- 50
  print(paste("Process ",d ," var_dummy: ",var_dummy, " numerosity_dummy: ",numerosity_dummy))
  
  for(n in 1:N_rep){
    generated_data = generate_data(var_dummy, numerosity_dummy)
    # save_generated_data(d, n, var_dummy, numerosity_dummy, generated_data)
        
    output = run_sa(generated_data)
    save_output(d, n, var_dummy, numerosity_dummy, output)
    
    print(paste("Run ",n ," done for core ",d))
  }
}

# reset working directory to base
setwd(base) 