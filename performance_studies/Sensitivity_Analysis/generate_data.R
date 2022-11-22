library(MASS)

# Function to build the datasets used for the simulation. It defines 5 families 
# plus a 6th which is an outlier group
build_df <- function(mean_param, var_param, eig1, eig2, outlier_param, numerosity_vec) {
  # Initialize means
  mu1 = c(-5,5)
  mu2 = c(5,5)
  mu3 = c(5,-5)
  mu4 = c(0,0)
  mu5 = c(0,0)
  mu6 = c(-5,-5)
  
  # Initialize covariances
  Id = diag(2)
  ones = matrix(1,2,2)
  negones = ones
  negones[1,2]= -1
  negones[2,1]= -1
  eigval1_A = eig1
  eigval2_A = eig2
  eigval1_B = eig2
  eigval2_B = eig1
  
  S1 = Id
  S4 = eigval1_A*ones+eigval2_A*negones
  S5 = eigval1_B*ones+eigval2_B*negones
  
  # Sample from the distributions
  gr1 = mvrnorm(n = numerosity_vec[1], mean_param*mu1, var_param*S1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  gr2 = mvrnorm(n = numerosity_vec[1], mean_param*mu2, var_param*S1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  gr3 = mvrnorm(n = numerosity_vec[1], mean_param*mu3, var_param*S1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  gr4 = mvrnorm(n = numerosity_vec[1], mean_param*mu4, var_param*S4, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  gr5 = mvrnorm(n = numerosity_vec[1], mean_param*mu5, var_param*S5, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  gr6 = mvrnorm(n = numerosity_vec[2], mean_param*mu6, outlier_param*var_param*S1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  # Aggregate the samples in a single matrix
  tmp_data = rbind(gr1,gr2)
  tmp_data = rbind(tmp_data,gr3)
  tmp_data = rbind(tmp_data,gr4)
  tmp_data = rbind(tmp_data,gr5)
  all_data = rbind(tmp_data,gr6)
  
  # Build the labels vector
  idx1 = 1:numerosity_vec[1]                                              # known group 1
  idx2 = (numerosity_vec[1]+1):(2*numerosity_vec[1])                      # known group 2
  idx3 = (2*numerosity_vec[1]+1):(3*numerosity_vec[1])                    # known group 3
  idx4 = (3*numerosity_vec[1]+1):(4*numerosity_vec[1])                    # novelty 1
  idx5 = (4*numerosity_vec[1]+1):(5*numerosity_vec[1])                    # novelty 2
  idx6 = (5*numerosity_vec[1]+1):(5*numerosity_vec[1]+numerosity_vec[2])  # outliers
  
  labels = 1:(5*numerosity_vec[1]+numerosity_vec[2])
  labels[idx1] = 1
  labels[idx2] = 2
  labels[idx3] = 3
  labels[idx4] = 4
  labels[idx5] = 5
  labels[idx6] = 6
  
  # Choose the elements in the data for the training part
  idx1 = 1:50                                                                # known group 1
  idx2 = (numerosity_vec[1]+1):(numerosity_vec[1]+50)                        # known group 2
  idx3 = (2*numerosity_vec[1]+1):(2*numerosity_vec[1]+50)                    # known group 3
  
  train = rep(0,(5*numerosity_vec[1]+numerosity_vec[2]))
  train[idx1] = 1
  train[idx2] = 1
  train[idx3] = 1
  
  # Build the dataframe with points and labels
  df <- data.frame(all_data, labels, train)
  
  return(df)
}

# Function to handle the build_df function inside the for loops 
generate_data <- function(var_dummy, numerosity_dummy) {
  # Use build_df to actually generate the datasets
  mean_param = 1
  
  # We set two possible levels for the variance parameter, according to our previous studies
  if(var_dummy == 0) {
    var_param = 0.5
  } else if(var_dummy == 1) {
    var_param = 3
  } else {
    print("Wrong var input, the value should be either 0 or 1")
  }
  
  # We set the parameter for the outlier. In this way the variance parameter for 
  # the group of the outliers is 3*(var_param)^2
  outlier_param = 3*var_param
  
  # Those values are such that all the groups are nearly symmetric in all the directions
  eig1 = 0.9
  eig2 = 0.1
  
  # We set two possible levels for the numerosity parameters, according to our previous studies
  if(numerosity_dummy == 0) {
    numerosity_vec = c(100,5)
  } else if(numerosity_dummy == 1) {
    numerosity_vec = c(2000,100)
  } else {
    print("Wrong numerosity input, the value should be either 0 or 1")
  }
  
  df = build_df(mean_param, var_param, eig1, eig2, outlier_param, numerosity_vec)
  
  return(df)
}