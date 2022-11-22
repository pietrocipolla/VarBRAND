#Variables taken into account: gamma, lambda_0_dp, nu_0_dp, a_k_dir

model_eval <- function(Y_tot, labels_tot, p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                       lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                       PSI_0_MIX,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                       lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX) {
  
  # Evaluation of the model at the desidered point
  X <- foobar2::VarBrand(as.matrix(Y_tot), p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                         lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                         PSI_0_MIX,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                         lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX)
  labels_pred <- argmax(X$Phi_m_k, rows = TRUE)
  
  #Se il risultato è una matrice di Na restituisco Na altrimenti calcolo l'ari
  if(anyNA(X$Phi_m_k)) {
    ari <- NaN
    f1 <- NaN
  }
  else {
    ari <- ARI(labels_pred,labels_tot)
    f1 <- F1(labels_pred, labels_tot, J)
  }
  
  return(c(ari,f1))
}

explore_space <- function(Y_tot, Y_training, labels_tot, labels_training, p, J, Tr, n_iter, tol, gamma, 
                                 a_dir_k_in, mu_0_DP, PSI_0_DP, mu_0_MIX, nu_0_MIX, lambda_0_MIX, PSI_0_MIX, 
                                 PSI_0_MIX_FACTOR, M, Phi_m_k, eta_k, a_k_beta, b_k_beta, mu_var_DP, 
                                 PSI_var_DP, mu_VAR_MIX, nu_VAR_MIX, lambda_VAR_MIX, PSI_VAR_MIX, n_points) {
  
  #Generiamo i punti equidistanti nello spazio, al momento tra 1 e 10, in base all'input ricevuto
  gamma.space <- seq(1,10, length.out = n_points)
  lambda_0_DP.space <- seq(1,10, length.out = n_points)
  nu_0_DP.space <- seq(p+2,100, length.out = n_points)
  a_k_dir.space <- seq(0.1,1, length.out = n_points)
  
  #Costruiamo un dataset con tutti i valori da dare in input 
  X <- expand.grid(gamma.space, lambda_0_DP.space, nu_0_DP.space, a_k_dir.space)
  colnames(X) <- c("gamma", "lambda_0_DP", "nu_0_DP", "a_k_dir")
  n <- dim(X)[1]
  
  #Costruiamo la matrice degli output
  Y_out <- matrix(nrow = n,ncol = 2)
  
  #Per ogni ingresso del dataset estraiamo i valori e valutiamo il modello. Il valore in output è l'ari
  for(i in 1:n) {
    gamma <- X$gamma[i]
    lambda_0_DP <- X$lambda_0_DP[i]
    nu_0_DP <- X$nu_0_DP[i]
    a_dir_k <- X$a_k_dir[i] * a_dir_k_in
    
    b_k_beta = rep(gamma, Tr - 1)
    nu_var_DP = rep(nu_0_DP, Tr)
    lambda_var_DP = rep(lambda_0_DP, Tr)
    
    Y_out[i,] <- model_eval(as.matrix(Y_tot), labels_tot, p,J,Tr,n_iter,tol,gamma,a_dir_k, mu_0_DP, nu_0_DP,
                                 lambda_0_DP, PSI_0_DP,mu_0_MIX,nu_0_MIX,lambda_0_MIX,
                                 PSI_0_MIX,M,Phi_m_k,eta_k,a_k_beta,b_k_beta,mu_var_DP,nu_var_DP,
                                 lambda_var_DP,PSI_var_DP,mu_VAR_MIX,nu_VAR_MIX,lambda_VAR_MIX,PSI_VAR_MIX)
  }
  
  for(i in 1:dim(Y_out)[2]) X[,i+4] <- Y_out[,i]
  colnames(X) <- colnames(X) <- c("gamma", "lambda_0_DP", "nu_0_DP", "a_k_dir", "ARI", "F1")
  return(X)
}  

F1 <- function(labels_pred, labels_tot, J) {
  binary_true <- labels_tot < (J+1)
  binary_pred <- labels_pred < (J+1)
  
  true_pos <- sum(binary_true & binary_pred)
  true_neg <- sum(!binary_true & !binary_pred)
  false_pos <- sum(!binary_true & binary_pred)
  false_neg <- sum(binary_true & !binary_pred)
  
  f1 <- true_pos/(true_pos + 0.5*(false_pos + false_neg))
  return(f1)
}
