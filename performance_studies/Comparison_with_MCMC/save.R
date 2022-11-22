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











