## Save Generated Data ##

save_generated_data <- function(d, n, var_dummy, num_dummy, generated_data) {
  base_filename = paste('_',d,'_',var_dummy,'_',num_dummy,'_',n)
  
  write.table(generated_data,paste("dataset",base_filename,'.csv',  sep = ""), row.names = FALSE, sep=',')
}

## Save output ##
save_output <- function(d, n, var_dummy, num_dummy, output) {
  base_filename = paste('_',d,'_',var_dummy,'_',num_dummy,'_',n)
  
  write.table(output, file = paste("output", base_filename, ".csv", sep = ""), row.names = FALSE, sep=',')
}