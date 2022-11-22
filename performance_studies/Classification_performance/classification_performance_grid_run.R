wd <- "***//VarBRAND"
setwd(paste(wd, "//performance_studies//Classification_performance", sep =""))
working_directory <- getwd()

source(paste(working_directory, "//classification_performance_utils.R", sep = ""))
#define search grids
SEED_LIST = 17021890*1:50
p_list = c(2,3,5,7,10)
n_factor_list = c(0.5, 1, 2.5, 5, 10)

for(n_sim in n_factor_list){
  for(p_sim in p_list){
    for(seed_sim in SEED_LIST){
      print(paste("p=",p_sim,", n=",n_sim, ", sim = ",seed_sim/17021890,sep = ""))
      set.seed(seed_sim)
      output_df = generate_data(n_sim, p_sim)
      save_data(output_df$X,output_df$Y,output_df$p,output_df$cl_train_label_noise, output_df$G, output_df$cl_tot, seed_sim, n_sim)
      pico_varbrand(output_df$X,output_df$Y,output_df$p,output_df$cl_train_label_noise, output_df$G, output_df$cl_tot, seed_sim)
    }
  }
}

avg_results()
