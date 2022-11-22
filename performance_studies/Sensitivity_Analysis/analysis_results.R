library(tidyverse)
library(here)

res_files <- list.files(here("SA/results/results_SA/"))
res_files <- list.files(here("SA/results/2022-11-04 17:15:50 UTC/"))
# res_files <- list.files(here("SA/results/2022-11-05 17:15:48 UTC/"))
res_file_info <- str_split(res_files,pattern = "_|\\.",simplify = TRUE)
df_var_dummy_num_dummy_scenario <- (res_file_info[,3:5])

df_results <-
  map_dfr(
    1:length(res_files),
    ~ read_csv(file = paste0(
      here("SA/results/2022-11-04 17:15:50 UTC/"), res_files[.x]
      # here("SA/results/results_SA/"), res_files[.x]
    )) |>
      mutate(
        var_dummy = df_var_dummy_num_dummy_scenario[.x, 1],
        num_dummy = df_var_dummy_num_dummy_scenario[.x, 2],
        id_sim = df_var_dummy_num_dummy_scenario[.x, 3]
      )
  )


df_with_problems <- df_results |> 
  filter(is.na(ARI))

df_with_problems |> 
  print(n=Inf)
# Problems only with var dummy =1 and num_dummy=0
unique(df_with_problems$var_dummy)
unique(df_with_problems$num_dummy)
unique(df_with_problems$gamma)
unique(df_with_problems$lambda_0_DP)
unique(df_with_problems$nu_0_DP)

# Try to understand NA ----------------------------------------------------


