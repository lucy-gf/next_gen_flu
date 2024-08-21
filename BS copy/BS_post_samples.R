
## merging the posterior samples from inference into a list
## to be loaded in for BS_data_fcns.R

post_size <- 10000 
thinning_steps <- 50
burn_in <- 500000

post_samples_merged_simulation <- list()
count <- 1
for(sel_cntr in df_cntr_table$country){
  PSM_country_spec <- data.frame()
  country_code_input <- df_cntr_table$COUNTRY_CODE[df_cntr_table$country == sel_cntr]
  if(country_code_input == 'UK'){country_code_input <- 'GBR'}
  for(strain_index in c('INF_A','INF_B')){
    for(i in 1:df_cntr_table[,paste0(strain_index,'_NONSENT')][df_cntr_table$country==sel_cntr]){
      epid <- i
      if(file.exists(paste0("command_line_runs/", strain_index, '_bt/', sel_cntr,"_epid",
                             epid, "_",
                             post_size, "_", thinning_steps, "_", burn_in, ".rds"))){
      results_RDS <- readRDS(paste0("command_line_runs/", strain_index, '_bt/', sel_cntr,"_epid",
                                    epid, "_",
                                    post_size, "_", thinning_steps, "_", burn_in, ".rds"))
      
      post_samples <- data.frame(results_RDS)
      # colnames(post_samples) <- rep(c("reporting","trans","sus","infected","blank","blank"),length(epidemics_to_fit))
      colnames(post_samples) <- paste0(c("reporting","trans","sus","infected"))
      # make into real values
      post_samples$reporting <- exp(post_samples$reporting)
      post_samples$trans <- post_samples$trans/100
      post_samples$infected <- 10^(post_samples$infected)
      post_samples$timestep <- 1:nrow(post_samples)
      post_samples$incl <- (post_samples$timestep %% thinning_steps == 0 &
                              post_samples$timestep > burn_in)
      post_samples <- post_samples %>% filter(incl == T) %>% mutate(timestep=1:post_size)
      post_samples_m <- post_samples %>% pivot_longer(!c(timestep, incl)) %>% 
        rename(variable = name) %>% mutate(epidemic_id = epid,
                                           strain = strain_index)
      
      PSM_country_spec <- rbind(PSM_country_spec,post_samples_m)
      }
    }
  }
  post_samples_merged_simulation[[count]] <- PSM_country_spec
  count <- count + 1
}

names(post_samples_merged_simulation) <- c(df_cntr_table$COUNTRY_CODE[1:6], 'GBR')

saveRDS(post_samples_merged_simulation, file = "data_for_BS/post_samples_merged_simulation.rds")







