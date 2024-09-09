#### RUN THE FLU MODEL ####

#setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)
library(fluEvidenceSynthesis)
library(parallel)

source('vacc_types.R') # key outputs: vaccine_programs, vacc_type_list
source('transmission_model.R') # loads transmission model
source('demography.R') # will calculate weekly age- and vaccine-specific population
source('flu_sim.R') # runs the flu model

flu_parallel <- function(vaccine_type){
  many_flu(country,
           ageing, 
           risk_ratios = c(0,0,0,0), 
           epid_inputs,  
           vaccine_program = vaccine_programs[[vaccine_type]]
  )
}

cases_rds_list <- mclapply(1:6, flu_parallel, mc.cores=5)

saveRDS(cases_rds_list, file = paste0("data/vacc_", cluster_code, "_", scenario_name,".rds"))
























