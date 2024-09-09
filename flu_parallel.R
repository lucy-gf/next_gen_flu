#### RUN THE FLU MODEL ####

#setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)
library(fluEvidenceSynthesis)
library(parallel)

source('vacc_types.R') # key outputs: vaccine_programs, vacc_type_list
source('transmission_model.R') # loads transmission model
source('demography.R') # will calculate weekly age- and vaccine-specific population
source('flu_sim.R') # runs the flu model

## vaccine programs ##
cov_main <- 0.5 # coverage level
age_targeting <- 0:17
vaccine_programs <- c(
  fcn_vacc_prog(NA, 0, T),
  fcn_vacc_prog(age_targeting, cov_main)
)

## age groups - fixed for now ##
model_age_groups <- c(0,5,20,65)

## function to run ##
flu_parallel <- function(vaccine_type){
  
  epid_starts <- seq.Date(from = as.Date(paste0("01-11-2025"), '%d-%m-%Y'),
                          to = as.Date(paste0("01-11-2029"), '%d-%m-%Y'),
                          by = 365)
  
  epid_dt <- data.table(
    susceptibility = c(0.6,0.59,0.52,0.6,0.61),
    transmissibility = c(0.06,0.07,0.065,0.075,0.08),
    initial_infected = NULL,
    period_start_date = as.Date('01-01-2025',format='%d-%m-%Y'),
    epid_start_date = epid_starts,
    end_date = as.Date('01-01-2030',format='%d-%m-%Y'),
    vacc_calendar_start = '01-10',
    vacc_calendar_weeks = 12 # annual
  )
  
  many_flu(country = 'GBR',
           ageing = F, 
           risk_ratios = c(0,0,0,0), 
           epid_inputs = epid_dt,  
           vaccine_program = vaccine_programs[[vaccine_type]]
  )
}

## RUN OUTPUTS (parallelised across all vaccine types) ## 

infs_rds_list <- mclapply(1:6, flu_parallel, mc.cores=5)

## SAVE OUTPUTS ## 

saveRDS(infs_rds_list, file = paste0("outputs/vacc_", cluster_code, "_", scenario_name,".rds"))
























