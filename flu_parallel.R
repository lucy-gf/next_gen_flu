#### RUN THE FLU MODEL ####

# setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)
library(fluEvidenceSynthesis)
library(parallel)

# key outputs: vaccine_programs, vacc_type_list
source('vacc_types.R') 
# will calculate weekly age- and vaccine-specific population, also loads transmission model
source('functions/demography.R') 
# runs the flu model
source('functions/flu_sim.R') 

#### VACCINE PROGRAM DETAILS ####
cov_main <- 0.5 # coverage level in targeted age groups
age_targeting <- 0:17 # ages being targeted
vacc_calendar_start <- '01-10'
vacc_calendar_weeks <- 12

vaccine_programs <- c(
  fcn_vacc_prog(NA, 0,
                vacc_calendar_start, 
                vacc_calendar_weeks,
                T),
  fcn_vacc_prog(age_targeting, cov_main,
                vacc_calendar_start, 
                vacc_calendar_weeks)
)

## MODEL AGE GROUPS - fixed ##
model_age_groups <- c(0,5,20,65)
age_group_names <- paste0(model_age_groups,"-", c(model_age_groups[2:length(model_age_groups)],99))
# could generalise across more/different age groups, not done though.

#### EPIDEMIC DATA ####
epid_starts <- c(as.Date(paste0('01-11-', 2025:2029), format='%d-%m-%Y'))
epid_dt <- data.table(
  susceptibility = c(0.6,0.59,0.7,0.65,0.61),
  transmissibility = c(0.09,0.08,0.075,0.075,0.08),
  match = c(T,F,F,T,F),
  initial_infected = c(list(c(100,100,100,100)),list(c(100,100,100,100)),list(c(100,100,100,100)),
                       list(c(100,100,100,100)),list(c(100,100,100,100))),
  period_start_date = as.Date('01-01-2025',format='%d-%m-%Y'),
  epid_start_date = epid_starts,
  end_date = as.Date('01-01-2031',format='%d-%m-%Y'),
  r0_to_scale = NA
)

ageing_date <- '01-04'
iso3c <- 'GBR'

#### FUNCTION TO RUN ####
## only input is vaccine type, to parallelise over vt ##
flu_parallel <- function(vaccine_type){

  mf_output <- many_flu(country = iso3c,
                        ageing = T, 
                        ageing_date,
                        epid_inputs = epid_dt,  
                        vaccine_program = vaccine_programs[[vaccine_type]],
                        model_age_groups
                        )
  
  mf_output[, vacc_type := names(vaccine_programs)[vaccine_type]]

}

#### RUN OUTPUTS #### 
## (parallelised across all vaccine types) ##

infs_rds_list <- mclapply(1:length(vaccine_programs), flu_parallel, mc.cores=5)

#### SAVE OUTPUTS ####

saveRDS(infs_rds_list, file = paste0("outputs/vacc_", iso3c,".rds"))








