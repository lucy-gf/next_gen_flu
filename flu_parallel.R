#### RUN THE FLU MODEL ####

#setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)
library(fluEvidenceSynthesis)
library(parallel)

# key outputs: vaccine_programs, vacc_type_list
source('vacc_types.R') 
# will calculate weekly age- and vaccine-specific population, also loads transmission model
source('functions/demography.R') 
# runs the flu model
source('functions/flu_sim.R') 

#### FRAMEWORK ####

# **vacc_types.R** contains functions to define vaccine programs 
# for no vaccinations and 5 NGIVs (can be edited)
# (inputs: coverage level, targeted ages)
# (outputs: VE, mean immunity length, coverage across model age groups)

# **functions/demography.R** calculates weekly age- and vaccination-status specific population over
# the relevant time period (can remove ageing if needed), for a given country and vaccine program,
# and contains the function to calculate contact matrices from a given country and age-specific population

# **functions/transmission_model.R** contains the ODE model builder, 
# epidemic simulation function, and a function to calculate vaccination status-specific demography

# **functions/flu_sim.R** contains:
# one_flu(), which runs an epidemic,
# many_flu(), which takes a data-table of epidemic data and combines many epidemics, 
# dfn_vaccine_calendar(), which converts the vaccine program and epidemic dates into a vaccine calendar
# flu_doses(), which calculates how many vaccines were given (before wastage) in the same epidemics as many_flu()

# **flu_parallel.R** sets vaccine programs, 
# then runs many_flu() for some epidemic data, parallelised across each vaccine type


#### VACCINE PROGRAMS ####
cov_main <- 0.5 # coverage level in targeted age groups
age_targeting <- 0:17 # ages being targeted
vaccine_programs <- c(
  fcn_vacc_prog(NA, 0, T),
  fcn_vacc_prog(age_targeting, cov_main)
)

## age groups - fixed ##
model_age_groups <- c(0,5,20,65)
age_group_names <- paste0(model_age_groups,"-", c(model_age_groups[2:length(model_age_groups)],99))

#### FUNCTION TO RUN ####
## only input is vaccine type, to parallelise over vt ##
flu_parallel <- function(vaccine_type){
  
  epid_starts <- seq.Date(from = as.Date(paste0("01-11-2025"), '%d-%m-%Y'),
                          to = as.Date(paste0("01-11-2029"), '%d-%m-%Y'),
                          by = 365)
  
  epid_dt <- data.table(
    susceptibility = c(0.6,0.59,0.52,0.6,0.61),
    transmissibility = c(0.06,0.07,0.065,0.075,0.08),
    match = c(T,F,F,T,F),
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

#### RUN OUTPUTS #### 
## (parallelised across all vaccine types) ##

infs_rds_list <- mclapply(1:6, flu_parallel, mc.cores=5)

#### SAVE OUTPUTS ####

saveRDS(infs_rds_list, file = paste0("outputs/vacc_", cluster_code, "_", scenario_name,".rds"))
























