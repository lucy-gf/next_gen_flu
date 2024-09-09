#### MAIN FLU SIMULATION FUNCTIONS ####

## function to simulate one epidemic
one_flu <- function(
    demography = c(100,200,300,200), # 4-vector of unvaccinated pop sizes
    contact_matrix, # 4x4 contact matrix
    susceptibility, # in [0,1]
    transmissibility, 
    initial_infected, # 4-vector
    period_start_date = as.Date('01-01-2025',format='%d-%m-%Y'),
    epid_start_date = as.Date('01-11-2025',format='%d-%m-%Y'),
    end_date = as.Date('01-01-2030',format='%d-%m-%Y'), # end of period
    init_vacc = rep(0,4), # 4-vector 
    vacc_calendar_start, # annual
    vacc_calendar_weeks, # annual
    vacc_cov, # 4-vector
    vaccine_program # from list to be made
    ){
  
  # define vaccine calendar
  calendar_input <- as_vaccination_calendar(
    efficacy = c(0.4, 0.4, 0.4, 0.4), #c(rep(vacc_details$VE[1 + 2*(1 - match)], 3), vacc_details$VE[2 + 2*(1 - match)]),
    dates = seq.Date(epid_start_date, end_date, 7),
    coverage = matrix(0, nrow = length(seq.Date(epid_start_date, end_date, 7)), ncol = 4),
    no_age_groups = length(demography),
    no_risk_groups = 1
  )
    
  waning_rate <- vacc_type_list[[vaccine_program$vacc_type]]$imm_duration
  
  incidence_VS(
    demography,
    susceptibility,
    transmissibility,
    initial_infected,
    calendar_input,
    contact_matrix,
    waning_rate,
    vaccination_ratio_input,
    begin_date, 
    end_date, 
    age_groups_model,
    kappa_input = 1
  )
  
}

## function to combine multiple epidemics
many_flu <- function(
    country,
    ageing, # = T/F
    risk_ratios = c(0,0,0,0), # 4-vector
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program 
    ){
  
  ## dates, country, ageing made consistent across all epidemics
  ## need some monday checking on dates too
  
  ## make epid_inputs into a list
  
  output_df <- data.table()
  for(epidemic_i in 1:nrow(df_inputs)){ ## change from loop to speed up!
    output_df <- output_df + one_flu(df_inputs[epidemic_i,])
  }
  
  output_df # return epidemics
}

## function to calculate vaccine doses over same period
flu_doses <- function(){
  
}




















