#### MAIN FLU SIMULATION FUNCTIONS ####

model_age_groups <- c(0,5,20,65)

## function to simulate one epidemic
one_flu <- function(
    demography, # c(4-vector of unvaccinated pop sizes, 4-vector of vaccinated pop sizes)
    contact_matrix, # 4x4 contact matrix
    susceptibility, # in [0,1]
    transmissibility, 
    initial_infected, # 4-vector
    period_start_date = as.Date('01-01-2025',format='%d-%m-%Y'),
    epid_start_date,
    end_date, # end of period
    init_vacc = rep(0,4), # 4-vector 
    vacc_cov, # 4-vector
    vacc_dates, # annual
    vaccine_program # from list to be made
    ){
  
  # define vaccine calendar
  calendar_input <- 
    
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
  
  output_df <- data.table()
  for(epidemic_i in 1:nrow(df_inputs)){ ## change from loop to speed up!
    output_df <- output_df + one_flu(df_inputs[epidemic_i,])
  }
  
  output_df # return epidemics
}

## function to calculate vaccine doses over same period
flu_doses <- function(){
  
}




















