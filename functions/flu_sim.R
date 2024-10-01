#### MAIN FLU SIMULATION FUNCTIONS ####

library(countrycode)
country_itzs_names <- data.table(read_csv('data/country_itzs_names.csv', show_col_types=F))

## function to simulate one epidemic
one_flu <- function(
    country_code,
    demography, # 4-vector of pop sizes
    demography_dt,
    init_vacc, # 4-vector of proportions
    susceptibility, # in [0,1]
    transmissibility, 
    matching,
    initial_infected, # 4-vector
    period_start_date,
    epid_start_date,
    end_date, # end of period
    intended_r0 = NULL,
    vaccine_program # from vaccine_programs list
    ){
  
  country_name <- countrycode(country_code, origin='iso3c', destination='country.name')

  ## making sure dates are mondays!
  period_start_date <- last_monday(period_start_date)
  epid_start_date <- last_monday(epid_start_date)
  end_date <- last_monday(end_date)
  
  vacc_details <- vacc_type_list[[vaccine_program$vacc_type]]
  efficacy_input <- c(rep(vacc_details$VE[1 + 2*(1 - matching)], 3), vacc_details$VE[2 + 2*(1 - matching)])
  vacc_rel_inf <- vacc_details$rel_inf
  vacc_calendar_start <- vaccine_program$start
  vacc_calendar_weeks <- vaccine_program$weeks
  
  waning_rate <- 1/(365*vacc_details$imm_duration)
  
  poss_dates <- as.Date(unlist(lapply(as.Date(paste0(vacc_calendar_start, '-', (year(epid_start_date) - 1):(year(end_date) + 1)), format = '%d-%m-%Y'),
                                                last_monday)))
 
  week1 <- poss_dates[epid_start_date - poss_dates < 7*vacc_calendar_weeks &
                         epid_start_date - poss_dates >= 0]
  week2 <- poss_dates[poss_dates - epid_start_date <= (365+6) &
                        poss_dates - epid_start_date > 0]
  existing_cov_list <- c(
    demography_dt[week %in% c(week1,week2) & V==T,]$value/demography_dt[week %in% c(week1,week2) & V==T,]$total_as
  )
    
  # define vaccine calendar
  calendar_input <- dfn_vaccine_calendar(
    vacc_cov = vaccine_program$pop_coverage,
    existing_cov = existing_cov_list,
    dates_to_run = seq.Date(epid_start_date, end_date, 7),
    efficacy = efficacy_input,
    no_age_groups = length(demography),
    no_risk_groups = 1,
    vacc_calendar_start = vacc_calendar_start,
    vacc_calendar_weeks = vacc_calendar_weeks
  )
  
  contact_matrix <- fcn_contact_matrix(
    country_name, 
    country_name_altern = country_itzs_names[country %like% country_name]$country_altern,
    country_name_altern_2 = country_itzs_names[country %like% country_name]$country_altern_2,
    pop_model = demography
  )
  ## transform to per capita contacts
  contact_matrix_small <- t(t(contact_matrix)/demography) 
  
  ## scale to match a given r0 if necessary
  current_r0 <- fluEvidenceSynthesis::as_R0(
    transmission_rate = transmissibility,
    contact_matrix = contact_matrix_small,
    age_groups = demography*c((0.2*1 + 0.8*susceptibility),
                              rep(susceptibility,3))
  )
  # print(paste0('R0 = ', current_r0))
  if(is.na(intended_r0)){ intended_r0 <- NULL }
  if(!is.null(intended_r0)){
    contact_matrix_small <- contact_matrix_small*(intended_r0/current_r0)
  }
  
  dt <- incidence_VS(
    demography_input = demography,
    susceptibility,
    transmissibility,
    initial_infected,
    calendar_input,
    contacts = contact_matrix_small,
    waning_rate,
    init_vaccinated = init_vacc,
    begin_date = epid_start_date, 
    end_date, 
    age_groups_model,
    vacc_rel_inf
  )
  
  ## merging into total time period ##
  time_before <- seq.Date(period_start_date, (epid_start_date-1), by=7)
  output <- rbind(data.table(time = time_before, 
                       I1 = 0, I2 = 0, I3 = 0, I4 = 0,
                       IU1 = 0, IU2 = 0, IU3 = 0, IU4 = 0,
                       IV1 = 0, IV2 = 0, IV3 = 0, IV4 = 0), dt)
  
  output
}

## function to combine multiple epidemics

many_flu <- function(
    country,
    ageing, # = T/F
    ageing_date = NULL, # e.g. '01-04'
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program,
    model_age_groups
    ){
  
  dates_many_flu <- seq.Date(last_monday(min(epid_inputs$period_start_date)), 
                    last_monday(max(epid_inputs$end_date)), 
                    by=7)
  
  ## vaccination and ageing
  demography_dt <- fcn_weekly_demog(
    country,
    ageing,
    ageing_date,
    dates_in = dates_many_flu,
    demographic_start_year = year(min(epid_inputs$period_start_date)),
    vaccine_program,
    init_vaccinated = c(0,0,0,0),
    model_age_groups
  )
  
  output_dt <- data.table()
  for(epidemic_i in 1:nrow(epid_inputs)){ 
    
    epid_data <- epid_inputs[epidemic_i,]
    
    demog_flu <- demography_dt[week == last_monday(epid_data$epid_start_date)]
    demog_flu_pop <- demog_flu[V==T,]$total_as
    demog_flu_vacc_prop <- demog_flu[V==T,]$value/demog_flu[V==T,]$total_as
    
    flu_epid_output <- one_flu(
      country_code = country, # iso3c string
      demography = demog_flu_pop, # 4-vector of pop sizes
      demography_dt = demography_dt,
      init_vacc = demog_flu_vacc_prop, 
      susceptibility = epid_data$susceptibility, # in [0,1]
      transmissibility = epid_data$transmissibility, 
      matching = epid_data$match,
      initial_infected = unlist(epid_data$initial_infected), # 4-vector
      period_start_date = epid_data$period_start_date,
      epid_start_date = epid_data$epid_start_date,
      end_date = epid_data$end_date, # end of period
      intended_r0 = epid_data$r0_to_scale,
      vaccine_program = vaccine_program
    )
    
    if(epidemic_i == 1){
      output_dt <- flu_epid_output
    }else{
      output_dt[,2:13] <- output_dt[,2:13] + flu_epid_output[,2:13]
    }
    
  }
  
  output_dt ## return epidemics
  
}


#### ADDITIONAL FUNCTIONS ####

## function to produce uniformly increasing vaccination coverage
## see: https://blackedder.github.io/flu-evidence-synthesis/vaccination.html

change_coverage <- function(data, final_uptake, existing_uptake) {
  
  for(i in 1:length(final_uptake)){
    if(final_uptake[i] > 0){
      data[,i] <- seq(existing_uptake[i], final_uptake[i], by = (final_uptake[i]-existing_uptake[i])/(nrow(data)))[2:(nrow(data)+1)]
    }
  }  
   
  data
}

## function to define vaccine calendar for a given vaccine program

dfn_vaccine_calendar <- function(
    vacc_cov, # intended coverage, e.g. c(0.5 0.43, 0, 0)
    existing_cov,
    dates_to_run, # seq.Date between start of epidemic and end of period
    efficacy, # e.g. c(0.42, 0.42, 0.42, 0.28)
    no_age_groups, # 4
    no_risk_groups, # 1
    vacc_calendar_start, # e.g. '01-10'
    vacc_calendar_weeks, # e.g. 12
    next_cal = T
    )
  {
  
  dates_to_run <- as.Date(unlist(lapply(dates_to_run, last_monday)))
  
  coverage_matrix <- matrix(0, nrow = length(dates_to_run), ncol = no_age_groups)
  
  vc <- as_vaccination_calendar(
    dates = dates_to_run,
    efficacy = efficacy,
    coverage = coverage_matrix,
    no_age_groups = no_age_groups,
    no_risk_groups = no_risk_groups
  )
  
  ## will do current vaccination program *if* epidemic begins in the middle of a vaccination period,
  ## and then the *next* vaccination program
  
  ## *current* vaccine calendar (if applicable):
  possible_start_dates <- as.Date(unlist(lapply(as.Date(paste0(vacc_calendar_start, '-', (year(dates_to_run[1]) - 1):(year(dates_to_run[length(dates_to_run)]) + 1)), format = '%d-%m-%Y'),
                                                last_monday)))
  curr_start <- possible_start_dates[dates_to_run[1] - possible_start_dates < 7*vacc_calendar_weeks &
                                       dates_to_run[1] - possible_start_dates >= 0]
  if(length(curr_start) > 0){
    dates_to_vacc_curr <- seq.Date(curr_start, curr_start + 7*vacc_calendar_weeks - 1, 7)
    rows_to_vacc_curr <- which(dates_to_run %in% dates_to_vacc_curr) ## MONDAYS DON'T ALIGN
  }else{
    rows_to_vacc_curr <- 0
  }
  
  ## *next* vaccine calendar:
  next_start <- possible_start_dates[possible_start_dates - dates_to_run[1] <= (365+6) &
                                       possible_start_dates - dates_to_run[1] > 0]
  if(length(next_start)>1){warning('next_start includes more than one date')}
  dates_to_vacc_next <- seq.Date(next_start, next_start + 7*vacc_calendar_weeks - 1, 7)
  rows_to_vacc_next <- which(dates_to_run %in% dates_to_vacc_next) ## MONDAYS DON'T ALIGN
    
  ## changin vaccine coverage
  
  ## check number of weeks matches
  if((! (length(existing_cov)/no_age_groups) == (length(curr_start) + next_cal)) & 
     ((length(curr_start) + next_cal) > 0)){
    warning("Number of vaccination calendars isn't matching!")
  }
  existing_cov_curr <- unlist(ifelse((length(curr_start) > 0), 
                               list(existing_cov[1:4]), list(rep(0,4))))
  if(next_cal == T){
    existing_cov_next <- unlist(ifelse((length(curr_start) > 0), 
                                       list(existing_cov[5:8]), list(existing_cov[1:4])))
  }
  
  vc$calendar[rows_to_vacc_curr, 1:no_age_groups] <- change_coverage(matrix(0, nrow = length(rows_to_vacc_curr), ncol = no_age_groups), 
                                                                     vacc_cov, existing_cov_curr)
  if(next_cal == T){
    vc$calendar[rows_to_vacc_next, 1:no_age_groups] <- change_coverage(matrix(0, nrow = length(rows_to_vacc_next), ncol = no_age_groups), 
                                                                       vacc_cov, existing_cov_next)
  }
  
  vc
}


## function to match a given date to the most recent monday 
## (will leave any monday as is!)
last_monday <- function(date){
  
  move_back <- 7 - which(weekdays(seq.Date(from = date - 6,
                                           to = date,
                                           by = 1)) == 'Monday')
  
  monday <- date - move_back
                          
  if(!weekdays(monday) == 'Monday'){warning('Monday move wrong')}
  
  monday
  
}

## function to calculate vaccine doses over same period
flu_doses <- function(
    country,
    ageing, # = T/F
    ageing_date = NULL, # e.g. '01-04'
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program,
    model_age_groups){
  
  dates_many_flu <- seq.Date(last_monday(min(epid_inputs$period_start_date)), 
                             last_monday(max(epid_inputs$end_date)), 
                             by=7)
  
  doses <- fcn_annual_doses(
    country,
    ageing = T,
    ageing_date,
    dates_in = dates_many_flu,
    demographic_start_year = year(min(epid_inputs$period_start_date)),
    vaccine_program,
    init_vaccinated = c(0,0,0,0),
    model_age_groups
  )
  
  doses
}
















