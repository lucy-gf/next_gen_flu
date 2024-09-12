#### MAIN FLU SIMULATION FUNCTIONS ####

library(countrycode)
country_itzs_names <- data.table(read_csv('data/country_itzs_names.csv', show_col_types=F))

## function to simulate one epidemic
one_flu <- function(
    country_name,
    demography, # 4-vector of pop sizes
    init_vacc, # 4-vector of proportions
    susceptibility, # in [0,1]
    transmissibility, 
    initial_infected, # 4-vector
    period_start_date,
    epid_start_date,
    end_date, # end of period
    vacc_calendar_start, # annual
    vacc_calendar_weeks, # annual
    vaccine_program # from vaccine_programs list
    ){
  
  # demography <- c(1000,2000,3000,2000)
  # susceptibility <- 0.6
  # transmissibility <- 0.1
  # initial_infected <- c(20,20,20,20)
  # period_start_date <- as.Date('01-01-2025',format='%d-%m-%Y')
  # epid_start_date <- as.Date('01-11-2025',format='%d-%m-%Y')
  # end_date <- as.Date('01-01-2030',format='%d-%m-%Y')
  # init_vacc <- c(0.3,0.3,0,0)
  # age_groups_model = c(0,5,20,65)
  # vacc_rel_inf <- 1
  
  ## making sure dates are mondays!
  period_start_date <- last_monday(period_start_date)
  epid_start_date <- last_monday(epid_start_date)
  end_date <- last_monday(end_date)
  
  efficacy_input <- vacc_type_list[[vaccine_program$vacc_type]]$VE
  
  waning_rate <- 1/(365*vacc_type_list[[vaccine_program$vacc_type]]$imm_duration)
  
  # define vaccine calendar
  calendar_input <- dfn_vaccine_calendar(
    prop_vacc = init_vacc,
    vacc_cov = vaccine_program$pop_coverage,
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
  # transform to per capita contacts
  contact_matrix_small <- t(t(contact_matrix)/demography) 
  
  dt <- incidence_VS(
    demography,
    susceptibility,
    transmissibility,
    initial_infected,
    calendar_input,
    contact_matrix_small,
    waning_rate,
    init_vacc,
    begin_date, 
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
    ageing_date, # e.g. '01-04'
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program,
    model_age_groups
    ){
  
  ## dates, country, ageing made consistent across all epidemics
  ## need some monday checking on dates too
  
  ## make epid_inputs into a list
  # dates <- min(epid_input$date), max(min(epid_input$date))
  dates <- seq.Date(period_start_date, end_date, by=7)
  
  ## vaccination and ageing
  demography_dt <- fcn_weekly_demog(
    country,
    ageing,
    ageing_date,
    dates,
    demographic_start_year,
    vaccine_program,
    init_vaccinated,
    model_age_groups
  )
  
  output_df <- data.table()
  for(epidemic_i in 1:nrow(df_inputs)){ ## change from loop to speed up?
    
    one_flu(
      country_name = country, # string, need to choose whether to be country name or iso3c
      demography = c(1000,2000,3000,2000), # 4-vector of pop sizes
      init_vacc = c(0,0,0,0), 
      susceptibility = 0.6, # in [0,1]
      transmissibility = 0.1, 
      initial_infected = c(20,20,20,20), # 4-vector
      period_start_date = as.Date('04-01-2025',format='%d-%m-%Y'),
      epid_start_date = as.Date('01-11-2025',format='%d-%m-%Y'),
      end_date = as.Date('01-01-2027',format='%d-%m-%Y'), # end of period
      vacc_calendar_start = '01-10', 
      vacc_calendar_weeks = 12,  
      vaccine_program = vaccine_programs[[2]]
    )
    
    output_df <- output_df + one_flu(df_inputs[epidemic_i,])
  }
  
  output_df # return epidemics
}


#### ADDITIONAL FUNCTIONS ####

## function to produce uniformly increasing vaccination coverage
## see: https://blackedder.github.io/flu-evidence-synthesis/vaccination.html

change_coverage <- function(data, final_uptake) {
  sums <- data[nrow(data),]
  # If final uptake is zero in a group then we need to make some kind of assumption on uptake rate over time
  if (any(sums == 0)) {
    # warning("No prior information on uptake rate. Using constant uptake rate")
    col <- which(sums == 0)
    data[,col] <- seq(0, (nrow(data)-1))
    sums <- data[nrow(data),]    
  }
  for(i in 1:nrow(data)) {
    data[i,] <- data[i,]*final_uptake/sums
  }
  data
}

## function to define vaccine calendar for a given vaccine program

dfn_vaccine_calendar <- function(
    prop_vacc, # initial covergae, e.g. c(0.3, 0.3, 0, 0)
    vacc_cov, # intended coverage, e.g. c(0.5 0.43, 0, 0)
    dates_to_run, # seq.Date between start of epidemic and end of period
    efficacy, # e.g. c(0.42, 0.42, 0.42, 0.28)
    no_age_groups, # 4
    no_risk_groups, # 1
    vacc_calendar_start, # e.g. '01-10'
    vacc_calendar_weeks # e.g. 12
    )
  {
  
  if(sum(prop_vacc>1) + sum(vacc_cov>1)>0){stop('Init_vacc and prop_vacc must be proportions, not absolute numbers')}
  
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
  next_start <- possible_start_dates[possible_start_dates - dates_to_run[1] <= 365 &
                                       possible_start_dates - dates_to_run[1] > 0]
  dates_to_vacc_next <- seq.Date(next_start, next_start + 7*vacc_calendar_weeks - 1, 7)
  rows_to_vacc_next <- which(dates_to_run %in% dates_to_vacc_next) ## MONDAYS DON'T ALIGN
    
  ## changin vaccine coverage
  vc$calendar[rows_to_vacc_curr, 1:no_age_groups] <- change_coverage(matrix(0, nrow = length(rows_to_vacc_curr), ncol = no_age_groups), vacc_cov)
  vc$calendar[rows_to_vacc_next, 1:no_age_groups] <- change_coverage(matrix(0, nrow = length(rows_to_vacc_next), ncol = no_age_groups), vacc_cov)
  
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




# ## function to calculate vaccine doses over same period
# flu_doses <- function(){
#   
# }
















