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
    vacc_cov, # 4-vector
    vaccine_program # from vaccine_programs list
    ){
  
  efficacy_input <- vacc_type_list[[vaccine_program$vacc_type]]$VE
  
  waning_rate <- 1/(365*vacc_type_list[[vaccine_program$vacc_type]]$imm_duration)
  
  demography <- c(1000,2000,3000,2000)
  susceptibility <- 0.6
  transmissibility <- 0.1
  initial_infected <- c(20,20,20,20)
  period_start_date <- as.Date('01-01-2025',format='%d-%m-%Y')
  epid_start_date <- as.Date('01-11-2025',format='%d-%m-%Y')
  end_date <- as.Date('01-01-2030',format='%d-%m-%Y')
  init_vacc <- c(0.3,0.3,0,0)
  age_groups_model = c(0,5,20,65)
  vacc_rel_inf <- 1
  
  # define vaccine calendar
  calendar_input <- dfn_vaccine_calendar(
    prop_vacc = init_vacc,
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
    risk_ratios = c(0,0,0,0), # 4-vector
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program 
    ){
  
  ## dates, country, ageing made consistent across all epidemics
  ## need some monday checking on dates too
  
  ## make epid_inputs into a list
  
  one_flu(
    country_name = country, # string
    demography = c(1000,2000,3000,2000), # 4-vector of pop sizes
    init_vacc = c(0,0,0,0), 
    susceptibility = 0.6, # in [0,1]
    transmissibility = 0.1, 
    initial_infected = c(20,20,20,20), # 4-vector
    period_start_date = as.Date('04-01-2025',format='%d-%m-%Y'),
    epid_start_date = as.Date('01-11-2025',format='%d-%m-%Y'),
    end_date = as.Date('01-01-2030',format='%d-%m-%Y'), # end of period
    vacc_calendar_start = '01-10', 
    vacc_calendar_weeks = 12, 
    vacc_cov = c(0.3,0.3,0,0), 
    vaccine_program = vaccine_programs[[1]]
  )
  
  output_df <- data.table()
  for(epidemic_i in 1:nrow(df_inputs)){ ## change from loop to speed up!
    output_df <- output_df + one_flu(df_inputs[epidemic_i,])
  }
  
  output_df # return epidemics
}


dfn_vaccine_calendar(
  dates_to_run = seq.Date(epid_start_date, end_date, 7),
  efficacy = efficacy_input,
  no_age_groups = length(demography),
  no_risk_groups = 1,
  vacc_calendar_start = '01-10',
  vacc_calendar_weeks = 12
)

change_coverage <- function(data, final_uptake) {
  sums <- data[nrow(data),]
  # If final uptake is zero in a group then we need to make some kind of assumption on uptake rate over time
  if (any(sums == 0)) {
    warning("No prior information on uptake rate. Using constant uptake rate")
    col <- which(sums == 0)
    data[,col] <- seq(0, (nrow(data)-1))
    sums <- data[nrow(data),]    
  }
  for(i in 1:nrow(data)) {
    data[i,] <- data[i,]*final_uptake/sums
  }
  data
}

# function to define vaccine calendar from vaccine program
dfn_vaccine_calendar <- function(
    prop_vacc, # e.g. c(0.3, 0.3, 0, 0)
    dates_to_run, # seq.Date between start of epidemic and end of period
    efficacy, # e.g. c(0.42, 0.42, 0.42, 0.28)
    no_age_groups, # 4
    no_risk_groups, # 1
    vacc_calendar_start, # e.g. '01-10'
    vacc_calendar_weeks # e.g. 12
    )
  {
  
  if(sum(prop_vacc>1)>0){stop('Init_vacc must be proportions, not absolute numbers')}
  
  coverage_matrix <- matrix(0, nrow = length(dates_to_run), ncol = no_age_groups)
  coverage_matrix[1, ] <- prop_vacc
  
  ## vaccination calendar will encompass the first annual vaccination program which overlaps with or comes after the epidemic beginning
  
  year_to_vacc <- ifelse(dates_to_run[1] < (as.Date(paste0(vacc_calendar_start, '-', year(dates_to_run[1])), format = '%d-%m-%Y') + 7*vacc_calendar_weeks - 1),
                         year(dates_to_run[1]),
                         year(dates_to_run[1]) + 1)
    
  vaccination_program_dates_weekly <- seq.Date(as.Date(paste0(vacc_calendar_start, '-', year_to_vacc), format = '%d-%m-%Y'), 
                                               as.Date(paste0(vacc_calendar_start, '-', year_to_vacc), format = '%d-%m-%Y') + 7*vacc_calendar_weeks - 1, 
                                               7)
  vaccination_program_dates_daily <- seq.Date(as.Date(paste0(vacc_calendar_start, '-', year_to_vacc), format = '%d-%m-%Y'),
                                              as.Date(paste0(vacc_calendar_start, '-', year_to_vacc), format = '%d-%m-%Y') + 7*vacc_calendar_weeks - 1, 
                                              1)
  
  as_vaccination_calendar(
    dates = dates_to_run,
    efficacy = efficacy,
    coverage = coverage_matrix,
    no_age_groups = no_age_groups,
    no_risk_groups = no_risk_groups
  )
  
  
  
  coverage_matrix <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  new_coverage_matrix <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  # only implementing vaccinations which occur in that year (defined as between ageing dates)
  dates_numeric <- as.numeric(dates_to_run)
  year_of_first_vacc <- case_when(
    hemisphere == 'NH' & month(start_date) < as.numeric(substr(vaccine_program$SH_vacc_date, 4, 5)) ~
      year(start_date) - 1,
    hemisphere == 'NH' & month(start_date) >= as.numeric(substr(vaccine_program$SH_vacc_date, 4, 5)) ~
      year(start_date),
    hemisphere == 'SH' & month(start_date) < as.numeric(substr(vaccine_program$NH_vacc_date, 4, 5)) ~
      year(start_date),
    hemisphere == 'SH' & month(start_date) >= as.numeric(substr(vaccine_program$NH_vacc_date, 4, 5)) ~
      year(start_date) + 1
  )
  first_vacc <- c(as.Date(paste0(vacc_date, '-', year_of_first_vacc), format = '%d-%m-%Y'))
  dates_to_vacc <- first_vacc + 7*(0:(vaccine_program$weeks_vaccinating-1))
  nums_to_vacc <- as.numeric(dates_to_vacc)
  rows_to_vacc <- c()
  for(i in 1:length(nums_to_vacc)){
    rows_to_vacc <- c(rows_to_vacc, which((dates_numeric - nums_to_vacc[i]) >= 0 &
                                            (dates_numeric - nums_to_vacc[i]) < 7))
  }
  changed_cov <- change_coverage(coverage_matrix[1:(vaccine_program$weeks_vaccinating + 1),], vaccine_program$pop_coverage/vacc_details$coverage_pattern)
  if(year_of_first_vacc == start_year & vaccine_program$first_year_all == T){
    changed_cov <- change_coverage(coverage_matrix[1:(vaccine_program$weeks_vaccinating + 1),], vaccine_program$pop_coverage)
  }
  
  if(length(rows_to_vacc) == 0){
    coverage_matrix_minimal <- coverage_matrix
    dates <- dates_to_run
  }else{
    if(length(rows_to_vacc) == 1){
      new_coverage_matrix[c(rows_to_vacc, rows_to_vacc + 1),] <- changed_cov[1:(length(rows_to_vacc) + 1),]
      coverage_matrix_minimal <- new_coverage_matrix[1:(rows_to_vacc[length(rows_to_vacc)] + 1),]
      coverage_matrix_minimal <- rbind(coverage_matrix_minimal, coverage_matrix_minimal[nrow(coverage_matrix_minimal),])
      dates <- dates_to_run[1:(rows_to_vacc[length(rows_to_vacc)] + 1)]
    }else{
      new_coverage_matrix[c(rows_to_vacc, rows_to_vacc[length(rows_to_vacc)] + 1),] <- changed_cov[1:(length(rows_to_vacc) + 1),]
      coverage_matrix_minimal <- new_coverage_matrix[1:(rows_to_vacc[length(rows_to_vacc)] + 1),]
      coverage_matrix_minimal <- rbind(coverage_matrix_minimal, changed_cov[(length(rows_to_vacc) + 1),])
      dates <- dates_to_run[1:(rows_to_vacc[length(rows_to_vacc)] + 2)]
    }}
  vaccine_calendar <- as_vaccination_calendar(
    efficacy = c(rep(vacc_details$VE[1 + 2*(1 - match)], 3), vacc_details$VE[2 + 2*(1 - match)]),
    dates = dates,
    coverage = coverage_matrix_minimal,
    no_age_groups = length(age_group_names),
    no_risk_groups=1
  )
  
  if(nrow(unvaxed_pop)>0){
    unvaxed_pop[c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4')] <- as.numeric(unvaxed_pop[c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4')])
    prop_vacc <- c(unname(unvaxed_pop['V1']/(unvaxed_pop['U1'] + unvaxed_pop['V1'])), 
                   unname(unvaxed_pop['V2']/(unvaxed_pop['U2'] + unvaxed_pop['V2'])),
                   unname(unvaxed_pop['V3']/(unvaxed_pop['U3'] + unvaxed_pop['V3'])),
                   unname(unvaxed_pop['V4']/(unvaxed_pop['U4'] + unvaxed_pop['V4'])))
    for(i in 1:nrow(vaccine_calendar$calendar)){
      for(j in 1:4){
        if(!vaccine_calendar$calendar[i,j] == 0){
          vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
            (12 + i - rows_to_vacc[length(rows_to_vacc)])*(vaccine_program$pop_coverage[j] - unlist(prop_vacc)[j])/vaccine_program$weeks_vaccinating
        }
      }
    } 
  }
}







# ## function to calculate vaccine doses over same period
# flu_doses <- function(){
#   
# }
















