#### MAIN FLU SIMULATION FUNCTIONS ####

library(countrycode)
country_itzs_names <- data.table(read_csv(here::here('data','country_itzs_names.csv'), show_col_types=F))
country_itzs_names[codes=='TUR',country_altern_2 := 'Turkiye']

## function to simulate one epidemic
one_flu <- function(
    country_code,
    demography_dt,
    susceptibility, # in [0,1]
    transmissibility, 
    matching,
    initial_infected, # 4-vector
    period_start_date,
    epid_start_date,
    end_date, # end of period
    intended_r0 = NULL,
    vaccine_program,
    vaccine_used,
    vaccine_var,
    doses_dt,
    vacc_cov_vec
    ){
  
  country_name <- countrycode(country_code, origin='iso3c', destination='country.name')

  ## making sure dates are mondays!
  period_start_date <- last_monday(period_start_date)
  epid_start_date <- last_monday(epid_start_date)
  end_date <- last_monday(end_date)
  
  ## vaccine details, from vacc_type_list
  vacc_details <- vacc_type_list[[vaccine_program$vacc]]
  efficacy_input <- c(rep(vacc_details$VE[1 + 2*(1 - matching)], 3), vacc_details$VE[2 + 2*(1 - matching)])
  vacc_rel_inf <- vacc_details$rel_inf
  waning_rate <- vaccine_program$waning
  
  ## making vaccine calendar for epidemic simulation
  poss_dates <- as.Date(unlist(lapply(as.Date(paste0(vacc_calendar_start, '-', (year(epid_start_date) - 1):(year(end_date) + 1)), format = '%d-%m-%Y'),
                                                last_monday)))
 
  week1 <- poss_dates[epid_start_date - poss_dates < 7*vacc_calendar_weeks &
                         epid_start_date - poss_dates >= 0]
  week2 <- poss_dates[poss_dates - epid_start_date <= (365+6) &
                        poss_dates - epid_start_date > 0]
  existing_cov_list <- c(
    demography_dt[week %in% c(week1,week2) & V==T,]$value/demography_dt[week %in% c(week1,week2) & V==T,]$total_as
  )
  
  ## which year to vaccinate population?
  key_vacc_date <- case_when(month(epid_start_date) < as.numeric(substr(ageing_date,4,5)) & hemisphere_input == 'NH' ~ year(epid_start_date) - 1,
                             month(epid_start_date) >= as.numeric(substr(ageing_date,4,5)) & hemisphere_input == 'NH' ~ year(epid_start_date),
                             month(epid_start_date) < as.numeric(substr(ageing_date,4,5)) & hemisphere_input == 'SH' ~ year(epid_start_date),
                             month(epid_start_date) >= as.numeric(substr(ageing_date,4,5)) & hemisphere_input == 'SH' ~ year(epid_start_date) + 1)
  # make sure key_vacc_date is within range
  # (period often technically starts at the very end of previous year due to last_monday, so add 7 days)
  if(key_vacc_date < year(period_start_date + 7)){key_vacc_date <- year(period_start_date + 7)}
  if(key_vacc_date > year(end_date)){key_vacc_date <- year(end_date)}
  
  key_vacc_date_full <- as.Date(paste0(vacc_calendar_start,'-',key_vacc_date), format='%d-%m-%Y')
  demog_flu <- demography_dt[week == last_monday(key_vacc_date_full)]
  # vaccinate next year if epidemic starts in the middle of a vaccination period
  demog_flu_next_yr <- if(epid_start_date %in% seq.Date(last_monday(key_vacc_date_full), 
                                                        last_monday(key_vacc_date_full) + 7*vacc_calendar_weeks + 1,
                                                        by = 1)){
    demography_dt[week %in% last_monday(key_vacc_date_full + 365)]
  }else{data.table()}
  
  # prop unvaccinated
  prop_unv <- demog_flu[U==T]$value/demog_flu[U==T]$total_as
  tot_pop <- demog_flu[U==T]$total_as
  
  if(nrow(demog_flu_next_yr)>0){
    
    prop_unv_n <- demog_flu_next_yr[U==T]$value/demog_flu_next_yr[U==T]$total_as
    tot_pop_n <- demog_flu_next_yr[U==T]$total_as
    
    # define vaccine calendar
    calendar_input <- if(vaccine_var == 'doses'){
      dfn_vaccine_calendar_doses(
        vacc_cov = prop_unv*unname(unlist(doses_dt[year==key_vacc_date & vacc_scenario==vaccine_used[length(vaccine_used)]]$doses))/tot_pop,
        dates_to_run = seq.Date(epid_start_date, end_date, 7),
        key_vacc_date_full,
        efficacy = efficacy_input,
        no_age_groups = length(tot_pop),
        no_risk_groups = 1,
        vacc_calendar_start = vacc_calendar_start,
        vacc_calendar_weeks = vacc_calendar_weeks,
        next_cal = T,
        vacc_cov_next = prop_unv_n*unname(unlist(doses_dt[year==key_vacc_date+1 & vacc_scenario==vaccine_used[length(vaccine_used)]]$doses))/tot_pop_n
      )
    }else{
      dfn_vaccine_calendar_cov(
        vacc_cov = vacc_cov_vec,
        existing_cov = c((1 - prop_unv), (1 - prop_unv_n)), # vector of existing coverage in year one, then year two (if using)
        dates_to_run = seq.Date(epid_start_date, end_date, 7),
        efficacy = efficacy_input,
        no_age_groups = length(tot_pop),
        no_risk_groups = 1,
        vacc_calendar_start = vacc_calendar_start,
        vacc_calendar_weeks = vacc_calendar_weeks,
        next_cal = T
      )
    }
    
  }else{
    calendar_input <- if(vaccine_var == 'doses'){
      dfn_vaccine_calendar_doses(
        vacc_cov = prop_unv*unname(unlist(doses_dt[year==key_vacc_date & vacc_scenario==vaccine_used[length(vaccine_used)]]$doses))/tot_pop,
        dates_to_run = seq.Date(epid_start_date, end_date, 7),
        key_vacc_date_full,
        efficacy = efficacy_input,
        no_age_groups = length(tot_pop),
        no_risk_groups = 1,
        vacc_calendar_start = vacc_calendar_start,
        vacc_calendar_weeks = vacc_calendar_weeks,
        next_cal = T,
        vacc_cov_next = prop_unv*unname(unlist(doses_dt[year==key_vacc_date & vacc_scenario==vaccine_used[length(vaccine_used)]]$doses))/tot_pop
      )}else{
        dfn_vaccine_calendar_cov(
          vacc_cov = vacc_cov_vec,
          existing_cov = 1 - prop_unv,
          dates_to_run = seq.Date(epid_start_date, end_date, 7),
          efficacy = efficacy_input,
          no_age_groups = length(tot_pop),
          no_risk_groups = 1,
          vacc_calendar_start = vacc_calendar_start,
          vacc_calendar_weeks = vacc_calendar_weeks,
          next_cal = T
        )
      }
  }
  
  ## making contact matrix wrt aged population
  contact_matrix <- fcn_contact_matrix(
    country_name, 
    country_name_altern = country_itzs_names[codes %like% country_code]$country_altern,
    country_name_altern_2 = country_itzs_names[codes %like% country_code]$country_altern_2,
    pop_model = tot_pop
  )
  ## transform to per capita contacts
  contact_matrix_small <- t(t(contact_matrix)/tot_pop) 
  
  ## scale to match a given r0 if necessary
  current_r0 <- fluEvidenceSynthesis::as_R0(
    transmission_rate = transmissibility,
    contact_matrix = contact_matrix_small,
    age_groups = tot_pop
  )
  # print(paste0('R0 = ', current_r0))
  if(is.na(intended_r0)){ intended_r0 <- NULL }
  if(!is.null(intended_r0)){
    contact_matrix_small <- contact_matrix_small*(intended_r0/current_r0)
  }
  
  ## run transmission model
  dt <- incidence_VS(
    demography_input = tot_pop,
    susceptibility,
    transmissibility,
    initial_infected,
    calendar_input,
    contacts = contact_matrix_small,
    waning_rate,
    init_vaccinated = 1 - prop_unv,
    begin_date = epid_start_date, 
    end_date = min(end_date, epid_start_date + 548), # allowed to run for 1.5 years 
    model_age_groups,
    vacc_rel_inf,
    vaccine_var_in = vaccine_var
  )
  
  ## merging into total time period ##
  if(!epid_start_date==period_start_date & !max(dt$time)==end_date){
    time_before <- seq.Date(period_start_date, (epid_start_date-1), by=7)
    time_after <- seq.Date(max(dt$time)+7, end_date, by=7)
    output <- rbind(data.table(time = time_before, 
                               I1 = 0, I2 = 0, I3 = 0, I4 = 0,
                               IU1 = 0, IU2 = 0, IU3 = 0, IU4 = 0,
                               IV1 = 0, IV2 = 0, IV3 = 0, IV4 = 0), 
                    dt,
                    data.table(time = time_after, 
                               I1 = 0, I2 = 0, I3 = 0, I4 = 0,
                               IU1 = 0, IU2 = 0, IU3 = 0, IU4 = 0,
                               IV1 = 0, IV2 = 0, IV3 = 0, IV4 = 0))
  }
  if(!epid_start_date==period_start_date & max(dt$time)==end_date){
    time_before <- seq.Date(period_start_date, (epid_start_date-1), by=7)
    output <- rbind(data.table(time = time_before, 
                               I1 = 0, I2 = 0, I3 = 0, I4 = 0,
                               IU1 = 0, IU2 = 0, IU3 = 0, IU4 = 0,
                               IV1 = 0, IV2 = 0, IV3 = 0, IV4 = 0), dt)
  }
  if(epid_start_date==period_start_date & !max(dt$time)==end_date){
    time_after <- seq.Date(max(dt$time)+7, end_date, by=7)
    output <- rbind(dt,
                    data.table(time = time_after, 
                               I1 = 0, I2 = 0, I3 = 0, I4 = 0,
                               IU1 = 0, IU2 = 0, IU3 = 0, IU4 = 0,
                               IV1 = 0, IV2 = 0, IV3 = 0, IV4 = 0))
  }
  if(!epid_start_date==period_start_date & max(dt$time)==end_date){
    output <- dt
  }

  output
}

## function to combine multiple epidemics

many_flu <- function(
    country,
    ageing, # = T/F
    ageing_date = NULL, # e.g. '01-04'
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_used,
    vaccine_var,
    doses_dt = NULL,
    vacc_cov_vec = NULL,
    model_age_groups,
    demography_dt
    ){
  
  dates_many_flu <- seq.Date(last_monday(min(epid_inputs$period_start_date)), 
                    last_monday(max(epid_inputs$end_date)), 
                    by=7)
  
  ## which vaccine is used in each year/epidemic? 
  imm_dur_vec <- vaccine_used
  for(i in 1:length(imm_dur_vec)){
    imm_dur_vec[i] <- vacc_type_list[[imm_dur_vec[i]]]$imm_duration
  }
  imm_dur_vec <- as.numeric(imm_dur_vec)
  waning_vec <- 1/(365*imm_dur_vec)
  waning_dt <- data.table(year = year(min(epid_inputs$period_start_date)):(year(min(epid_inputs$period_start_date)) + length(waning_vec) - 1),
                          waning = waning_vec,
                          vacc = vaccine_used)
  
  ## loop over number of epidemics in simulation
  output_dt <- data.table()
  for(epidemic_i in 1:nrow(epid_inputs)){ 
    
    epid_data <- epid_inputs[epidemic_i,]
    
    row <- max(year(epid_data$period_start_date), 
               year(epid_data$epid_start_date) - 1 + 
                 (as.numeric(substr(vacc_calendar_start,4,5)) <
                    as.numeric(substr(ageing_date,4,5))))
    vaccine_used_row <- waning_dt[year %in% row]
    
    flu_epid_output <- one_flu(
      country_code = country, # iso3c string
      demography_dt = demography_dt,
      susceptibility = epid_data$susceptibility, # in [0,1]
      transmissibility = epid_data$transmissibility, 
      matching = epid_data$match,
      initial_infected = rep(epid_data$initial_infected, 4), # 4-vector
      period_start_date = epid_data$period_start_date,
      epid_start_date = epid_data$epid_start_date,
      end_date = epid_data$end_date,
      intended_r0 = epid_data$r0_to_scale,
      vaccine_program = vaccine_used_row,
      vaccine_used = vaccine_used,
      vaccine_var,
      doses_dt,
      vacc_cov_vec
    )
    
    # if(is.na(sum(rowSums(flu_epid_output %>% select(starts_with('I')))))){print(epidemic_i)}
    # plot(flu_epid_output$time, rowSums(flu_epid_output[,c('I1','I2','I3','I4')]),type='l',ylim=c(0,10e6))
    # print(epidemic_i)
    # print(paste0('AR: ',round(unname(unlist(colSums(flu_epid_output[,2:5])/demography_dt[U==T & week==min(demography_dt$week)]$total_as)),2)))
    
    ## merge outputs
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
## no longer using though

change_coverage <- function(data, final_uptake, existing_uptake) {
  
  for(i in 1:length(final_uptake)){
    if(final_uptake[i] > 0){
      data[,i] <- seq(existing_uptake[i], final_uptake[i], by = (final_uptake[i]-existing_uptake[i])/(nrow(data)))[2:(nrow(data)+1)]
    }
  }  
   
  data
}

## function to define vaccine calendar for a given vaccine program

dfn_vaccine_calendar_doses <- function(
    vacc_cov, # intended doses, e.g. c(5e6, 8e7, 0, 0)
    dates_to_run, # seq.Date between start of epidemic and end of period
    key_vacc_date_full,
    efficacy, # e.g. c(0.42, 0.42, 0.42, 0.28)
    no_age_groups, # 4
    no_risk_groups, # 1
    vacc_calendar_start, # e.g. '01-10'
    vacc_calendar_weeks, # e.g. 12
    next_cal = T,
    vacc_cov_next = NULL
){
  
  if(next_cal == T & is.null(vacc_cov_next)){
    stop('Need vacc_cov_next entry')
  }
  
  dates_to_run <- last_monday(dates_to_run)
  
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
  # possible_start_dates <- as.Date(unlist(lapply(as.Date(paste0(vacc_calendar_start, '-', (year(dates_to_run[1]) - 1):(year(dates_to_run[length(dates_to_run)]) + 1)), format = '%d-%m-%Y'),
  #                                               last_monday)))
  # curr_start <- possible_start_dates[dates_to_run[1] - possible_start_dates < 7*vacc_calendar_weeks &
  #                                      dates_to_run[1] - possible_start_dates >= 0 |
  #                                      possible_start_dates - dates_to_run[1] < 365 - 7*vacc_calendar_weeks &
  #                                      possible_start_dates - dates_to_run[1] > 0 ]
  curr_start <- last_monday(key_vacc_date_full)
  
  if(length(curr_start) > 0){
    dates_to_vacc_curr <- seq.Date(curr_start, curr_start + 7*vacc_calendar_weeks - 1, 7)
    rows_to_vacc_curr <- which(dates_to_run %in% dates_to_vacc_curr) 
  }else{
    rows_to_vacc_curr <- 0
  }
  
  if(length(curr_start) > 0 & year(curr_start[1]) == start_year_of_analysis + years_of_analysis - 1){
    next_cal <- F
  }
  
  ## *next* vaccine calendar:
  next_start <- last_monday(key_vacc_date_full + 365)
  if(length(next_start)>1){next_start <- next_start[1]}
  dates_to_vacc_next <- seq.Date(next_start, next_start + 7*vacc_calendar_weeks - 1, 7)
  rows_to_vacc_next <- which(dates_to_run %in% dates_to_vacc_next) 
  
  # ## check number of weeks matches
  # if((! (length(existing_cov)/no_age_groups) == (length(curr_start) + next_cal)) & 
  #    ((length(curr_start) + next_cal) > 0)){
  #   warning("Number of vaccination calendars isn't matching!")
  # }
  # existing_cov_curr <- unlist(ifelse((length(curr_start) > 0), 
  #                                    list(existing_cov[1:4]), list(rep(0,4))))
  # if(next_cal == T){
  #   existing_cov_next <- unlist(ifelse((length(curr_start) > 0), 
  #                                      list(existing_cov[5:8]), list(existing_cov[1:4])))
  # }
  
  vc$calendar[rows_to_vacc_curr, 1:no_age_groups] <- suppressWarnings(t(matrix(vacc_cov/vacc_calendar_weeks, 
                                                              ncol = length(rows_to_vacc_curr), nrow = no_age_groups)))
  # errors arise if the epidemic starts after that year's vaccine program - all good!
                                                                   
  if(next_cal == T){
    vc$calendar[rows_to_vacc_next, 1:no_age_groups] <- t(matrix(vacc_cov_next/vacc_calendar_weeks, 
                                                                ncol = length(rows_to_vacc_next), nrow = no_age_groups))
  }
  
  vc
}

## function to define vaccine calendar for a given vaccine program 

dfn_vaccine_calendar_cov <- function(
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
  
  dates_to_run <- last_monday(dates_to_run)
  
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
  
  if(length(curr_start) > 0 & year(curr_start[1]) == start_year_of_analysis + years_of_analysis - 1){
    next_cal <- F
  }
  
  ## *next* vaccine calendar:
  next_start <- possible_start_dates[possible_start_dates - dates_to_run[1] <= (365+6) &
                                       possible_start_dates - dates_to_run[1] > 0]
  if(length(next_start)>1){warning('next_start includes more than one date')}
  dates_to_vacc_next <- seq.Date(next_start, next_start + 7*vacc_calendar_weeks - 1, 7)
  rows_to_vacc_next <- which(dates_to_run %in% dates_to_vacc_next) 
    
  ## changing vaccine coverage
  
  ## check number of weeks matches
  if((! (length(existing_cov)/no_age_groups) == (length(curr_start) + next_cal)) & 
     ((length(curr_start) + next_cal) > 0)){
    stop("Number of vaccination calendars isn't matching!")
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
last_monday <- function(dates){
  
  output <- c()
  
  for(date in dates){
    move_back <- 7 - which(weekdays(seq.Date(from = as.Date(date) - 6,
                                             to = as.Date(date),
                                             by = 1)) == 'Monday')
    
    monday <- as.Date(date) - move_back
    
    if(!weekdays(monday) == 'Monday'){warning('Monday move wrong')}
    
    output <- c(output, monday)
  }
  
  as.Date(output)
  
}

## function to calculate vaccine doses over same period
flu_doses <- function(
    country,
    ageing, # = T/F
    ageing_date = NULL, # e.g. '01-04'
    epid_inputs, # data.table of vectors of inputs for one_flu
    vaccine_program,
    m_a_g){
  
  dates_many_flu <- seq.Date(last_monday(min(epid_inputs$period_start_date)), 
                             last_monday(max(epid_inputs$end_date)), 
                             by=7)
  
  vacc_name <- names(vacc_type_list)[vaccine_type]
  
  vaccine_used_vec <- if(vaccine_variable == 'doses'){
    # if using doses, NGIVs are often introduced years after the start of the epidemic period
    doses[vacc_scenario == vacc_name & model_age_group==1]$vacc_used
  }else{
    # if using coverage, assuming NGIVs available each year (adjust here if not!)
    rep(vacc_name, (year(epid_dt$end_date[1]) - year(epid_dt$period_start_date[1])))
  }
  
  doses <- fcn_annual_doses(
    country,
    ageing,
    ageing_date,
    dates_in = dates_many_flu,
    demographic_start_year = year(min(epid_inputs$period_start_date)),
    vaccine_used = vaccine_used_vec,
    vaccine_var = vaccine_variable,
    doses_dt = if(vaccine_variable == 'doses'){doses}else{NULL},
    vacc_cov_vec = if(vaccine_variable == 'coverage'){cov_vec}else{NULL},
    init_vaccinated = c(0,0,0,0),
    model_age_groups = m_a_g
  )
  
  doses
}
















