#### NATIONAL DEMOGRAPHY ####

library(data.table)
library(readr)
library(wpp2022)
data(popF); data(popM)

#### model functions ####
source(here::here('next_gen_flu','functions','transmission_model.R'))
source(here::here('next_gen_flu','functions','contact_matr_fcns.R'))

#### POPULATION SIZE FUNCTIONS  ####
pop_hist_WPP_data <- data.table(read_csv(here::here('next_gen_flu','data','pop_hist_WPP_data.csv'), show_col_types = F))
pop_proj_WPP_data <- data.table(read_csv(here::here('next_gen_flu','data','pop_proj_WPP_data.csv'), show_col_types = F))

fcn_pop_all <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data[name %in% country & Year == year_demog]
  }else{
    data_in <- pop_proj_WPP_data[name %in% country & Year == year_demog]
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in[,!c('name', 'Type', 'Year')])))
  unlist(lapply(pop_in/5, function(x) rep(x,5)))
}

fcn_pop_model <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data[name %in% country & Year == year_demog]
  }else{
    data_in <- pop_proj_WPP_data[name %in% country & Year == year_demog]
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in[,!c('name', 'Type', 'Year')])))
  c(pop_in[1], sum(pop_in[2:4]), sum(pop_in[5:13]), sum(pop_in[14:21]))
}

fcn_pop_5_to_75 <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data[name %in% country & Year == year_demog]
  }else{
    data_in <- pop_proj_WPP_data[name %in% country & Year == year_demog]
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in[,!c('name', 'Type', 'Year')])))
  c(pop_in[1:15], sum(pop_in[16:21]))
}

fcn_yr_res_pop <- function(pop){
  # c(rep((pop$U1+pop$V1)/5, 5), rep((pop$U2+pop$V2)/15, 15),
  #   rep((pop$U3+pop$V3)/45, 45), rep((pop$U4+pop$V4)/5, 5))
  c((pop$U1+pop$V1), (pop$U2+pop$V2),
    (pop$U3+pop$V3), (pop$U4+pop$V4))
}

fcn_vri <- function(pop){
  c(unlist(c(unname(pop['V1']/(pop['U1'] + pop['V1'])), 
             unname(pop['V2']/(pop['U2'] + pop['V2'])),
             unname(pop['V3']/(pop['U3'] + pop['V3'])),
             unname(pop['V4']/(pop['U4'] + pop['V4'])))))
}

#### VACCINATING AND AGEING ####

## demographic data calculation ##
LT_WPP_data <- data.table(read_csv(here::here('next_gen_flu','data','LT_WPP_data.csv'), show_col_types = F))
CBR_WPP_data <- data.table(read_csv(here::here('next_gen_flu','data','CBR_WPP_data.csv'), show_col_types = F))
demog_data <- data.frame(country = unique(CBR_WPP_data$Name),
                         CBR = NA, M1 = NA, M2 = NA, M3 = NA, M4 = NA)

# assuming rectangular age distributions within each age group so equal
# weighting on each nqx in the group:
for(i in 1:nrow(demog_data)){
  country_input <- LT_WPP_data[Name == demog_data$country[i]]
  country_input[nrow(country_input), n := 1]
  country_input <- country_input[rep(1:nrow(country_input), country_input$n),]
  country_input$x <- 0:100
  
  demog_data$CBR[i] <- (CBR_WPP_data[Name == demog_data$country[i]])$CBR_indiv
  demog_data$M1[i] <- mean(country_input[country_input$x %in% model_age_groups[1]:(model_age_groups[2]-1),]$dying_prob_annual)
  demog_data$M2[i] <- mean(country_input[country_input$x %in% model_age_groups[2]:(model_age_groups[3]-1),]$dying_prob_annual)
  demog_data$M3[i] <- mean(country_input[country_input$x %in% model_age_groups[3]:(model_age_groups[4]-1),]$dying_prob_annual)
  demog_data$M4[i] <- 1/(country_input[x == model_age_groups[4]])$LEx
}

## function ##
fcn_weekly_demog <- function(country,
                             ageing,
                             ageing_date = NULL, 
                             dates_in,
                             demographic_start_year,
                             vaccine_used,
                             doses_dt,
                             init_vaccinated,
                             model_age_groups){
  
  if((ageing == T) & (is.null(ageing_date))){
    stop('Needs an ageing date.')
  }
  if(ageing == F){
    ageing_date <- '14-02'
  } # setting arbitrary date
  
  ## SETTING DEMOGRAPHIC PARAMETERS
  country_name <- ifelse(country=='XKX', 'Kosovo', countrycode(country, origin='iso3c', destination='country.name'))
  country_names <- unname(unlist(country_itzs_names[country==country_name, c('country','country_altern','country_altern_2')]))
  if(length(country_names)==0){
    country_names <- unname(unlist(country_itzs_names[country_altern==country_name, c('country','country_altern','country_altern_2')]))
  }
  start_pop <- fcn_pop_model(country_names, year(dates_in[1])) 
  pars <- demog_data[demog_data$country %in% country_names,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/(model_age_groups[2:4] - model_age_groups[1:3]), 0) 
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  if(ageing == F){
    CBR <- 0
    RH_matrix <- diag(1, nrow=4)
  } # enforcing no demographic changes
  
  ## VE - arbitrary
  efficacy_input <- c(0,0,0,0)
  imm_dur_vec <- vaccine_used
  for(i in 1:length(imm_dur_vec)){
    imm_dur_vec[i] <- vacc_type_list[[imm_dur_vec[i]]]$imm_duration
  }
  imm_dur_vec <- as.numeric(imm_dur_vec)
  waning_vec <- 1/(365*imm_dur_vec)
  waning_dt <- data.table(year = demographic_start_year:(max(doses$year)),
                          waning = waning_vec,
                          vacc = vaccine_used)
  
  # arbitrary, doesn't matter what it is
  contact_matrix_input <- fcn_contact_matrix(country_name = 'France',
                                             country_name_altern = NA,
                                             country_name_altern_2 = NA,
                                             pop_model = c(10,10,10,10)) 
  
  ## MAKING OUTPUT MATRIX
  output <- data.frame(country = country[1],
                       week = dates_in,
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4') # useful later
  
  ## INIT CONDITIONS
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  ## SETTING ACTION DATES
  start_year <- year(dates_in[1])
  end_year <- year(dates_in[length(dates_in)])
  vacc_date <- as.Date(paste0(vacc_calendar_start, '-', start_year), format ='%d-%m-%Y')
  
  dates <- sort(c(as.Date(paste0(vacc_calendar_start, '-', start_year:end_year), format ='%d-%m-%Y'),
                  as.Date(paste0(ageing_date, '-', start_year:end_year), format ='%d-%m-%Y')))
  vacc_first <- isTRUE(month(dates[1]) == as.numeric(substr(vacc_calendar_start,4,5)))
  dates <- as.Date(unlist(lapply(dates, last_monday)))
  dates <- dates[dates %in% seq.Date(dates_in[1], dates_in[length(dates_in)], by=1)]
  
  ## FILLING FIRST SECTION 
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1) # days 1-by-1
  dates_to_run <- seq.Date(from = days1[1], 
                           to = days1[length(days1)],
                           by = 7) # days by week
  # definitely no vaccination in first section
  coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  vaccine_calendar0 <- as_vaccination_calendar(
    efficacy = rep(0.5, 4),
    dates = dates_to_run,
    coverage = coverage_matrix0,
    no_age_groups = 4, no_risk_groups=1
  )
  
  pop1 <- output[1,]
  input1 <- fcn_vaccinated_demography(
    demography_input = fcn_yr_res_pop(pop1),
    calendar_input = vaccine_calendar0,
    contacts = contact_matrix_input, 
    waning_rate = waning_dt[year == max(demographic_start_year, 
                                        year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
    vaccination_ratio_input = fcn_vri(pop1),
    begin_date = dates_to_run[1], 
    end_date = dates_to_run[length(dates_to_run)],  
    age_groups_model = model_age_groups
  ) %>% mutate(t = as.Date(t))
  output[output$week %in% dates_to_run,columns] <- input1 %>% select(!c(t))
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == F){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      # need to use waning rate of that year's vaccine, UNLESS it's ageing and that comes first, UNLESS it's year 1 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(length(dates)>3){
  for(loop_index in 3:(length(dates) - 1)){
    if((vacc_first + loop_index) %% 2 == 1){
      births <- CBR*sum(output[output$week == action_week,columns])
      output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
      output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) 
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      
      calendar_input <- dfn_vaccine_calendar_doses(
        vacc_cov = c(0,0,0,0),
        dates_to_run = dates_to_run,
        key_vacc_date_full = action_week,
        efficacy = efficacy_input,
        no_age_groups = length(start_pop),
        no_risk_groups = 1,
        vacc_calendar_start = vacc_calendar_start,
        vacc_calendar_weeks = vacc_calendar_weeks,
        next_cal = F
      )
      
      input2 <- fcn_vaccinated_demography(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = calendar_input,
        contacts = contact_matrix_input, 
        waning_rate = waning_dt[year == max(demographic_start_year, 
                                            year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
        vaccination_ratio_input = fcn_vri(pop2),
        begin_date = dates_to_run[1], 
        end_date = dates_to_run[length(dates_to_run)],  
        age_groups_model = model_age_groups
      ) %>% mutate(t = as.Date(t))
      
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
    }else{
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) 
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      
      calendar_input <- dfn_vaccine_calendar_doses(
        vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                            vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
        dates_to_run = dates_to_run,
        key_vacc_date_full = action_week,
        efficacy = efficacy_input,
        no_age_groups = length(start_pop),
        no_risk_groups = 1,
        vacc_calendar_start = vacc_calendar_start,
        vacc_calendar_weeks = vacc_calendar_weeks,
        next_cal = F
      )
      
      input2 <- fcn_vaccinated_demography(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = calendar_input,
        contacts = contact_matrix_input, 
        waning_rate = waning_dt[year == max(demographic_start_year, 
                                            year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
        vaccination_ratio_input = fcn_vri(pop2),
        begin_date = dates_to_run[1], 
        end_date = dates_to_run[length(dates_to_run)],  
        age_groups_model = model_age_groups
      ) %>% mutate(t = as.Date(t))
      
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
    }
    action_week <- dates_to_run[length(dates_to_run)]
    }
  }
  
  ## FILLING IN FINAL SECTION
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  
  output <- output %>% pivot_longer(!c(country, year, week)) %>% 
    mutate(U = grepl('U', name),
           V = grepl('V', name),
           age_grp = case_when(grepl('1', name) ~ paste0(model_age_groups[1], '-', model_age_groups[2]-1),
                               grepl('2', name) ~ paste0(model_age_groups[2], '-', model_age_groups[3]-1),
                               grepl('3', name) ~ paste0(model_age_groups[3], '-', model_age_groups[4]-1),
                               grepl('4', name) ~ paste0(model_age_groups[4], '+'))) %>% 
    group_by(week, age_grp) %>% mutate(total_as = sum(value)) %>% ungroup()
  
  output$age_grp <- factor(output$age_grp, levels=unique(output$age_grp))
  
  output <- data.table(output)
  
  return(output)
  
}


#### CONTACT MATRICES ####
load(here::here('next_gen_flu','data','contact_all.rdata'))

fcn_contact_matrix <- function(country_name, country_name_altern,
                               country_name_altern_2, pop_model){
  age_low_vals <- seq(0,75,5)
  low_vec <- which(age_low_vals %in% (model_age_groups + (-model_age_groups %% 5)))# rounding up to nearest 5-years
  model_age_groups_fcn <- data.frame(agegroup_name=age_group_names, duration=diff(c(model_age_groups,120)),
                                  wpp_agegroup_low=low_vec, wpp_agegroup_high=c(low_vec[2:4]-1, length(age_low_vals)),
                                  popul=pop_model)
  # age groups corresponding to Prem et al. 2021 matrices
  standard_age_groups <- fun_cntr_agestr(i_cntr = c(country_name, country_name_altern, country_name_altern_2),
                                         i_year="2020",
                                         age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120))
  
  # modify contact matrix to correspond to our age groups
  if(country_name=='Kosovo'){
    sel_cntr_code <- 'XKX'
  }else{
    sel_cntr_code <- countrycode(country_name, origin = 'country.name', destination = 'iso3c')
  }
  # using 'similar' countries for those not in contact_all from Prem et al.:
  if(sel_cntr_code %in% setdiff(country_itzs_names$codes, names(contact_all))){
    sel_cntr_code <- case_when(sel_cntr_code == 'AUS' ~ 'NZL',
                               sel_cntr_code == 'SOM' ~ 'ETH',
                               sel_cntr_code == 'LBN' ~ 'SYR',
                               sel_cntr_code == 'JPN' ~ 'KOR',
                               sel_cntr_code == 'TWN' ~ 'CHN',
                               sel_cntr_code == 'XKX' ~ 'SRB',
                               sel_cntr_code == 'NCL' ~ 'VUT',
                               sel_cntr_code == 'GUF' ~ 'SUR',
                               sel_cntr_code == 'HTI' ~ 'DOM',)
    standard_age_groups <- case_when(sel_cntr_code == 'NZL' ~ fun_cntr_agestr(i_cntr = 'New Zealand',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'ETH' ~ fun_cntr_agestr(i_cntr = 'Ethiopia',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SYR' ~ fun_cntr_agestr(i_cntr = 'Syrian Arab Republic',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'KOR' ~ fun_cntr_agestr(i_cntr = 'Republic of Korea',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'CHN' ~ fun_cntr_agestr(i_cntr = 'China',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SRB' ~ fun_cntr_agestr(i_cntr = 'Serbia',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'VUT' ~ fun_cntr_agestr(i_cntr = 'Vanuatu',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SUR' ~ fun_cntr_agestr(i_cntr = 'Suriname',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'DOM' ~ fun_cntr_agestr(i_cntr = 'Dominican Republic',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),)
  }
  C_m_merged_nonrecipr <- fun_create_red_C_m(C_m_full=contact_all[[sel_cntr_code]],
                                             model_agegroups=model_age_groups_fcn,
                                             orig_age_groups_duration=standard_age_groups$duration,
                                             orig_age_groups_sizes=standard_age_groups$values)
  # make it reciprocal for the larger group
  contact_matrix_output <- fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                                               age_group_sizes = model_age_groups_fcn$popul)
  return(contact_matrix_output)
}


#### ANNUAL VACCINE DOSES ####

## function ##
fcn_annual_doses <- function(country,
                             ageing,
                             ageing_date = NULL, 
                             dates_in,
                             demographic_start_year,
                             vaccine_used,
                             doses_dt,
                             init_vaccinated,
                             model_age_groups){
  
  if((ageing == T) & (is.null(ageing_date))){
    stop('Needs an ageing date.')
  }
  if(ageing == F){
    ageing_date <- '14-02'
  } # setting arbitrary date
  
  ## SETTING DEMOGRAPHIC PARAMETERS
  country_name <- ifelse(country=='XKX', 'Kosovo', countrycode(country, origin='iso3c', destination='country.name'))
  country_names <- unname(unlist(country_itzs_names[country==country_name, c('country','country_altern','country_altern_2')]))
  if(length(country_names)==0){
    country_names <- unname(unlist(country_itzs_names[country_altern==country_name, c('country','country_altern','country_altern_2')]))
  }
  start_pop <- fcn_pop_model(country_names, year(dates_in[1])) 
  pars <- demog_data[demog_data$country %in% country_names,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/(model_age_groups[2:4] - model_age_groups[1:3]), 0) 
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  if(ageing == F){
    CBR <- 0
    RH_matrix <- diag(1, nrow=4)
  } # enforcing no demographic changes
  
  ## VE - arbitrary
  efficacy_input <- c(0,0,0,0)
  imm_dur_vec <- vaccine_used
  for(i in 1:length(imm_dur_vec)){
    imm_dur_vec[i] <- vacc_type_list[[imm_dur_vec[i]]]$imm_duration
  }
  imm_dur_vec <- as.numeric(imm_dur_vec)
  waning_vec <- 1/(365*imm_dur_vec)
  waning_dt <- data.table(year = demographic_start_year:(max(doses$year)),
                          waning = waning_vec,
                          vacc = vaccine_used)
  
  # arbitrary, doesn't matter what it is
  contact_matrix_input <- fcn_contact_matrix(country_name = 'France',
                                             country_name_altern = NA,
                                             country_name_altern_2 = NA,
                                             pop_model = c(10,10,10,10)) 
  
  ## MAKING OUTPUT MATRIX
  output <- data.frame(country = country[1],
                       week = dates_in,
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA,
                       vaccs1 = NA, vaccs2 = NA, vaccs3 = NA, vaccs4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4,
                                                            vaccs1, vaccs2, vaccs3, vaccs4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4',
               'vaccs1', 'vaccs2', 'vaccs3', 'vaccs4')
  
  ## INIT CONDITIONS
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  ## SETTING ACTION DATES
  start_year <- year(dates_in[1])
  end_year <- year(dates_in[length(dates_in)])
  vacc_date <- as.Date(paste0(vacc_calendar_start, '-', start_year), format ='%d-%m-%Y')
  
  dates <- sort(c(as.Date(paste0(vacc_calendar_start, '-', start_year:end_year), format ='%d-%m-%Y'),
                  as.Date(paste0(ageing_date, '-', start_year:end_year), format ='%d-%m-%Y')))
  vacc_first <- isTRUE(month(dates[1]) == as.numeric(substr(vacc_calendar_start,4,5)))
  dates <- as.Date(unlist(lapply(dates, last_monday)))
  dates <- dates[dates %in% seq.Date(dates_in[1], dates_in[length(dates_in)], by=1)]
  
  ## FILLING FIRST SECTION 
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1) # days 1-by-1
  dates_to_run <- seq.Date(from = days1[1], 
                           to = days1[length(days1)],
                           by = 7) # days by week
  # definitely no vaccination in first section
  coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  vaccine_calendar0 <- as_vaccination_calendar(
    efficacy = rep(0.5, 4),
    dates = dates_to_run,
    coverage = coverage_matrix0,
    no_age_groups = 4, no_risk_groups=1
  )
  
  pop1 <- output[1,]
  input1 <- fcn_vaccinated_demography_doses(
    demography_input = fcn_yr_res_pop(pop1),
    calendar_input = vaccine_calendar0,
    contacts = contact_matrix_input, 
    waning_rate = waning_dt[year == max(demographic_start_year, 
                                        year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
    vaccination_ratio_input = fcn_vri(pop1),
    begin_date = dates_to_run[1], 
    end_date = dates_to_run[length(dates_to_run)],  
    age_groups_model = model_age_groups
  ) %>% mutate(t = as.Date(t))
  output[output$week %in% dates_to_run,columns] <- input1 %>% select(!c(t))
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == F){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      # need to use waning rate of that year's vaccine, UNLESS it's ageing and that comes first, UNLESS it's year 1 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(length(dates)>3){
    for(loop_index in 3:(length(dates) - 1)){
      if((vacc_first + loop_index) %% 2 == 1){
        births <- CBR*sum(output[output$week == action_week,columns])
        output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
        output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
        days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) 
        dates_to_run <- seq.Date(from = days2[1], 
                                 to = days2[length(days2)],
                                 by = 7) # days by week
        pop2 <- output[output$week == action_week,]
        
        calendar_input <- dfn_vaccine_calendar_doses(
          vacc_cov = c(0,0,0,0),
          dates_to_run = dates_to_run,
          efficacy = efficacy_input,
          no_age_groups = length(start_pop),
          no_risk_groups = 1,
          vacc_calendar_start = vacc_calendar_start,
          vacc_calendar_weeks = vacc_calendar_weeks,
          next_cal = F
        )
        
        input2 <- fcn_vaccinated_demography_doses(
          demography_input = fcn_yr_res_pop(pop2),
          calendar_input = calendar_input,
          contacts = contact_matrix_input, 
          waning_rate = waning_dt[year == max(demographic_start_year, 
                                              year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
          vaccination_ratio_input = fcn_vri(pop2),
          begin_date = dates_to_run[1], 
          end_date = dates_to_run[length(dates_to_run)],  
          age_groups_model = model_age_groups
        ) %>% mutate(t = as.Date(t))
        
        output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      }else{
        days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) 
        dates_to_run <- seq.Date(from = days2[1], 
                                 to = days2[length(days2)],
                                 by = 7) # days by week
        pop2 <- output[output$week == action_week,]
        
        calendar_input <- dfn_vaccine_calendar_doses(
          vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                              vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
          dates_to_run = dates_to_run,
          key_vacc_date_full = action_week,
          efficacy = efficacy_input,
          no_age_groups = length(start_pop),
          no_risk_groups = 1,
          vacc_calendar_start = vacc_calendar_start,
          vacc_calendar_weeks = vacc_calendar_weeks,
          next_cal = F
        )
        
        input2 <- fcn_vaccinated_demography_doses(
          demography_input = fcn_yr_res_pop(pop2),
          calendar_input = calendar_input,
          contacts = contact_matrix_input, 
          waning_rate = waning_dt[year == max(demographic_start_year, 
                                              year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
          vaccination_ratio_input = fcn_vri(pop2),
          begin_date = dates_to_run[1], 
          end_date = dates_to_run[length(dates_to_run)],  
          age_groups_model = model_age_groups
        ) %>% mutate(t = as.Date(t))
        
        output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      }
      action_week <- dates_to_run[length(dates_to_run)]
    }
  }
  
  ## FILLING IN FINAL SECTION
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov = c(0,0,0,0),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) 
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    
    calendar_input <- dfn_vaccine_calendar_doses(
      vacc_cov =  (1 - fcn_vri(pop2))*unname(unlist(doses_dt[year==year(dates_to_run[1]) &
                                          vacc_scenario == vaccine_used[length(vaccine_used)]]$doses))/fcn_yr_res_pop(pop2),
      dates_to_run = dates_to_run,
      key_vacc_date_full = action_week,
      efficacy = efficacy_input,
      no_age_groups = length(start_pop),
      no_risk_groups = 1,
      vacc_calendar_start = vacc_calendar_start,
      vacc_calendar_weeks = vacc_calendar_weeks,
      next_cal = F
    )
    
    input2 <- fcn_vaccinated_demography_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = calendar_input,
      contacts = contact_matrix_input, 
      waning_rate = waning_dt[year == max(demographic_start_year, 
                                          year(dates_to_run[1]) - 1 + vacc_first)]$waning, 
      vaccination_ratio_input = fcn_vri(pop2),
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
      age_groups_model = model_age_groups
    ) %>% mutate(t = as.Date(t))
    
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  
  output_v <- output %>% select(country, year, starts_with('vacc')) %>% 
      pivot_longer(!c(country,year)) %>% mutate(age_grp = substr(name,6,6)) %>% 
      group_by(country, year, age_grp) %>% summarise(vaccs = max(value))
    
  output_v$age_grp <- factor(output_v$age_grp, levels=unique(output_v$age_grp))
    
  output_v <- data.table(output_v)
    
  return(output_v)
  
}




