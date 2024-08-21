
### DATA ON THE COUNTRY-SPECIFIC DEMOGRAPHY IN EACH PROJECTION COUNTRY

#setwd("~/Desktop/research asst/Global Code")
# library(ggplot2)
# library(readr)
# library(dplyr)
# library(socialmixr)
# library(fluEvidenceSynthesis)
# library(data.table)
# library(odin)
# library(parallel)
# library(countrycode)
# library(readxl)
# library(tidyverse)

#source('BS/BS_data_fcns.R')
source('fcns/inference_function.R')

## 2025 AGE STRUCTURE, using WPP 2022

## functions to calculate country-specific age structure in a given year 
## (historical estimates for year < 2025, projections for year >= 2025)
## either 1-year or in the model age groups

# pop_hist_WPP_data <- read_xlsx("data/WPP2022_POP_F02_1_POPULATION_5-YEAR_AGE_GROUPS_BOTH_SEXES.xlsx",
#                                sheet = 1, skip = 16) %>% 
#   select(!c(Index, Variant, Notes, `Location code`, `ISO3 Alpha-code`, `ISO2 Alpha-code`,
#             `SDMX code**`, `Parent code`)) %>% 
#   filter(Type == 'Country/Area') %>% rename(name = `Region, subregion, country or area *`)
# write_csv(pop_hist_WPP_data, file = "data_for_BS/pop_hist_WPP_data.csv")
pop_hist_WPP_data <- read_csv('data_for_BS/pop_hist_WPP_data.csv')
# pop_proj_WPP_data <- read_xlsx("data/WPP2022_POP_F02_1_POPULATION_5-YEAR_AGE_GROUPS_BOTH_SEXES.xlsx",
#                                sheet = 2, skip = 16) %>% 
#   select(!c(Index, Variant, Notes, `Location code`, `ISO3 Alpha-code`, `ISO2 Alpha-code`,
#             `SDMX code**`, `Parent code`)) %>% 
#   filter(Type == 'Country/Area') %>% rename(name = `Region, subregion, country or area *`)
# write_csv(pop_proj_WPP_data, file = "data_for_BS/pop_proj_WPP_data.csv")
pop_proj_WPP_data <- read_csv('data_for_BS/pop_proj_WPP_data.csv')

fcn_pop_all <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data %>% filter(name %in% country, Year == year_demog)
  }else{
    data_in <- pop_proj_WPP_data %>% filter(name %in% country, Year == year_demog)
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in %>% select(!c(name, Type, Year)))))
  unlist(lapply(pop_in/5, function(x) rep(x,5)))
}
fcn_pop_model <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data %>% filter(name %in% country, Year == year_demog)
  }else{
    data_in <- pop_proj_WPP_data %>% filter(name %in% country, Year == year_demog)
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in %>% select(!c(name, Type, Year)))))
  c(pop_in[1], sum(pop_in[2:4]), sum(pop_in[5:13]), sum(pop_in[14:21]))
}

### FIXED DEMOGRAPHIC PARAMETERS
# using WPP 2022 PROJECTIONS FOR 2025

## NEED TO CHECK THAT WPP 2022 (AND 2017) HAVE THE RIGHT COUNTRIES TO MATCH PREM ET AL. 

start_year <- 2025

# LT_WPP_data_large <- read_xlsx("data/WPP2022_MORT_F07_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx",
#                                sheet = 2, skip = 16)
# LT_WPP_data <- LT_WPP_data_large %>% select(`Region, subregion, country or area *`, Year, `Age (x)`,
#                                             `Age interval (n)`, `Probability of dying q(x,n)`,
#                                             `Expectation of life e(x)`) %>% 
#   rename(Name = `Region, subregion, country or area *`,
#          x = `Age (x)`, n = `Age interval (n)`, dying_prob = `Probability of dying q(x,n)`,
#          LEx = `Expectation of life e(x)`) %>% 
#   filter(Year == start_year) %>% mutate(dying_prob_annual = dying_prob/n)
# write_csv(LT_WPP_data, file = "data_for_bs/LT_WPP_data.csv")
LT_WPP_data <- read_csv('data_for_BS/LT_WPP_data.csv')

# CBR_WPP_data_large <- read_xlsx("data/UN_PPP2022_Output_CBR.xlsx",
#                                sheet = 2, skip = 16)
# CBR_WPP_data <- CBR_WPP_data_large %>% select(`Region, subregion, country or area *`, Type,
#                                               `2025`) %>%
#   filter(Type == 'Country/Area') %>% 
#   rename(Name = `Region, subregion, country or area *`, CBR = `2025`) %>% 
#   mutate(CBR_indiv = as.numeric(CBR)/1000)
# write_csv(CBR_WPP_data, file = "data_for_bs/CBR_WPP_data.csv")
CBR_WPP_data <- read_csv('data_for_BS/CBR_WPP_data.csv')

# data-frame

demog_data <- data.frame(country = unique(CBR_WPP_data$Name),
                         CBR = NA, M1 = NA, M2 = NA, M3 = NA, M4 = NA)
# remove countries not in the final list when list finalised

# assuming rectangular age distributions within each age group so equal
# weighting on each nqx in the group:
for(i in 1:nrow(demog_data)){
  country_input <- LT_WPP_data %>% filter(Name == demog_data$country[i])
  demog_data$CBR[i] <- (CBR_WPP_data %>% filter(Name == demog_data$country[i]))$CBR_indiv
  demog_data$M1[i] <- (country_input[country_input$x == 0,]$dying_prob_annual +
    4*country_input[country_input$x == 1,]$dying_prob_annual)/5
  demog_data$M2[i] <- mean(country_input[country_input$x %in% 5:19,]$dying_prob_annual)
  demog_data$M3[i] <- mean(country_input[country_input$x %in% 20:60,]$dying_prob_annual)
  demog_data$M4[i] <- 1/(country_input %>% filter(x == 65))$LEx
}


# function to take a vector of c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4') 
# to a one-year stratified population
fcn_yr_res_pop <- function(pop){
  c(rep((pop$U1+pop$V1)/5, 5), rep((pop$U2+pop$V2)/15, 15),
    rep((pop$U3+pop$V3)/45, 45), rep((pop$U4+pop$V4)/5, 5))
}

change_coverage <- function(data, final_uptake){
  sums <- data[nrow(data),]
  # If final uptake is zero in a group then we need to make some kind of assumption on uptake rate over time
  if (any(sums == 0)) {
    col <- which(sums == 0)
    data[,col] <- seq(0, (nrow(data)-1))
    sums <- data[nrow(data),]    
  }
  for(i in 1:nrow(data)) {
    data[i,] <- data[i,]*final_uptake/sums
  }
  data
}

## NEW FUNCTION, VACCINATES USING VACCINE_CALENDAR
fcn_weekly_demog <- function(country, demographic_pars = demog_data,
                                  start_year = 2025, years = 30,
                                  hemisphere,
                                  pop_coverage,
                                  weeks_vaccinating,
                                  imm_duration, 
                                  coverage_pattern,
                                  first_year_all = T,
                                  NH_vacc_date = '01-10',
                                  SH_vacc_date = '01-04',
                                  init_vaccinated = c(0,0,0,0)){
  ## SETTING DEMOGRAPHIC PARAMETERS
  start_pop <- fcn_pop_model(country, start_year)
  pars <- demographic_pars[demographic_pars$country %in% country,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/5, 1/15, 1/45, 0)
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  #RH_matrix <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4, ncol=4) ## CHANGE BACK PLS
  contact_matrix_input <- fcn_contact_matrix(country_name = 'France',
                     country_name_altern = NA,
                     country_name_altern_2 = NA,
                     pop_model = c(10,10,10,10)) # doesn't matter what it is
  
  ## MAKING OUTPUT MATRIX
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  output <- data.frame(country = country[1],
                       week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                       to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                       by = 7),
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4') # useful later
  
  ## INIT CONDITIONS
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  ## SETTING ACTION DATES
  if(hemisphere == 'NH'){
    vacc_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y')
    ageing_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y')
  }else{
    vacc_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
    ageing_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
  }
  dates <- sort(c(as.Date(paste0(SH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y'),
                  as.Date(paste0(NH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y')))
  vacc_first <- isTRUE(vacc_date == dates[1]) # equiv to hemisphere == 'SH'
  
  ## FILLING FIRST SECTION 
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1) # days 1-by-1
  dates_to_run <- seq.Date(from = days1[1], 
                           to = days1[length(days1)],
                           by = 7) # days by week
  # def no vaccination in first section
  coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  vaccine_calendar0 <- as_vaccination_calendar(
    efficacy = rep(0.5, 4),
    dates = dates_to_run,
    coverage = coverage_matrix0,
    no_age_groups = 4, no_risk_groups=1
  )
  
  pop1 <- output[1,]
  input1 <- incidence_function_fit_demog(
    demography_input = fcn_yr_res_pop(pop1),
    calendar_input = vaccine_calendar0,
    contacts = contact_matrix_input, 
    waning_rate = 1/(365*imm_duration),
    vaccination_ratio_input = list(
      # proportion of the population who have *been* vaccinated
      prop_vaccine_compartments = c(unname(pop1['V1']/(pop1['U1'] + pop1['V1'])), 
                                    unname(pop1['V2']/(pop1['U2'] + pop1['V2'])),
                                    unname(pop1['V3']/(pop1['U3'] + pop1['V3'])),
                                    unname(pop1['V4']/(pop1['U4'] + pop1['V4'])),
                                    rep(0,8)), 
      # proportion of those *vaccinated* who were vaccinated *effectively*
      prop_R_vaccinated = rep(0,12), # not relevant
      # proportion of those *unvaccinated* who are immune
      prop_R = rep(0,12)), # not relevant
    begin_date = dates_to_run[1], # CHANGE
    end_date = dates_to_run[length(dates_to_run)],  # CHANGE
    age_groups_model = c(5,20,65)
  ) %>% mutate(t = as.Date(t))
  output[output$week %in% dates_to_run,columns] <- input1 %>% select(!c(t))
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == F){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    # coverage_new = rbind(matrix(c(rep(pop_coverage, each = 12), rep(0,8*12)), nrow=12), rep(0,12))
    # vaccine_calendar <- list(
    #   efficacy = rep(0.5, 12),
    #   dates = dates_to_run[1:(weeks_vaccinating + 1)],
    #   calendar = coverage_new,
    #   no_age_groups = 4, no_risk_groups=1
    # )
    # 
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  for(loop_index in 3:(length(dates) - 1)){
    if((vacc_first + loop_index) %% 2 == 1){
      births <- CBR*sum(output[output$week == action_week,columns])
      output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
      output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
      vaccine_calendar0 <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run,
        coverage = coverage_matrix0,
        no_age_groups = 4, no_risk_groups=1
      )
      
      input2 <-  incidence_function_fit_demog(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar0,
        contacts = contact_matrix_input, 
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }else{
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
      #coverage_matrix <- change_coverage(coverage_matrix, 1/(1 + (1 - pop_coverage)*coverage_pattern/pop_coverage))
      if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
      vaccine_calendar <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run[1:(weeks_vaccinating + 1)],
        coverage = coverage_matrix,
        no_age_groups = 4, no_risk_groups=1
      )
      vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
      
      prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
      for(i in 1:nrow(vaccine_calendar$calendar)){
        for(j in 1:4){
          if(!vaccine_calendar$calendar[i,j] == 0){
            vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
              i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
          }
        }
      }
        
      input2 <-  incidence_function_fit_demog(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar,
        contacts = contact_matrix_input,
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }
  }
  
  ## FILLING IN FINAL SECTION
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1)  # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                   unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                   unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                   unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
    for(i in 1:nrow(vaccine_calendar$calendar)){
      for(j in 1:4){
        if(!vaccine_calendar$calendar[i,j] == 0){
          vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
            i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
        }
      }
    }
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }
  
  output <- output %>% pivot_longer(!c(country, year, week)) %>% 
    mutate(U = grepl('U', name),
           V = grepl('V', name),
           age_grp = case_when(grepl('1', name) ~ '0-5',
                               grepl('2', name) ~ '5-20',
                               grepl('3', name) ~ '20-65',
                               grepl('4', name) ~ '65+')) %>% 
    group_by(week, age_grp) %>% mutate(total_as = sum(value)) %>% ungroup()
  
  output$age_grp <- factor(output$age_grp, levels=unique(output$age_grp))
  
  return(output)
}

# x <- 0.6

# wd_ex <- rbind(fcn_weekly_demog(country = 'Australia',
#                           hemisphere = 'SH',
#                           pop_coverage = c(x, x, 0, x),
#                           imm_duration = 0.5, # in years 
#                           coverage_pattern = 1,
#                           weeks_vaccinating = 12) %>% mutate(vp=1),
#                fcn_weekly_demog(country = 'Australia',
#                                 hemisphere = 'SH',
#                                 pop_coverage = c(x, x, 0, x),
#                                 imm_duration = 1, # in years 
#                                 coverage_pattern = 1,
#                                 weeks_vaccinating = 12) %>% mutate(vp=2),
#                fcn_weekly_demog(country = 'Australia',
#                                 hemisphere = 'SH',
#                                 pop_coverage = c(x, x, 0, x),
#                                 imm_duration = 2, # in years 
#                                 coverage_pattern = 1,
#                                 weeks_vaccinating = 12) %>% mutate(vp=3),
#                fcn_weekly_demog(country = 'Australia',
#                                 hemisphere = 'SH',
#                                 pop_coverage = c(x, x, 0, x),
#                                 imm_duration = 3, # in years 
#                                 coverage_pattern = 1,
#                                 weeks_vaccinating = 12) %>% mutate(vp=4),
#                fcn_weekly_demog(country = 'Australia',
#                                 hemisphere = 'SH',
#                                 pop_coverage = c(x, x, 0, x),
#                                 imm_duration = 5, # in years 
#                                 coverage_pattern = 1,
#                                 weeks_vaccinating = 12) %>% mutate(vp=5))
# 
# # wd_ex %>% #filter(age_grp %in% c('0-5','65+')) %>%
# #   ggplot() +
# #   geom_line(aes(x=week, y=value/1000000, group=name, col=V), lwd=0.8) +
# #   geom_line(aes(x=week, y=total_as/1000000), lwd=0.8) +
# #   theme_minimal() + ylab('Population size, millions') + xlab('Year') +
# #   ggtitle(paste0('Population projections by vaccination status')) +
# #   facet_grid(age_grp~vp, scales='free') +
# #   theme(text=element_text(size=14))
# 
# wd_ex %>% filter(V==T, year %in% 2025:2030) %>% 
#   ggplot() +
#   #geom_hline(yintercept = x, lty=2, alpha=0.8, lwd=1) +
#   geom_line(aes(x=week, y=value/total_as, group=vp, col=as.factor(vp)), lwd=0.8) +
#   theme_bw() + ylab('Vaccinated proportion of population') + xlab('Year') +
#   ggtitle(paste0('Population projections')) +
#   facet_grid(age_grp~., scales='free') + ylim(c(0,NA)) + 
#   theme(text=element_text(size=14)) + 
#   labs(color = "Vaccine type") +
#   scale_color_manual(labels = c("Current", "Improved \nminimal",
#                                 "Improved \nefficacy", "Improved \nbreadth",
#                                 "Universal"), values = 2:6)
# 
# wd_ex %>% filter(V==T, year %in% 2025:2028,
#                  age_grp == '0-5') %>% 
#   ggplot() +
#   geom_hline(yintercept = x, lty=2, alpha=0.8, lwd=1) +
#   geom_line(aes(x=week, y=value/total_as, group=vp, col=as.factor(vp)), lwd=0.8) +
#   theme_bw() + ylab('Vaccinated proportion of population') + xlab('Year') +
#   ggtitle(paste0('Population projections, age group 20-65')) +
#   facet_grid(age_grp~., scales='free') + ylim(c(0,NA)) + 
#   theme(text=element_text(size=14)) + 
#   labs(color = "Vaccine type") +
#   scale_color_manual(labels = c("Current", "Improved \nminimal",
#                                 "Improved \nefficacy", "Improved \nbreadth",
#                                 "Universal"), values = 2:6)
# 
# wd_ex %>% #filter(age_grp %in% c('0-5','65+'), V==T) %>% 
#   ggplot() + 
#   geom_line(aes(x=week, y=value/total_as, group=name, col=V), lwd=0.8) +
#   # geom_line(aes(x=week, y=total_as/1000000), lwd=0.8) +
#   theme_minimal() + xlab('Year') +
#   ggtitle(paste0('Population projections by vaccination status')) + 
#   facet_grid(age_grp~vp, scales='free') + 
#   theme(text=element_text(size=14))
# 
# wd_ex %>% group_by(week, V) %>% mutate(total_V = sum(value)) %>% 
#   ungroup() %>% group_by(week) %>% mutate(total = sum(value)) %>% 
#   ungroup() %>% 
#   ggplot() + 
#   geom_line(aes(x=week, y=total_V/1000000, group=V, col=V), lwd=0.8) +
#   geom_line(aes(x=week, y=total/1000000), lwd=0.8) +
#   theme_minimal() + ylab('Population size, millions') + xlab('Year') +
#   ggtitle(paste0('Population projections by vaccination status')) + 
#   theme(text=element_text(size=14))
# 
# wd_ex %>% 
#   ggplot() + 
#   geom_line(aes(x=week, y=value/total_as, group=name, col=V), lwd=0.8) +
#   geom_line(aes(x=week, y=total_as/total_as), lwd=0.8) +
#   theme_minimal() + ylab('Population size, millions') + xlab('Year') +
#   ggtitle(paste0('Population projections by vaccination status')) + 
#   facet_grid(.~age_grp) + 
#   theme(text=element_text(size=14))


## FUNCTION TO ONLY EXTRACT UNVACCINATED POPULATION ON EACH VACCINATION WEEK

fcn_weekly_demog2 <- function(country, demographic_pars = demog_data,
                             start_year = 2025, years = 30,
                             hemisphere,
                             pop_coverage,
                             weeks_vaccinating,
                             imm_duration, 
                             coverage_pattern,
                             first_year_all = T,
                             NH_vacc_date = '01-10',
                             SH_vacc_date = '01-04',
                             init_vaccinated = c(0,0,0,0)){
  
  output_pop_unvax <- data.frame()
  
  ## SETTING DEMOGRAPHIC PARAMETERS
  start_pop <- fcn_pop_model(country, start_year)
  pars <- demographic_pars[demographic_pars$country %in% country,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/5, 1/15, 1/45, 0)
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  #RH_matrix <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4, ncol=4) ## CHANGE BACK PLS
  contact_matrix_input <- fcn_contact_matrix(country_name = 'France',
                                             country_name_altern = NA,
                                             country_name_altern_2 = NA,
                                             pop_model = c(10,10,10,10)) # doesn't matter what it is
  
  ## MAKING OUTPUT MATRIX
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  output <- data.frame(country = country[1],
                       week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                       to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                       by = 7),
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4') # useful later
  
  ## INIT CONDITIONS
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  ## SETTING ACTION DATES
  if(hemisphere == 'NH'){
    vacc_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y')
    ageing_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y')
  }else{
    vacc_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
    ageing_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
  }
  dates <- sort(c(as.Date(paste0(SH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y'),
                  as.Date(paste0(NH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y')))
  vacc_first <- isTRUE(vacc_date == dates[1]) # equiv to hemisphere == 'SH'
  
  ## FILLING FIRST SECTION 
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1) # days 1-by-1
  dates_to_run <- seq.Date(from = days1[1], 
                           to = days1[length(days1)],
                           by = 7) # days by week
  # def no vaccination in first section
  coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  vaccine_calendar0 <- as_vaccination_calendar(
    efficacy = rep(0.5, 4),
    dates = dates_to_run,
    coverage = coverage_matrix0,
    no_age_groups = 4, no_risk_groups=1
  )
  
  pop1 <- output[1,]
  input1 <- incidence_function_fit_demog(
    demography_input = fcn_yr_res_pop(pop1),
    calendar_input = vaccine_calendar0,
    contacts = contact_matrix_input, 
    waning_rate = 1/(365*imm_duration),
    vaccination_ratio_input = list(
      # proportion of the population who have *been* vaccinated
      prop_vaccine_compartments = c(unname(pop1['V1']/(pop1['U1'] + pop1['V1'])), 
                                    unname(pop1['V2']/(pop1['U2'] + pop1['V2'])),
                                    unname(pop1['V3']/(pop1['U3'] + pop1['V3'])),
                                    unname(pop1['V4']/(pop1['U4'] + pop1['V4'])),
                                    rep(0,8)), 
      # proportion of those *vaccinated* who were vaccinated *effectively*
      prop_R_vaccinated = rep(0,12), # not relevant
      # proportion of those *unvaccinated* who are immune
      prop_R = rep(0,12)), # not relevant
    begin_date = dates_to_run[1], # CHANGE
    end_date = dates_to_run[length(dates_to_run)],  # CHANGE
    age_groups_model = c(5,20,65)
  ) %>% mutate(t = as.Date(t))
  output[output$week %in% dates_to_run,columns] <- input1 %>% select(!c(t))
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == F){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    output_pop_unvax <- rbind(output_pop_unvax, (c(pop2, year(action_week))))
    
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    # coverage_new = rbind(matrix(c(rep(pop_coverage, each = 12), rep(0,8*12)), nrow=12), rep(0,12))
    # vaccine_calendar <- list(
    #   efficacy = rep(0.5, 12),
    #   dates = dates_to_run[1:(weeks_vaccinating + 1)],
    #   calendar = coverage_new,
    #   no_age_groups = 4, no_risk_groups=1
    # )
    # 
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    output_pop_unvax <- rbind(output_pop_unvax, (c(pop2, year(action_week))))
    
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  for(loop_index in 3:(length(dates) - 1)){
    if((vacc_first + loop_index) %% 2 == 1){
      births <- CBR*sum(output[output$week == action_week,columns])
      output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
      output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
      vaccine_calendar0 <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run,
        coverage = coverage_matrix0,
        no_age_groups = 4, no_risk_groups=1
      )
      
      input2 <-  incidence_function_fit_demog(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar0,
        contacts = contact_matrix_input, 
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }else{
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      output_pop_unvax <- rbind(output_pop_unvax, unlist(c(pop2, year(action_week))))
      
      coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
      #coverage_matrix <- change_coverage(coverage_matrix, 1/(1 + (1 - pop_coverage)*coverage_pattern/pop_coverage))
      if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
      vaccine_calendar <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run[1:(weeks_vaccinating + 1)],
        coverage = coverage_matrix,
        no_age_groups = 4, no_risk_groups=1
      )
      vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
      
      prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                     unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                     unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                     unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
      for(i in 1:nrow(vaccine_calendar$calendar)){
        for(j in 1:4){
          if(!vaccine_calendar$calendar[i,j] == 0){
            vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
              i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
          }
        }
      }
      
      input2 <-  incidence_function_fit_demog(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar,
        contacts = contact_matrix_input,
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }
  }
  
  ## FILLING IN FINAL SECTION
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1)  # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    output_pop_unvax <- rbind(output_pop_unvax, unlist(c(pop2, year(action_week))))
    
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                   unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                   unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                   unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
    for(i in 1:nrow(vaccine_calendar$calendar)){
      for(j in 1:4){
        if(!vaccine_calendar$calendar[i,j] == 0){
          vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
            i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
        }
      }
    }
    
    input2 <-  incidence_function_fit_demog(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }
  
  output <- output %>% pivot_longer(!c(country, year, week)) %>% 
    mutate(U = grepl('U', name),
           V = grepl('V', name),
           age_grp = case_when(grepl('1', name) ~ '0-5',
                               grepl('2', name) ~ '5-20',
                               grepl('3', name) ~ '20-65',
                               grepl('4', name) ~ '65+')) %>% 
    group_by(week, age_grp) %>% mutate(total_as = sum(value)) %>% ungroup()
  
  output$age_grp <- factor(output$age_grp, levels=unique(output$age_grp))
  
  return(output_pop_unvax)
}



##_______________________________________________________________________________________________
## FUNCTION TRACKING DOSES GIVEN:
fcn_vaccine_doses <- function(country, demographic_pars = demog_data,
                             start_year = 2025, years = 30,
                             hemisphere,
                             pop_coverage,
                             weeks_vaccinating,
                             imm_duration, 
                             coverage_pattern,
                             first_year_all = T,
                             NH_vacc_date = '01-10',
                             SH_vacc_date = '01-04',
                             init_vaccinated = c(0,0,0,0)){
  ## SETTING DEMOGRAPHIC PARAMETERS
  start_pop <- fcn_pop_model(country, start_year)
  pars <- demographic_pars[demographic_pars$country %in% country,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/5, 1/15, 1/45, 0)
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  #RH_matrix <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4, ncol=4) ## CHANGE BACK PLS
  contact_matrix_input <- fcn_contact_matrix(country_name = 'France',
                                             country_name_altern = NA,
                                             country_name_altern_2 = NA,
                                             pop_model = c(10,10,10,10)) # doesn't matter what it is
  
  ## MAKING OUTPUT MATRIX
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  output <- data.frame(country = country[1],
                       week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                       to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                       by = 7),
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA,
                       vaccs1 = NA, vaccs2 = NA, vaccs3 = NA, vaccs4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4,
                                                            vaccs1, vaccs2, vaccs3, vaccs4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4', 
               'vaccs1', 'vaccs2', 'vaccs3', 'vaccs4') # useful later
  
  ## INIT CONDITIONS
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  ## SETTING ACTION DATES
  if(hemisphere == 'NH'){
    vacc_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y')
    ageing_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y')
  }else{
    vacc_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
    ageing_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
  }
  dates <- sort(c(as.Date(paste0(SH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y'),
                  as.Date(paste0(NH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y')))
  vacc_first <- isTRUE(vacc_date == dates[1]) # equiv to hemisphere == 'SH'
  
  ## FILLING FIRST SECTION 
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1) # days 1-by-1
  dates_to_run <- seq.Date(from = days1[1], 
                           to = days1[length(days1)],
                           by = 7) # days by week
  # def no vaccination in first section
  coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
  vaccine_calendar0 <- as_vaccination_calendar(
    efficacy = rep(0.5, 4),
    dates = dates_to_run,
    coverage = coverage_matrix0,
    no_age_groups = 4, no_risk_groups=1
  )
  
  pop1 <- output[1,]
  input1 <- incidence_function_fit_doses(
    demography_input = fcn_yr_res_pop(pop1),
    calendar_input = vaccine_calendar0,
    contacts = contact_matrix_input, 
    waning_rate = 1/(365*imm_duration),
    vaccination_ratio_input = list(
      # proportion of the population who have *been* vaccinated
      prop_vaccine_compartments = c(unname(pop1['V1']/(pop1['U1'] + pop1['V1'])), 
                                    unname(pop1['V2']/(pop1['U2'] + pop1['V2'])),
                                    unname(pop1['V3']/(pop1['U3'] + pop1['V3'])),
                                    unname(pop1['V4']/(pop1['U4'] + pop1['V4'])),
                                    rep(0,8)), 
      # proportion of those *vaccinated* who were vaccinated *effectively*
      prop_R_vaccinated = rep(0,12), # not relevant
      # proportion of those *unvaccinated* who are immune
      prop_R = rep(0,12)), # not relevant
    begin_date = dates_to_run[1], # CHANGE
    end_date = dates_to_run[length(dates_to_run)],  # CHANGE
    age_groups_model = c(5,20,65)
  ) %>% mutate(t = as.Date(t))
  output[output$week %in% dates_to_run,columns] <- input1 %>% select(!c(t))
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == F){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[2] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    # coverage_new = rbind(matrix(c(rep(pop_coverage, each = 12), rep(0,8*12)), nrow=12), rep(0,12))
    # vaccine_calendar <- list(
    #   efficacy = rep(0.5, 12),
    #   dates = dates_to_run[1:(weeks_vaccinating + 1)],
    #   calendar = coverage_new,
    #   no_age_groups = 4, no_risk_groups=1
    # )
    # 
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = dates[3] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    if(first_year_all == T){
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage)
    }else{
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    }
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
  }
  action_week <- dates_to_run[length(dates_to_run)]
  
  for(loop_index in 3:(length(dates) - 1)){
    if((vacc_first + loop_index) %% 2 == 1){
      births <- CBR*sum(output[output$week == action_week,columns])
      output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
      output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
      vaccine_calendar0 <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run,
        coverage = coverage_matrix0,
        no_age_groups = 4, no_risk_groups=1
      )
      
      input2 <-  incidence_function_fit_doses(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar0,
        contacts = contact_matrix_input, 
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }else{
      days2 <- seq.Date(from = action_week, to = dates[loop_index + 1] + 6, by = 1) # days 1-by-1
      dates_to_run <- seq.Date(from = days2[1], 
                               to = days2[length(days2)],
                               by = 7) # days by week
      pop2 <- output[output$week == action_week,]
      coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
      coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
      #coverage_matrix <- change_coverage(coverage_matrix, 1/(1 + (1 - pop_coverage)*coverage_pattern/pop_coverage))
      if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
      vaccine_calendar <- as_vaccination_calendar(
        efficacy = rep(0.5, 4),
        dates = dates_to_run[1:(weeks_vaccinating + 1)],
        coverage = coverage_matrix,
        no_age_groups = 4, no_risk_groups=1
      )
      vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
      
      prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                     unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                     unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                     unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
      for(i in 1:nrow(vaccine_calendar$calendar)){
        for(j in 1:4){
          if(!vaccine_calendar$calendar[i,j] == 0){
            vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
              i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
          }
        }
      }
      
      input2 <-  incidence_function_fit_doses(
        demography_input = fcn_yr_res_pop(pop2),
        calendar_input = vaccine_calendar,
        contacts = contact_matrix_input,
        waning_rate = 1/(365*imm_duration),
        vaccination_ratio_input = list(
          # proportion of the population who have *been* vaccinated
          prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                        unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                        unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                        unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                        rep(0,8)), 
          # proportion of those *vaccinated* who were vaccinated *effectively*
          prop_R_vaccinated = rep(0,12), # not relevant
          # proportion of those *unvaccinated* who are immune
          prop_R = rep(0,12)), # not relevant
        begin_date = dates_to_run[1], # CHANGE
        end_date = dates_to_run[length(dates_to_run)],  # CHANGE
        age_groups_model = c(5,20,65)
      ) %>% mutate(t = as.Date(t))
      output[output$week %in% dates_to_run,columns] <- input2 %>% select(!c(t))
      action_week <- dates_to_run[length(dates_to_run)]
    }
  }
  
  ## FILLING IN FINAL SECTION
  if(vacc_first == T){
    births <- CBR*sum(output[output$week == action_week,columns])
    output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == action_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == action_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1) # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix0 <- matrix(0, nrow = length(dates_to_run), ncol = 4)
    vaccine_calendar0 <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run,
      coverage = coverage_matrix0,
      no_age_groups = 4, no_risk_groups=1
    )
    
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar0,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }else{
    days2 <- seq.Date(from = action_week, to = output$week[nrow(output)] + 6, by = 1)  # days 1-by-1
    dates_to_run <- seq.Date(from = days2[1], 
                             to = days2[length(days2)],
                             by = 7) # days by week
    pop2 <- output[output$week == action_week,]
    coverage_matrix <- matrix(0, nrow = (weeks_vaccinating + 1), ncol = 4)
    coverage_matrix <- change_coverage(coverage_matrix, pop_coverage/coverage_pattern)
    if(weeks_vaccinating == 1){coverage_matrix <- rbind(coverage_matrix, coverage_matrix[2,])}
    vaccine_calendar <- as_vaccination_calendar(
      efficacy = rep(0.5, 4),
      dates = dates_to_run[1:(weeks_vaccinating + 1)],
      coverage = coverage_matrix,
      no_age_groups = 4, no_risk_groups=1
    )
    vaccine_calendar$calendar[,1:4] <- rbind(coverage_matrix[2:nrow(coverage_matrix),1:4], rep(0,4))
    
    prop_vacc <- c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                   unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                   unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                   unname(pop2['V4']/(pop2['U4'] + pop2['V4'])))
    for(i in 1:nrow(vaccine_calendar$calendar)){
      for(j in 1:4){
        if(!vaccine_calendar$calendar[i,j] == 0){
          vaccine_calendar$calendar[i,j] <- unlist(prop_vacc)[j] + 
            i*(pop_coverage[j] - unlist(prop_vacc)[j])/weeks_vaccinating
        }
      }
    }
    
    input2 <-  incidence_function_fit_doses(
      demography_input = fcn_yr_res_pop(pop2),
      calendar_input = vaccine_calendar,
      contacts = contact_matrix_input, 
      waning_rate = 1/(365*imm_duration),
      vaccination_ratio_input = list(
        # proportion of the population who have *been* vaccinated
        prop_vaccine_compartments = c(unname(pop2['V1']/(pop2['U1'] + pop2['V1'])), 
                                      unname(pop2['V2']/(pop2['U2'] + pop2['V2'])),
                                      unname(pop2['V3']/(pop2['U3'] + pop2['V3'])),
                                      unname(pop2['V4']/(pop2['U4'] + pop2['V4'])),
                                      rep(0,8)), 
        # proportion of those *vaccinated* who were vaccinated *effectively*
        prop_R_vaccinated = rep(0,12), # not relevant
        # proportion of those *unvaccinated* who are immune
        prop_R = rep(0,12)), # not relevant
      begin_date = dates_to_run[1], # CHANGE
      end_date = dates_to_run[length(dates_to_run)],  # CHANGE
      age_groups_model = c(5,20,65)
    ) %>% mutate(t = as.Date(t))
    output[output$week %in% dates_to_run,columns] <- input2 %>% filter(t %in% output$week) %>% select(!c(t))
  }
  
  
  # output <- output %>% pivot_longer(!c(country, year, week)) %>% 
  #   mutate(U = grepl('U', name),
  #          V = grepl('V', name),
  #          age_grp = case_when(grepl('1', name) ~ '0-5',
  #                              grepl('2', name) ~ '5-20',
  #                              grepl('3', name) ~ '20-65',
  #                              grepl('4', name) ~ '65+')) %>% 
  #   group_by(week, age_grp) %>% mutate(total_as = sum(value)) %>% ungroup()
  # 
  # output$age_grp <- factor(output$age_grp, levels=unique(output$age_grp))
  # v_output <- output %>% mutate(
  #   vaccs1nocum = vaccs1 - lag(vaccs1),
  #   vaccs2nocum = vaccs2 - lag(vaccs2),
  #   vaccs3nocum = vaccs3 - lag(vaccs3),
  #   vaccs4nocum = vaccs4 - lag(vaccs4)
  # ) %>% 
  #   select(year, week, vaccs1nocum, vaccs2nocum, vaccs3nocum, vaccs4nocum) %>% 
  #   filter(vaccs1nocum >= 0, vaccs2nocum >= 0, vaccs3nocum >= 0, vaccs4nocum >= 0) %>% 
  #   group_by(year) %>% summarise(vaccs1 = sum(vaccs1nocum), vaccs2 = sum(vaccs2nocum),
  #                                vaccs3 = sum(vaccs3nocum), vaccs4 = sum(vaccs4nocum))
  
  wasted <- output %>% mutate(
    vaccs1nocum = vaccs1 - lag(vaccs1),
    vaccs2nocum = vaccs2 - lag(vaccs2),
    vaccs3nocum = vaccs3 - lag(vaccs3),
    vaccs4nocum = vaccs4 - lag(vaccs4),
    uvaccs1nocum = vaccs1 - lead(vaccs1),
    uvaccs2nocum = vaccs2 - lead(vaccs2),
    uvaccs3nocum = vaccs3 - lead(vaccs3),
    uvaccs4nocum = vaccs4 - lead(vaccs4)
  ) %>% filter((vaccs1nocum + uvaccs1nocum +
                  vaccs2nocum + uvaccs2nocum +
                  vaccs3nocum + uvaccs3nocum +
                  vaccs4nocum + uvaccs4nocum) < 0 & 
                 (vaccs1nocum + vaccs2nocum + vaccs3nocum + vaccs4nocum >= 0)) %>% 
    group_by(year) %>% mutate(
      totv1 = sum(vaccs1nocum), totv2 = sum(vaccs2nocum), 
      totv3 = sum(vaccs3nocum), totv4 = sum(vaccs4nocum),
    ) %>% ungroup() %>% mutate(
      year_lag = year - lag(year)
    ) %>% filter(is.na(year_lag) | year_lag == 1) %>% mutate(
      wasted1 = V1*totv1/U1, wasted2 = V2*totv2/U2,
      wasted3 = V3*totv3/U3, wasted4 = V4*totv4/U4
    ) %>% select(year, totv1, totv2, totv3, totv4, wasted1, wasted2, wasted3, wasted4)
  
  return(wasted)
}


# vd <- fcn_vaccine_doses(country = 'Australia',
#                  hemisphere = 'SH',
#                  pop_coverage = rep(0.6, 4),
#                  imm_duration = 1, # in years
#                  coverage_pattern = 1,
#                  weeks_vaccinating = 12)



## _______________________________________________________________

## DEMOGRAPHY TRACKER (ONLY AGEING)

fcn_ageing <- function(country, demographic_pars = demog_data,
                       start_year = 2025, years = 30){
  # demographic_pars the data-frame of CBR, M1, M2, M3, M4
  start_pop <- fcn_pop_model(country, start_year)
  pars <- demographic_pars[demographic_pars$country %in% country,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  demog_df <- data.frame(country = country[1], 
                         year = rep(start_year:(start_year+years-1), each = 4),
                         age_grp = rep(c('0-5','5-20','20-65','65+')),
                         pop = NA)
  demog_df[demog_df$year == start_year,]$pop <- start_pop
  
  ageing_vec <- c(1/5, 1/15, 1/45, 0)
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  
  for(year_index in (start_year+1):(start_year+years-1)){
    yr_before <- demog_df$pop[demog_df$year==year_index - 1]
    births <- CBR*sum(yr_before)
    demog_df[demog_df$year==year_index,]$pop <- c(births,0,0,0) + 
      c(yr_before%*%RH_matrix)
  }
  demog_df$age_grp <- factor(demog_df$age_grp, levels=unique(demog_df$age_grp))
  
  return(demog_df)
}

dem_example <- fcn_ageing('Italy', demog_data)

### WEEKLY TRACKER OF VACCINATIONS, AGEING AND IMMUNITY WANING

require(deSolve)

waning_de <- function(pop, imm_duration_input, days){
  waning <- function(times, init, parameters) {
    with(as.list(c(init, parameters)), {
      # Unvaccinated population
      U1 = init[1]
      U2 = init[3]
      U3 = init[5]
      U4 = init[7]
      # Vaccinated population
      V1 = init[2]
      V2 = init[4]
      V3 = init[6]
      V4 = init[8]
      # Differential equations
      # dU1 <- (V1/(imm_duration*365))
      # dU2 <- (V2/(imm_duration*365))
      # dU3 <- (V3/(imm_duration*365))
      # dU4 <- (V4/(imm_duration*365))
      # dV1 <- (- V1/(imm_duration*365))
      # dV2 <- (- V2/(imm_duration*365))
      # dV3 <- (- V3/(imm_duration*365))
      # dV4 <- (- V4/(imm_duration*365))
      dU1 <- 0#(V1/(imm_duration*365))
      dU2 <- 0#(V2/(imm_duration*365))
      dU3 <- 0#(V3/(imm_duration*365))
      dU4 <- 0#(V4/(imm_duration*365))
      dV1 <- 0#(- V1/(imm_duration*365))
      dV2 <- 0#(- V2/(imm_duration*365))
      dV3 <- 0#(- V3/(imm_duration*365))
      dV4 <- 0#(- V4/(imm_duration*365))
      return(list(c(dU1, dV1, dU2, dV2, dU3, dV3, dU4, dV4)))
    })
  }
  parameters <- c(imm_duration = imm_duration_input) 
  
  init <- pop 
  times <- 0:length(days)
  
  model <- ode(y=init, times=times, func=waning, parms=parameters)
  data <- data.frame(model)
  weekly_data <- data[seq(1, nrow(data), 7), ]
  return(weekly_data) 
}

## ORIGINAL FUNCTION, WHICH VACCINATES BY MATRIX MULTIPLICATION
fcn_weekly_demog_orig <- function(country, demographic_pars = demog_data,
                                  start_year = 2025, years = 30,
                                  hemisphere = 'NH',
                                  pop_coverage,
                                  imm_duration, 
                                  coverage_pattern,
                                  first_year_all = T,
                                  NH_vacc_date = '01-10',
                                  SH_vacc_date = '01-04',
                                  init_vaccinated = c(0,0,0,0)){
  # demographic_pars the data-frame of CBR, M1, M2, M3, M4
  start_pop <- fcn_pop_model(country, start_year)
  pars <- demographic_pars[demographic_pars$country %in% country,]
  CBR <- pars$CBR; M1 <- pars$M1; M2 <- pars$M2; M3 <- pars$M3; M4 <- pars$M4
  ageing_vec <- c(1/5, 1/15, 1/45, 0)
  mort_vec <- c(M1, M2, M3, M4)
  RH_matrix <- matrix(c((1-ageing_vec[1])*(1-mort_vec[1]), 0, 0, 0,
                        (1-mort_vec[1])*ageing_vec[1], (1-ageing_vec[2])*(1-mort_vec[2]), 0, 0,
                        0, (1-mort_vec[2])*ageing_vec[2], (1-ageing_vec[3])*(1-mort_vec[3]), 0,
                        0, 0, (1-mort_vec[3])*ageing_vec[3], (1-ageing_vec[4])*(1-mort_vec[4])), nrow=4)
  
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  output <- data.frame(country = country[1],
                       week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                       to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                       by = 7),
                       U1 = NA, U2 = NA, U3 = NA, U4 = NA,
                       V1 = NA, V2 = NA, V3 = NA, V4 = NA)
  output <- output %>% mutate(year = year(week)) %>% select(country, year, week,
                                                            U1, V1, U2, V2, U3, V3, U4, V4)
  columns <- c('U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4')
  
  # setting 'action dates' 
  if(hemisphere == 'NH'){
    vacc_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y')
    ageing_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y')
  }else{
    vacc_date <- as.Date(paste0(SH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
    ageing_date <- as.Date(paste0(NH_vacc_date, '-', start_year), format ='%d-%m-%Y') 
  }
  
  # initial conditions
  output[1,c('V1','V2','V3','V4')] <- start_pop*init_vaccinated
  output[1,c('U1','U2','U3','U4')] <- start_pop - output[1,c('V1','V2','V3','V4')] 
  
  dates <- sort(c(as.Date(paste0(SH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y'),
                  as.Date(paste0(NH_vacc_date, '-', start_year:(start_year + years - 1)), format ='%d-%m-%Y')))
  
  days1 <- seq.Date(from = output$week[1], to = dates[1] + 6, by = 1)
  input1 <- waning_de(pop = unname(unlist(output[1,columns])),
                      imm_duration_input = imm_duration,
                      days = days1) %>% select(!time)
  output[output$week %in% days1[1:(length(days1) - 7)],columns] <- input1[2:(nrow(input1)),]
  last_week <- output[output$week %in% days1,]$week[length(output[output$week %in% days1,]$week) - 1]
  next_week <- last_week + 7
  
  if(month(dates[1]) == month(vacc_date)){
    vacc_first <- T
    if(first_year_all == T){
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*(1 - pop_coverage)
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')] + output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*pop_coverage
    }else{
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*(1 - pop_coverage/coverage_pattern)
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')] + output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*pop_coverage/coverage_pattern
    }
  }else{
    vacc_first <- F
    births <- CBR*sum(output[output$week == last_week,columns])
    output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
  }
  
  days2 <- seq.Date(from = next_week, to = dates[2] + 7, by = 1)
  input2 <- waning_de(pop = unname(unlist(output[output$week == next_week,columns])),
                      imm_duration_input = imm_duration,
                      days = days2) %>% select(!time)
  output[output$week %in% days2[1:(length(days2) - 7)],columns] <- input2[2:nrow(input2),]
  
  last_week <- output[output$week %in% days2,]$week[length(output[output$week %in% days2,]$week) - 1]
  next_week <- last_week + 7
  
  if(vacc_first == F){
    if(first_year_all == T){
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*(1 - pop_coverage)
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')] + output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*pop_coverage
    }else{
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*(1 - pop_coverage/coverage_pattern)
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')] + output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*pop_coverage/coverage_pattern
    }
  }else{
    births <- CBR*sum(output[output$week == last_week,columns])
    output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
    output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
  }
  
  for(loop_index in 3:(length(dates))){
    days <- seq.Date(from = next_week, to = dates[loop_index] + 7, by = 1)
    input <- waning_de(pop = unname(unlist(output[output$week == next_week,columns])),
                       imm_duration_input = imm_duration,
                       days = days) %>% select(!time)
    output[output$week %in% days[1:(length(days) - 7)],columns] <- input[2:(nrow(output[output$week %in% days,])),]
    
    last_week <- output[output$week %in% days,]$week[length(output[output$week %in% days,]$week) - 1]
    next_week <- last_week + 7
    
    if((vacc_first + loop_index) %% 2 == 0){
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*(1 - pop_coverage)
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')] + output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]*pop_coverage
    }else{
      births <- CBR*sum(output[output$week == last_week,columns])
      output[output$week == next_week,c('U1', 'U2', 'U3', 'U4')] <- c(births,0,0,0) + unname(unlist(output[output$week == last_week,c('U1', 'U2', 'U3', 'U4')]))%*%RH_matrix
      output[output$week == next_week,c('V1', 'V2', 'V3', 'V4')] <- unname(unlist(output[output$week == last_week,c('V1', 'V2', 'V3', 'V4')]))%*%RH_matrix
    }
  }
  
  ## then add in a final differential equation section
  days3 <- seq.Date(from = next_week, to = output$week[nrow(output)] + 7, by = 1)
  input3 <- waning_de(pop = unname(unlist(output[output$week == next_week,columns])),
                      imm_duration_input = imm_duration,
                      days = days3) %>% select(!time)
  output[output$week %in% days3[2:(length(days3))],columns] <- input3[2:(nrow(output[output$week %in% days3,])),]
  
  output <- output %>% pivot_longer(!c(country, year, week)) %>% 
    mutate(U = grepl('U', name),
           V = grepl('V', name),
           age_grp = case_when(grepl('1', name) ~ '0-5',
                               grepl('2', name) ~ '5-20',
                               grepl('3', name) ~ '20-65',
                               grepl('4', name) ~ '65+')) %>% 
    group_by(week, age_grp) %>% mutate(total_as = sum(value)) %>% ungroup()
  
  output$age_grp <- factor(output$age_grp, levels=unique(output$age_grp))
  
  return(output)
}





