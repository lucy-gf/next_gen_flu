#### NATIONAL DEMOGRAPHY ####

library(data.table)
library(readr)

source('functions/transmission_model.R')
source('functions/contact_matr_fcns.R')

model_age_groups <- c(0,5,20,65)

#### POPULATION SIZES ####
pop_hist_WPP_data <- data.table(read_csv('data/pop_hist_WPP_data.csv', show_col_types = F))
pop_proj_WPP_data <- data.table(read_csv('data/pop_proj_WPP_data.csv', show_col_types = F))

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

#### VACCINATING AND AGEING ####

## demographic data ##
# LT_WPP_data <- data.table(read_csv('data/LT_WPP_data.csv', show_col_types = F))
# CBR_WPP_data <- data.table(read_csv('data/CBR_WPP_data.csv', show_col_types = F))
# demog_data <- data.frame(country = unique(CBR_WPP_data$Name),
#                          CBR = NA, M1 = NA, M2 = NA, M3 = NA, M4 = NA)
# 
# # assuming rectangular age distributions within each age group so equal
# # weighting on each nqx in the group:
# for(i in 1:nrow(demog_data)){
#   country_input <- LT_WPP_data[Name == demog_data$country[i]]
#   demog_data$CBR[i] <- (CBR_WPP_data[Name == demog_data$country[i]])$CBR_indiv
#   demog_data$M1[i] <- (country_input[country_input$x == 0,]$dying_prob_annual +
#                          4*country_input[country_input$x == 1,]$dying_prob_annual)/5
#   demog_data$M2[i] <- mean(country_input[country_input$x %in% 5:19,]$dying_prob_annual)
#   demog_data$M3[i] <- mean(country_input[country_input$x %in% 20:60,]$dying_prob_annual)
#   demog_data$M4[i] <- 1/(country_input[x == 65])$LEx
# }
# write_csv(demog_data, 'data/demog_data.csv')
demog_data <- data.table(read_csv('data/demog_data.csv', show_col_types=F))

## function ##
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
    
    input2 <- incidence_function_fit_demog(
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
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],   
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
      begin_date = dates_to_run[1], 
      end_date = dates_to_run[length(dates_to_run)],  
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




#### CONTACT MATRICES ####
load("data/contact_all.rdata")

fcn_contact_matrix <- function(country_name, country_name_altern,
                               country_name_altern_2, pop_model){
  model_age_groups_fcn <- data.frame(agegroup_name=age_group_names, duration=diff(c(model_age_groups,120)),
                                  wpp_agegroup_low=c(1,2,5,14), wpp_agegroup_high=c(1,4,13,16),
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

## function to scale to match a given r0?
# fcn_scale_r0 <- function(){
#   
# }


















