
## SIMULATING X-YEAR TIME SERIES FOR EACH CLUSTER/COUNTRY IN INF A+B ##

## DATA AND KEY FUNCTIONS ##
# setwd("~/Desktop/research asst/Global Code")
library(readr)
library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
#library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(zoo)

## ** LOADING IN DATA **

df_cntr_table <- read.csv("data_for_cluster/df_cntr_table.csv") # saved from epid_identif_cont_matrs.R
list_contact_matr <- readRDS("data_for_cluster/list_contact_matr.RDS") # saved from output
df_epid_threshold <- read.csv("data_for_cluster/df_epid_threshold.csv") # saved from output
df_epid_lims <- read.csv("data_for_cluster/df_epid_lims.csv") # saved from epid_identif_cont_matrs.R
matches <- read.csv("data_for_cluster/matches.csv") # saved from vaccine_inputs.R
cov_data <- read.csv("data_for_cluster/cov_data.csv") # saved from vaccine_inputs.R
## Brazil coverage from https://doi.org/10.1186/s12889-020-09409-7 
## (filling in the 2018-2019 values but we aren't using those years' data)
braz_cov <- data.frame(country = rep('Brazil',10),
                       year = 2010:2019,
                       X0.5 = rep(0,10)*(0.5/5) + c(0, 0.29, 0.30, 0.36, 0.79, 0.79, 0.86, 0.75, rep(0.75,2))*(4.5/5),
                       X5.20 = c(0, 0.01, 0.01, 0.06, 0.07, 0.06, 0.08, 0.10, rep(0.10, 2)),
                       X20.65 = c(0, 0.03, 0.03, 0.09, 0.11, 0.11, 0.12, 0.15, rep(0.15, 2)),
                       X65. = rep(0.8, 10))
cov_data <- rbind(cov_data, braz_cov)

infection_delays <- c(0.8, 1.8) # 0.8 and 1.8 day

## need to also load in demog_data and fcn_ageing

## dataframe of all epidemics in exemplar countries ##

unique_epids <- data.frame(country = c(rep(df_cntr_table$country, 
                                         times = df_cntr_table$INF_A_NONSENT),
                                       rep(df_cntr_table$country, 
                                           times = df_cntr_table$INF_B_NONSENT)),
                           code = c(rep(df_cntr_table$COUNTRY_CODE, 
                                        times = df_cntr_table$INF_A_NONSENT),
                                    rep(df_cntr_table$COUNTRY_CODE, 
                                        times = df_cntr_table$INF_B_NONSENT)),
                           strain = c(rep('INF_A', sum(df_cntr_table$INF_A_NONSENT)),
                                      rep('INF_B', sum(df_cntr_table$INF_B_NONSENT))),
                           year = NA,
                           month = NA,
                           day = NA)
unique_epids <- unique_epids %>% mutate(code = case_when(!code=='UK' ~ code,
                                                         code=='UK' ~ 'GBR'))

for(strain_index in c('INF_A','INF_B')){
  for(sel_cntr in df_cntr_table$country){
    if(sel_cntr == 'United Kingdom'){sel_cntr <- c(sel_cntr,'UK')}
    date_inputs <- as.Date(unname(unlist(df_epid_lims %>%
      filter(country %in% sel_cntr, STRAIN == strain_index, metasource == 'NONSENT') %>%
      ungroup() %>% dplyr::select(start_date))))
    unique_epids$year[unique_epids$country %in% sel_cntr &
                        unique_epids$strain == strain_index] <- year(date_inputs) 
    unique_epids$month[unique_epids$country %in% sel_cntr &
                        unique_epids$strain == strain_index] <- month(date_inputs) 
    unique_epids$day[unique_epids$country %in% sel_cntr &
                         unique_epids$strain == strain_index] <- day(date_inputs) 
  }
}

# making all month & days two digits and each country code 3 characters long
unique_epids <- unique_epids %>% mutate(month = ifelse(month < 10, paste0(0, month), month),
                                        day = ifelse(day < 10, paste0(0, day), day), 
                                        epidID = paste0(code, substr(strain, 5, 5),
                                            substr(year, 3, 4), month, day)) %>% arrange(epidID) 

# adding epidemic id to match up with post_samples_merged_simulation
vec <- c()
for(i in 1:length(rle(unique_epids$strain)$lengths)){
  vec <- c(vec, 1:rle(unique_epids$strain)$lengths[i])
}
unique_epids$epidemic_index <- vec

# write_csv(unique_epids, file = "output/unique_epids.csv")

## data-frame with all years ##

allyrs <- data.frame(country = c(rep(df_cntr_table$country, 
                              each = 10),
                          rep(df_cntr_table$country, 
                              each = 10)),
              code = c(rep(df_cntr_table$COUNTRY_CODE, 
                           each = 10),
                       rep(df_cntr_table$COUNTRY_CODE, 
                           each = 10)),
              strain = c(rep('INF_A', 70),
                         rep('INF_B', 70)),
              year = rep(2010:2019, 14))

allyrs <- allyrs %>% mutate(code = case_when(
                                    code != 'UK' ~ code,
                                    code == 'UK' ~ 'GBR'),
                            yearID = paste0(code, substr(strain, 5, 5),
                                            substr(year, 3, 4)),
                            num_epis = NA,
                            epidID = NA)

for(i in 1:nrow(allyrs)){
  allyrs$num_epis[i] = nrow(unique_epids %>% filter(grepl(allyrs$yearID[i], epidID)))
  allyrs$epidID[i] = paste0(unname(unique_epids %>% filter(grepl(allyrs$yearID[i], epidID)) %>% 
                                     dplyr::select(epidID)))
}                            

# write_csv(allyrs, file = "output/allyrs.csv")

## loading in clusters, matching to exemplar countries ##
flu_ITZ_clusters <- read_csv("data_for_BS/flu_ITZ_clusters.csv") 

cluster_codes <- data.frame(cluster = 1:7, 
                            cluster_name = unique(flu_ITZ_clusters$cluster_name),
                            code = c('ARG','CAN','AUS','CHN','GHA','GBR','TUR'))

## new_clusters defined in cluster_expansion.R
clusters <- read_csv("data_for_BS/new_clustering.csv")

clusters <- clusters %>% arrange(cluster_name) %>% mutate(cluster_code = NA)

for(i in 1:nrow(clusters)){
  clusters$cluster_code[i] = cluster_codes$code[cluster_codes$cluster_name == clusters$cluster_name[i]]
}

## loading in the inference parameter samples: ##

#post_samples_merged_simulation <- readRDS("data_for_BS/post_samples_merged_simulation.rds")

## function to bootstrap the parameters: ##

fcn_sample_epidemics <- function(country_code, sampled_years){ # k in 1:7, cluster number
  #sampled_years <- sample(2010:2019, size = years, replace = T)
  #sampled_years <- 2010:2019 # to fix at reality
  years <- nrow(sampled_years)
  sampled_epis <- data.frame(exemplar_code = country_code,
                             simulation_cal_year = 1:years,
                             year = sampled_years$years,
                             N_A_match = sampled_years$N_A_match,
                             N_B_match = sampled_years$N_B_match,
                             S_A_match = sampled_years$S_A_match,
                             S_B_match = sampled_years$S_B_match,
                             num_epis = NA, epidID = NA)
  for(i in 1:nrow(sampled_epis)){
    sampled_epis$num_epis[i] <- sum(unname(unlist(allyrs %>% 
                                                filter(code == country_code,
                                                       year %in% sampled_epis$year[i]) %>% 
                                                select(num_epis))))
  }
  sampled_epis <- sampled_epis[rep(1:years, sampled_epis$num_epis),]
  rownames(sampled_epis) <- 1:nrow(sampled_epis)
  for(i in 1:years){
    sampled_epis$epidID[sampled_epis$simulation_cal_year == i] <- c(unlist(unname(unique_epids %>% 
                                       filter(code == country_code,
                                              year %in% sampled_epis$year[sampled_epis$simulation_cal_year == i]) %>% 
                                       select(epidID))))
  }
  post_samples_merged_country_specific <- post_samples_merged_simulation[[country_code]]
  post_size <- max(post_samples_merged_country_specific$timestep)
  sampled_epis <- sampled_epis %>% mutate(strain = paste0('INF_', substr(epidID, 4,4)),
                                          timestep = sample(1:post_size, 
                                                            size=nrow(sampled_epis), replace=T),
                                          month = substr(epidID, 7, 8),
                                          day = substr(epidID, 9, 10),
                                          sus = NA,
                                          trans = NA, 
                                          infected = NA,
                                          reporting = NA)
  for(i in 1:nrow(sampled_epis)){
    parameters <- post_samples_merged_country_specific %>% filter(strain == paste0('INF_',substr(sampled_epis$epidID[i], 4, 4)),
                        epidemic_id == unique_epids$epidemic_index[unique_epids$epidID == sampled_epis$epidID[i]],
                        timestep == sampled_epis$timestep[i])
    sampled_epis$sus[i] <- unlist(unname(parameters %>% filter(variable == 'sus') %>% select(value)))
    sampled_epis$trans[i] <- unlist(unname(parameters %>% filter(variable == 'trans') %>% select(value)))
    sampled_epis$infected[i] <- unlist(unname(parameters %>% filter(variable == 'infected') %>% select(value)))
    sampled_epis$reporting[i] <- unlist(unname(parameters %>% filter(variable == 'reporting') %>% select(value)))
  }
  return(sampled_epis)
}

## function to calculate the contact matrix for a given country ##

# set model age limits, same as inference
age_limits <- c(0,5,20,65)
age_group_names <- paste0(age_limits,"-", c(age_limits[2:length(age_limits)],99))

load("data_for_BS/contact_all.rdata")
source("fcns/fcns.R")

## function to re-weight Prem et al. contact matrices to a specific demography
## (1-year age groups in pop_all, in the model age groups in pop_model)

fcn_contact_matrix <- function(country_name, country_name_altern,
                               country_name_altern_2, pop_model){
  model_age_groups <<- data.frame(agegroup_name=age_group_names, duration=diff(c(age_limits,120)),
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
  if(sel_cntr_code %in% setdiff(clusters$codes, names(contact_all))){
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
                                             model_agegroups=model_age_groups,
                                             orig_age_groups_duration=standard_age_groups$duration,
                                             orig_age_groups_sizes=standard_age_groups$values)
  # make it reciprocal for the larger group
  contact_matrix_output <- fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                                                       age_group_sizes = model_age_groups$popul)
  return(contact_matrix_output)
}

source('BS/BS_demography.R')

## function to run one epidemic ##

source("fcns/inference_function.R")
source("fcns/2_1b_model_epidemic_yearcross.R")

fcn_run_epidemic <- function(sus, trans, strain,
                             sel_cntr,
                             demography, contact_matrix_input,
                             start_date, end_of_period,
                             vaccine_program,
                             match, vacc_spec_pop, 
                             vacc_date, hemisphere,
                             init = 1, unvaxed_pop,
                             direct = F){ 
  
  dates_to_run <- seq.Date(from = as.Date(start_date), 
                           to = as.Date(end_of_period) + 7,
                           by = 7) # end after the final week
  
  vacc_details <- vacc_type_list[[vaccine_program$vacc_type]]
  
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
  
  prop_vacc_start <- list(
    # proportion of the population who have *been* vaccinated
    prop_vaccine_compartments = c(unname(vacc_spec_pop['V1']/(vacc_spec_pop['U1'] + vacc_spec_pop['V1'])), 
                                  unname(vacc_spec_pop['V2']/(vacc_spec_pop['U2'] + vacc_spec_pop['V2'])),
                                  unname(vacc_spec_pop['V3']/(vacc_spec_pop['U3'] + vacc_spec_pop['V3'])),
                                  unname(vacc_spec_pop['V4']/(vacc_spec_pop['U4'] + vacc_spec_pop['V4'])),
                                  rep(0,8)), 
    # proportion of those *vaccinated* who were vaccinated *effectively*
    prop_R_vaccinated = c(rep(vacc_details$VE[1 + 2*(1 - match)], 3), vacc_details$VE[2 + 2*(1 - match)], 
                          rep(0,8)),
    # proportion of those *unvaccinated* who are immune
    prop_R = rep(0,12)
  )
  
  # pars: c(reporting, transmissibility, susceptibility, initial infected)
  pars <- c(NA, trans, sus, init) 
  
  if(direct == F){
    fit <- incidence_function_fit_VS(
      demography_input = demography,
      parameters = pars,
      calendar_input = vaccine_calendar,
      contact_ids_sample = NA, 
      contacts = contact_matrix_input,
      waning_rate = 1/(365*vacc_details$imm_duration),
      vaccination_ratio_input = prop_vacc_start,
      begin_date = start_date, 
      end_date = end_of_period,  
      year_to_run = year(vaccine_calendar$dates[1]), 
      efficacy_now =rep(0,12), 
      efficacy_next=rep(0,12),
      efficacy_next2 =rep(0,12), 
      previous_summary = NA, 
      age_groups_model = c(5,20,65),
      kappa_input = vaccine_program$relative_vacc_infec
    )
  }
  if(direct == T){
    fit <- incidence_function_fit_VS_direct(
      demography_input = demography,
      parameters = pars,
      calendar_input = vaccine_calendar,
      contact_ids_sample = NA, 
      contacts = contact_matrix_input,
      waning_rate = 1/(365*vacc_details$imm_duration),
      vaccination_ratio_input = prop_vacc_start,
      begin_date = start_date, 
      end_date = end_of_period,  
      year_to_run = year(vaccine_calendar$dates[1]), 
      efficacy_now =rep(0,12), 
      efficacy_next=rep(0,12),
      efficacy_next2 =rep(0,12), 
      previous_summary = NA, 
      age_groups_model = c(5,20,65),
      kappa_input = vaccine_program$relative_vacc_infec
    )
  }
  
  fit <- data.table(fit)
  return(fit)
}

fcn_run_original_epidemic <- function(sus, trans, strain,
                             sel_cntr,
                             demography, contact_matrix_input,
                             start_date, end_of_period,
                             vaccine_program,
                             match, 
                             vacc_date, hemisphere,
                             init = 1){ 
  
  dates_to_run <- c(as.Date(start_date), as.Date(end_of_period) + 7) # end after the final week
  
  vaccine_calendar <- as_vaccination_calendar(
    efficacy = rep(0,length(age_group_names)),
    dates = dates_to_run,
    coverage = matrix(
      0,
      nrow = 3, 
      ncol = length(age_group_names)
    ),
    no_age_groups = length(age_group_names),
    no_risk_groups=1
  )
  year_for_vacc = year(start_date)
  hemisphere_input <- hemisphere
  year_index = year_for_vacc
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(sel_cntr %in% c('Canada' ,'United Kingdom') & 
     month(vaccine_calendar$dates[1]) %in% 1:5){
    year_index <- year_index - 1
  }
  if(year_index < 2010){year_index <- 2010} # quick fix here for if the pushback is too far
  if(sel_cntr %in% cov_data$country){
    coverage = as.numeric(c(unname(cov_data %>% filter(country==sel_cntr,
                                                       year==year_index) %>% select(!country:year)))) 
  }else{coverage <- rep(0,4)}
  
  efficacy = unlist(unname(matches %>% filter(year == year_index,
                                              strain_match == strain,
                                              hemisphere == hemisphere_input
  ) %>% select(!hemisphere:match)))
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(coverage, 3), #proportion of individuals vaccinated
    prop_R_vaccinated = rep(efficacy, 3), #proportion of vaccinated individuals in Rv compartment (as opposed to Sv)
    prop_R = rep(0,12) # proportion of individuals in R (unvaccinated)
    # set to 0 as immunity is not assumed to be related to the previous season
  )
  # print(coverage)
  # print(efficacy)
  
  # pars: c(reporting, transmissibility, susceptibility, initial infected)
  pars <- c(NA, trans, sus, init) 
  
  fit <- incidence_function_fit_VS(
    demography_input = demography,
    parameters = pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = contact_matrix_input,
    waning_rate = 0,
    vaccination_ratio_input = prop_vacc_start,
    begin_date = start_date, 
    end_date = end_of_period,  
    year_to_run = year(vaccine_calendar$dates[1]), 
    efficacy_now =rep(0,12), 
    efficacy_next=rep(0,12),
    efficacy_next2 =rep(0,12), 
    previous_summary = NA, 
    age_groups_model = c(5,20,65),
    kappa_input = vaccine_program$relative_vacc_infec
  )
  fit <- data.table(fit)
  return(fit)
}

# change_coverage <- function(data, final_uptake) {
#   sums <- data[nrow(data),]
#   # If final uptake is zero in a group then we need to make some kind of assumption on uptake rate over time
#   if (any(sums == 0)) {
#     col <- which(sums == 0)
#     data[,col] <- seq(0, (nrow(data)-1))
#     sums <- data[nrow(data),]    
#   }
#   for(i in 1:nrow(data)) {
#     data[i,] <- data[i,]*final_uptake/sums
#   }
#   data
# }
# 
# new_coverage <- change_coverage(coverage_matrix_minimal, rep(0.8, ncol(coverage_matrix_minimal)))

# function to check vaccination status-specific population 6 months into epidemic
fcn_pop_check <- function(sus, trans, strain,
                             sel_cntr,
                             demography, contact_matrix_input,
                             start_date, end_of_period,
                             vaccine_program,
                             match, vacc_spec_pop, 
                             vacc_date, hemisphere,
                             init = 1, forward, unvaxed_pop){ 
  
  dates_to_run <- seq.Date(from = as.Date(start_date), 
                           to = as.Date(end_of_period) + 7,
                           by = 7) # end after the final week
  
  vacc_details <- vacc_type_list[[vaccine_program$vacc_type]]
  
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
  
  prop_vacc_start <- list(
    # proportion of the population who have *been* vaccinated
    prop_vaccine_compartments = c(unname(vacc_spec_pop['V1']/(vacc_spec_pop['U1'] + vacc_spec_pop['V1'])), 
                                  unname(vacc_spec_pop['V2']/(vacc_spec_pop['U2'] + vacc_spec_pop['V2'])),
                                  unname(vacc_spec_pop['V3']/(vacc_spec_pop['U3'] + vacc_spec_pop['V3'])),
                                  unname(vacc_spec_pop['V4']/(vacc_spec_pop['U4'] + vacc_spec_pop['V4'])),
                                  rep(0,8)), 
    # proportion of those *vaccinated* who were vaccinated *effectively*
    prop_R_vaccinated = c(rep(vacc_details$VE[1 + 2*(1 - match)], 3), vacc_details$VE[2 + 2*(1 - match)], 
                          rep(0,8)),
    # proportion of those *unvaccinated* who are immune
    prop_R = rep(0,12)
  )
  
  # pars: c(reporting, transmissibility, susceptibility, initial infected)
  pars <- c(NA, trans, sus, init) 
  
  fit <- incidence_function_fit_demog(
    demography_input = demography,
    calendar_input = vaccine_calendar,
    contacts = contact_matrix_input,
    waning_rate = 1/(365*vacc_details$imm_duration),
    vaccination_ratio_input = prop_vacc_start,
    begin_date = start_date, 
    end_date = end_of_period,  
    age_groups_model = c(5,20,65)
  )
  
  fit <- fit[(forward + 1),]
  return(fit)
}

fcn_run_epidemic_inference <- function(sus, trans, strain,
                             sel_cntr,
                             demography, contact_matrix_input,
                             start_date, year_for_vacc = NA,
                             end_of_period,
                             vaccine_program = NA, init = 1){ 
  
  dates_to_run <- c(as.Date(start_date), as.Date(end_of_period) + 7) # end after the final week
  
  vaccine_calendar <- as_vaccination_calendar(
    efficacy = rep(0,length(age_group_names)),
    dates = dates_to_run,
    coverage = matrix(
      0,
      nrow = 3, 
      ncol = length(age_group_names)
    ),
    no_age_groups = length(age_group_names),
    no_risk_groups=1
  )
  
  hemisphere_input <- df_cntr_table[df_cntr_table$country == sel_cntr,]$hemisphere
  year_index = year_for_vacc
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(hemisphere_input == 'NH' & 
     sel_cntr %in% c('Canada', 'United Kingdom') & # ignore countries with no vacc
     month(vaccine_calendar$dates[1]) %in% 1:7){
     year_index <- year_index - 1
  }
  if(year_index < 2010){year_index <- 2010} # quick fix here for if the pushback is too far
  coverage = as.numeric(c(unname(cov_data %>% filter(country==sel_cntr,
                                                     year==year_index) %>% select(!country:year))))
  efficacy = unlist(unname(matches %>% filter(year == year_index,
                                              strain_match == strain,
                                              hemisphere == hemisphere_input
  ) %>% select(!hemisphere:match)))
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(coverage, 3), #proportion of individuals vaccinated
    prop_R_vaccinated = rep(efficacy, 3), #proportion of vaccinated individuals in Rv compartment (as opposed to Sv)
    prop_R = rep(0,12) # proportion of individuals in R (unvaccinated)
    # set to 0 as immunity is not assumed to be related to the previous season
  )
  
  # pars: c(reporting, transmissibility, susceptibility, initial infected)
  pars <- c(NA, trans, sus, init) 
  
  fit <- incidence_function_fit(
    demography_input = demography,
    parameters = pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = contact_matrix_input,
    waning_rate = 0,
    vaccination_ratio_input = prop_vacc_start,
    begin_date = vaccine_calendar$dates[1], 
    end_date = vaccine_calendar$dates[length(vaccine_calendar$dates)],  
    year_to_run = year(vaccine_calendar$dates[1]), 
    efficacy_now =rep(0,12), 
    efficacy_next=rep(0,12),
    efficacy_next2 =rep(0,12), 
    previous_summary = NA, 
    age_groups_model = c(5,20,65)
  )
  fit <- data.table(fit)
  weekly_cases <- fit[, sum(.SD, na.rm=TRUE), by="time"]
  return(weekly_cases)
}

fcn_run_epidemic_total_inf <- function(sus, trans, strain,
                                       sel_cntr,
                                       demography, contact_matrix_input,
                                       start_date, year_for_vacc = NA,
                                       end_of_period,
                                       vaccine_program = NA, init = 1){ 
  
  dates_to_run <- c(as.Date(start_date), as.Date(end_of_period) + 7) # end after the final week
  
  vaccine_calendar <- as_vaccination_calendar(
    efficacy = rep(0,length(age_group_names)),
    dates = dates_to_run,
    coverage = matrix(
      0,
      nrow = 3, 
      ncol = length(age_group_names)
    ),
    no_age_groups = length(age_group_names),
    no_risk_groups=1
  )
  
  hemisphere_input <- df_cntr_table[df_cntr_table$country == sel_cntr,]$hemisphere
  year_index = year_for_vacc
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(hemisphere_input == 'NH' & 
     sel_cntr %in% c('Canada', 'United Kingdom') & # ignore countries with no vacc
     month(vaccine_calendar$dates[1]) %in% 1:7){
    year_index <- year_index - 1
  }
  if(year_index < 2010){year_index <- 2010} 
  coverage = as.numeric(c(unname(cov_data %>% filter(country==sel_cntr,
                                                     year==year_index) %>% select(!country:year))))
  efficacy = unlist(unname(matches %>% filter(year == year_index,
                                              strain_match == strain,
                                              hemisphere == hemisphere_input
  ) %>% select(!hemisphere:match)))
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(coverage, 3), #proportion of individuals vaccinated
    prop_R_vaccinated = rep(efficacy, 3), #proportion of vaccinated individuals in Rv compartment (as opposed to Sv)
    prop_R = rep(0,12) # proportion of individuals in R (unvaccinated)
    # set to 0 as immunity is not assumed to be related to the previous season
  )
  
  # pars: c(reporting, transmissibility, susceptibility, initial infected)
  pars <- c(NA, trans, sus, init) 
  
  fit <- incidence_function_fit_total_inf(
    demography_input = demography,
    parameters = pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = contact_matrix_input,
    waning_rate = 0,
    vaccination_ratio_input = prop_vacc_start,
    begin_date = vaccine_calendar$dates[1], 
    end_date = vaccine_calendar$dates[length(vaccine_calendar$dates)],  
    year_to_run = year(vaccine_calendar$dates[1]), 
    efficacy_now =rep(0,12), 
    efficacy_next=rep(0,12),
    efficacy_next2 =rep(0,12), 
    previous_summary = NA, 
    age_groups_model = c(5,20,65)
  )
  fit <- data.table(fit)

  return(fit)
}

## function to identify how much further back than the identified start date to ##
## start the epidemic from 10^1 initial cases

fcn_identify_start <- function(sus_input, trans_input, infected, date, country, 
                               strain_input, yr_res_pop, contact_matrix_small,
                               ageing_dates, hemisphere, simulation_year,
                               start_year = 2025, years = 10){
 
  start_date_input <- as.Date(date)
  
  # what year does the most recent ageing occur?
  if(hemisphere=='NH'){
    ageing_date_no_year <- ageing_dates[1]
  }; if(hemisphere=='SH'){
    ageing_date_no_year <- ageing_dates[2]
  }
  ageing_year <- year(start_date_input) - 
    isTRUE(month(start_date_input) < as.numeric(substr(ageing_date_no_year, 4, 5)))
  
  # what date does the most recent ageing occur on?
  if(hemisphere=='NH'){
    ageing_date <- as.Date(paste0(ageing_dates[1], '-', ageing_year), format ='%d-%m-%Y')
  }; if(hemisphere=='SH'){
    ageing_date <- as.Date(paste0(ageing_dates[2], '-', ageing_year), format ='%d-%m-%Y')
  }

  # end date pushed back a year to pre-empt cutoffs
  end_of_period_input <- start_date_input + 365*4
  output_epid <- fcn_run_epidemic_inference(sus = sus_input, 
                                  trans = trans_input, 
                                  strain = strain_input,
                                  sel_cntr = country,
                                  demography = yr_res_pop, 
                                  contact_matrix_input = contact_matrix_small,
                                  start_date = start_date_input, 
                                  year_for_vacc = year(start_date_input), 
                                  end_of_period = end_of_period_input,
                                  init = 1)
  # need the number of cases produced in the first week of an epidemic with the initial condition:
  # threshold_w_inits <- fcn_run_epidemic_inference(sus = sus_input, 
  #                                       trans = trans_input, 
  #                                       strain = strain_input,
  #                                       sel_cntr = country,
  #                                       demography = yr_res_pop, 
  #                                       contact_matrix_input = contact_matrix_small,
  #                                       start_date = start_date_input,
  #                                       year_for_vacc = year(start_date_input), 
  #                                       end_of_period = as.Date(start_date_input + 14),
  #                                       init = log10(unlist(infected)))$V1[1]
  # if the threshold is never crossed, return the time difference in *peaks* instead 
  #if(sum(output_epid$V1 >= threshold_w_inits) == 0){
    epid_w_inits <- fcn_run_epidemic_inference(sus = sus_input, 
                                          trans = trans_input, 
                                          strain = strain_input,
                                          sel_cntr = country,
                                          demography = yr_res_pop, 
                                          contact_matrix_input = contact_matrix_small,
                                          start_date = start_date_input,
                                          year_for_vacc = year(start_date_input), 
                                          end_of_period = end_of_period_input, 
                                          init = log10(unlist(infected)))
    peak_size_1 <- max(epid_w_inits$V1)
    peak_1 <- epid_w_inits$time[match(peak_size_1, epid_w_inits$V1)]
    peak_size_2 <- max(output_epid$V1)
    peak_2 <- output_epid$time[match(peak_size_2, output_epid$V1)]
    peak_diff <- peak_2 - peak_1
    
    # has this pushed the start back *before* the most recent ageing date?
    # OR the start of the inference period
    pushed_back_start <- start_date_input - peak_diff
    if((year(pushed_back_start) < year(start_date_input)) & (simulation_year == 1)){
      epid_find_inits <- fcn_run_epidemic_total_inf(sus = sus_input, 
                                                    trans = trans_input, 
                                                    strain = strain_input,
                                                    sel_cntr = country,
                                                    demography = yr_res_pop, 
                                                    contact_matrix_input = contact_matrix_small,
                                                    start_date = pushed_back_start,
                                                    year_for_vacc = year(start_date_input), 
                                                    end_of_period = end_of_period_input, 
                                                    init = 1)
      nye = as.Date(paste0('31-12-', year(pushed_back_start)), format='%d-%m-%Y') 
      init_nye <- epid_find_inits[epid_find_inits$time <= nye,]$total_inf[nrow(epid_find_inits[epid_find_inits$time <= nye,])]
      init_nye <- init_nye/4 # should be initial infected *in each age group*
      return(c(init_nye, 'nye'))
    }
   else{if(pushed_back_start <= ageing_date){
     epid_find_inits <- fcn_run_epidemic_total_inf(sus = sus_input, 
                                                   trans = trans_input, 
                                                   strain = strain_input,
                                                   sel_cntr = country,
                                                   demography = yr_res_pop, 
                                                   contact_matrix_input = contact_matrix_small,
                                                   start_date = pushed_back_start,
                                                   year_for_vacc = year(start_date_input), 
                                                   end_of_period = end_of_period_input, 
                                                   init = 1)
     init_ageing_date <- epid_find_inits[epid_find_inits$time <= ageing_date,]$total_inf[nrow(epid_find_inits[epid_find_inits$time <= ageing_date,])]
     init_ageing_date <- init_ageing_date/4 # should be initial infected *in each age group*
     return(c(init_ageing_date, 'ageing_date'))
   }
    else{
      return(peak_diff)
    }}
}

## function to scale contact matrices to match r0
fcn_scale_contacts <- function(trans,
                               contacts_small, 
                               age_groups,
                               r0){
  r0_here <- fluEvidenceSynthesis::as_R0(
    transmission_rate = trans,
    contact_matrix = contacts_small, 
    age_groups = age_groups
  )
  return(contacts_small*r0/r0_here)
}

## function to run the bootstrapped epidemics ##

fcn_simulate_epidemics <- function(country_code_input, 
                                   hemisphere, 
                                   sampled_epidemics_input, 
                                   vacc_prog,
                                   start_year = 2025,
                                   years = 30,
                                   r0_scale = T,
                                   direct_input = F){ 
  # name the country
  if(country_code_input == 'XKX'){
    country <- 'Kosovo'
    country_altern <- 'Kosovo (under UNSC res. 1244)'
    country_altern_2 <- 'Kosovo'
  }else{
    country <- countrycode(country_code_input, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == country_code_input]
    country_altern_2 <- clusters$country_altern_2[clusters$codes == country_code_input]
  }
  exemplar_code_input <- clusters$cluster_code[clusters$codes == country_code_input]
  # demography data-frame
  vacc_date_input <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  vacc_details <- vacc_type_list[[vacc_prog$vacc_type]]
  age_structure <- fcn_weekly_demog(country = c(country, country_altern, country_altern_2),
                                     pop_coverage = vacc_prog$pop_coverage,
                                     weeks_vaccinating = vacc_prog$weeks_vaccinating,
                                     first_year_all = vacc_prog$first_year_all,
                                     NH_vacc_date = vacc_prog$NH_vacc_date,
                                     SH_vacc_date = vacc_prog$SH_vacc_date,
                                     init_vaccinated = vacc_prog$init_vaccinated,
                                     imm_duration = vacc_details$imm_duration, # in years 
                                     coverage_pattern = vacc_details$coverage_pattern,
                                     hemisphere = hemisphere,
                                     start_year = start_year, years = years) 
  age_structure2 <- fcn_weekly_demog2(country = c(country, country_altern, country_altern_2),
                                    pop_coverage = vacc_prog$pop_coverage,
                                    weeks_vaccinating = vacc_prog$weeks_vaccinating,
                                    first_year_all = vacc_prog$first_year_all,
                                    NH_vacc_date = vacc_prog$NH_vacc_date,
                                    SH_vacc_date = vacc_prog$SH_vacc_date,
                                    init_vaccinated = vacc_prog$init_vaccinated,
                                    imm_duration = vacc_details$imm_duration, # in years 
                                    coverage_pattern = vacc_details$coverage_pattern,
                                    hemisphere = hemisphere,
                                    start_year = start_year, years = years) 
  if(sum(grepl('X2025L', colnames(age_structure2)))==0){age_structure2 <- rename(age_structure2,c('X2025L'='X2025'))}
  epidemics <- sampled_epidemics_input %>% filter(exemplar_code == exemplar_code_input)
  # each week needs to start on a Monday
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  
  cases_df <- data.frame(week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                         to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                         by = 7),
                         IU1A = 0, IU2A = 0, IU3A = 0, IU4A = 0,
                         IV1A = 0, IV2A = 0, IV3A = 0, IV4A = 0,
                         IU1B = 0, IU2B = 0, IU3B = 0, IU4B = 0,
                         IV1B = 0, IV2B = 0, IV3B = 0, IV4B = 0)
  
  # removing any epidemics which start over a month before the first vaccination date,
  # as these will not be relevant to vaccine effectiveness
  vacc_date <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  epidemics <- epidemics %>% filter(! (simulation_cal_year == 1 & 
                                      (as.numeric(substr(vacc_date, 4, 5)) - as.numeric(month)) > 1))
  
  pop_check_random <- sample(1:nrow(epidemics), 3)
    
  for(j in 1:nrow(epidemics)){
    data_sample <<- epidemics[j,]
    # may need to shift the start date forward to the next monday to match the weeks
    start_date_late <- as.Date(paste0(as.numeric(data_sample$day), '-', as.numeric(data_sample$month), '-', 
                                       (start_year + data_sample$simulation_cal_year - 1)), '%d-%m-%Y')
    original_date <- as.Date(paste0(as.numeric(data_sample$day), '-', data_sample$month, '-', 
                                    data_sample$year), '%d-%m-%Y')
    if(!is.na(data_sample$pushback)){
      pushback <- data_sample$pushback
      start_date_input <- start_date_late - pushback
      init_input <- 1
    }
    if(is.na(data_sample$pushback)){
      if(is.na(data_sample$init_ageing_date)){
        start_date_input <- as.Date(paste0('01-01-', start_year), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_nye)
      }else{
        ageing_day <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 1, 2))
        ageing_month <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 4, 5))
        ageing_year_start <- year(start_date_late) 
        if(month(start_date_late) < ageing_month){
          ageing_year_start <- ageing_year_start - 1
        }
        if((month(start_date_late) = ageing_month) & (day(start_date_late) < ageing_day)){
          ageing_year_start <- ageing_year_start - 1
        }
        start_date_input <- as.Date(paste0(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 
                                           '-', ageing_year_start), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_ageing_date)
      }
    }
    monday_input <- which(weekdays(seq.Date(from = start_date_input, 
                                            to = start_date_input + 6, 
                                            by = 1)) == 'Monday')
    start_date_input <- start_date_input + (monday_input - 1)
    end_of_period_input <- start_date_input + 365*2
    
    # define demography
    year_index <- (start_year + data_sample$simulation_cal_year - 1)
    group_pop_vs <- (age_structure %>% filter(week == start_date_input))$value
    names(group_pop_vs) <- (age_structure %>% filter(week == start_date_input))$name
    group_pop <- (age_structure %>% filter(week == start_date_input, U == T))$total_as 
    contact_matrix <- fcn_contact_matrix(country_name = country,
                                         country_name_altern = country_altern,
                                         country_name_altern_2 = country_altern_2,
                                         pop_model = group_pop)
    yr_res_pop <- c(rep(group_pop[1]/5, 5), rep(group_pop[2]/15, 15),
                    rep(group_pop[3]/45, 45), rep(group_pop[4]/5, 5))
    year_of_first_vacc <- case_when(
      hemisphere == 'NH' & month(start_date_input) < as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input) - 1,
      hemisphere == 'NH' & month(start_date_input) >= as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) < as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) >= as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input) + 1
    )
   
    match_epid <- unname(unlist(data_sample[,paste0(substr(hemisphere, 1, 1), '_', substr(data_sample$strain, 5, 5), '_match')]))
    contact_matrix_small <- t(t(contact_matrix)/group_pop) 
    
    if(r0_scale == T){
      contact_matrix_small <- fcn_scale_contacts(
        trans = data_sample$trans,
        contacts_small = contact_matrix_small, 
        age_groups = group_pop,
        r0 = data_sample$r0)
    }
    
    output_epid <- fcn_run_epidemic(sus = unlist(data_sample$sus), 
                                      trans = unlist(data_sample$trans), 
                                      strain = data_sample$strain,
                                      sel_cntr = country,
                                      demography = yr_res_pop, 
                                      contact_matrix_input = contact_matrix_small,
                                      start_date = start_date_input, 
                                      end_of_period = end_of_period_input,
                                      vaccine_program = vacc_prog,
                                      match = match_epid,
                                      init = init_input,
                                      vacc_spec_pop = group_pop_vs,
                                      vacc_date = vacc_date_input,
                                      hemisphere = hemisphere,
                                      unvaxed_pop = (age_structure2 %>% filter(X2025L == year_of_first_vacc)),
                                      direct = direct_input) 
    ## RANDOMLY CHECKING IN THREE EPIDEMICS THAT AGE- AND VACC STATUS-SPECIFIC POPULATION
    ## SIZES ARE THE SAME TEN WEEKS IN
    if(j %in% pop_check_random){
      k <- 10
      # k <- 1:50
      pop_check_demog <- (age_structure %>% filter(week %in% (start_date_input + 7*(k))))
      pop_check_epid <- fcn_pop_check(sus = unlist(data_sample$sus),
                                      trans = unlist(data_sample$trans),
                                      strain = data_sample$strain,
                                      sel_cntr = country,
                                      demography = yr_res_pop,
                                      contact_matrix_input = contact_matrix_small,
                                      start_date = start_date_input,
                                      end_of_period = end_of_period_input,
                                      vaccine_program = vacc_prog,
                                      match = match_epid,
                                      init = init_input,
                                      vacc_spec_pop = group_pop_vs,
                                      vacc_date = vacc_date_input,
                                      forward = k, hemisphere = hemisphere,
                                      unvaxed_pop = (age_structure2 %>% filter(X2025L == year_of_first_vacc)))
      # rbind((c(pop_check_demog$value, pop_check_demog$week[1])),
      # as.numeric(pop_check_epid))
      # plot(pop_check_demog[pop_check_demog$name=='V4',]$value, type='l', ylim=c(0, 2000000))
      # par(new=T)
      # plot(pop_check_epid$V4, type='l', ylim=c(0, 2000000), col=2)

      if(sum(abs(as.numeric(c(pop_check_demog$value, pop_check_demog$week[1])) -
                 as.numeric(pop_check_epid))) > 100){
        print(paste0('Population sizes differing, ', country, ', epid ', j, ', by = ',
                     sum(abs(as.numeric(c(pop_check_demog$value, pop_check_demog$week[1])) -
                               as.numeric(pop_check_epid)))))
      }
    }
    output_epid <- output_epid %>% filter(!year(time) < start_year,
                                          !year(time) > start_year + years - 1) %>% 
        select(!c(I1, I2, I3, I4))
    if(sum(as.matrix(output_epid %>% select(!time), ncol=8) < 0) > 0){stop('<0')}
    cases_df[cases_df$week %in% output_epid$time, 
               (2:9 + 8*(data_sample$strain=='INF_B'))] <- cases_df[cases_df$week %in% output_epid$time, 
                                                                    (2:9 + 8*(data_sample$strain=='INF_B'))] + 
        as.matrix(output_epid %>% select(!time), ncol=8)
  }
  return(cases_df)
}

## function to run the original epidemics 
## in exemplar countries, with vaccination fixed 
## at observed coverage ##
pop_hist_WPP_data <- read_csv('data_for_BS/pop_hist_WPP_data.csv')
pop_hist_WPP_data2015 <- data.table(pop_hist_WPP_data)[Year==2015]

fcn_simulate_original_epidemics <- function(country_code_input, 
                                   hemisphere, 
                                   sampled_epidemics_input, 
                                   vacc_prog,
                                   start_year = 2010,
                                   years = 10){ 
  # name the country
  if(country_code_input == 'XKX'){
    country <- 'Kosovo'
    country_altern <- 'Kosovo (under UNSC res. 1244)'
    country_altern_2 <- 'Kosovo'
  }else{
    country <- countrycode(country_code_input, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == country_code_input]
    country_altern_2 <- clusters$country_altern_2[clusters$codes == country_code_input]
  }
  exemplar_code_input <- clusters$cluster_code[clusters$codes == country_code_input]
  # demography data-frame
  vacc_date_input <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  vacc_details <- vacc_type_list[[vacc_prog$vacc_type]]
  
  epidemics <- sampled_epidemics_input %>% filter(exemplar_code == exemplar_code_input)
  # each week needs to start on a Monday
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  
  cases_df <- data.frame(week = seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                         to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                         by = 7),
                         IU1A = 0, IU2A = 0, IU3A = 0, IU4A = 0,
                         IV1A = 0, IV2A = 0, IV3A = 0, IV4A = 0,
                         IU1B = 0, IU2B = 0, IU3B = 0, IU4B = 0,
                         IV1B = 0, IV2B = 0, IV3B = 0, IV4B = 0)
  
  # removing any epidemics which start over a month before the first vaccination date,
  # as these will not be relevant to vaccine effectiveness
  vacc_date <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  epidemics <- epidemics %>% filter(! (simulation_cal_year == 1 & 
                                         (as.numeric(substr(vacc_date, 4, 5)) - as.numeric(month)) > 1))
  
  pop_check_random <- sample(1:nrow(epidemics), 3)
  
  for(j in 1:nrow(epidemics)){
    data_sample <<- epidemics[j,]
    # may need to shift the start date forward to the next monday to match the weeks
    start_date_late <- as.Date(paste0(as.numeric(data_sample$day), '-', as.numeric(data_sample$month), '-', 
                                      (start_year + data_sample$simulation_cal_year - 1)), '%d-%m-%Y')
    original_date <- as.Date(paste0(as.numeric(data_sample$day), '-', data_sample$month, '-', 
                                    data_sample$year), '%d-%m-%Y')
    if(!is.na(data_sample$pushback)){
      pushback <- data_sample$pushback
      start_date_input <- start_date_late - pushback
      init_input <- 1
    }
    if(is.na(data_sample$pushback)){
      if(is.na(data_sample$init_ageing_date)){
        start_date_input <- as.Date(paste0('01-01-', start_year), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_nye)
      }else{
        ageing_day <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 1, 2))
        ageing_month <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 4, 5))
        ageing_year_start <- year(start_date_late) 
        if(month(start_date_late) < ageing_month){
          ageing_year_start <- ageing_year_start - 1
        }
        if((month(start_date_late) = ageing_month) & (day(start_date_late) < ageing_day)){
          ageing_year_start <- ageing_year_start - 1
        }
        start_date_input <- as.Date(paste0(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 
                                           '-', ageing_year_start), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_ageing_date)
      }
    }
    monday_input <- which(weekdays(seq.Date(from = start_date_input, 
                                            to = start_date_input + 6, 
                                            by = 1)) == 'Monday')
    start_date_input <- start_date_input + (monday_input - 1)
    end_of_period_input <- start_date_input + 365*2
    
    # define demography
    year_index <- (start_year + data_sample$simulation_cal_year - 1)
    # group_pop <- pop_age(wpp_age(country, 2015), age.limits=c(0,5,20,65))$population 
    long_pop <- 1000*unlist(unname(pop_hist_WPP_data2015[name %in% c(country, country_altern, country_altern_2),4:24]))
    group_pop <- c(long_pop[1], sum(long_pop[2:4]), sum(long_pop[5:13]), sum(long_pop[14:21]))
    contact_matrix <- fcn_contact_matrix(country_name = country,
                                         country_name_altern = country_altern,
                                         country_name_altern_2 = country_altern_2,
                                         pop_model = group_pop)
    yr_res_pop <- c(rep(group_pop[1]/5, 5), rep(group_pop[2]/15, 15),
                    rep(group_pop[3]/45, 45), rep(group_pop[4]/5, 5))
    
    year_of_first_vacc <- case_when(
      hemisphere == 'NH' & month(start_date_input) < as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input) - 1,
      hemisphere == 'NH' & month(start_date_input) >= as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) < as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) >= as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input) + 1
    )
    
    match_epid <- unname(unlist(data_sample[,paste0(substr(hemisphere, 1, 1), '_', substr(data_sample$strain, 5, 5), '_match')]))
    contact_matrix_small <- t(t(contact_matrix)/group_pop) 
    
    contact_matrix_scaled <- fcn_scale_contacts(
      trans = data_sample$trans,
      contacts_small = contact_matrix_small,
      age_groups = group_pop,
      r0 = data_sample$r0)
    
    output_epid <- fcn_run_original_epidemic(sus = unlist(data_sample$sus), 
                                    trans = unlist(data_sample$trans), 
                                    strain = data_sample$strain,
                                    sel_cntr = country,
                                    demography = yr_res_pop, 
                                    contact_matrix_input = contact_matrix_scaled,
                                    start_date = start_date_input, 
                                    end_of_period = end_of_period_input,
                                    vaccine_program = vacc_prog,
                                    match = match_epid,
                                    init = init_input,
                                    vacc_date = vacc_date_input,
                                    hemisphere = hemisphere) 
    output_epid <- output_epid %>% filter(!year(time) < start_year,
                                          !year(time) > start_year + years - 1) %>% 
      select(!c(I1, I2, I3, I4))
    
    # print(paste0('Epidemic ', j, ', attack rate: ', sum(output_epid[,2:8])/sum(yr_res_pop)))
    
    cases_df[cases_df$week %in% output_epid$time, 
             (2:9 + 8*(data_sample$strain=='INF_B'))] <- cases_df[cases_df$week %in% output_epid$time, 
                                                                  (2:9 + 8*(data_sample$strain=='INF_B'))] + 
      as.matrix(output_epid %>% select(!time), ncol=8)
  }
  return(cases_df)
}


fcn_r0_ratio <- function(country_code_input, 
                                   hemisphere, 
                                   sampled_epidemics_input, 
                                   vacc_prog,
                                   start_year = 2025,
                                   years = 30){ 
  # name the country
  if(country_code_input == 'XKX'){
    country <- 'Kosovo'
    country_altern <- 'Kosovo (under UNSC res. 1244)'
    country_altern_2 <- 'Kosovo'
  }else{
    country <- countrycode(country_code_input, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == country_code_input]
    country_altern_2 <- clusters$country_altern_2[clusters$codes == country_code_input]
  }
  exemplar_code_input <- clusters$cluster_code[clusters$codes == country_code_input]
  # demography data-frame
  vacc_date_input <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  vacc_details <- vacc_type_list[[vacc_prog$vacc_type]]
  age_structure <- fcn_weekly_demog(country = c(country, country_altern, country_altern_2),
                                    pop_coverage = vacc_prog$pop_coverage,
                                    weeks_vaccinating = vacc_prog$weeks_vaccinating,
                                    first_year_all = vacc_prog$first_year_all,
                                    NH_vacc_date = vacc_prog$NH_vacc_date,
                                    SH_vacc_date = vacc_prog$SH_vacc_date,
                                    init_vaccinated = vacc_prog$init_vaccinated,
                                    imm_duration = vacc_details$imm_duration, # in years 
                                    coverage_pattern = vacc_details$coverage_pattern,
                                    hemisphere = hemisphere,
                                    start_year = start_year, years = years) 
 
  epidemics <- sampled_epidemics_input %>% filter(exemplar_code == exemplar_code_input)
  # each week needs to start on a Monday
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  
  # removing any epidemics which start over a month before the first vaccination date,
  # as these will not be relevant to vaccine effectiveness
  vacc_date <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
  epidemics <- epidemics %>% filter(! (simulation_cal_year == 1 & 
                                         (as.numeric(substr(vacc_date, 4, 5)) - as.numeric(month)) > 1))
  
  r0_ratios <- c()
  
  for(j in 1:nrow(epidemics)){
    data_sample <<- epidemics[j,]
    # may need to shift the start date forward to the next monday to match the weeks
    start_date_late <- as.Date(paste0(as.numeric(data_sample$day), '-', as.numeric(data_sample$month), '-', 
                                      (start_year + data_sample$simulation_cal_year - 1)), '%d-%m-%Y')
    original_date <- as.Date(paste0(as.numeric(data_sample$day), '-', data_sample$month, '-', 
                                    data_sample$year), '%d-%m-%Y')
    if(!is.na(data_sample$pushback)){
      pushback <- data_sample$pushback
      start_date_input <- start_date_late - pushback
      init_input <- 1
    }
    if(is.na(data_sample$pushback)){
      if(is.na(data_sample$init_ageing_date)){
        start_date_input <- as.Date(paste0('01-01-', start_year), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_nye)
      }else{
        ageing_day <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 1, 2))
        ageing_month <- as.numeric(substr(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 4, 5))
        ageing_year_start <- year(start_date_late) 
        if(month(start_date_late) < ageing_month){
          ageing_year_start <- ageing_year_start - 1
        }
        if((month(start_date_late) = ageing_month) & (day(start_date_late) < ageing_day)){
          ageing_year_start <- ageing_year_start - 1
        }
        start_date_input <- as.Date(paste0(ifelse(hemisphere=='NH', vacc_prog$SH_vacc_date, vacc_prog$NH_vacc_date), 
                                           '-', ageing_year_start), format = '%d-%m-%Y')
        init_input <- log10(data_sample$init_ageing_date)
      }
    }
    monday_input <- which(weekdays(seq.Date(from = start_date_input, 
                                            to = start_date_input + 6, 
                                            by = 1)) == 'Monday')
    start_date_input <- start_date_input + (monday_input - 1)
    end_of_period_input <- start_date_input + 365*2
    
    # define demography
    year_index <- (start_year + data_sample$simulation_cal_year - 1)
    group_pop_vs <- (age_structure %>% filter(week == start_date_input))$value
    names(group_pop_vs) <- (age_structure %>% filter(week == start_date_input))$name
    group_pop <- (age_structure %>% filter(week == start_date_input, U == T))$total_as 
    contact_matrix <- fcn_contact_matrix(country_name = country,
                                         country_name_altern = country_altern,
                                         country_name_altern_2 = country_altern_2,
                                         pop_model = group_pop)
    yr_res_pop <- c(rep(group_pop[1]/5, 5), rep(group_pop[2]/15, 15),
                    rep(group_pop[3]/45, 45), rep(group_pop[4]/5, 5))
    year_of_first_vacc <- case_when(
      hemisphere == 'NH' & month(start_date_input) < as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input) - 1,
      hemisphere == 'NH' & month(start_date_input) >= as.numeric(substr(vacc_prog$SH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) < as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input),
      hemisphere == 'SH' & month(start_date_input) >= as.numeric(substr(vacc_prog$NH_vacc_date, 4, 5)) ~
        year(start_date_input) + 1
    )
    
    match_epid <- unname(unlist(data_sample[,paste0(substr(hemisphere, 1, 1), '_', substr(data_sample$strain, 5, 5), '_match')]))
    contact_matrix_small <- t(t(contact_matrix)/group_pop) 
    r0_here <- fluEvidenceSynthesis::as_R0(
      transmission_rate = data_sample$trans,
      contact_matrix = contact_matrix_small, 
      age_groups = group_pop
    )
    
  r0_ratios <- c(r0_ratios, data_sample$r0/r0_here)
  
  }
  
  return(data.frame(code = country_code_input, epidID = epidemics$epidID,
                    r0 = epidemics$r0, r0_ratio = r0_ratios))
}






