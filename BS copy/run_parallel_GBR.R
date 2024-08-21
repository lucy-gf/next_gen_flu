
## SIMULATING X-YEAR TIME SERIES FOR EACH CLUSTER/COUNTRY IN INF A+B ##

## RUNNING WITH VARIOUS VACCINATION SCENARIOS ##
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")
length_simulation <- 30
n_simulations <- 100
start_year <- 2025
years <- 30

## change c_number loop for each ITZ 
c_number <- 4
print(c_number)

c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
            "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
            "Southern America")[c_number]
cluster_code <- c('GHA','TUR','CHN','GBR','CAN','AUS','ARG')[c_number]
hemisphere_vacc <- c('NH','NH','NH','NH','NH','SH','SH')[c_number]

scenario_name <- c('none', 'base', 'low_cov', 'rel_inf', 'depth', 'breadth', 'direct',
                   'BASE50','LOW20','BASE50_REL_INF','BASE50_BREADTH','BASE50_DEPTH')[10]

run_in_parallel <- function(k){
  targeting <- paste0('ct_', k)
  if(scenario_name == 'none'){
    scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]]
  }else{
    scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]][1:5 + (k-1)*5]
  }
  
  sampled_epidemics <- rbind(read_csv(paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                             "_", n_simulations,"_",cluster_code, "_wr0.csv"),show_col_types = F))
  
  ITZ <- clusters %>% filter(cluster_name == c_name)
  n <- nrow(ITZ)
  
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                             to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                             by = 7))
  
  cases_rds <- list()
  
  for(vacc_prog_num in 1:length(scenarios)){
    
    cases_30_100 <- data.table(country = rep(rep(unique(ITZ$country), each=n_weeks), n_simulations),
                               country_code = rep(rep(unique(ITZ$codes), each=n_weeks), n_simulations),
                               simulation_index = rep(1:n_simulations, each = n*n_weeks),
                               week = rep(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                                   to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                                   by = 7), length(unique(ITZ$country))*n_simulations),
                               IU1A = 0, IU2A = 0, IU3A = 0, IU4A = 0,
                               IV1A = 0, IV2A = 0, IV3A = 0, IV4A = 0,
                               IU1B = 0, IU2B = 0, IU3B = 0, IU4B = 0,
                               IV1B = 0, IV2B = 0, IV3B = 0, IV4B = 0)
    
    for(sims in 1:n_simulations){
      start_time <- Sys.time()
      for(country_code_index in unique(cases_30_100$country_code)){
        country_time <- Sys.time()
        cases_30_100[cases_30_100$simulation_index == sims &
                       cases_30_100$country_code == country_code_index, 5:20] <-
          fcn_simulate_epidemics(country_code_input = country_code_index,
                                 hemisphere = hemisphere_vacc,
                                 sampled_epidemics_input =
                                   sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                 vacc_prog = scenarios[[vacc_prog_num]]) %>% select(!week)
        # print(paste0(countrycode(country_code_index, origin='iso3c', destination = 'country.name')))
        # print(Sys.time() - country_time)
        write.table(paste0(
          'Vaccination program: ', vacc_prog_num,
          ', Simulations complete: ', sims,
          ', Country = ', country_code_index,
          ', time taken: ', round(Sys.time() - start_time, digits = 2)
        ), file = paste0(cluster_code,'_', scenario_name,'_output_',k,'.txt'))
      }
      print(paste0('Vaccination program: ', vacc_prog_num,
                   ', Simulations complete: ', sims,
                   ', time taken: ', round(Sys.time() - start_time, digits = 2)))
    }
    cases_30_100_short <- copy(cases_30_100)
    cases_30_100_short[, year := year(week)][, week := NULL]
    cases_30_100_short <- cases_30_100_short[,lapply(.SD,sum), by=c('country','country_code','simulation_index','year')]
    cases_rds[[vacc_prog_num]] <- cases_30_100_short
  }
  names(cases_rds) <- names(scenarios)
  
  cases_rds
}

library(parallel)

cases_rds_list <- mclapply(1:5, run_in_parallel, mc.cores=5)

saveRDS(cases_rds_list, file = paste0("data/vacc_output_",scenario_name,
                                      "/vacc_", cluster_code, 
                                      "_", scenario_name,".rds"))

x <- 1
saveRDS(x, file = paste0("data/vacc_output_",scenario_name,"/dud.rds"))








