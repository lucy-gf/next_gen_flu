
## CALCULATING AGE-SPECIFIC ANNUAL DOSES GIVEN IN BS_EPID_SIMULATIONS ##

## RUNNING WITH VARIOUS VACCINATION SCENARIOS ##
#setwd("~/Desktop/research asst/Global Code")
source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")
length_simulation <- 30
n_simulations <- 100

start_year <- 2025 
years <- 30

## CHANGE THIS EACH TIME
scenario_name <- c('base', 'low_cov', 'BASE50', 'LOW20')[4]
scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]]

# loops through each ITZ and vaccination program:
for(c_number in 1:7){
    c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                "Southern America")[c_number]
    cluster_code <- c('GHA','TUR','CHN','GBR','CAN','AUS','ARG')[c_number]
    hemisphere <- c('NH','NH','NH','NH','NH','SH','SH')[c_number]

    sampled_epidemics <- rbind(read_csv(paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                               "_", n_simulations,"_",cluster_code, ".csv")))

    ITZ <- clusters %>% filter(cluster_name == c_name)
    n <- nrow(ITZ)

    monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                            to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                            by = 1)) == 'Monday')
    n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                               to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                               by = 7))
    
    n_vacc_programs <- length(scenarios)
      
    vacc_doses <- data.frame(country = rep(rep(unique(ITZ$country), each=years), n_vacc_programs),
                             country_code = rep(rep(unique(ITZ$codes), each=years), n_vacc_programs),
                             vacc_program = rep(1:n_vacc_programs, each = n*years),
                             year = start_year:(start_year + years - 1),
                             v1 = NA, v2 = NA, v3 = NA, v4 = NA,
                             w1 = NA, w2 = NA, w3 = NA, w4 = NA)
    
    for(vacc_prog_num in 1:n_vacc_programs){
      
    vacc_prog <- scenarios[[vacc_prog_num]]
    vacc_details <- vacc_type_list[[vacc_prog$vacc_type]]

      for(country_code_index in unique(vacc_doses$country_code)){
        country_time <- Sys.time()
        
        if(country_code_index == 'XKX'){
          country <- 'Kosovo'
          country_altern <- 'Kosovo (under UNSC res. 1244)'
          country_altern_2 <- 'Kosovo'
        }else{
          country <- countrycode(country_code_index, origin = 'iso3c', destination = 'country.name')
          country_altern <- clusters$country_altern[clusters$codes == country_code_index]
          country_altern_2 <- clusters$country_altern_2[clusters$codes == country_code_index]
        }
        exemplar_code_input <- clusters$cluster_code[clusters$codes == country_code_index]
        vacc_date_input <- unlist(unname(vacc_prog[paste0(hemisphere, '_vacc_date')]))
        
        vacc_doses[vacc_doses$country_code == country_code_index &
                   vacc_doses$vacc_program == vacc_prog_num, 5:12] <-
          fcn_vaccine_doses(country = c(country, country_altern, country_altern_2),
                           pop_coverage = vacc_prog$pop_coverage,
                           weeks_vaccinating = vacc_prog$weeks_vaccinating,
                           first_year_all = vacc_prog$first_year_all,
                           NH_vacc_date = vacc_prog$NH_vacc_date,
                           SH_vacc_date = vacc_prog$SH_vacc_date,
                           init_vaccinated = vacc_prog$init_vaccinated,
                           imm_duration = vacc_details$imm_duration, # in years 
                           coverage_pattern = vacc_details$coverage_pattern,
                           hemisphere = hemisphere,
                           start_year = start_year, years = years) %>% select(!year)
      
      }
    print(paste0('Vaccination program: ', vacc_prog_num,
                 ', cluster: ', c_name))
    write_csv(vacc_doses, file = paste0("data/vacc_doses_", scenario_name,"/vacc_doses_", 
                                        cluster_code, "_",  scenario_name, ".csv"))
    }
}





