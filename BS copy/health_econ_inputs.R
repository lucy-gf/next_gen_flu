## PLOTTING PROJECTION OUTPUTS
setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")

library(readr)
library(dplyr)
library(data.table)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(bayestestR)

scenario_name <- c('none', 'base', 'low_cov', 'rel_inf', 'depth', 'breadth', 'direct')[2]

start_time <- Sys.time()

c_number <- 2
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                        "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                        "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]
print(paste0(c_code, ', ', scenario_name))

# load no vacc cases
itz_cases_no_vacc <- data.table(readRDS(paste0("data/vacc_output_base/vacc_", c_code, '_none_ct_1.rds'))[[1]])

itz_cases_no_vacc_y <- copy(itz_cases_no_vacc)
itz_cases_no_vacc_y <- itz_cases_no_vacc_y[, year := year(week)]
itz_cases_no_vacc_y[, c("country","week") := NULL]  
no_vacc <- itz_cases_no_vacc_y[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'simulation_index', 'year')]
no_vacc[, scenario := 'no_vacc']

econ_inp <- copy(no_vacc)

print('no_vacc done')
print(Sys.time() - start_time)

for(ct in 1:5){
  list <-  readRDS(paste0("data/vacc_output_", scenario_name, "/no-sync.nosync/vacc_", c_code, '_', scenario_name, '_ct_', ct,'.rds'))
  for(vt in 1:5){
    itz_cases <- data.table(list[[vt]])
    
    itz_cases_y <- copy(itz_cases)
    itz_cases_y <- itz_cases_y[, year := year(week)]
    itz_cases_y[, c("country","week") := NULL]  
    dt <- itz_cases_y[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'simulation_index', 'year')]
    dt[, scenario := paste0('ct_',ct,'_vt_',vt)]
    econ_inp <- rbind(econ_inp, dt)

    print(paste0(ct, ' ', vt, ' done'))
    print(Sys.time() - start_time)
  }
}

save(econ_inp, file = paste0('data/vacc_output_',scenario_name,'/econ_inp_', c_code, '.Rdata'))

print(Sys.time() - start_time)








