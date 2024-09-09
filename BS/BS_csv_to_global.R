## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")

library(readr)
library(dplyr)
library(data.table)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)

scenario_name <- c('none', 'base', 'low_cov', 'rel_inf', 'depth', 'breadth', 'direct')[5]

if(exists("nat_ann")){rm(nat_ann)}
if(exists("global_weekly")){rm(global_weekly)}

for(c_number in 1:7){
  c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/nat_ann_', c_code, '.Rdata'))){
    print(paste0(c_code, ' exists'))
    load(paste0('data/vacc_output_', scenario_name,'/nat_ann_', c_code, '.Rdata'))
    if(exists('nat_ann')){
      nat_ann <- rbind(nat_ann, nat_ann_c)
    }else{
      nat_ann <- copy(nat_ann_c)
    }
  }
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/global_weekly_', c_code, '.Rdata'))){
    load(paste0('data/vacc_output_', scenario_name,'/global_weekly_', c_code, '.Rdata'))
    if(exists('global_weekly')){
      global_weekly <- rbind(global_weekly, global_weekly_c)
      global_weekly <- global_weekly[, lapply(.SD, sum, na.rm=T), by = c('simulation_index', 'week', 'scenario')]
    }else{
      global_weekly <- copy(global_weekly_c)
    }
  }
  print(c_code)
}

save(nat_ann, file = paste0('data/vacc_output_',scenario_name,'/nat_ann.Rdata'))
save(global_weekly, file = paste0('data/vacc_output_',scenario_name,'/global_weekly.Rdata'))
























