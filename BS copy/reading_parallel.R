## READING BASE50 ETC. RUNS FROM HPC ##
setwd("~/Desktop/research asst/Global Code")

library(data.table)
library(readr)

scenario_name <- c('none', 'base', 'low_cov', 'rel_inf', 'depth', 'breadth', 'direct',
                   'BASE50','LOW20','BASE50_REL_INF','BASE50_BREADTH','BASE50_DEPTH')[12]
print(scenario_name)

for(i in 1:7){
  country_code <- c('GHA','TUR','CHN','GBR','CAN','AUS','ARG')[i]
  print(country_code)

  data_list <- readRDS(paste0('data/vacc_output_', scenario_name, '/vacc_', country_code, '_', scenario_name, '.rds'))

  itz_cases_no_vacc <- data.table(readRDS(paste0("data/vacc_output_base/vacc_", country_code, '_none_ct_1.rds'))[[1]])
  itz_cases_no_vacc <- itz_cases_no_vacc[, year := year(week)]
  itz_cases_no_vacc[, c("country","week") := NULL]
  no_vacc <- itz_cases_no_vacc[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'simulation_index', 'year')]
  no_vacc[, scenario := 'no_vacc']

  econ_inp <- copy(no_vacc)

  for(ct in 1:5){
    for(vt in 1:5){
      itz_cases_ec <- data.table((unname(data_list[[ct]][vt]))[[1]])
      itz_cases_ec[, c("country") := NULL]
      itz_cases_ec[, scenario := paste0('ct_',ct,'_vt_',vt)]
      # print(paste0('ct_',ct,'_vt_',vt))
      econ_inp <- rbind(econ_inp, itz_cases_ec)
    }
  }

  save(econ_inp, file = paste0('data/vacc_output_',scenario_name,'/econ_inp_', country_code, '.Rdata'))
}

# econ_inp[, tot := rowSums(econ_inp[,4:19])]
# avs <- econ_inp[,c('country_code','year','simulation_index','scenario','tot')]
# avs <- avs[, lapply(.SD, mean), by=c('country_code','year','scenario')]
# source('BS/BS_colors.R')
#
# ggplot(avs[country_code=='GHA']) +
#   geom_line(aes(x=year, y=tot/1000000, group=scenario,
#                 col=substr(scenario,9,9))) +
#   theme_bw() + facet_grid(.~substr(scenario,4,4)) +
#   scale_color_manual(values = vt_colors)




