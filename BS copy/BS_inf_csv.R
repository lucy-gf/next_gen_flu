
## 
setwd("~/Desktop/research asst/Global Code")
library(readr)
library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(bayestestR)
library(viridis)
library(patchwork)
library(countrycode)

source("BS/BS_vaccine_programs.R")
source("BS/BS_colors.R")

folder <- c('vacc_output_base','vacc_output_no_scaling','vacc_output_low_cov',
            'vacc_output_breadth','vacc_output_direct','vacc_output_depth',
            'vacc_output_rel_inf')[1]
print(folder)

clusters <- data.table(read_csv("data/new_clustering.csv"))

eti50L <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[25] + 0.75*sort(x)[26])
  }else{stop('length(x) != 100')}
}
eti50U <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[76] + 0.75*sort(x)[75])
  }else{stop('length(x) != 100')}
}
eti95L <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[3] + 0.475*sort(x)[4])
  }else{stop('length(x) != 100')}
}
eti95U <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[98] + 0.475*sort(x)[97])
  }else{stop('length(x) != 100')}
}

load(paste0('data/',folder,'/nat_ann.Rdata')) # nrow 186*30*26*5 = 725400
load(paste0('data/',folder,'/global_weekly.Rdata')) # nrow 1565*100*26 = 4069000

# nat_ann_m <- data.table(melt(nat_ann, id.vars = c("country_code", "year", "scenario", "measure")))
# nat_ann_m[, vacc_status := substr(variable,2,2)]
# nat_ann_m[, strain := substr(variable,4,4)]
# nat_ann_m[grepl('1', variable), age_grp := "0-5"]
# nat_ann_m[grepl('2', variable), age_grp := "5-20"]
# nat_ann_m[grepl('3', variable), age_grp := "20-65"]
# nat_ann_m[grepl('4', variable), age_grp := "65+"]
# nat_ann_m[, ct := substr(scenario,4,4)]
# nat_ann_m[, vt := substr(scenario,9,9)]
# nat_ann_m[grepl('no_vacc', scenario), ct := "0"]
# nat_ann_m[grepl('no_vacc', scenario), vt := "0"]
# 
# nat_ann_novs <- nat_ann_m[measure == 'median',]
# nat_ann_novs[, c('ct','vt','variable','vacc_status','measure','strain','age_grp','year') := NULL]
# nat_ann_novs <- nat_ann_novs[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'scenario')]
# nat_ann_averted <- copy(nat_ann_novs)
# nat_ann_nv <- nat_ann_novs[scenario == 'no_vacc',]
# nat_ann_averted[nat_ann_nv, on = c("country_code"), nv := i.value]
# nat_ann_averted[, averted := nv - value]
# nat_ann_max <- nat_ann_averted[scenario == 'ct_5_vt_5',]
# nat_ann_averted[nat_ann_max, on = c("country_code"), max_av := i.averted]
# nat_ann_averted[, prop_of_max := averted/max_av]
# nat_ann_averted[, ct := substr(scenario, 4,4)]
# nat_ann_averted[, vt := substr(scenario, 9,9)]
# 
# clusters <- data.table(clusters)
# clusters[,country_code := codes]
# nat_ann_averted[clusters, on = c("country_code"), itz := i.cluster_name]

global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]

global_weekly_cum <- copy(global_weekly_m)
global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
global_weekly_cum[, value:=NULL]

global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
global_weekly_cum[global_weekly_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
global_weekly_cum[,averted := nv_cum-cum.sum]
# global_weekly_cum[,nv_cum:=NULL]

global_weekly_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
  dt[, measure := meas]
  global_weekly_meas <- rbind(global_weekly_meas, dt)
  print(meas)
}
global_weekly_meas[,simulation_index:=NULL]

global_cumulative <- copy(global_weekly_meas)
global_cumulative[,averted:=NULL]
global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                            week+scenario~measure,
                                            value.var = "cum.sum")
global_cumulative_wide[, ct := substr(scenario, 4,4)]
global_cumulative_wide[, vt := substr(scenario, 9,9)]

global_averted <- copy(global_weekly_meas)
global_averted[,cum.sum:=NULL]
global_averted_wide <- dcast.data.table(global_averted, 
                                           week+scenario~measure,
                                           value.var = "averted")
global_averted_wide[, ct := substr(scenario, 4,4)]
global_averted_wide[, vt := substr(scenario, 9,9)]

global_av_save <- global_averted_wide[week==max(global_averted_wide$week)]
global_av_save <- global_av_save[!scenario=='no_vacc']
global_av_save[, c('scenario','week','eti50U', 'eti50L'):=NULL]
global_av_save[vt==1, vt_n:='Current']
global_av_save[vt==2, vt_n:='Improved (minimal)']
global_av_save[vt==3, vt_n:='Improved (efficacy)']
global_av_save[vt==4, vt_n:='Improved (breadth)']
global_av_save[vt==5, vt_n:='Universal']
global_av_save[ct==1, ct_n:='0-4']
global_av_save[ct==2, ct_n:='0-10']
global_av_save[ct==3, ct_n:='0-17']
global_av_save[ct==4, ct_n:='65+']
global_av_save[ct==5, ct_n:='0-17, 65+']
setnames(global_av_save, 'ct_n', 'Age-targeting strategy')
setnames(global_av_save, 'vt_n', 'Vaccine type')
global_av_save[, median := round(median/(30*1000000),2)]
global_av_save[, eti95L := round(eti95L/(30*1000000),2)]
global_av_save[, eti95U := round(eti95U/(30*1000000),2)]
global_av_save <- global_av_save[,c('Age-targeting strategy','Vaccine type','median','eti95L','eti95U')]
setnames(global_av_save, 'median','Median infections averted (millions)')
global_av_save[, '95% CI' := paste0('(',eti95L,', ', eti95U,')')]
global_av_save[, c('eti95L','eti95U'):=NULL]
write_csv(global_av_save, paste0("output/plots/",folder,"/averted.csv"))

global_c_save <- global_cumulative_wide[week==max(global_cumulative_wide$week)]
global_c_save[, c('scenario','week','eti50U', 'eti50L'):=NULL]
global_c_save[vt==1, vt_n:='Current']
global_c_save[vt==2, vt_n:='Improved (minimal)']
global_c_save[vt==3, vt_n:='Improved (efficacy)']
global_c_save[vt==4, vt_n:='Improved (breadth)']
global_c_save[vt==5, vt_n:='Universal']
global_c_save[ct==1, ct_n:='0-4']
global_c_save[ct==2, ct_n:='0-10']
global_c_save[ct==3, ct_n:='0-17']
global_c_save[ct==4, ct_n:='65+']
global_c_save[ct==5, ct_n:='0-17, 65+']
setnames(global_c_save, 'ct_n', 'Age-targeting strategy')
setnames(global_c_save, 'vt_n', 'Vaccine type')
global_c_save[, median := round(median/(30*1000000),2)]
global_c_save[, eti95L := round(eti95L/(30*1000000),2)]
global_c_save[, eti95U := round(eti95U/(30*1000000),2)]
global_c_save <- global_c_save[,c('Age-targeting strategy','Vaccine type','median','eti95L','eti95U')]
setnames(global_c_save, 'median','Median cumulative infections (millions)')
global_c_save[, '95% CI' := paste0('(',eti95L,', ', eti95U,')')]
global_c_save[, c('eti95L','eti95U'):=NULL]
write_csv(global_c_save, paste0("output/plots/",folder,"/cumulative.csv"))

scenario_name <- c('base')
if(folder == 'vacc_output_low_cov'){scenario_name <- 'low_cov'}
vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  if(folder == 'vacc_output_low_cov'){
    vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_low_cov/vacc_doses_", c_code, "_",
                                                    scenario_name, ".csv"),show_col_type=F) %>% 
                          mutate(cluster_code = c_code))
  }else{
    vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_base/vacc_doses_", c_code, "_",
                                                    scenario_name, ".csv"),show_col_type=F) %>% 
                          mutate(cluster_code = c_code))
  }
}

vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  group_by(cluster_code, year, name, vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(wasted = grepl('w', name), age_grp = substr(name, 2, 2)) %>%
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))
if(folder=='vacc_output_breadth'){
  vacc_doses_g <- rbind(vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,])
  vacc_doses_g$vacc_type <- rep(1:5, each = nrow(vacc_doses_g)/5)
}

vacc_doses_g <- data.table(vacc_doses_g) 
vacc_save <- vacc_doses_g[,c('year','vacc_type','age_cov','value')]
vacc_save <- vacc_save[,lapply(.SD,sum),by=c('year','vacc_type','age_cov')]
vacc_save[, year:=NULL]
vacc_save <- vacc_save[,lapply(.SD,mean),by=c('vacc_type','age_cov')]
vacc_save[, value := round(value/1000000, 2)]
setnames(vacc_save, 'value','Vaccine doses')
setnames(vacc_save, 'vacc_type','vt')
setnames(vacc_save, 'age_cov','ct')
vacc_save[vt==1, vt_n:='Current']
vacc_save[vt==2, vt_n:='Improved (minimal)']
vacc_save[vt==3, vt_n:='Improved (efficacy)']
vacc_save[vt==4, vt_n:='Improved (breadth)']
vacc_save[vt==5, vt_n:='Universal']
vacc_save[ct==1, ct_n:='0-4']
vacc_save[ct==2, ct_n:='0-10']
vacc_save[ct==3, ct_n:='0-17']
vacc_save[ct==4, ct_n:='65+']
vacc_save[ct==5, ct_n:='0-17, 65+']
setnames(vacc_save, 'vt_n','Vaccine type')
setnames(vacc_save, 'ct_n','Age-targeting strategy')
vacc_save[, c('ct','vt'):=NULL]
vacc_save <- vacc_save[, c('Age-targeting strategy','Vaccine type','Vaccine doses')]
write_csv(vacc_save, paste0("output/plots/",folder,"/doses.csv"))

















