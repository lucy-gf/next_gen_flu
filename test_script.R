#### TEST SCRIPT ####

#### load relevant packages ####
library(here)
source(here::here('setup','packages.R'))

#### colour schemes etc. ####
# same as https://github.com/lucy-gf/flu_model_LG
source(here::here('setup','aesthetics.R'))

################################################
############## set key parameters ##############
################################################

model_age_groups <- c(0,5,18,65)
age_group_names <- paste0(model_age_groups,"-", c(model_age_groups[2:length(model_age_groups)],99))

start_year_of_analysis <- 2025
years_of_analysis <- 5

simulations <- 100
ageing <- T # are the populations being aged in the simulations?
key_dates <- c('01-04', '01-10') # vaccination and ageing dates (hemisphere-dependent)
vacc_calendar_weeks <- 12 # number of weeks in vaccination program

################################################
################################################
################################################

#### load flu functions ####
source(here::here('functions/flu_parallel.R'))

#### read in test epidemics ####
iso3c_input <- 'GBR'; itz_input <- 'GBR'
hemisphere_input <- 'NH' # North because using GBR data here
epid_dt <- data.table(read_csv(here::here('data','test_epids.csv'), show_col_types=F))

#### choose vaccine variable ####
vaccine_variable <- c('doses','coverage')[2] # using MMGH doses or % coverage?

#### define coverage if using ####
if(vaccine_variable == 'coverage'){
  
  # define percentage coverage intended
  cov_val <- 0.5
  
  # define age groups targeted, e.g. here <10yos and 65+yos
  cov_ages <- c(0:10, 65:101)
    
  # what % coverage in each model age group?
  cov_vec <- coverage_vector(cov_ages, cov_val, model_age_groups)
  
}else{
  
  #### read in test doses ####
  doses <- data.table(read_csv(here::here('data','test_doses.csv'), show_col_types=F))
  
}

ageing_date <<- ifelse(hemisphere_input=='NH', key_dates[1], key_dates[2])
ageing_day <<- as.numeric(substr(ageing_date, 1, 2))
ageing_month <<- as.numeric(substr(ageing_date, 4, 5))
vacc_calendar_start <<- ifelse(hemisphere_input=='NH', key_dates[2], key_dates[1])

infs_rds_list <- mclapply(1:length(vacc_type_list), flu_parallel, mc.cores=length(vacc_type_list))
infs_dt <- rbindlist(infs_rds_list)

#### SAVE OUTPUTS ####

saveRDS(infs_dt, file = here::here('outputs','data','epi',paste0('vacc_',itz_input,'.rds')))
# infs_dt <- readRDS(here::here('outputs','data','epi','vacc_GBR.rds'))

#### PLOT OUTPUTS ####

infs_dt$tot <- rowSums(infs_dt[,2:5])
infs_out <- infs_dt[,c('vacc_type','tot','simulation_index')][, lapply(.SD, cumsum), by=c('vacc_type','simulation_index')]
infs_out$time <- infs_dt$time

infs_out %>% 
  group_by(time,vacc_type) %>% 
  summarise(med = median(tot),
            eti95L = quantile(tot, 0.025),
            eti95U = quantile(tot, 0.975)) %>% 
  ggplot() +
  geom_ribbon(aes(time, ymin=eti95L/1e6, ymax=eti95U/1e6,
                  fill = vacc_type), alpha=0.4) +
  geom_line(aes(time, med/1e6, col=vacc_type),lwd=0.8) +
  theme_bw() + scale_color_manual(values = vtn_colors) +
  scale_fill_manual(values = vtn_colors) +
  ylab('Infections (millions)') +
  facet_wrap(vacc_type~.) + theme(legend.position = 'none')

infs_out %>% 
  group_by(time,vacc_type) %>% 
  summarise(med = median(tot),
            eti95L = quantile(tot, 0.025),
            eti95U = quantile(tot, 0.975)) %>% 
  ggplot() +
  geom_ribbon(aes(time, ymin=eti95L/1e6, ymax=eti95U/1e6,
                  fill = vacc_type), alpha=0.4) +
  geom_line(aes(time, med/1e6, col=vacc_type),lwd=0.8) +
  theme_bw() + scale_color_manual(values = vtn_colors) +
  scale_fill_manual(values = vtn_colors) +
  ylab('Infections (millions)') + theme(legend.position = 'none')



## doses

doses_check <- flu_doses(
    country = 'GBR',
    ageing = T, 
    ageing_date = '01-04',
    epid_inputs = epid_dt[simulation_index==1],
    vaccine_program = 1,
    m_a_g = model_age_groups
    )
ggplot(doses_check) + geom_line(aes(year, vaccs, col = age_grp))






  
  
  