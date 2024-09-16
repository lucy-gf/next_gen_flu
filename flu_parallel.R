#### RUN THE FLU MODEL ####

# setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)
library(fluEvidenceSynthesis)
library(parallel)
library(patchwork)

# key outputs: vaccine_programs, vacc_type_list
source('vacc_types.R') 
# will calculate weekly age- and vaccine-specific population, also loads transmission model
source('functions/demography.R') 
# runs the flu model
source('functions/flu_sim.R') 

#### VACCINE PROGRAM DETAILS ####
cov_main <- 0.5 # coverage level in targeted age groups
age_targeting <- 0:17 # ages being targeted
vacc_calendar_start <- '01-10'
vacc_calendar_weeks <- 12

vaccine_programs <- c(
  fcn_vacc_prog(NA, 0,
                vacc_calendar_start, 
                vacc_calendar_weeks,
                T),
  fcn_vacc_prog(age_targeting, cov_main,
                vacc_calendar_start, 
                vacc_calendar_weeks)
)

## MODEL AGE GROUPS - fixed ##
model_age_groups <- c(0,5,20,65)
age_group_names <- paste0(model_age_groups,"-", c(model_age_groups[2:length(model_age_groups)],99))
# could generalise across more/different age groups, not done though.

#### EPIDEMIC DATA ####
epid_starts <- c(as.Date(paste0('01-11-', 2025:2029), format='%d-%m-%Y'))
epid_dt <- data.table(
  susceptibility = c(0.6,0.59,0.7,0.65,0.61),
  transmissibility = c(0.09,0.08,0.075,0.075,0.08),
  match = c(T,F,F,T,F),
  initial_infected = c(list(c(100,100,100,100)),list(c(100,100,100,100)),list(c(100,100,100,100)),
                       list(c(100,100,100,100)),list(c(100,100,100,100))),
  period_start_date = as.Date('01-01-2025',format='%d-%m-%Y'),
  epid_start_date = epid_starts,
  end_date = as.Date('01-01-2031',format='%d-%m-%Y'),
  vacc_calendar_start = '01-10',
  vacc_calendar_weeks = 12 # annual
)

ageing_date <- '01-04'
iso3c <- 'GBR'

#### FUNCTION TO RUN ####
## only input is vaccine type, to parallelise over vt ##
flu_parallel <- function(vaccine_type){

  mf_output <- many_flu(country = iso3c,
                        ageing = T, 
                        ageing_date,
                        epid_inputs = epid_dt,  
                        vaccine_program = vaccine_programs[[vaccine_type]],
                        model_age_groups
                        )
  
  mf_output[, vacc_type := names(vaccine_programs)[vaccine_type]]

}

#### RUN OUTPUTS #### 
## (parallelised across all vaccine types) ##

infs_rds_list <- mclapply(1:length(vaccine_programs), flu_parallel, mc.cores=5)

#### SAVE OUTPUTS ####

saveRDS(infs_rds_list, file = paste0("outputs/vacc_", iso3c,".rds"))

#### PLOTS ####

infs_out <- data.table()
for(i in 1:length(infs_rds_list)){
  add_on <- infs_rds_list[[i]]
  add_on <- add_on[, c('time','vacc_type','I1','I2','I3','I4')]
  add_on[, tot := I1 + I2 + I3 + I4]
  infs_out <- rbind(infs_out, add_on)
}

infs_cum <- infs_out[, c('vacc_type','tot')][, lapply(.SD, cumsum), by=c('vacc_type')]
infs_cum[, time:=infs_out$time]

# colors
vt_colors <- c('no_vacc' = '#000000', 'current' = '#d91818', 'improved_minimal' = '#e2790e', 
               'improved_efficacy' = '#eacb2c', 'improved_breadth' = '#62adc1', 'universal' = '#324da0')

incidence <- ggplot(infs_out) + 
  geom_line(aes(x=time, y=tot/1000000, col=vacc_type), lwd=0.8) +
  ylab('Infections (millions)') + xlab('') + 
  scale_color_manual(values = vt_colors) + 
  labs(col='Vaccine type') +
  theme_minimal() + theme(text = element_text(size = 14))

cumulative <- ggplot(infs_cum) + 
  geom_line(aes(x=time, y=tot/1000000, col=vacc_type), lwd=0.8) +
  ylab('Cumulative infections (millions)') + xlab('Time') + 
  labs(col='Vaccine type') + 
  scale_color_manual(values = vt_colors) + 
  theme_minimal() + theme(text = element_text(size = 14))
  
incidence + cumulative + plot_layout(guides='collect',
                                     nrow=2)

# side_by_side <- dcast(infs_cum, time ~ vacc_type, value.var = 'tot'); View(side_by_side)












