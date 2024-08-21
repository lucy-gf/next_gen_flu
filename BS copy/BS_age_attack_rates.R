
## ANNUAL AGE-SPECIFIC ATTCAK RATES IN EACH COUNTRY
#setwd("~/Desktop/research asst/Global Code")
library(readr)
library(countrycode)

source('BS/BS_data_fcns.R')
source("BS/BS_colors.R")

clusters <- read_csv("data/new_clustering.csv")
annual_country <- read_csv(paste0('data/vacc_output/annual_country.csv'))

annual_cases <- annual_country %>% mutate(age1cases = IU1A + IU1B + IV1A + IV1B,
                                          age2cases = IU2A + IU2B + IV2A + IV2B,
                                          age3cases = IU3A + IU3B + IV3A + IV3B,
                                          age4cases = IU4A + IU4B + IV4A + IV4B) %>% 
  select(scenario, country_code, year, age1cases, age2cases, age3cases, age4cases)

annual_demog <- annual_cases %>% mutate(age1pop = NA, age2pop=NA,
                                        age3pop = NA, age4pop=NA)

annual_country$tot <- rowSums(annual_country[,4:19])
gg <- annual_country %>% group_by(year, scenario) %>% 
  mutate(tot_agg = sum(tot)) %>% ungroup() %>% 
  filter(country_code=='BHS') %>% 
  group_by(scenario) %>% mutate(cum = cumsum(tot_agg)) %>% 
  ungroup() %>% mutate(vt = substr(scenario, 4,4), ct = substr(scenario, 9,9)) %>% 
  filter(!scenario=='no_vacc')

gg %>% ggplot() + 
  geom_line(aes(x=year, y=tot_agg, col=vt), lwd=0.8) + 
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                        ct = supp.labs.age)) +
  scale_color_manual(values = vt_colors) +
  theme_bw() + ylim(c(0,NA))

gg %>% ggplot() + 
  geom_line(aes(x=year, y=cum, col=vt), lwd=0.8) + 
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_color_manual(values = vt_colors) +
  theme_bw() + ylim(c(0,NA))

weekly_agg <- read_csv(paste0('data/vacc_output/weekly_agg.csv')) 
weekly_agg$tot <- rowSums(weekly_agg[,4:19])
gg2 <- weekly_agg %>% group_by(scenario, week) %>% 
  mutate(med = median(tot)) %>% ungroup() %>% filter(simulation_index==1) %>% 
  group_by(scenario) %>% mutate(cum = cumsum(med)) %>% 
  ungroup() %>% mutate(vt = substr(scenario, 4,4), ct = substr(scenario, 9,9)) %>% 
  filter(!scenario=='no_vacc')

gg2 %>% ggplot() + 
  geom_line(aes(x=week, y=med/1000000, col=vt), lwd=0.8) + 
  facet_grid(ct~vt, labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_color_manual(values = vt_colors) +
  theme_bw() + ylim(c(0,NA))


for(cc in unique(annual_demog$country_code)){
  if(cc == 'XKX'){
    country <- 'Kosovo'
    country_altern <- 'Kosovo (under UNSC res. 1244)'
    country_altern_2 <- 'Kosovo'
  }else{
    country <- countrycode(cc, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == cc]
    country_altern_2 <- clusters$country_altern_2[clusters$codes == cc]
  }
  
  annual_demog[annual_demog$country_code == cc, 8:11] <- fcn_weekly_demog(
                   country = c(country, country_altern, country_altern_2),
                   pop_coverage = c(0,0,0,0),
                   weeks_vaccinating = 12,
                   first_year_all = T,
                   NH_vacc_date = '01-10',
                   SH_vacc_date = '01-04',
                   init_vaccinated = c(0,0,0,0),
                   imm_duration = 1, # in years 
                   coverage_pattern = 1,
                   hemisphere = 'NH',
                   start_year = 2025, years = 30) %>% 
    filter(month(week) == 5, day(week) %in% 1:7, U == T) %>% 
    pivot_wider(id_cols=c(year), names_from = age_grp, values_from = total_as) %>% 
    slice(rep(row_number(), 26)) %>% 
    select(!year) 
}

attack_rates <- annual_demog %>% mutate(age1 = age1cases/age1pop,
                                        age2 = age2cases/age2pop,
                                        age3 = age3cases/age3pop,
                                        age4 = age4cases/age4pop) %>% 
  select(scenario, country_code, year, age1, age2, age3, age4) %>% 
  pivot_longer(!c(scenario, country_code, year)) %>% 
  mutate(vt = substr(scenario, 4,4), ct = substr(scenario, 9,9)) %>% 
  mutate(age = substr(name, 4,4))

attack_rates %>% filter(country_code == 'CAN',
                        !scenario == 'no_vacc') %>% ggplot() +
  geom_line(aes(x=year, y=value, group=age, colour=age), lwd=0.8) +
  theme_bw() + 
  xlab('Year') + ylab('Annual attack rate') +
  labs(col = 'Age group') + 
  scale_color_manual(values = age_colors1,
                     labels = c('0-5','5-20','20-65','65+')) +
  theme(text=element_text(size=14)) + ylim(c(0,NA)) +
  facet_grid(ct~vt, labeller = labeller(vt = supp.labs,
                                        ct = supp.labs.age)) 

attack_rates %>% filter(country_code == 'CAN',
                        !scenario == 'no_vacc') %>% ggplot() +
  geom_line(aes(x=year, y=value, group=vt, colour=vt), lwd=0.8) +
  theme_bw() + 
  xlab('Year') + ylab('Annual attack rate') +
  labs(col = 'Vaccine type') + 
  scale_color_manual(values = vt_colors,
                     labels = supp.labs) +
  theme(text=element_text(size=14)) + ylim(c(0,NA)) +
  facet_grid(age~ct, labeller = labeller(vt = supp.labs,
                                        ct = supp.labs.age,
                                        age = supp.labs.agegrps)) 

attack_rates %>% filter(country_code == 'CAN',
                        scenario == 'no_vacc') %>% ggplot() +
  geom_line(aes(x=year, y=value, group=age, colour=age), lwd=0.8) +
  theme_bw() + 
  xlab('Year') + ylab('Annual attack rate') + ylim(c(0,NA)) +
  labs(col = 'Age group') + 
  scale_color_manual(values = age_colors1,
                     labels = c('0-5','5-20','20-65','65+')) +
  theme(text=element_text(size=14)) 


  
  
  
  







