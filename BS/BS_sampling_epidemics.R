
## SELECTING 30 YEARS, MATCHING VALUES, THE EPIDEMICS, AND THEIR PARAMETERS, 100 TIMES

source("BS/BS_data_fcns.R")
source("BS/BS_post_samples.R")

length_simulation <- 30
n_simulations <- 100

# matches <- read.csv("data_for_cluster/matches.csv") # saved from vaccine_inputs.R
# match_probability <- data.frame(`N_A` = NA,`N_B` = NA,`S_A` = NA,`S_B` = NA)
# for(hemisphere in c('N','S')){
#   for(strain in c('A','B')){
#     mini_df <- matches[grepl(hemisphere, matches$hemisphere) &
#                                grepl(strain, matches$strain_match),]
#     match_probability[paste0(hemisphere,'_',strain)] <-
#       sum(mini_df$match)/nrow(mini_df)
#   }
# }
#
# sampled_years_df <- data.frame() # will be 30*100=3000 years
# 
# sampled_years_vec <- sample(2010:2019, size = length_simulation*n_simulations, replace = T)
# sampled_years_df <- data.frame(index = rep(1:n_simulations, each=length_simulation),
#                                simulation_year = rep(1:30, n_simulations),
#                                years = sampled_years_vec, 
#                                N_A_match = sample(c(T,F), size = length_simulation*n_simulations, 
#                                                replace = T, prob=c(match_probability$N_A, 1-match_probability$N_A)),
#                                N_B_match = sample(c(T,F), size = length_simulation*n_simulations, 
#                                                replace = T, prob=c(match_probability$N_A, 1-match_probability$N_B)),
#                                 S_A_match = sample(c(T,F), size = length_simulation*n_simulations, 
#                                                replace = T, prob=c(match_probability$N_A, 1-match_probability$S_A)),
#                                S_B_match = sample(c(T,F), size = length_simulation*n_simulations, 
#                                                replace = T, prob=c(match_probability$N_A, 1-match_probability$S_B)))
# write_csv(sampled_years_df, file = paste0("data_for_BS/sampled_years_", length_simulation,
#                                           "_", n_simulations, ".csv"))

sampled_years_df <- read_csv("data_for_BS/sampled_years_30_100.csv")

for(cntr_num in 1:7){
## "ARG" "AUS" "CAN" "CHN" "GBR" "GHA" "TUR"
country_code_index <- unique((clusters %>% arrange(cluster_code))$cluster_code)[cntr_num]
hemisphere_input <- c('SH', 'SH', 'NH', 'NH', 'NH', 'NH', 'NH')[cntr_num]

sampled_epidemics <- data.frame() # will be around 45,000 epidemics
print("country_code_index") 
for(index in (1:100)){
    epids_to_add <- fcn_sample_epidemics(country_code = country_code_index,
                                         sampled_years = sampled_years_df[sampled_years_df$index == index,]) %>%
      mutate(simulation_index = index, pushback = NA, init_ageing_date = NA, init_nye = NA, pushback_T = F)
    code <- epids_to_add$exemplar_code[1]
    country_name <- countrycode(code, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == code]
    contact_matrix <- as.matrix(unname(list_contact_matr[[country_name]])) 
    yr_res_pop <- unlist(lapply(pop_age(wpp_age(country_name, 2015))$population/5, function(x) rep(x,5)))
    group_pop <- pop_age(wpp_age(country_name, 2015), age.limits=c(0,5,20,65))$population # size of each age group
    contact_matrix_small <- t(t(contact_matrix)/group_pop)
    for(i in 1:nrow(epids_to_add)){
      original_date <- as.Date(paste0(as.numeric(epids_to_add$day[i]), '-', epids_to_add$month[i], '-',
                                      epids_to_add$year[i]), '%d-%m-%Y')
      pb_input <- fcn_identify_start(sus_input = unlist(epids_to_add$sus[i]),
                                                     trans_input = unlist(epids_to_add$trans[i]),
                                                     infected = unlist(epids_to_add$infected[i]),
                                                     date = original_date,
                                                     country = country_name,
                                                     strain_input = epids_to_add$strain[i],
                                                     yr_res_pop = yr_res_pop,
                                                     contact_matrix_small = contact_matrix_small,
                                                     ageing_dates = c('01-04','01-10'),
                                                     hemisphere = hemisphere_input,
                                                     simulation_year = unlist(epids_to_add$simulation_cal_year[i]))
      if(length(pb_input) == 1){
        epids_to_add$pushback[i] <- pb_input
      }else{
        epids_to_add$pushback_T[i] <- F
        if(pb_input[2] == 'ageing_date'){
        epids_to_add$init_ageing_date[i] <- pb_input[1]
      }else{
        epids_to_add$init_nye[i] <- pb_input[1]
      }}
    }
    sampled_epidemics <- rbind(sampled_epidemics,epids_to_add)
    print(paste0('# simulation = ', index))
    write_csv(sampled_epidemics, file = paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                             "_", n_simulations,'_',country_code_index,".csv"))
}


#sampled_epidemics_input <- read.csv("data/sampled_epidemics_30_100.csv")

}





