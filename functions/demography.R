#### NATIONAL DEMOGRAPHY ####

#setwd("~/Desktop/research asst/next_gen_flu")

library(data.table)

model_age_groups <- c(0,5,20,65)

#### POPULATION SIZES ####
pop_hist_WPP_data <- data.table(read_csv('data/pop_hist_WPP_data.csv', show_col_types = F))
pop_proj_WPP_data <- data.table(read_csv('data/pop_proj_WPP_data.csv', show_col_types = F))

fcn_pop_all <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data[name %in% country & Year == year_demog]
  }else{
    data_in <- pop_proj_WPP_data[name %in% country & Year == year_demog]
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in[,!c('name', 'Type', 'Year')])))
  unlist(lapply(pop_in/5, function(x) rep(x,5)))
}
fcn_pop_model <- function(country, year_demog = 2025){
  if(year_demog < 2022){
    data_in <- pop_hist_WPP_data[name %in% country & Year == year_demog]
  }else{
    data_in <- pop_proj_WPP_data[name %in% country & Year == year_demog]
  }
  pop_in <- 1000*as.numeric(unname(unlist(data_in[,!c('name', 'Type', 'Year')])))
  c(pop_in[1], sum(pop_in[2:4]), sum(pop_in[5:13]), sum(pop_in[14:21]))
}

#### AGEING ####
LT_WPP_data <- data.table(read_csv('data/LT_WPP_data.csv', show_col_types = F))
CBR_WPP_data <- data.table(read_csv('data/CBR_WPP_data.csv', show_col_types = F))

demog_data <- data.frame(country = unique(CBR_WPP_data$Name),
                         CBR = NA, M1 = NA, M2 = NA, M3 = NA, M4 = NA)

# assuming rectangular age distributions within each age group so equal
# weighting on each nqx in the group:
for(i in 1:nrow(demog_data)){
  country_input <- LT_WPP_data[Name == demog_data$country[i]]
  demog_data$CBR[i] <- (CBR_WPP_data[Name == demog_data$country[i]])$CBR_indiv
  demog_data$M1[i] <- (country_input[country_input$x == 0,]$dying_prob_annual +
                         4*country_input[country_input$x == 1,]$dying_prob_annual)/5
  demog_data$M2[i] <- mean(country_input[country_input$x %in% 5:19,]$dying_prob_annual)
  demog_data$M3[i] <- mean(country_input[country_input$x %in% 20:60,]$dying_prob_annual)
  demog_data$M4[i] <- 1/(country_input[x == 65])$LEx
}

#### CONTACT MATRICES ####
load("data/contact_all.rdata")

fcn_contact_matrix <- function(country_name, country_name_altern,
                               country_name_altern_2, pop_model){
  model_age_groups <<- data.frame(agegroup_name=age_group_names, duration=diff(c(age_limits,120)),
                                  wpp_agegroup_low=c(1,2,5,14), wpp_agegroup_high=c(1,4,13,16),
                                  popul=pop_model)
  # age groups corresponding to Prem et al. 2021 matrices
  standard_age_groups <- fun_cntr_agestr(i_cntr = c(country_name, country_name_altern, country_name_altern_2),
                                         i_year="2020",
                                         age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120))
  # modify contact matrix to correspond to our age groups
  if(country_name=='Kosovo'){
    sel_cntr_code <- 'XKX'
  }else{
    sel_cntr_code <- countrycode(country_name, origin = 'country.name', destination = 'iso3c')
  }
  # using 'similar' countries for those not in contact_all from Prem et al.:
  if(sel_cntr_code %in% setdiff(clusters$codes, names(contact_all))){
    sel_cntr_code <- case_when(sel_cntr_code == 'AUS' ~ 'NZL',
                               sel_cntr_code == 'SOM' ~ 'ETH',
                               sel_cntr_code == 'LBN' ~ 'SYR',
                               sel_cntr_code == 'JPN' ~ 'KOR',
                               sel_cntr_code == 'TWN' ~ 'CHN',
                               sel_cntr_code == 'XKX' ~ 'SRB',
                               sel_cntr_code == 'NCL' ~ 'VUT',
                               sel_cntr_code == 'GUF' ~ 'SUR',
                               sel_cntr_code == 'HTI' ~ 'DOM',)
    standard_age_groups <- case_when(sel_cntr_code == 'NZL' ~ fun_cntr_agestr(i_cntr = 'New Zealand',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'ETH' ~ fun_cntr_agestr(i_cntr = 'Ethiopia',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SYR' ~ fun_cntr_agestr(i_cntr = 'Syrian Arab Republic',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'KOR' ~ fun_cntr_agestr(i_cntr = 'Republic of Korea',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'CHN' ~ fun_cntr_agestr(i_cntr = 'China',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SRB' ~ fun_cntr_agestr(i_cntr = 'Serbia',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'VUT' ~ fun_cntr_agestr(i_cntr = 'Vanuatu',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'SUR' ~ fun_cntr_agestr(i_cntr = 'Suriname',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),
                                     sel_cntr_code == 'DOM' ~ fun_cntr_agestr(i_cntr = 'Dominican Republic',
                                                                              i_year="2020",
                                                                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120)),)
  }
  C_m_merged_nonrecipr <- fun_create_red_C_m(C_m_full=contact_all[[sel_cntr_code]],
                                             model_agegroups=model_age_groups,
                                             orig_age_groups_duration=standard_age_groups$duration,
                                             orig_age_groups_sizes=standard_age_groups$values)
  # make it reciprocal for the larger group
  contact_matrix_output <- fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                                               age_group_sizes = model_age_groups$popul)
  return(contact_matrix_output)
}

## function to scale to match a given r0?



















