#### RUN THE FLU MODEL ####

# key outputs: vaccine_programs, vacc_type_list
source(here::here('functions','vacc_types.R'))
# will calculate weekly age- and vaccine-specific population, also loads transmission model
source(here::here('functions','demography.R'))
# runs the flu model
source(here::here('functions','flu_sim.R'))

#### FUNCTION TO RUN ####
## only input is vaccine type, to parallelise over vt ##
flu_parallel <- function(vaccine_type){
  
  # total_start_time <- Sys.time()
  
  dates_many_flu <- seq.Date(last_monday(min(epid_dt$period_start_date)), 
           last_monday(max(epid_dt$end_date)), 
           by=7)
  
  vacc_name <- names(vacc_type_list)[vaccine_type]
  
  vaccine_used_vec <- if(vaccine_variable == 'doses'){
    # if using doses, NGIVs are often introduced years after the start of the epidemic period
    doses[vacc_scenario == vacc_name & model_age_group==1]$vacc_used
  }else{
    # if using coverage, assuming NGIVs available each year (adjust here if not!)
    rep(vacc_name, (year(epid_dt$end_date[1]) - year(epid_dt$period_start_date[1])))
  }
  
  ## vaccination and ageing
  demography_dt <- fcn_weekly_demog(
    country = iso3c_input,
    ageing,
    ageing_date,
    dates_in = dates_many_flu,
    demographic_start_year = start_year_of_analysis,
    vaccine_used = vaccine_used_vec,
    vaccine_var = vaccine_variable,
    doses_dt = if(vaccine_variable == 'doses'){doses}else{NULL},
    vacc_cov_vec = if(vaccine_variable == 'coverage'){cov_vec}else{NULL},
    init_vaccinated = c(0,0,0,0),
    model_age_groups
  )
  # demography_dt %>% mutate(prop = value/total_as) %>% filter(V==T) %>% ggplot() + geom_line(aes(week, prop, col=age_grp))
  if(min(demography_dt$value) < 0){ # quick fix if any vaccination issues (there shouldn't be)
    print(paste0('Negative values in demography_dt, iso3c = ', iso3c_input,', vaccine type = ', vaccine_type))
  }
  
  mf_output <- data.table()
  
  # loop over 1:100 simulations
  for(sim_index in unique(epid_dt$simulation_index)){
    start_time <- Sys.time()
    
    # run flu simulations
    mf_output_si <- many_flu(country = iso3c_input,
                             ageing, 
                             ageing_date,
                             epid_inputs = epid_dt[simulation_index==sim_index],  
                             vaccine_used = vaccine_used_vec,
                             vaccine_var = vaccine_variable,
                             doses_dt = if(vaccine_variable == 'doses'){doses}else{NULL},
                             vacc_cov_vec = if(vaccine_variable == 'coverage'){cov_vec}else{NULL},
                             model_age_groups,
                             demography_dt
    )
    mf_output_si <- mf_output_si[year(time) >= start_year_of_analysis] # in case epidemic started pre-2025 
    mf_output_si[, vacc_type := names(vacc_type_list)[vaccine_type]] # add vaccine name
    mf_output_si[, simulation_index := sim_index] # add simulation number
    
    # printing if there is an NA error (shouldn't happen)
    if(is.na(sum(rowSums(mf_output_si %>% select(starts_with('I')))))){
      print(paste0('vt = ', vaccine_type, ', sim_index = ', sim_index, ' - is.na'))
    }
    
    # merge output 
    if(nrow(mf_output)==0){
      mf_output <- mf_output_si
    }else{
      mf_output <- rbind(mf_output, mf_output_si)
    }
    
    # if(!file.exists(here::here('output','data','epi',paste0(itz_input,'_text')))){
    #   dir.create(file.path(here::here('output','data','epi',paste0(itz_input,'_text'))))
    # }
    
    # print txt file to keep track of simulations
    # if(sim_index == 1 | sim_index %% 10 == 0){
    #   writeLines(paste0(iso3c_input, ', simulation ', sim_index, ', time taken = ', round(Sys.time() - start_time,2),
    #                     ', number of epids = ', nrow(epid_dt[simulation_index==sim_index]),
    #                     ', total time on country = ', round(Sys.time() - total_start_time,2)),
    #              paste0('output/data/epi/',paste0(itz_input,'_text'),'/',paste0(vaccine_type, '_text.txt')))  
    # }
    
  }
  
  mf_output
}







