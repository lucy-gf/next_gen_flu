#### RUN THE FLU MODEL ####

# key outputs: vaccine_programs, vacc_type_list
source(here::here('next_gen_flu','vacc_types.R'))
# will calculate weekly age- and vaccine-specific population, also loads transmission model
source(here::here('next_gen_flu','functions','demography.R'))
# runs the flu model
source(here::here('next_gen_flu','functions','flu_sim.R'))

#### FUNCTION TO RUN ####
## only input is vaccine type, to parallelise over vt ##
flu_parallel <- function(vaccine_type){
  
  total_start_time <- Sys.time()
  
  dates_many_flu <- seq.Date(last_monday(min(epid_dt$period_start_date)), 
           last_monday(max(epid_dt$end_date)), 
           by=7)
  
  vacc_name <- names(vacc_type_list)[vaccine_type]
  vaccine_used_vec <- doses[vacc_scenario == vacc_name & model_age_group==1]$vacc_used
  
  ## vaccination and ageing
  demography_dt <- fcn_weekly_demog(
    country = iso3c_input,
    ageing,
    ageing_date,
    dates_in = dates_many_flu,
    demographic_start_year = start_year_of_analysis,
    vaccine_used = vaccine_used_vec,
    doses_dt = doses,
    init_vaccinated = c(0,0,0,0),
    model_age_groups
  )
  
  mf_output <- data.table()
  
  for(sim_index in unique(epid_dt$simulation_index)){
    start_time <- Sys.time()
    
    mf_output_si <- many_flu(country = iso3c_input,
                          ageing = T, 
                          ageing_date,
                          epid_inputs = epid_dt[simulation_index==sim_index],  
                          vaccine_used = vaccine_used_vec,
                          doses_dt = doses,
                          model_age_groups,
                          demography_dt
    )
    mf_output_si <- mf_output_si[year(time) >= start_year_of_analysis]
    mf_output_si[, vacc_type := names(vacc_type_list)[vaccine_type]]
    mf_output_si[, simulation_index := sim_index]
    if(nrow(mf_output)==0){
      mf_output <- mf_output_si
    }else{
      mf_output <- rbind(mf_output, mf_output_si)
    }
    
    write.table(paste0(iso3c_input, ', simulation ', sim_index, ', time taken = ', round(Sys.time() - start_time,2),
                 ', number of epids = ', nrow(epid_dt[simulation_index==sim_index]),
                ', total time on country = ', round(Sys.time() - total_start_time,2)),
                file = here::here('output','data','epi',paste0(itz_input),paste0(vaccine_type, '_text.txt')))
  }
  
  mf_output
}







