#### VACCINE TYPES ####

# vaccine types:
vacc_type_list <- list(
  `0` = list(
    imm_duration = 0.5, # immunity duration
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    rel_inf = 1 # relative infectiousness of vaccinated individuals
  ),
  A.1 = list(
    imm_duration = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    rel_inf = 1
  ),
  A.2 = list(
    imm_duration = 0.5,
    VE = c(0.9, 0.7, 0.42, 0.28),
    rel_inf = 1
  ),
  B.1 = list(
    imm_duration = 2,
    VE = c(0.9, 0.7, 0.7, 0.4),
    rel_inf = 1
  ),
  B.2 = list(
    imm_duration = 3,
    VE = c(0.7, 0.46, 0.7, 0.46),
    rel_inf = 1
  ),
  C = list(
    imm_duration = 5,
    VE = c(0.9, 0.7, 0.9, 0.7),
    rel_inf = 1
  )
)

## function to map targeted ages to proportions of each model age group
coverage_vector <- function(
    target, # vector of targeted ages
    cov, # scalar
    age_groups # model age group cutoffs
    ){
  # adding 0 if not in
  if(! 0 %in% age_groups){
    age_groups <- c(0,age_groups)
  }
  
  coverage_vec <- c()
  
  for(group in 1:length(age_groups)){
    ages <- unlist(ifelse(group<length(age_groups), list(age_groups[group]:(age_groups[group+1]-1)), list(age_groups[group]:101)))
    prop_coverage <- length(intersect(ages,target))/length(ages) # proportion of age group targeted
    coverage_vec <- c(coverage_vec, prop_coverage*cov)
  }
  
  coverage_vec
}

## function to map age targeting, coverage input, vaccination program details to a list
fcn_vacc_prog <- function(target_ages, cov_input, 
                          vacc_calendar_start,
                          vacc_calendar_weeks, 
                          none = F){
  if(none == T){
    output_list <- list(
      c(list(pop_coverage = c(0,0,0,0),
             vacc_type = 'current',
             start = vacc_calendar_start,
             weeks = vacc_calendar_weeks)))
    names(output_list) <- paste0('no_vacc')
    return(output_list)
  }
  
  output_list <- list()
  
  for(vt in names(vacc_type_list)){
    output_list <- c(output_list,
                        list(list(pop_coverage = coverage_vector(target_ages, 
                                                                 cov_input,
                                                                 model_age_groups),
                               vacc_type = vt,
                               start = vacc_calendar_start,
                               weeks = vacc_calendar_weeks
                               )
                          )
                        )
  }
  
  names(output_list) <- names(vacc_type_list)
  output_list 
}

make_cov_matrix <- function(cov_input){
  data.table(cov_input %>% select(year,age_grp,cov) %>% 
    pivot_wider(id_cols=year,names_from=age_grp,names_prefix='age_grp_',values_from = cov))
}

fcn_vacc_prog_mmgh <- function(cov_input, # matrix
                               vacc_calendar_start,
                               vacc_calendar_weeks
                               ){
  
  cov1_matrix <- make_cov_matrix(cov_input[vaccine==1])
  cov2_matrix <- make_cov_matrix(cov_input[vaccine==2])
  cov0_matrix <- cov1_matrix
  cov0_matrix[, c('age_grp_1','age_grp_2','age_grp_3','age_grp_4') := 0]
  
  output_list <- list(
      c(list(pop_coverage = cov0_matrix,
             vacc_type = 'current',
             start = vacc_calendar_start,
             weeks = vacc_calendar_weeks)))
  
  for(vt in names(vacc_type_list)){
    output_list <- c(output_list,
                     list(list(pop_coverage = case_when(vt=='current' ~ cov1_matrix,
                                                        !vt=='current' ~ cov2_matrix),
                               vacc_type = vt,
                               start = vacc_calendar_start,
                               weeks = vacc_calendar_weeks
                     )
                     )
    )
  }
  
  names(output_list) <- c('no_vacc', names(vacc_type_list))
  output_list 
}




















