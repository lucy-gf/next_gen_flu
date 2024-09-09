#### VACCINE TYPES ####

# vaccine types:
vacc_type_list <- list(
  current = list(
    imm_duration = 0.5, # immunity duration
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    rel_inf = 1 # relative infectiousness of vaccinated individuals
  ),
  improved_minimal = list(
    imm_duration = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    rel_inf = 1
  ),
  improved_efficacy = list(
    imm_duration = 2,
    VE = c(0.9, 0.7, 0.7, 0.4),
    rel_inf = 1
  ),
  improved_breadth = list(
    imm_duration = 3,
    VE = c(0.7, 0.46, 0.7, 0.46),
    rel_inf = 1
  ),
  universal = list(
    imm_duration = 5,
    VE = c(0.9, 0.7, 0.9, 0.7),
    rel_inf = 1
  )
)

# targeted populations
coverage_vector <- function(
    target, # vector of targeted ages
    cov, # scalar
    age_groups = c(5,20,65) # model age group cutoffs
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

fcn_vacc_prog <- function(target_ages, cov_input, none = F){
  if(none == T){
    output_list <- list(
      c(list(pop_coverage = c(0,0,0,0),
             vacc_type = 'current')))
    names(output_list) <- paste0('no_vacc')
    return(output_list)
  }
  output_list <- list(
    c(list(pop_coverage = coverage_vector(target_ages, cov_input),
           vacc_type = 'current')),
    c(list(pop_coverage = coverage_vector(target_ages, cov_input),
           vacc_type = 'improved_minimal')),
    c(list(pop_coverage = coverage_vector(target_ages, cov_input),
           vacc_type = 'improved_efficacy')),
    c(list(pop_coverage = coverage_vector(target_ages, cov_input),
           vacc_type = 'improved_breadth')),
    c(list(pop_coverage = coverage_vector(target_ages, cov_input),
           vacc_type = 'universal'))
  )
  names(output_list) <- c('current','improved_minimal','improved_efficacy','improved_breadth','universal')
  output_list 
}






















