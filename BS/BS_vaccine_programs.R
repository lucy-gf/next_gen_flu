
## SIMULATING X-YEAR TIME SERIES FOR EACH CLUSTER/COUNTRY IN INF A+B ##

## RUNNING WITH VARIOUS VACCINATION SCENARIOS ##
library(dplyr)

# vaccine types:
vacc_type_list <- list(
  current = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    mismatch = T
  ),
  improved_minimal = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_efficacy = list(
    imm_duration = 2,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.7, 0.4),
    mismatch = T
  ),
  improved_breadth = list(
    imm_duration = 3,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.7, 0.46),
    mismatch = F
  ),
  universal = list(
    imm_duration = 5,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.9, 0.7),
    mismatch = F
  ),
  current_breadth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = T
  ),
  improved_minimal_breadth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_efficacy_breadth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.7, 0.4),
    mismatch = T
  ),
  improved_breadth_breadth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.7, 0.46),
    mismatch = F
  ),
  universal_breadth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.9, 0.7),
    mismatch = F
  ),
  current_depth = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = T
  ),
  improved_minimal_depth = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = T
  ),
  improved_efficacy_depth = list(
    imm_duration = 2,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = T
  ),
  improved_breadth_depth = list(
    imm_duration = 3,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = F
  ),
  universal_depth = list(
    imm_duration = 5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), 
    mismatch = F
  )
)

# default inputs for some vaccine program properties:
default_v <- list(
  weeks_vaccinating = 12,
  first_year_all = T,
  NH_vacc_date = '01-10',
  SH_vacc_date = '01-04',
  init_vaccinated = c(0, 0, 0, 0)
)

cov_main <- 0.7
cov_low <- 0.4

coverage_vector <- function(target, cov){
  if(! target %in% 1:5){
    stop("Target scheme must be in 1:5")
  }
  case_when(
    target == 1 ~ c(cov, 0, 0, 0),
    target == 2 ~ c(cov, cov*(6/15), 0, 0),
    target == 3 ~ c(cov, cov*(13/15), 0, 0),
    target == 4 ~ c(0, 0, 0, cov),
    target == 5 ~ c(cov, cov*(13/15), 0, cov)
  )
}

fcn_vacc_prog <- function(target_input, cov_input, relative_vacc_infec_input, none = F){
  if(none == T){
    output_list <- list(
      c(list(pop_coverage = c(0,0,0,0),
             vacc_type = 'current',
             relative_vacc_infec = 1),
        default_v))
    names(output_list) <- paste0('vt_0_ct_0')
    return(output_list)
  }
  output_list <- list(
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'current',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_minimal',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_efficacy',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_breadth',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'universal',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v)
  )
  names(output_list) <- paste0('vt_', 1:5, '_ct_', target_input)
  output_list
}

# vaccine programs:

vaccine_programs_none <- c(
  fcn_vacc_prog(1, 0, 1, T)
)

vaccine_programs_base <- c(
  fcn_vacc_prog(1, cov_main, 1),
  fcn_vacc_prog(2, cov_main, 1),
  fcn_vacc_prog(3, cov_main, 1),
  fcn_vacc_prog(4, cov_main, 1),
  fcn_vacc_prog(5, cov_main, 1)
)

vaccine_programs_low_cov <- c(
  fcn_vacc_prog(1, cov_low, 1),
  fcn_vacc_prog(2, cov_low, 1),
  fcn_vacc_prog(3, cov_low, 1),
  fcn_vacc_prog(4, cov_low, 1),
  fcn_vacc_prog(5, cov_low, 1)
)

vaccine_programs_rel_inf <- c(
  fcn_vacc_prog(1, cov_main, 0.5),
  fcn_vacc_prog(2, cov_main, 0.5),
  fcn_vacc_prog(3, cov_main, 0.5),
  fcn_vacc_prog(4, cov_main, 0.5),
  fcn_vacc_prog(5, cov_main, 0.5)
)

vaccine_programs_breadth <- c(
  fcn_vacc_prog(1, cov_main, 1),
  fcn_vacc_prog(2, cov_main, 1),
  fcn_vacc_prog(3, cov_main, 1),
  fcn_vacc_prog(4, cov_main, 1),
  fcn_vacc_prog(5, cov_main, 1)
)
for(i in 1:length(vaccine_programs_breadth)){
  vaccine_programs_breadth[[i]]$vacc_type <- paste0(vaccine_programs_breadth[[i]]$vacc_type,
                                                    "_breadth")
}

vaccine_programs_depth <- c(
  fcn_vacc_prog(1, cov_main, 1),
  fcn_vacc_prog(2, cov_main, 1),
  fcn_vacc_prog(3, cov_main, 1),
  fcn_vacc_prog(4, cov_main, 1),
  fcn_vacc_prog(5, cov_main, 1)
)
for(i in 1:length(vaccine_programs_depth)){
  vaccine_programs_depth[[i]]$vacc_type <- paste0(vaccine_programs_depth[[i]]$vacc_type,
                                                    "_depth")
}

new_cov_main <- 0.5
new_cov_low <- 0.2

vaccine_programs_BASE50 <- c(
  fcn_vacc_prog(1, new_cov_main, 1),
  fcn_vacc_prog(2, new_cov_main, 1),
  fcn_vacc_prog(3, new_cov_main, 1),
  fcn_vacc_prog(4, new_cov_main, 1),
  fcn_vacc_prog(5, new_cov_main, 1)
)

vaccine_programs_LOW20 <- c(
  fcn_vacc_prog(1, new_cov_low, 1),
  fcn_vacc_prog(2, new_cov_low, 1),
  fcn_vacc_prog(3, new_cov_low, 1),
  fcn_vacc_prog(4, new_cov_low, 1),
  fcn_vacc_prog(5, new_cov_low, 1)
)

vaccine_programs_BASE50_REL_INF <- c(
  fcn_vacc_prog(1, new_cov_main, 0.5),
  fcn_vacc_prog(2, new_cov_main, 0.5),
  fcn_vacc_prog(3, new_cov_main, 0.5),
  fcn_vacc_prog(4, new_cov_main, 0.5),
  fcn_vacc_prog(5, new_cov_main, 0.5)
)

vaccine_programs_BASE50_BREADTH <- c(
  fcn_vacc_prog(1, new_cov_main, 1),
  fcn_vacc_prog(2, new_cov_main, 1),
  fcn_vacc_prog(3, new_cov_main, 1),
  fcn_vacc_prog(4, new_cov_main, 1),
  fcn_vacc_prog(5, new_cov_main, 1)
)
for(i in 1:length(vaccine_programs_BASE50_BREADTH)){
  vaccine_programs_BASE50_BREADTH[[i]]$vacc_type <- paste0(vaccine_programs_BASE50_BREADTH[[i]]$vacc_type,
                                                    "_breadth")
}

vaccine_programs_BASE50_DEPTH <- c(
  fcn_vacc_prog(1, new_cov_main, 1),
  fcn_vacc_prog(2, new_cov_main, 1),
  fcn_vacc_prog(3, new_cov_main, 1),
  fcn_vacc_prog(4, new_cov_main, 1),
  fcn_vacc_prog(5, new_cov_main, 1)
)
for(i in 1:length(vaccine_programs_BASE50_DEPTH)){
  vaccine_programs_BASE50_DEPTH[[i]]$vacc_type <- paste0(vaccine_programs_BASE50_DEPTH[[i]]$vacc_type,
                                                           "_depth")
}

vaccine_programs_merged <- list(
  vaccine_programs_none = vaccine_programs_none,
  vaccine_programs_base = vaccine_programs_base,
  vaccine_programs_low_cov = vaccine_programs_low_cov,
  vaccine_programs_rel_inf = vaccine_programs_rel_inf,
  vaccine_programs_breadth = vaccine_programs_breadth,
  vaccine_programs_depth = vaccine_programs_depth,
  vaccine_programs_BASE50 = vaccine_programs_BASE50,
  vaccine_programs_LOW20 = vaccine_programs_LOW20,
  vaccine_programs_BASE50_REL_INF = vaccine_programs_BASE50_REL_INF,
  vaccine_programs_BASE50_BREADTH = vaccine_programs_BASE50_BREADTH,
  vaccine_programs_BASE50_DEPTH = vaccine_programs_BASE50_DEPTH
)















