#### TRANSMISSION MODEL ####

incidence_VS <- function(
    demography_input,
    susceptibility,
    transmissibility,
    initial_infected,
    calendar_input,
    contacts,
    waning_rate,
    vaccination_ratio_input,
    begin_date, 
    end_date, 
    age_groups_model,
    kappa_input = 1
){
  
  risk_ratios_input <- matrix(c(rep(0,8)), ncol = 4 , byrow = T) # Not using risk groups so setting this here for now
  age_group_sizes <- stratify_by_age(demography_input, age_groups_model)
  population_stratified <- stratify_by_risk(age_group_sizes, risk_ratios_input)
  
  # define model timings
  interval = 7
  t <- as.numeric(seq(begin_date, end_date, interval))
  # define age group inputs
  no_groups <- length(population_stratified)
  
  initial_infected <- rep(10^parameters[4],length(age_group_sizes)) #6)
  initial_infected <- stratify_by_risk(initial_infected, risk_ratios_input)
  susceptibility <- c((0.2*1 + 0.8*parameters[3]), rep(parameters[3],3)) # 1 in [0,1), parameters[3] in [1,\infty)
  transmissibility = parameters[2]
  
  initial_vaccinated_prop = unlist(vaccination_ratio_input[[1]])
  initial_Rv_prop = unlist(vaccination_ratio_input[[2]])
  initial_R_prop = unlist(vaccination_ratio_input[[3]])
  #Assume that all R become susceptible again at the start of each posterior
  
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover risk groups
  new_cij <- matrix(rep(0,12*12), nrow = 12)
  for (k in 1:3) {
    for (l in 1:3) {
      lk <- (k - 1)*4 + 1
      ll <- (l - 1)*4 + 1
      new_cij[lk:(lk + 4 - 1), ll:(ll + 4 - 1)] <- contacts
    }
  }
  
  if(sum((population_stratified - initial_infected) < 0) > 0){
    parameters[4] <- min(log10(population_stratified[1:4]))
    initial_infected <- rep(10^parameters[4],length(age_group_sizes)) #6)
    initial_infected <- stratify_by_risk(initial_infected, risk_ratios_input)
    initial_infected <- initial_infected*(1 - initial_vaccinated_prop) - c(rep(1,4),rep(0,8))
  }else{
    initial_infected <- initial_infected*(1 - initial_vaccinated_prop) - c(rep(1,4),rep(0,8))
  }
  
  # specify the model
  mod <- gen_seeiir_ag_vacc_waning_updated$new(
    no_groups = no_groups,
    cij = new_cij,
    trans = transmissibility,
    pop = population_stratified,
    I0 = initial_infected,
    V0 = initial_vaccinated_prop,
    R0 = initial_R_prop,
    RV0 = initial_Rv_prop,
    susc = rep(susceptibility,3),
    alpha = calendar_input$efficacy[1:no_groups],
    omega = waning_rate,
    dates = calendar_input$dates,
    calendar = matrix(calendar_input$calendar, ncol = 4*3),
    gamma1 = 2/infection_delays[1],
    gamma2 = 2/infection_delays[2], 
    num_vac_start = rep(0,no_groups), 
    kappa = kappa_input
  )
  
  # run the model
  y_run <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  # calculate the cumulative values
  y <- mod$transform_variables(y_run)$cumI
  yU <- mod$transform_variables(y_run)$cumIU
  yV <- mod$transform_variables(y_run)$cumIV
  # Returning the differences in cumulative infections from one time-step to the other
  y_transform <- function(y_input){
    y_input <- data.table(y_input[2:(nrow(y_input)), ] - y_input[1:(nrow(y_input) - 1), ])
    if(ncol(y_input) == 1){y_input <- data.table(unname(t(y_input)))}
    y_input$time <- seq(begin_date,length.out = nrow(y_input), by = interval)
    y_input <- data.table(y_input)
  }
  y <- y_transform(y)
  yU <- y_transform(yU)
  yV <- y_transform(yV)
  
  output_y <- cbind(time = y$time, y[,V1:V4], yU[,V1:V4], yV[,V1:V4])
  colnames(output_y) <- c('time', 'I1', 'I2', 'I3', 'I4',
                          'IU1', 'IU2', 'IU3', 'IU4', 'IV1', 'IV2', 'IV3', 'IV4')
  
  return(output_y)
}


incidence_VS(
    demography_input,
    parameters,
    calendar_input,
    contacts,
    waning_rate,
    vaccination_ratio_input,
    begin_date, 
    end_date, 
    age_groups_model,
    kappa_input = 1
)



library(odin)

flu_odin <- odin::odin({
  # Number of groups
  no_groups <- user()
  
  # INITIAL CONDITIONS
  # Population size by age/risk group
  pop[] <- user()
  # Start vaccinated by age/risk group
  V0[] <- user()
  R0[] <- user()
  RV0[] <- user()
  # Initial infection by age/risk group
  I0[] <- user()
  num_vac_start[] <- user()
  
  # MODEL PARAMETERS
  # Susceptibility
  susc[] <- user()
  
  # Transmissibility
  trans <- user()
  
  # Latent periods
  gamma1 <- user()
  gamma2 <- user()
  
  # Vaccine related variables 
  dates[] <- user()
  calendar[,] <- user()
  # waning of vaccine immunity
  omega <- user()
  
  # kappa in [0,1] the relative (decreased) infectiousness of vaccinated indivs
  kappa <- user()
  
  # efficacy
  alpha[] <- user()
  
  # Contact matrix
  cij[,] <- user()
  
  # Force of infection
  lambda[] <- trans * susc[i] * (sum(sij[i,]))
  
  # Vaccination. The rate is a step function that changes at each date according
  # to the passed calendar
  vI[] <- interpolate(dates, calendar, "constant")
  
  sumN[] <- if (vI[i]>0) (S[i]+E1[i]+E2[i]+I1[i]+I2[i]+R[i]) else 0
  v[] <- if (sumN[i]>0) (1 - ((1 - vI[i])/(sumN[i]/pop[i]))) else 0
  
  # Transmission matrix
  sij[,] <- cij[i,j] * (I1[j] + I2[j] + kappa*(I1v[j] + I2v[j]))
  
  # Newly infected
  newInf[] <- lambda[i] * S[i]
  newInfv[] <- lambda[i] * Sv[i]
  
  # THE DERIVATIVES OF THE SEEIIR MODEL
  # Derivatives of the not vaccinated group
  deriv(S[])  <- + omega*(Sv[i] + Rev[i]) - newInf[i] - v[i] * S[i] 
  deriv(E1[]) <- + omega*E1v[i] + newInf[i] - gamma1 * E1[i] - v[i] * E1[i] 
  deriv(E2[]) <- + omega*E2v[i] + gamma1 * (E1[i] - E2[i]) - v[i] * E2[i]
  deriv(I1[]) <- + omega*I1v[i] + gamma1 * E2[i]  - gamma2 * I1[i] - v[i] * I1[i] 
  deriv(I2[]) <- + omega*I2v[i] + gamma2 * (I1[i] - I2[i]) - v[i] * I2[i] 
  deriv(R[])  <- + omega*Rnv[i] + gamma2 * I2[i] - v[i] * R[i] 
  
  # Derivatives vaccination group
  deriv(Sv[])  <- - omega*Sv[i]  - newInfv[i] + v[i] * (1-alpha[i]) * S[i] 
  deriv(E1v[]) <- - omega*E1v[i] + newInfv[i] - gamma1 * E1v[i] + v[i] * E1[i]
  deriv(E2v[]) <- - omega*E2v[i] + gamma1 * (E1v[i] - E2v[i]) + v[i] * E2[i]
  deriv(I1v[]) <- - omega*I1v[i] + gamma1 * E2v[i]  - gamma2 * I1v[i] + v[i] * I1[i]
  deriv(I2v[]) <- - omega*I2v[i] + gamma2 * (I1v[i] - I2v[i]) + v[i] * I2[i]
  deriv(Rv[])  <- - omega*Rv[i]  + gamma2 * I2v[i] + v[i] * (R[i] + alpha[i] * S[i])
  deriv(Rnv[])  <- - omega*Rnv[i]  + gamma2 * I2v[i] + v[i] * (R[i])
  deriv(Rev[])  <- - omega*Rev[i] + v[i] * (alpha[i] * S[i])
  deriv(VT[]) <- + v[i]*(S[i] + E1[i] + E2[i] + I1[i] + I2[i] + R[i])
  
  # Tracking the cumulative amount of infections over time for output of incidence
  deriv(cumI[]) <- newInf[i] + newInfv[i]
  deriv(cumIU[]) <- newInf[i] 
  deriv(cumIV[]) <- newInfv[i]
  
  # Initial value of the variables
  initial(S[1:no_groups]) <- pop[i]*(1-V0[i])*(1-R0[i]) - I0[i]
  initial(E1[1:no_groups]) <- 0
  initial(E2[1:no_groups]) <- 0
  initial(I1[1:no_groups]) <- I0[i]
  initial(I2[1:no_groups]) <- 0
  initial(R[1:no_groups]) <- pop[i]*(1-V0[i])*(R0[i])
  initial(cumI[1:no_groups]) <- 0
  initial(cumIU[1:no_groups]) <- 0
  initial(cumIV[1:no_groups]) <- 0
  
  initial(Sv[1:no_groups]) <- (pop[i]*V0[i])*(1-RV0[i])
  initial(E1v[1:no_groups]) <- 0
  initial(E2v[1:no_groups]) <- 0
  initial(I1v[1:no_groups]) <- 0
  initial(I2v[1:no_groups]) <- 0
  initial(Rv[1:no_groups]) <- (pop[i]*V0[i])*(RV0[i])
  initial(Rnv[1:no_groups]) <- 0
  initial(Rev[1:no_groups]) <- (pop[i]*V0[i])*(RV0[i])
  initial(VT[1:no_groups]) <- num_vac_start[i]
  
  # Set dimension of all variables/parameters
  dim(dates) <- user()
  dim(calendar) <- user()
  
  dim(pop) <- no_groups
  dim(I0) <- no_groups
  dim(V0) <- no_groups
  dim(R0) <- no_groups
  dim(RV0) <- no_groups
  dim(num_vac_start) <- no_groups
  dim(susc) <- no_groups
  dim(lambda) <- no_groups
  dim(v) <- no_groups
  dim(vI) <- no_groups
  dim(sumN) <- no_groups
  dim(alpha) <- no_groups
  dim(cij) <- c(no_groups, no_groups)
  dim(sij) <- c(no_groups, no_groups)
  
  dim(S) <- no_groups
  dim(E1) <- no_groups
  dim(E2) <- no_groups
  dim(I1) <- no_groups
  dim(I2) <- no_groups
  dim(R) <- no_groups
  dim(Sv) <- no_groups
  dim(E1v) <- no_groups
  dim(E2v) <- no_groups
  dim(I1v) <- no_groups
  dim(I2v) <- no_groups
  dim(Rv) <- no_groups
  dim(Rnv) <- no_groups
  dim(Rev) <- no_groups
  dim(cumI) <- no_groups
  dim(cumIU) <- no_groups
  dim(cumIV) <- no_groups
  dim(newInf) <- no_groups
  dim(newInfv) <- no_groups
  dim(VT) <- no_groups
})
























