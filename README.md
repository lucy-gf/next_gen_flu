# next_gen_flu

Reproducible flexible next-generation vaccine flu model. Takes some vaccine programs and epidemic data, and outputs weekly vaccination-status and age-specific influenza infections.

## **vacc_types.R** 
Contains functions to define vaccine programs for no vaccinations and 5 NGIVs (can be edited)
# (inputs: coverage level, targeted ages)
# (outputs: VE, mean immunity length, coverage across model age groups)

# **functions/demography.R** calculates weekly age- and vaccination-status specific population over
# the relevant time period (can remove ageing if needed), for a given country and vaccine program,
# and contains the function to calculate contact matrices from a given country and age-specific population

# **functions/transmission_model.R** contains the ODE model builder, 
# epidemic simulation function, and a function to calculate vaccination status-specific demography

# **functions/flu_sim.R** contains:
# one_flu(), which runs an epidemic,
# many_flu(), which takes a data-table of epidemic data and combines many epidemics, 
# dfn_vaccine_calendar(), which converts the vaccine program and epidemic dates into a vaccine calendar
# flu_doses(), which calculates how many vaccines were given (before wastage) in the same epidemics as many_flu()

# **flu_parallel.R** sets vaccine programs, 
# then runs many_flu() for some epidemic data, parallelised across each vaccine type