# next_gen_flu

Reproducible flexible next-generation vaccine flu model. Takes some vaccine programs and epidemic data, and outputs weekly vaccination-status and age-specific influenza infections.

Based on https://github.com/lucy-gf/flu_model_LG

## ```vacc_types.R```
Contains functions to define vaccine programs for no vaccinations and 5 NGIVs (can be edited).

- Inputs: coverage level, targeted ages
- Outputs: VE, mean immunity length, coverage across model age groups

##```functions/demography.R```
Calculates weekly age- and vaccination-status specific population over the relevant time period (can remove ageing if needed), for a given country and vaccine program, and contains the function to calculate contact matrices from a given country and age-specific population.

## ```functions/transmission_model.R```
Contains the ODE model builder, epidemic simulation function, a function to calculate vaccination status-specific demography.

## ```functions/flu_sim.R```
Contains:
- ```one_flu()```, which runs an epidemic,
- ```many_flu()```, which takes a data-table of epidemic data and combines many epidemics,
- ```dfn_vaccine_calendar()```, which converts the vaccine program and epidemic dates into a vaccine calendar
- ```flu_doses()```, which calculates how many vaccines were given (before wastage) in the same epidemics as ```many_flu()```

## ```flu_parallel.R```
Sets vaccine programs, then runs ```many_flu()``` for some epidemic data, parallelised across each vaccine type