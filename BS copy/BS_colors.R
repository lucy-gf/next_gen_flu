
## COLORS PALETTES AND LABELS FOR OUTPUTS
#setwd("~/Desktop/research asst/Global Code")
library(viridis)

strain_colors <- c('tot' = 'black', 'totA' = '#7d66ac', 'totB' = '#e483a4')
strain_colors1 <- c('A' = '#7d66ac', 'B' = '#e483a4')
strain_colors2 <- c('INF_A' = '#7d66ac', 'INF_B' = '#e483a4')
vt_colors <- c('1' = '#d91818', '2' = '#e2790e', '3' = '#eacb2c', '4' = '#62adc1', '5' = '#324da0')
age_colors <- c('0-5' = '#FDE725FF', '5-20' = '#5DC863FF', '20-65' = '#21908CFF', '65+' = '#3B528BFF', 'Total' = '#440154FF')
age_colors1 <- c('1' = '#FDE725FF', '2' = '#5DC863FF', '3' = '#21908CFF', '4' = '#3B528BFF', '5' = '#440154FF')
cluster_colors2 <- c("Asia-Europe"="#21908CFF", "Southern America" = '#65156EFF', "Oceania-Melanesia-Polynesia" = '#CC4678FF', 
                     "Europe" = '#ffeba4', "Eastern and Southern Asia" = '#51C56AFF', "Northern America" = '#F89441FF', 
                     "Africa" = '#31688EFF')  
exemplar_colors <- c("TUR"="#21908CFF", "ARG" = '#65156EFF', "AUS" = '#CC4678FF', 
                     "GBR" = '#ffeba4', "CHN" = '#51C56AFF', "CAN" = '#F89441FF', 
                     "GHA" = '#31688EFF')  
income_colors <- c('High-income countries' = '#FDE725FF', 'Upper-middle-income countries' = '#21908CFF', 
                'Lower-middle-income countries' = '#CC4678FF', 'Low-income countries' = '#440154FF')
WHO_colors <- c("EMR"="#21908CFF", "WPR" = '#CC4678FF', 
                     "EUR" = '#65156EFF', "SEAR" = '#51C56AFF', "AMR" = '#F89441FF', 
                     "AFR" = '#31688EFF')  
age_targ_colors <- c('1' = '#FDE725FF', '2' = '#5DC863FF', '3' = '#21908CFF', '4' = '#3B528BFF', '5' = '#440154FF')

supp.labs <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs) <- c(1:5)
supp.labs1 <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs1) <- c(0:4)
supp.labs.n <- c('Current','Improved\n(minimal)','Improved\n(efficacy)','Improved\n(breadth)','Universal')
names(supp.labs.n) <- c(1:5)

who_region_labs <- c('African\nRegion','Region of\nthe Americas','South-East\nAsian Region','European\nRegion',
                     'Eastern\nMediterranean Region', 'Western\nPacific Region')
names(who_region_labs) <- c('AFR','AMR','EMR','EUR','SEAR','WPR')
who_region_labs2 <- c('African Region','Region of the Americas','South-East Asian Region','European Region',
                     'Eastern Mediterranean Region', 'Western Pacific Region')
names(who_region_labs2) <- c('AFR','AMR','EMR','EUR','SEAR','WPR')

supp.labs.age <- c('0-4','0-10','0-17','65+','0-17, 65+','None')
names(supp.labs.age) <- c(1:5,'None')
supp.labs.cov <- supp.labs.age

var_labs <- c('Deaths','Hospitalisations')
names(var_labs) <- c('deaths_av','hospitalisations_av')

supp.labs.agegrps <- c('0-4','5-19','20-64','65+')
names(supp.labs.agegrps) <- c(1:4)

supp.labs.strain <- c('Total','Influenza A','Influenza B')
names(supp.labs.strain) <- c('tot','totA','totB')

supp.labs.ITZ <- c("Africa", "Asia-Europe", "Eastern and\nSouthern Asia",
                   "Europe", "Northern\nAmerica", "Oceania-\nMelanesia-\nPolynesia",
                   "Southern\nAmerica")
names(supp.labs.ITZ) <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")
supp.labs2 <- c("Asia-\nEurope", "Africa", "Europe", "Southern\nAmerica",           
                "Oceania-\nMelanesia-\nPolynesia", "Eastern and\nSouthern Asia",  
                "Northern\nAmerica")
names(supp.labs2) <-  c('TUR','GHA','GBR','ARG','AUS','CHN','CAN')

supp.labs.country <- c("Ghana", "Turkey", "China","United\nKingdom", "Canada", "Australia","Argentina")
names(supp.labs.country) <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")

supp.labs.ITZ2 <- c("Africa", "Asia-Europe", "Eastern and\nSouthern Asia",
                   "Europe", "Northern\nAmerica", "Oceania-\nMelanesia-\nPolynesia",
                   "Southern\nAmerica")
names(supp.labs.ITZ2) <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                           "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                           "Southern America")


## CI functions 

eti50L <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[25] + 0.75*sort(x)[26])
  }else{stop('length(x) != 100')}
}
eti50U <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[76] + 0.75*sort(x)[75])
  }else{stop('length(x) != 100')}
}
eti95L <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[3] + 0.475*sort(x)[4])
  }else{stop('length(x) != 100')}
}
eti95U <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[98] + 0.475*sort(x)[97])
  }else{stop('length(x) != 100')}
}








