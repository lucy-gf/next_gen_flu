
### ADDING INITIAL R0 TO SAMPLED EPIDEMICS ###
library(fluEvidenceSynthesis)
library(dplyr)
library(data.table)
library(ggplot2)
source("BS/BS_colors.R")

age.group.limits <- c(5,20,65) 

### SAMPLED EPIDEMICS ###

for(c_num in 1:7){
  cluster_code <- c('ARG','AUS','CAN','CHN','GHA','TUR','GBR')[c_num]
  sampled_epidemics <- read_csv(paste0("data_for_BS/sampled_epidemics_30_100_",cluster_code, ".csv"))
  sel_cntr <- c('Argentina','Australia','Canada','China','Ghana','Turkey','United Kingdom')[c_num]

  contact_matrix <- as.matrix(unname(list_contact_matr[[c_num]]))
  yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))
  group_ages <- stratify_by_age(yr_res_pop, limits = age.group.limits)
  contact_matrix_small <- t(t(contact_matrix)/group_ages)

  sampled_epidemics$r0 <- NA

  for(k_row in 1:nrow(sampled_epidemics)){
    sampled_epidemics$r0[k_row] <- fluEvidenceSynthesis::as_R0(
      transmission_rate = sampled_epidemics$trans[k_row],
      contact_matrix = contact_matrix_small,
      # age_groups = group_ages*c((0.2*1 + 0.8*sampled_epidemics$sus[k_row]),
      #                           rep(sampled_epidemics$sus[k_row],3))
      age_groups = group_ages
    )
    if(k_row %% 300 == 0){
      print(round(k_row/nrow(sampled_epidemics), digits=2))
    }
  }

  write_csv(sampled_epidemics, file = paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                             "_", n_simulations,'_',cluster_code,"_wr0.csv"))
}

sampled_epidemics <- data.frame()
for(cluster_code in c('ARG','AUS','CAN','CHN','GHA','TUR','GBR')){
  sampled_epidemics <- rbind(sampled_epidemics,
                             read_csv(paste0("data_for_BS/sampled_epidemics_30_100_",cluster_code, "_wr0.csv")))
}

sampled_epidemics %>% ggplot() +
  geom_point(aes(x=trans, y=r0, col=exemplar_code)) +
  scale_color_manual(values=exemplar_colors) +
  theme_minimal() + xlim(c(0,NA)) + ylim(c(0,NA)) +
  # xlab('Susceptibility') +
  ylab('R0') +
  theme_minimal() +
  # facet_grid(exemplar_code~strain, scales='fixed') +
  labs(col='Country') +
  theme(text=element_text(size=14))

sampled_epidemics %>% ggplot() +
  geom_histogram(aes(x=r0, fill=exemplar_code), bins=100) +
  scale_fill_manual(values=exemplar_colors) +
  theme_bw() +
  # facet_grid(exemplar_code~., scales='fixed') +
  xlab('R0') + ylab('Density') +
  theme(text=element_text(size=14))






