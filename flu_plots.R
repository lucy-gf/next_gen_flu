#### PLOTS ####

library(patchwork)
library(countrycode)

iso3c <- 'GBR'

infs_rds_list <- readRDS(infs_rds_list, file = paste0("outputs/vacc_", iso3c,".rds"))

infs_out <- data.table()
for(i in 1:length(infs_rds_list)){
  add_on <- infs_rds_list[[i]]
  add_on <- add_on[, c('time','vacc_type','I1','I2','I3','I4')]
  add_on[, tot := I1 + I2 + I3 + I4]
  infs_out <- rbind(infs_out, add_on)
}

infs_cum <- infs_out[, c('vacc_type','tot')][, lapply(.SD, cumsum), by=c('vacc_type')]
infs_cum[, time:=infs_out$time]

# colors
vt_colors <- c('no_vacc' = '#000000', 'current' = '#d91818', 'improved_minimal' = '#e2790e', 
               'improved_efficacy' = '#eacb2c', 'improved_breadth' = '#62adc1', 'universal' = '#324da0')

incidence <- ggplot(infs_out) + 
  geom_line(aes(x=time, y=tot/1000000, col=vacc_type), lwd=0.8) +
  ylab('Infections (millions)') + xlab('') + 
  scale_color_manual(values = vt_colors) + 
  labs(col='Vaccine type') +
  theme_minimal() + theme(text = element_text(size = 14))

cumulative <- ggplot(infs_cum) + 
  geom_line(aes(x=time, y=tot/1000000, col=vacc_type), lwd=0.8) +
  ylab('Cumulative infections (millions)') + xlab('Time') + 
  labs(col='Vaccine type') + 
  scale_color_manual(values = vt_colors) + 
  theme_minimal() + theme(text = element_text(size = 14))
  
incidence + ggtitle(paste0('Influenza, ', countrycode(iso3c,
                                            destination = 'country.name',
                                            origin = 'iso3c'))) + 
  cumulative + plot_layout(guides='collect', nrow=2)

ggsave(paste0("plots/influenza_", iso3c, '.png'),
              width=32,height=24,units="cm")






