#### FUNCTIONS TO PRODUCE/REDUCE CONTACT MATRICES ####

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# create reduced contact matrix ----------------------------------------------------------
fun_create_red_C_m <- function(C_m_full,model_agegroups,orig_age_groups_duration,orig_age_groups_sizes){
  n_age <- nrow(model_agegroups)
  C_m=matrix(0,nrow=n_age, ncol=n_age)
  rownames(C_m)=model_agegroups$agegroup_name; colnames(C_m)=model_agegroups$agegroup_name
  for (i_row in 1:n_age){
    for (j_col in 1:n_age){
      # we are merging or splitting age groups, there are 3 possibilities for an age group in the NEW matrix:
      # 1) same 2) smaller than the group in original CM 3) larger than group in the original CM
      # (it is an *average* of contacts per person, and we have no resolution within age bands)
      #
      # if the 'i' group (C[i,j]) is the *same* as original or *smaller* contact rate UNCHANGED
      if (model_agegroups$wpp_agegroup_low[i_row]==model_agegroups$wpp_agegroup_high[i_row]) {
        # if contact (j) group same or smaller as original
        if (model_agegroups$wpp_agegroup_low[j_col]==model_agegroups$wpp_agegroup_high[j_col]) {
          # proportionality by 'time window' size (age x to x+n)
          f_dur=model_agegroups$duration[j_col]/orig_age_groups_duration[model_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=(C_m_full[model_agegroups$wpp_agegroup_low[i_row],
                                     model_agegroups$wpp_agegroup_low[j_col]])*f_dur
          # print("'i' group (C[i,j]) is the same/smaller as original & 'j' group same/smaller as original"); 
        } else { 
          # contact numbers BY contact groups (j) we need to SUM
          group_span=model_agegroups$wpp_agegroup_low[j_col]:model_agegroups$wpp_agegroup_high[j_col]
          agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
          C_m[i_row,j_col]=sum(C_m_full[i_row,group_span]) # agegroup_weights*
          # print("'i' group (C[i,j]) is the *same* as original or *smaller* & 
          #       'j' is larger than original group"); print(c(i_row,j_col))
        } # end of 'i' smaller or same as original
      } else { 
        # if 'i' in C[i,j] is a bigger age band -> weighted average of the contact rates of participant groups
        group_span=model_agegroups$wpp_agegroup_low[i_row]:model_agegroups$wpp_agegroup_high[i_row]
        agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
        # if 'j' is same/smaller -> contact rate with original group proportionally divided
        if (model_agegroups$wpp_agegroup_low[j_col]==model_agegroups$wpp_agegroup_high[j_col]) {
          f_dur=model_agegroups$duration[j_col]/orig_age_groups_duration[model_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=
            sum((orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span]))*C_m_full[group_span,j_col])*f_dur
          # print("'i' group (C[i,j]) is bigger age band than original & 
          #       'j' is same as original group"); print(c(i_row,j_col))
          
        } else {
          # if 'j' larger -> SUM of contacts by contact groups
          j_group_span=model_agegroups$wpp_agegroup_low[j_col]:model_agegroups$wpp_agegroup_high[j_col]
          # contact numbers BY participant groups we need to AVERAGE (proport to pop size) over
          # contact numbers BY contact groups we need to SUM
          C_m[i_row,j_col]=sum(agegroup_weights*unlist(lapply(group_span, function(x) sum(C_m_full[x,j_group_span]))))
          # print("'i' group (C[i,j]) is bigger age band than original & 
          #       'j' is bigger than original group"); print(c(i_row,j_col))
          
        }
      }
      
    } 
  }
  C_m 
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# funcn create reciprocal matrix ----------------------------------------------------------
fun_recipr_contmatr <- function(C_m_full,age_group_sizes){
  all_perms=merge(1:nrow(C_m_full),1:nrow(C_m_full)) # permutations(n=nrow(C_m_full),r=2,repeats.allowed=T)
  N_tot=sum(age_group_sizes)
  C_m_full_symm=matrix(0,nrow=nrow(C_m_full),ncol=nrow(C_m_full))
  for (k in 1:nrow(all_perms)) { 
    i=all_perms[k,1]; j=all_perms[k,2]
    C_m_full_symm[i,j]=(C_m_full[i,j] + C_m_full[j,i]*(age_group_sizes[j]/age_group_sizes[i]))/2
  }
  colnames(C_m_full_symm)=colnames(C_m_full); rownames(C_m_full_symm)=rownames(C_m_full) 
  C_m_full_symm
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# age structure of country ----------------------------------------------------------
fun_cntr_agestr <- function(i_cntr,i_year,age_low_vals,age_high_vals){
  age_groups=data.frame(age_low=seq(0,75,5), age_high=c(seq(4,74,5),100))
  if (!any((.packages()) %in% "wpp2019")) {library(wpp2019)}; if (!exists("popF")) {data("pop")}
  if(i_cntr[1] == 'Kosovo'){
    pop_hist_WPP_data <- data.frame(pop_hist_WPP_data)
    vals <- unname(unlist(pop_hist_WPP_data %>% filter(grepl('Kos', name), Year == 2020) %>% 
                            select(!c(name, Type, Year))))
    cntr_agestr = data.frame(agegroups=popF[popF$name %in% 'France',"age"],
                             values=vals)
  }else{
    cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],
                           values=popF[popF$name %in% i_cntr,i_year] + popM[popM$name %in% i_cntr,i_year])
  }
  agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
  N_tot=cntr_agestr$values[agegr_truthvals]
  N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals])
  N_tot=N_tot*1e3; # N_tot
  data.frame(age_low=age_low_vals, age_high=age_high_vals,values=N_tot, 
             duration=(age_high_vals-age_low_vals)+1) %>% mutate(proportion=values/sum(values))
}




