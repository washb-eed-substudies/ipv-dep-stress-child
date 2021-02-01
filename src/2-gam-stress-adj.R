rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-cesd-pss-covariates-stresslab.RDS"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex", "birthord", "momage","momheight","momedu", 
         "hfiacat", "Ncomp", "watmin", "walls", "floor", "tr", 
         "HHwealth")

Wvars[!(Wvars %in% colnames(d))]

#Add in time varying covariates:
Wvars_mhle_t3<-c("mhle_month_t3")
Wvars_cesd_t2<-c("cesd_month_t2")
Wvars_pss_dad_t3<-c("pss_dad_month_t3")

Wvars2_F2<-c("ageday_ut2", "month_ut2") 
Wvars3_vital<-c("ageday_t3_vital", "month_vt3") 
Wvars3_salimetrics<-c("ageday_t3_salimetrics", "month_lt3") 
Wvars3_oragene<-c("ageday_t3_oragene", "month_ot3") 

#maternal covariates
W_t3mat_F2 <- c(Wvars, Wvars2_F2, Wvars_mhle_t3) %>% unique(.)
W_t3mat_vitals <- c(Wvars, Wvars3_vital, Wvars_mhle_t3) %>% unique(.)
W_t3mat_salimetrics <- c(Wvars, Wvars3_salimetrics, Wvars_mhle_t3) %>% unique(.)
W_t3mat_oragene <- c(Wvars, Wvars3_oragene, Wvars_mhle_t3) %>% unique(.)

W_t2mat_F2 <- c(Wvars, Wvars2_F2, Wvars_cesd_t2) %>% unique(.)
W_t2mat_vitals <- c(Wvars, Wvars3_vital, Wvars_cesd_t2) %>% unique(.)
W_t2mat_salimetrics <- c(Wvars, Wvars3_salimetrics, Wvars_cesd_t2) %>% unique(.)
W_t2mat_oragene <- c(Wvars, Wvars3_oragene, Wvars_cesd_t2) %>% unique(.)

#paternal covariates
W_t3pat_vitals <- c(Wvars, Wvars3_vital, Wvars_pss_dad_t3) %>% unique(.)
W_t3pat_salimetrics <- c(Wvars, Wvars3_salimetrics, Wvars_pss_dad_t3) %>% unique(.)
W_t3pat_oragene <- c(Wvars, Wvars3_oragene, Wvars_pss_dad_t3) %>% unique(.)



#Loop over exposure-outcome pairs
pick_covariates_H1 <- function(j){
  if(grepl("t2_", j)){Wset = W_t3mat_F2}
  if(grepl("t3_gcr", j)){Wset = W_t3mat_oragene}
  if(grepl("map|hr", j)){Wset = W_t3mat_vitals}
  if(grepl("saa|cort", j)){
    if(grepl("residual", j)){Wset = W_t3mat_salimetrics}
    else{Wset = c(W_t3mat_salimetrics, "t3_col_time_z01_cont")}
  }
  return(Wset)
}

#### Hypothesis 1 ####
# Maternal exposure to cumulative lifetime IPV measured at Year 2 is negatively associated with child telomere length measured at Year 2
Xvars <- c("life_viol_any_t3")            
Yvars <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca",
           "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa", 
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12") 


#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H1(j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}


#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H1_adj_models, here("models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/H1_adj_spline_data.RDS"))




#### Hypothesis 2 ####
pick_covariates_H2 <- function(i, j){
  if(grepl("mom", i)){
    if(grepl("t3_gcr", j)){Wset = W_t3mat_oragene}
    if(grepl("map|hr", j)){Wset = W_t3mat_vitals}
    if(grepl("saa|cort", j)){
      if(grepl("residual", j)){Wset = W_t3mat_salimetrics}
      else{Wset = c(W_t3mat_salimetrics, "t3_col_time_z01_cont")}
    }
  }else{
    if(grepl("t3_gcr", j)){Wset = W_t3pat_oragene}
    if(grepl("map|hr", j)){Wset = W_t3pat_vitals}
    if(grepl("saa|cort", j)){
      if(grepl("residual", j)){Wset = W_t3pat_salimetrics}
      else{Wset = c(W_t3pat_salimetrics, "t3_col_time_z01_cont")}
    }
  }
  return(Wset)
}

Xvars <- c("pss_sum_mom_t3", "pss_sum_dad_t3")            
Yvars <- c("t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H2(i, j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_adj_plot_list <- NULL
H2_adj_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_adj_plot_list[[i]] <-  simul_plot$p
  H2_adj_plot_data <-  rbind(H2_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H2_adj_models, here("models/H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_splines.RDS"))

#Save plot data
saveRDS(H2_adj_plot_data, here("figure-data/H2_adj_splint_data.RDS"))




#### Hypothesis 3 ####
pick_covariates_H3 <- function(i, j){
  if(grepl("t2", i)){
    if(grepl("t2_", j)){Wset = W_t2mat_F2}
    if(grepl("t3_gcr", j)){Wset = W_t2mat_oragene}
    if(grepl("map|hr", j)){Wset = W_t2mat_vitals}
    if(grepl("saa|cort", j)){
      if(grepl("residual", j)){Wset = W_t2mat_salimetrics}
      else{Wset = c(W_t2mat_salimetrics, "t3_col_time_z01_cont")}
    }
  }else{
    if(grepl("t3_gcr", j)){Wset = W_t3mat_oragene}
    if(grepl("map|hr", j)){Wset = W_t3mat_vitals}
    if(grepl("saa|cort", j)){
      if(grepl("residual", j)){Wset = W_t3mat_salimetrics}
      else{Wset = c(W_t3mat_salimetrics, "t3_col_time_z01_cont")}
    }
  }
  return(Wset)
}

Xvars <- c("cesd_sum_t2", "cesd_sum_t2_binary")            
Yvars <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca",
           "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H3(i, j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

Xvars <- c("cesd_sum_ee_t3", "cesd_sum_ee_t3_binary")            
Yvars <- c("t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H3(i, j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}


#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  if(grepl("binary", H3_adj_models$X[i])){
    preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_adj_plot_list <- NULL
H3_adj_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_adj_plot_list[[i]] <-  simul_plot$p
  H3_adj_plot_data <-  rbind(H3_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H3_adj_models, here("models/H3_adj_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
saveRDS(H3_adj_plot_data, here("figure-data/H3_adj_spline_data.RDS"))

