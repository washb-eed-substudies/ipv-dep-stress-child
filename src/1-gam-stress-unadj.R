rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS("/Users/sophiatan/Downloads/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_stress==1)
names(d)

#Loop over exposure-outcome pairs

#### Hypothesis 1a ####
# Maternal exposure to cumulative lifetime IPV measured at Year 2 is associated with child stress biomarkers
Xvars <- c("life_viol_any_t3")            
Yvars <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca",
           "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  preds <- predict_gam_diff(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  if (grepl("t2", res$Y)) {
    preds$res$time <- "t2"
  } else {
    preds$res$time <- "t3"
  }
  H1_res <-  bind_rows(H1_res , preds$res)
}
H1_res$H <- 1

#Make list of plots
H1_plot_list <- NULL
H1_plot_data <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  simul_plot <- gam_simul_CI(H1_models$fit[i][[1]], H1_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_plot_list[[i]] <-  simul_plot$p
  H1_plot_data <-  rbind(H1_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H1_models, here("models/H1_models.RDS"))

#Save results
saveRDS(H1_res, here("results/unadjusted/H1_res.RDS"))


#Save plots
#saveRDS(H1_plot_list, here("figure-objects/H1_unadj_splines.RDS"))

#Save plot data
saveRDS(H1_plot_data, here("figure-data/H1_unadj_spline_data.RDS"))



#### Hypothesis 2 ####
#Parental perceived stress at Year 2 is associated with child stress biomarkers. 
Xvars <- c("pss_sum_mom_t3", "pss_sum_dad_t3")            
Yvars <- c("t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H2_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H2_models <- bind_rows(H2_models, res)
  }
}

#Get primary contrasts
H2_res <- NULL
for(i in 1:nrow(H2_models)){
  res <- data.frame(X=H2_models$X[i], Y=H2_models$Y[i])
  preds <- predict_gam_diff(fit=H2_models$fit[i][[1]], d=H2_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_res <-  bind_rows(H2_res , preds$res)
}
H2_res$H <- 2
H2_res$time <- "t3C"


#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_models)){
  res <- data.frame(X=H2_models$X[i], Y=H2_models$Y[i])
  simul_plot <- gam_simul_CI(H2_models$fit[i][[1]], H2_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H2_models, here("models/H2_models.RDS"))

#Save results
saveRDS(H2_res, here("results/unadjusted/H2_res.RDS"))


#Save plots
#saveRDS(H2_plot_list, here("figure-objects/H2_unadj_splines.RDS"))

#Save plot data
saveRDS(H2_plot_data, here("figure-data/H2_unadj_spline_data.RDS"))



#### Hypothesis 3 ####
#Maternal depression at Years 1 and 2, is associated with child stress biomarkers. 
Xvars <- c("cesd_sum_t2", "cesd_sum_t2_binary")            
Yvars <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca",
           "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H3_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H3_models <- bind_rows(H3_models, res)
  }
}

Xvars <- c("cesd_sum_ee_t3", "cesd_sum_ee_t3_binary")            
Yvars <- c("t3_saa_slope", "t3_saa_z01", "t3_saa_z02",
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03",
           "t3_map", "t3_hr_mean", 
           "t3_gcr_mean", "t3_gcr_cpg12") 

for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H3_models <- bind_rows(H3_models, res)
  }
}

#Get primary contrasts
H3_res <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  if(grepl("binary", H3_models$X[i])){
    preds <- predict_gam_diff(fit=H3_models$fit[i][[1]], d=H3_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H3_models$fit[i][[1]], d=H3_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  if(grepl("t2", res$Y)){
    preds$res$time <- "t2"
  }else{
    if(grepl("t3", res$X)){
      preds$res$time <- "t3C"
    }else{preds$res$time <- "t3S"}
  }
  H3_res <-  bind_rows(H3_res , preds$res)
}
H3_res$H <- 3

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  simul_plot <- gam_simul_CI(H3_models$fit[i][[1]], H3_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H3_models, here("models/H3_models.RDS"))

#Save results
saveRDS(H3_res, here("results/unadjusted/H3_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, here("figure-objects/H3_unadj_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("figure-data/H3_unadj_spline_data.RDS"))

