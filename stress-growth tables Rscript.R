rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R")) 

# load enrollment characteristics and results
# d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
# H4 <- readRDS(here('results/unadjusted/H4_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))
# H4adj <- readRDS(here('results/adjusted/H4_adj_res.RDS'))


#### Functions for growth tables ####
growth_tbl <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as a csv
  
  tbl <- data.table(" " = character(), " " = character(), " " = character(), " " = character(),
                    " Outcome, 75th Percentile v. 25th Percentile" = character(), " " = character(), " " = character(), " " = character(), 
                    " " = character(), " " = character(), " " = character(), " " = character())
  tbl <- rbind(tbl, list(" ", " ", " ", " ", "Unadjusted", " ", " ", " ", "Fully adjusted", " ", " ", " "))
  tbl <- rbind(tbl, list(" ", "Outcome", "25th Percentile", "75th Percentile", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value"))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list("", out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","","","","",""))
    }
  }
  tbl
}

growth_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  tbl <- data.table(matrix(nrow=0, ncol=12))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list(" ", out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","", "","","",""))
    }
  }
  
  # format for export
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = name, "V2" = "Outcome", "V3" = "25th Percentile", "V4" = "75th Percentile",
                                             "V5" = "Predicted Outcome at 25th Percentile", "V6" = "Predicted Outcome at 75th Percentile", "V7" = "Coefficient (95% CI)", "V8" = "P-value",
                                             "V9" = "Predicted Outcome at 25th Percentile", "V10" = "Predicted Outcome at 75th Percentile", "V11" = "Coefficient (95% CI)", "V12" = "P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","", "Unadjusted", "Fully adjusted"), colwidths=c(1,1,1,1,4,4))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- add_header_row(flextbl, values = c("","","","", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,8))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  flextbl
}



#### MAIN TABLES ####
#### Table 1 ####
# Characteristics of participants
nperc <- function(vector){
  n <- sum(vector==1, na.rm=T)
  perc <- round(n/sum(!is.na(vector))*100)
  paste(n, " (", perc, "%)", sep="")
}

mediqr <- function(vector){
  quantiles <- round(quantile(vector, na.rm=T), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
               mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
               mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
               mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
               mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
               nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
               mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
               nperc(d$life_viol_any_t3))

tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
                   "C2" = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
                           "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
                           "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
                           "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
                           "Intimate partner violence"),
                   "C3" = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
                           "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
                           "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
                           "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
                           "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                           "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                           "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
                           "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
                   "C4" = n_med_col)

tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
tbl1flex <- set_header_labels(tbl1flex,
                        values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
tbl1flex <- autofit(tbl1flex, part = "all")
tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
tbl1flex <- fit_to_width(tbl1flex, max_width=8)
names(tbl1)<- c("","","","n (%) or median (IQR)")


#### Table 2 #### HYPOTHESIS 1 

exposure <- c("life_viol_any_t3")
outcome <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")
expo_var <- c("Maternal Exposure to IPV at any time")
out_var <- c("IPF(2a)-III", "2,3-dinor-iPF(a2)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI",  "Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Perssure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

tbl2 <- growth_tbl("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl2flex <- growth_tbl_flex("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)

#### Table 3 #### HYPOTHESIS 2

exposure <- c("pss_sum_mom_t3", "pss_sum_dad_t3")
outcome <- c("t3_saa_slope",  "t3_saa_z01", "t3_saa_z02", "t3_cort_slope","t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12" )
expo_var <- c("Maternal Percieved Stress", "Paternal Percieved Stress")
out_var <- c("Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Perssure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

tbl3 <- growth_tbl("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl3flex <- growth_tbl_flex("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)

#### Table 4 #### HYPOTHESIS 3

exposure <- c("cesd_sum_t2", "cesd_sum_ee_t3")
outcome <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")
expo_var <- c( "CES-D score at Year 1", "CES-D score at Year 2")
out_var <- c("IPF(2a)-III", "2,3-dinor-iPF(a2)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI",  "Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Perssure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")
 
tbl4 <- growth_tbl("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl4flex <- growth_tbl_flex("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)



#### SAVE TABLES #### 

# write.csv(tbl1, file=here("tables/main/stress-growth-table1.csv")) --> in save_as_dox "Table 1" = tbl1flex,
write.csv(tbl2, file("tables/main/stress-growth-table2.csv"))
write.csv(tbl3, here('tables/main/stress-growth-table3.csv'))
write.csv(tbl4, here('tables/main/stress-growth-table4.csv'))

save_as_docx( "Table 2" = tbl2flex, "Table 3" = tbl3flex, "Table 4" = tbl4flex, path=here("tables/main/IPV-Stress-Dep_main_tables.docx"))

