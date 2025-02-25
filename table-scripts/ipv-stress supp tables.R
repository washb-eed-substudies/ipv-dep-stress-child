rm(list=ls())

library('flextable')
library('officer')
library('data.table')
source(here::here("0-config.R")) 
source(here::here("table-functions.R"))
here::here()

# load enrollment characteristics and results
# d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))



#### MAIN TABLES ####
#### Table 1 ####
# Characteristics of participants
# nperc <- function(vector){
#   n <- sum(vector==1, na.rm=T)
#   perc <- round(n/sum(!is.na(vector))*100)
#   paste(n, " (", perc, "%)", sep="")
# }
# 
# mediqr <- function(vector){
#   quantiles <- round(quantile(vector, na.rm=T), 2)
#   paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
# }
# 
# n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
#                mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
#                mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
#                mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
#                mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
#                nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
#                mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
#                nperc(d$life_viol_any_t3))
# 
# tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
#                    "C2" = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
#                            "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
#                            "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
#                            "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
#                            "Intimate partner violence"),
#                    "C3" = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
#                            "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
#                            "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
#                            "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
#                            "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
#                    "C4" = n_med_col)
# 
# tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
# tbl1flex <- set_header_labels(tbl1flex,
#                         values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
# tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
# tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
# tbl1flex <- autofit(tbl1flex, part = "all")
# tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
# tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
# tbl1flex <- fit_to_width(tbl1flex, max_width=8)
# names(tbl1)<- c("","","","n (%) or median (IQR)")


#### Table 2 #### HYPOTHESIS 1: IPV ~ CHILD STRESS BIOMARKERS

exposure <- c("life_viol_any_t3")
outcome <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")
expo_var <- c("Maternal Exposure to IPV at any time")
out_var <- c("IPF(2a)-III", "2,3-dinor-iPF(a2)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI", "Compiled F2-isoprostanes Score",  "Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Pressure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

tbl2 <- growth_tbl("Exposure to IPV", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl2flex <- growth_tbl_flex("Exposure to IPV", expo_var, out_var, exposure, outcome, H1, H1adj)

#### Table 3 #### HYPOTHESIS 2: PERCEIVED STRESS ~ CHILD STRESS BIOMARKERS 

exposure <- c("pss_sum_mom_t3", "pss_sum_dad_t3")
outcome <- c("t3_saa_slope",  "t3_saa_z01", "t3_saa_z02", "t3_cort_slope","t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12" )
expo_var <- c("Maternal Perceived Stress", "Paternal Perceived Stress")
out_var <- c("Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Pressure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

tbl3 <- growth_tbl("Parental Stress", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl3flex <- growth_tbl_flex("Parental Stress", expo_var, out_var, exposure, outcome, H2, H2adj)

#### Table 4 #### HYPOTHESIS 3: DEPRESSION ~ CHILD STRESS BIOMARKERS

exposure <- c("cesd_sum_t2", "cesd_sum_ee_t3", "cesd_sum_t2_binary", "cesd_sum_ee_t3_binary")
outcome <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")
expo_var <- c( "CES-D score at Year 1", "CES-D score at Year 2",  "Binary CES-D at Year 1", "Binary CES-D at  Year  2")
out_var <- c("IPF(2a)-III", "2,3-dinor-iPF(a2)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI", "Compiled F2-isoprostanes Score", "Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Pressure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

tbl4 <- growth_tbl("Maternal Depression", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl4flex <- growth_tbl_flex("Maternal Depression", expo_var, out_var, exposure, outcome, H3, H3adj)


#### SAVE TABLES #### 

# write.csv(tbl1, file=here("tables/main/stress-growth-table1.csv")) --> in save_as_dox "Table 1" = tbl1flex,
write.csv(tbl2, file("tables/supplementary/stress-growth-table1.csv"))
write.csv(tbl3, here('tables/supplementary/stress-growth-table2.csv'))
write.csv(tbl4, here('tables/supplementary/stress-growth-table3.csv'))

save_as_docx( "Table 1" = tbl2flex, "Table 2" = tbl3flex, "Table 3" = tbl4flex, 
              path=here("tables/supplementary/IPV-Stress-Dep_supp_tables082222.docx"),
              pr_section = sect_properties)

