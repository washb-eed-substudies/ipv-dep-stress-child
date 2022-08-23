library(data.table)
source(here::here("0-config.R"))
library(tidyverse)
library(flextable)
library(officer)

d <- readRDS("/Users/sophiatan/Downloads/bangladesh-cleaned-master-data.RDS") 
d <- d %>% filter(ipv_stress == 1)
#d <- readRDS("/Users/alexissilvera/Desktop/?./WASHB + RISE/IPV, Dep, Stress - Child SAP/ipv-dep-stress-child/ipv-cesd-pss-covariates-stresslab.RDS")

writeqntle<-function(vector) {
  quantiles<-round(quantile(vector, na.rm=TRUE), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

nperc <- function(vector){
  n <- sum(vector==1, na.rm=T)
  perc <- round(n/sum(!is.na(vector))*100)
  paste(n, " (", perc, "%)", sep="")
}

#PARENTAL TABLES 
mom_lab <-c(" ", "Maternal Perceived Stress", "Paternal Perceived Stress", "CES-D score at Year 1", "CES-D score at Year 2", " ", " ", "Maternal Exposure to IPV at any time")
mom_d <- d %>% distinct(dataid, .keep_all = T)
mom <-c("Median (25th, 75th percentile)", 
        writeqntle(mom_d$pss_sum_mom_t3),  writeqntle(mom_d$pss_sum_dad_t3), writeqntle(mom_d$cesd_sum_t2),  writeqntle(mom_d$cesd_sum_ee_t3),  
        " ", "N (Percent Exposed)",
        nperc(mom_d$life_viol_any_t3))



#CHILD TABLES 
child_lab_t2 <-c(" ", "IPF(2a)-III", "2,3-dinor-iPF(a2)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI", "Compiled F2-isoprostanes Score")  
                 
child_lab_t3 <-c(" ", "Pre to post-stress change in slope of sAA", "Pre-Stressor SAA", "Post-Stressor SAA", "Change in Slope of Cortisol", "Pre-Stressor Cortisol", "Post-Stressor Cortisol", "Mean Arterial Perssure", "Resting Heart Rate", "Entire promoter region (39 assayed CpG sites)","NGFI-A transcription factor binding site (CpG site #12)")

child_t2 <- c("Median (25th, 75th percentile)", writeqntle(d$t2_f2_8ip), writeqntle(d$t2_f2_23d), writeqntle(d$t2_f2_VI), writeqntle(d$t2_f2_12i), writeqntle(d$t2_iso_pca))

child_t3 <- c("Median (25th, 75th percentile)", writeqntle(d$t3_saa_slope), writeqntle(d$t3_saa_z01), writeqntle(d$t3_saa_z02), writeqntle(d$t3_cort_slope), writeqntle(d$t3_cort_z01), writeqntle(d$t3_cort_z03), writeqntle(d$t3_map), writeqntle(d$t3_hr_mean), writeqntle(d$t3_gcr_mean), writeqntle(d$t3_gcr_cpg12))

mom_tbl<-data.table(" "= mom_lab,
                    "At Enrollment"=mom)

child_tbl_t2<-data.table(" "= child_lab_t2,
                      "Age 14 Months"=child_t2)
                      
child_tbl_t3<-data.table(" " = child_lab_t3, 
                      "Age 28 Months"=child_t3)

view(mom_tbl)
view(child_tbl_t2)
view(child_tbl_t3)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)

save_as_docx("Maternal Exposure to IPV, Stress, & Depression" = flextable(mom_tbl), path="/Users/alexissilvera/Desktop/WASHB + RISE/WASHB Manuscript/mom_tbl1.docx", 
             pr_section = sect_properties) 

save_as_docx("Child Biomarkers 14 Months" = flextable(child_tbl), path="/Users/alexissilvera/Desktop/WASHB + RISE/WASHB Manuscript/child_tbl1_t2.docx", 
             pr_section = sect_properties) 
save_as_docx("Child Biomarkers 28 Months" = flextable(child_tbl), path="//Users/alexissilvera/Desktop/WASHB + RISE/WASHB Manuscript/child_tbl1_t3.docx", 
             pr_section = sect_properties) 



characteristics <- function(d, child_char = NULL, child_char_names = NULL, mom_char = NULL, mom_char_names = NULL) {
  nperc <- function(vector){
    n <- sum(vector==1, na.rm=T)
    perc <- round(n/sum(!is.na(vector))*100)
    paste(n, " (", perc, "%)", sep="")
  }
  
  mediqr <- function(vector){
    quantiles <- round(quantile(vector, na.rm=T), 2)
    paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
  }
  
  child <- c('sex', child_char,'laz_t1','waz_t1','whz_t1','hcz_t1',
             'laz_t2','waz_t2','whz_t2','hcz_t2',
             'laz_t3','waz_t3','whz_t3','hcz_t3','diar7d_t2','diar7d_t3')
  
  mom <- c('momage', 'momheight', 'momeduy', mom_char, 'cesd_sum_t2', 'cesd_sum_ee_t3', 'pss_sum_mom_t3', 'life_viol_any_t3')
  
  n_med_col <- NULL
  for (var in c(child, mom)) {
    if (var %in% c('sex', 'diar7d_t2', 'diar7d_t3', 'life_viol_any_t3', 'hfiacat_ind') | is.factor(d[[var]])) {
      if (var == 'sex') {
        n <- sum(d$sex=='female', na.rm=T)
        perc <- round(n/sum(!is.na(d$sex))*100)
        n_med_col <- c(n_med_col, paste(n, " (", perc, "%)", sep=""))
      }else {
        d[[var]] <- na_if(d[[var]], "Missing")
        n_med_col <- c(n_med_col, nperc(d[[var]]))
      }
    }else {
      n_med_col <- c(n_med_col, mediqr(d[[var]]))
    }
  }
  
  tbl1 <- data.table("C1" = c("Child", rep("", length(child)-1),"Mother", rep("",length(mom)-1)),
                     "C2" = c("", rep("", length(child_char)), "Anthropometry (3 months)","","","",
                              "Anthropometry (14 months)","","","",
                              "Anthropometry (28 months)","","","", "Diarrhea (14 months)", "Diarrhea (28 months)","",
                              "Anthropometry at enrollment", "Education", rep("", length(mom_char)), "Depression 14 Months", "Depression 28 Months", "Perceived stress 28 Months", 
                              "Intimate partner violence"),
                     "C3" = c("Female", child_char_names,
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                              "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
                              mom_char_names, "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
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
  tbl1flex
}

enroll <- characteristics(mom_d)
enroll2 <- characteristics(d)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)
save_as_docx("Maternal characteristics" = enroll, "Child characteristics"=enroll2, path="tables/main/enrollment.docx", 
             pr_section = sect_properties) 


#write.csv(tbl, file=here('tables/biomarkers.csv'))
#print(xtable(tbl), type="html", file=here("tables/biomarkers.html"))