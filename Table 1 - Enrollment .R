library(data.table)
source(here::here("0-config.R"))
library(tidyverse)
library(flextable)
library(officer)

d <- readRDS("/Users/alexissilvera/Desktop/?./WASHB + RISE/IPV, Dep, Stress - Child SAP/ipv-dep-stress-child/ipv-cesd-pss-covariates-stresslab.RDS")

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
mom <-c("Median (25th, 75th percentile)", 
        writeqntle(d$pss_sum_mom_t3),  writeqntle(d$pss_sum_dad_t3), writeqntle(d$cesd_sum_t2),  writeqntle(d$cesd_sum_ee_t3),  
        " ", "N (Percent Exposed)",
        nperc(d$life_viol_any_t3))



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

#write.csv(tbl, file=here('tables/biomarkers.csv'))
#print(xtable(tbl), type="html", file=here("tables/biomarkers.html"))