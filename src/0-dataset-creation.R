rm(list=ls())

source(here::here("0-config.R"))


d <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-dm-ee-ipv-cesd-pss-covariates-stresslab.csv"))


## add household wealth
d_hhwealth <- read.csv("C:/Users/Sophia/Documents/ee-secondary/sophia scripts/hhwealth.csv")
dfull <- left_join(d, d_hhwealth, by="dataid")

# convert time of day of pre-stressor measurement of cortisol and sAA into continuous variable
time_day <- dfull$t3_col_time_z01
time_split <- str_split(time_day, ":")
cont_time <- function(list_hr_min){
  # takes in list of time
  # first number is hour of the day
  # second number in list is minute of the hour
  num_time <- as.numeric(unlist(list_hr_min))
  num_time[1]+num_time[2]/60
}
dfull$t3_col_time_z01_cont <- sapply(time_split, cont_time)

# add variables to turn cesd into binary variables
# classify top 25% of mothers in sample as experiencing high depressive symptoms
cesd_t2_q<-quantile(dfull$cesd_sum_t2, na.rm=T)[4]
cesd_t3_q<-quantile(dfull$cesd_sum_ee_t3, na.rm=T)[4]
dfull$cesd_sum_t2_binary<-ifelse(dfull$cesd_sum_t2 >= cesd_t2_q, 
                                               1, 0)
dfull$cesd_sum_ee_t3_binary<-ifelse(dfull$cesd_sum_ee_t3 >= cesd_t3_q, 
                                                  1, 0)


saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-cesd-pss-covariates-stresslab.RDS"))
