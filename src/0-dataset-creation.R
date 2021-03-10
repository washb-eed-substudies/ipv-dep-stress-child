rm(list=ls())

source(here::here("0-config.R"))


d <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-dm-ee-ipv-cesd-pss-covariates-stresslab (1).csv"))

# transform outcome distributions
d <- d %>% 
  mutate(
    t3_saa_z01_raw=t3_saa_z01, 
    t3_saa_z02_raw=t3_saa_z02, 
    t3_cort_z01_raw=t3_cort_z01, 
    t3_cort_z03_raw=t3_cort_z03, 
    t2_f2_8ip_raw=t2_f2_8ip, 
    t2_f2_23d_raw=t2_f2_23d, 
    t2_f2_VI_raw=t2_f2_VI,
    t2_f2_12i_raw=t2_f2_12i, 
    t3_gcr_mean_raw=t3_gcr_mean, 
    t3_gcr_cpg12_raw=t3_gcr_cpg12,
    t3_saa_z01=log(t3_saa_z01), 
    t3_saa_z02=log(t3_saa_z02), 
    t3_cort_z01=log(t3_cort_z01), 
    t3_cort_z03=log(t3_cort_z03), 
    t2_f2_8ip=log(t2_f2_8ip), 
    t2_f2_23d=log(t2_f2_23d), 
    t2_f2_VI=log(t2_f2_VI),
    t2_f2_12i=log(t2_f2_12i), 
    t3_gcr_mean2=logit(t3_gcr_mean/100), 
    t3_gcr_cpg12=logit(t3_gcr_cpg12/100))


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



#---------------------------------------------------------------------------------------------
# calc combined iso variable
#---------------------------------------------------------------------------------------------

# Combined exposures: To determine overall oxidative stress, we will combine our four measures of urinary F2-
# isoprostanes: iPF(2a)-III, 2,3-dinor-iPF(2a)-III, iPF(2a)-VI, and 8,12-iso-iPF(2a)-VI that are
# consistent with prior oxidative stress operationalization into a single score.29 We will use
# the first principal component of a principal components analysis of the four measures of
# urinary F2-isoprostanes as the oxidative stress score if all measures are correlated with each other (P-value < 0.2), 
# otherwise we will analyze the urinary F2-isoprostanes separately.

#Get correlation of isoprostanes
iso <- d %>% select(c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i"))       

cor(iso, use="pairwise.complete.obs")
cor.test(iso[,1], iso[,2])$p.value < 0.2
cor.test(iso[,1], iso[,3])$p.value < 0.2
cor.test(iso[,1], iso[,4])$p.value < 0.2
cor.test(iso[,2], iso[,3])$p.value < 0.2
cor.test(iso[,2], iso[,4])$p.value < 0.2
cor.test(iso[,3], iso[,4])$p.value < 0.2


#isoprostanes are significantly correlated, so calculate 1st principal component

df <-  dfull %>% select(c("childid","t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i")) %>% as.data.frame()
dim(df)
df <- df[complete.cases(df),]
dim(df)

# #Select assets and separate out ID
id<-subset(df, select=c("childid")) #drop subjectid
df<-df[,which(!(colnames(df) %in% c("childid")))]

##Computing the principal component using eigenvalue decomposition ##
princ.return <- princomp(df) 

## To get the first principal component in a variable ##
load <- loadings(princ.return)[,1]   

pr.cp <- as.matrix(df) %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 

df$t2_iso_pca <- as.numeric(pr.cp) ## Gives us the 1st PC in numeric form in pr.

#examine combined score
hist(df$t2_f2_8ip)
hist(df$t2_f2_23d)
hist(df$t2_f2_VI)
hist(df$t2_f2_12i)
hist(df$t2_iso_pca)

#check direction between individual isoprostanes and the combined score
ggplot(df, aes(x=t2_f2_8ip, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_23d, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_VI, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_12i, y=t2_iso_pca)) + geom_point() + geom_smooth()

#merge combined score back into main dataset
df.pca <- data.frame(childid=id, t2_iso_pca=df$t2_iso_pca)

dfull <- left_join(dfull, df.pca, by="childid")
hist(dfull$t2_iso_pca)



############# Check covariate missingness ###################
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


for(exposure in c("life_viol_any_t3", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "pss_sum_dad_t3")){
  d_sub <- subset(dfull, !is.na(dfull[,exposure]))
  if(exposure %in% c("life_viol_any_t3", "cesd_sum_ee_t3", "pss_sum_mom_t3")){
    vars <- c(Wvars, Wvars_mhle_t3)
  }else if(exposure %in% c("cesd_sum_t2")){
    vars <- c(Wvars, Wvars_cesd_t2)
  }else{
    vars <- c(Wvars, Wvars_pss_dad_t3)
  }
  W_sub <- d_sub %>% select(all_of(vars))  
  
  miss_sub <- data.frame(name = names(W_sub), missing = colSums(is.na(W_sub)), missing_prop = colSums(is.na(W_sub))/nrow(d_sub), row.names = c(1:ncol(W_sub)))
  for (i in 1:nrow(miss_sub)) {
    miss_sub$class[i] <- class(W_sub[,which(colnames(W_sub) == miss_sub[i, 1])])
  }
  print(exposure)
  print(miss_sub)
}

for(outcome in c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca",
                 "t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_residual_saa", 
                 "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_residual_cort",
                 "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")){
  d_sub <- subset(dfull, !is.na(dfull[,outcome]))
  if(grepl("t2_", outcome)){vars = Wvars2_F2}
  if(grepl("t3_gcr", outcome)){vars = Wvars3_oragene}
  if(grepl("map|hr", outcome)){vars = Wvars3_vital}
  if(grepl("saa|cort", outcome)){
    if(grepl("residual", outcome)){vars = Wvars3_salimetrics}
    else{vars = c(Wvars3_salimetrics, "t3_col_time_z01_cont")}
  }
  W_sub <- d_sub %>% select(all_of(c(Wvars, vars)))  
  
  miss_sub <- data.frame(name = names(W_sub), missing = colSums(is.na(W_sub)), missing_prop = colSums(is.na(W_sub))/nrow(d_sub), row.names = c(1:ncol(W_sub)))
  for (i in 1:nrow(miss_sub)) {
    miss_sub$class[i] <- class(W_sub[,which(colnames(W_sub) == miss_sub[i, 1])])
  }
  print(outcome)
  print(miss_sub)
}

# roof has low variability
mean(dfull$roof, na.rm=T)
sd(dfull$roof, na.rm=T)
# remove roof from covariates

saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-cesd-pss-covariates-stresslab.RDS"))
