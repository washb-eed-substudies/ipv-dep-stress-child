rm(list=ls())
source(here::here("0-config.R"))
library(tibble)
data <- tibble(x = -10:100, y= -10:100)
head(data)
here::here()

d <- readRDS("/Users/alexissilvera/ipv-dep-stress-child/figure-data/figures/ipv-cesd-pss-covariates-stresslab.RDS")
exposures_y1 <- c("cesd_sum_t2","cesd_sum_t2_binary")
outcomes_y1 <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca")
exposures_y2 <- c("cesd_sum_ee_t3", "cesd_sum_ee_t3_binary","life_viol_any_t3", "pss_sum_mom_t3", "pss_sum_dad_t3")
outcomes_y2 <- c("t3_saa_slope", "t3_saa_z01", "t3_saa_z02", "t3_cort_slope", "t3_cort_z01", "t3_cort_z03", "t3_map", "t3_hr_mean", "t3_gcr_mean", "t3_gcr_cpg12")


#function for filtering for only participants with at least one  outcome
filtering <- function(row){
  any(!is.na(row))}
    
y1_has_exposures<-d[apply(select(d, all_of(exposures_y1)), 1, filtering),]
y1_has_both<-y1_has_exposures[apply(select(y1_has_exposures, all_of(c(outcomes_y1, outcomes_y2))), 1, filtering),]
y1_clusters<-length(unique(y1_has_both$clusterid))
y1_n<-nrow(y1_has_both)

if(is.null(exposures_y2)){
  y2_has_exposures <- y1_has_exposures
}else{
  y2_has_exposures <- d[apply(select(d, all_of(exposures_y2)), 1, filtering),]
}
y2_has_both<-y2_has_exposures[apply(select(y2_has_exposures, all_of(outcomes_y2)), 1, filtering),]
y2_clusters<-length(unique(y2_has_both$clusterid))
y2_n<-nrow(y2_has_both)

data %>% 
  ggplot(aes(x, y)) +
  scale_x_continuous(minor_breaks = seq(10, 100, 10)) +
  scale_y_continuous(minor_breaks = seq(-10, 100, 10)) +
  theme_void() ->
  p

p +
  geom_rect(xmin = 25, xmax=75, ymin=96, ymax=100, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=98,label= '13,279 compounds assessed for eligibility', size=3) ->
  p

p +
  geom_rect(xmin = 54, xmax=92, ymin=84, ymax=94, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 73, y=89, label= 'Excluded: 7,728 compounds \n 7,429 compounds excluded to create buffer zones\n 219 compounds did not meet enrollment criteria\n 80 compounds declined to participate', size=3) +
  annotate('text', x= 10, y=89,label= 'Enrollment', size=4) +
  
  
  geom_rect(xmin = 20, xmax=80, ymin=76, ymax=82, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=79,label= '720 clusters created and randomly allocated across 7 arms \n 5,551 compounds randomly allocated across 7 arms \n 2 of 7 arms selected into substudy', size=3)  +
  annotate('text', x= 10, y=79,label= 'Allocation', size=4) +
  
  geom_rect(xmin = 27, xmax=73, ymin=64, ymax=74, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=70,label= "paste(bold('                            Control arm and\nNutrition + Water + Sanitation + Handwashing arm'))", parse=TRUE, size=3) +
  annotate('text', x= 50, y=68.5,label= '\n270 clusters \n2,068 households', size=3) +
  annotate('text', x= 10, y=63,label= 'Subsample Target', size=4) +
  
  
  geom_rect(xmin = 76, xmax=96, ymin=58, ymax=68, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 86, y=64.8,label= "paste(bold('Number of clusters not \n selected into substudy'))", parse=TRUE, size=3) + 
  annotate('text', x= 86, y=63.4,label= "paste(bold('Year 1 '))", parse=TRUE, size=3) + 
  annotate('text', x= 86, y=62.1,label= '139 clusters', size=3) + 
  annotate('text', x= 86, y=60.9,label= "paste(bold('Year 2'))", parse=TRUE, size=3) + 
  annotate('text', x= 86, y=59.7,label= '135 clusters', size=3) +
  
  geom_rect(xmin = 42, xmax=58, ymin=52, ymax=62, color='black',
            fill='white', size=0.25) +
  annotate('text', x=50, y=61.2,label= "paste(bold('Year 1'))", parse=TRUE, size=3) + 
  annotate('text', x=50, y=60.2,label= '\n\n131 clusters \n996 children', size=3) +
  annotate('text', x=50, y=56.2,label= "paste(bold('Year 2'))", parse=TRUE, size=3) + 
  annotate('text', x=50, y=55.2,label= '\n\n135 clusters \n1,021 children', size=3) +
  
  
  geom_rect(xmin = 37, xmax=63, ymin=26, ymax=50, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=48.9,label= "paste(bold('Year 1'))", parse=TRUE, size=3) + 
  annotate('text', x= 50, y=45.1,label= '\n\n240 children lost to follow-up \n23 moved \n45 absent \n76 withdrew \n66 no live birth \n30 child death ', size=3) + 
  annotate('text', x= 50, y=37.7,label= "paste(bold('Year 2'))", parse=TRUE, size=3) + 
  annotate('text', x= 50, y=32.3,label= '\n  25 new children measured \n262 children lost to follow-up \n63 moved \n5 absent \n90 withdrew \n67 no live birth \n37 child death ', size=3) + 
  annotate('text', x= 10, y=38.1,label= 'Follow-up', size=4) +
  
  geom_rect(xmin = 42, xmax=58, ymin=14, ymax=24, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=23.2,label= "paste(bold('Year 1'))", parse=TRUE, size=3) +
  annotate('text', x= 50, y=22.2,label= '\n\n131 clusters \n756 children', size=3) + 
  annotate('text', x= 50, y=18.2,label= "paste(bold('Year 2'))", parse=TRUE, size=3) +
  annotate('text', x= 50, y=17.2,label= '\n\n135 clusters \n759 children  ', size=3) + 
  annotate('text', x= 10, y=19,label= 'Subsample Enrollment', size=4) +
  
  
  geom_rect(xmin = 36, xmax=64, ymin=4, ymax=12, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=10.2,label= "paste(bold('Year 1'))", parse=TRUE, size=3) + 
  annotate('text', x= 50, y=9.6,label= paste("\n", 756-y1_n, '  missing exposure or outcome', sep=""), size=3) + 
  annotate('text', x= 50, y=7.2,label= "paste(bold('Year 2'))", parse=TRUE, size=3) + 
  annotate('text', x= 50, y=6.6,label= paste("\n", 759-y2_n, '  missing exposure or outcome', sep=""), size=3) + 
  annotate('text', x= 10, y=8,label= 'Specimen Collection', size=4) +
  
  
  geom_rect(xmin = 42, xmax=58, ymin=-8, ymax=2, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=1.2,label= "paste(bold('Year 1'), sep='')", parse=T, size=3) +
  annotate('text', x= 50, y=.2,label= paste("\n\n", y1_clusters, ' clusters \n ', y1_n, ' children', sep=''), size=3) +
  annotate('text', x= 50, y=-3.8,label= "paste(bold('Year 2'), sep='')", parse=T, size=3) +
  annotate('text', x= 50, y=-4.8,label= paste("\n\n", y2_clusters, ' clusters \n ', y2_n, ' children', sep=''), size=3) +
  annotate('text', x= 10, y=-3,label= 'Analysis', size=4) ->
  p


p +
  geom_segment(
    x=50, xend=50, y=96, yend=82, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=54, y=89, yend=89, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=76, yend=74, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=64, yend=62, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=76, y=63, yend=63, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=52, yend=50, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=26, yend=24, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=14, yend=12, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=50, y=4, yend=2, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) ->
  p
p

# YOU MAY NEED TO CHANGE THE FILE PATHS HERE
ggsave(p, file = here("figures/enrollment_figure1.jpg"), height=14, width=9)


