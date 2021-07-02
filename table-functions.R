library('flextable')
library('officer')

### page settings ###
sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)

#### Functions for gam tables ####
growth_tbl <- function(name, expo_var, out_var, exposure, outcome, results, results_adj, adj_only=F){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data frame with adjusted results
  ### adj_only: T or F if T will produce table with only the adjusted results, otherwise will display all results together
  
  ### this function produces a table that can be saved as a csv
  if(adj_only){
    tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                      "Outcome, 75th Percentile v. 25th Percentile" = character(),
                      " " = character(), " " = character(), " " = character())
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Adjusted", " ", " ", " "))
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                           "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value"))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
        if (nrow(filtered_adj)==0){
          skipped<-T
          next
        }
        
        sigadj <- ifelse(filtered_adj$BH.Pval < 0.2, "*", "")
        pvaladj <- paste(round(filtered_adj$Pval, 2), sigadj, sep="")
        
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, pvaladj))
          skipped<-F
        }else {
          tbl <- rbind(tbl, list("", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, pvaladj))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","",""))
      }
    }
  }else{
    tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                      " Outcome, 75th Percentile v. 25th Percentile" = character(), " " = character(), " " = character(), " " = character(), " " = character(),
                      " " = character(), " " = character(), " " = character(), " " = character(), " " = character())
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Unadjusted", " ", " ", " ", " ", "Adjusted", " ", " ", " ", " "))
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                           "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", "FDR adjusted P-value", 
                           "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", "FDR adjusted P-value"))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_res <- results[results$Y==out & results$X==exp,]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
        adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
        if (nrow(filtered_res)==0){
          skipped<-T
          next
        }
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped<-F
        }else {
          tbl <- rbind(tbl, list("", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","","","","","","","",""))
      }
    }
  }
  tbl
}

growth_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj, adj_only=F, exp_col_size = 1, out_col_size = 1){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  ### adj_only: T or F if T will produce table with only the adjusted results, otherwise will display all results together
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  if(adj_only){
    tbl <- data.table(matrix(nrow=0, ncol=9))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
        if (nrow(filtered_adj)==0){
          skipped<-T
          next
        }
        
        sigadj <- ifelse(filtered_adj$BH.Pval < 0.2, "*", "")
        pvaladj <- paste(round(filtered_adj$Pval, 2), sigadj, sep="")
        
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, pvaladj))
          skipped=F
        }else {
          tbl <- rbind(tbl, list(" ", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, pvaladj))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","",""))
      }
    }
    flextbl<-flextable(tbl, col_keys=names(tbl))
    flextbl <- set_header_labels(flextbl,
                                 values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                               "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", 
                                               "V8" = "Coefficient (95% CI)", "V9" = "P-value"))
    flextbl <- add_header_row(flextbl, values = c("","","","","", "Adjusted"), colwidths=c(1,1,1,1,1,4))
    flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,4))

  }else{
    tbl <- data.table(matrix(nrow=0, ncol=15))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_res <- results[results$Y==out & results$X==exp,]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
        adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
        if (nrow(filtered_res)==0){
          skipped<-T
          next
        }
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped=F
        }else {
          tbl <- rbind(tbl, list(" ", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","", "","","","","","",""))
      }
    }
    
    flextbl<-flextable(tbl, col_keys=names(tbl))
    flextbl <- set_header_labels(flextbl,
                                 values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                               "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", "V8" = "Coefficient (95% CI)", "V9" = "P-value", "V10" = "FDR Corrected P-value",
                                               "V11" = "Predicted Outcome at 25th Percentile", "V12" = "Predicted Outcome at 75th Percentile", "V13" = "Coefficient (95% CI)", "V14" = "P-value", "V15" = "FDR Corrected P-value"))
    flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Adjusted"), colwidths=c(1,1,1,1,1,5,5))
    flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,10))
  }
  
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  
  if(adj_only){
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the adjusted analyses", colwidths = 9)
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "* P-value < 0.2 after adjusting for multiple comparisons using the Benjamini-Hochberg procedure", colwidths = 9)
    flextbl <- fontsize(flextbl, part = "all", size = 7)
    flextbl <- width(flextbl, 1:9, width=c(exp_col_size, out_col_size, .3, .55, .55, .8, .8, 1, .5))
    
  }else{
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the unadjusted analyses", colwidths = 15)
    flextbl <- fontsize(flextbl, part = "all", size = 6)
    flextbl <- width(flextbl, 1:15, width=c(exp_col_size, out_col_size, .3, .5, .5, .5, .5, 1, .3, .5, .5, .5, 1, .3, .5))
    
  }
  
  flextbl
}

subgroup_tbl <- function(name, expo_var, out_var, sub_var, exposure, outcome, subgroup, results){
  # build table
  tbl <- data.table(matrix(nrow=0, ncol=10))
  skippedexp<-F
  for (k in 1:length(subgroup)){
    num.sub <- 0
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        sub <- subgroup[k]
        exp <- exposure[i]
        out <- outcome[j]
        
        filtered_adj <- results[results$Y==out & results$X==exp & results$V==sub,]
        v1 <- paste(round(filtered_adj$`point.diff`[1], 2), " (", round(filtered_adj$`lb.diff`[1], 2), ", ", round(filtered_adj$`ub.diff`[1], 2), ")", sep="")
        v2 <- paste(round(filtered_adj$`point.diff`[2], 2), " (", round(filtered_adj$`lb.diff`[2], 2), ", ", round(filtered_adj$`ub.diff`[2], 2), ")", sep="")
        
        if (nrow(filtered_adj)==0){
          skippedexp<-T
          next
        }
        
        if((i==1 & j==1)|num.sub==0){
          s_name <- sub_var[k]
          num.sub <- num.sub+1
        }else{
          s_name <- " "
        }
        if(j==1|skippedexp==T){
          e_name <- expo_var[i]
          skippedexp <- F
        }else{
          e_name <- " "
        }
        
        tbl <- rbind(tbl, list(s_name, e_name, out_var[j], filtered_adj$N[1], round(filtered_adj$Vlevel[1], 2),
                               v1, round(filtered_adj$Pval[1], 2), round(filtered_adj$BH.Pval[1], 2), "", ""))
        tbl <- rbind(tbl, list(" ", " ", " ", " ", round(filtered_adj$Vlevel[2], 2), 
                               v2, round(filtered_adj$Pval[2], 2), round(filtered_adj$BH.Pval[2], 2), 
                               round(filtered_adj$int.Pval[2], 2), round(filtered_adj$BH.int.Pval[2], 2)))
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, as.list(rep("",10)))
      }
    }
    if (k != length(subgroup)){
      tbl <- rbind(tbl, as.list(rep("",10)))
    }
  }
  
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ",
                                             "V5" = " ", "V6" = "Coefficient (95% CI)", "V7" = "P-value", "V8" = "FDR Corrected P-value",
                                             "V9" = "Interaction P-value", "V10" = "FDR Corrected Interaction P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","","","Adjusted"), colwidths=c(1,1,1,1,1,5))
  flextbl <- add_header_row(flextbl, values = c("Effect Modifier", name, "Outcome", "N", "Modifier value", 
                                                "Outcome, 75th Percentile v. 25th Percentile of Exposure"), colwidths=c(1,1,1,1,1,5))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2, 3), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  flextbl
}

hr_tbl <- function(name, expo_var, out_var, exposure, outcome, results, results_adj, adj_only=F){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data frame with adjusted results
  ### adj_only: T or F if T will produce table with only the adjusted results, otherwise will display all results together
  
  ### this function produces a table that can be saved as a csv
  if(adj_only){
    tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                      "Outcome, 75th Percentile v. 25th Percentile" = character(),
                      " " = character())
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Adjusted", " "))
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Hazard Ratio (95% CI)", "P-value"))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        adj <- paste(round(filtered_adj$`point.HR`, 2), " (", round(filtered_adj$`lb.HR`, 2), ", ", round(filtered_adj$`ub.HR`, 2), ")", sep="")
        if (nrow(filtered_adj)==0){
          skipped<-T
          next
        }
        
        sigadj <- ifelse(filtered_adj$BH.Pval < 0.2, "*", "")
        pvaladj <- paste(round(filtered_adj$Pval, 2), sigadj, sep="")
        
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 adj, pvaladj))
          skipped<-F
        }else {
          tbl <- rbind(tbl, list("", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 adj, pvaladj))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","",""))
      }
    }
  }else{
    tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                      "Outcome, 75th Percentile v. 25th Percentile" = character(), " " = character(), " " = character(),
                      " " = character(), " " = character(), " " = character())
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Unadjusted", " ", " ", "Adjusted", " ", " "))
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                           "Hazard Ratio (95% CI)", "P-value", "FDR adjusted P-value", 
                           "Hazard Ratio (95% CI)", "P-value", "FDR adjusted P-value"))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_res <- results[results$Y==out & results$X==exp,]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        unadj <- paste(round(filtered_res$`point.HR`, 2), " (", round(filtered_res$`lb.HR`, 2), ", ", round(filtered_res$`ub.HR`, 2), ")", sep="")
        adj <- paste(round(filtered_adj$`point.HR`, 2), " (", round(filtered_adj$`lb.HR`, 2), ", ", round(filtered_adj$`ub.HR`, 2), ")", sep="")
        if (nrow(filtered_res)==0){
          skipped<-T
          next
        }
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped<-F
        }else {
          tbl <- rbind(tbl, list("", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","","","",""))
      }
    }
  }
  tbl
}

hr_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj, adj_only=F, exp_col_size = 1, out_col_size = 1){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  ### adj_only: T or F if T will produce table with only the adjusted results, otherwise will display all results together
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  if(adj_only){
    tbl <- data.table(matrix(nrow=0, ncol=7))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        adj <- paste(round(filtered_adj$`point.HR`, 2), " (", round(filtered_adj$`lb.HR`, 2), ", ", round(filtered_adj$`ub.HR`, 2), ")", sep="")
        if (nrow(filtered_adj)==0){
          skipped<-T
          next
        }
        
        sigadj <- ifelse(filtered_adj$BH.Pval < 0.2, "*", "")
        pvaladj <- paste(round(filtered_adj$Pval, 2), sigadj, sep="")
        
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 adj, pvaladj))
          skipped=F
        }else {
          tbl <- rbind(tbl, list(" ", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 adj, pvaladj))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","",""))
      }
    }
    flextbl<-flextable(tbl, col_keys=names(tbl))
    flextbl <- set_header_labels(flextbl,
                                 values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                               "V6" = "Hazard Ratio (95% CI)", "V7" = "P-value"))
    flextbl <- add_header_row(flextbl, values = c("","","","","", "Adjusted"), colwidths=c(1,1,1,1,1,2))
    flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,2))
    
  }else{
    tbl <- data.table(matrix(nrow=0, ncol=11))
    skipped<-F
    for (i in 1:length(exposure)) {
      for (j in 1:length(outcome)) {
        exp <- exposure[i]
        out <- outcome[j]
        filtered_res <- results[results$Y==out & results$X==exp,]
        filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
        unadj <- paste(round(filtered_res$`point.HR`, 2), " (", round(filtered_res$`lb.HR`, 2), ", ", round(filtered_res$`ub.HR`, 2), ")", sep="")
        adj <- paste(round(filtered_adj$`point.HR`, 2), " (", round(filtered_adj$`lb.HR`, 2), ", ", round(filtered_adj$`ub.HR`, 2), ")", sep="")
        if (nrow(filtered_res)==0){
          skipped<-T
          next
        }
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped=F
        }else {
          tbl <- rbind(tbl, list(" ", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                                 unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                                 adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","","","",""))
      }
    }
    
    flextbl<-flextable(tbl, col_keys=names(tbl))
    flextbl <- set_header_labels(flextbl,
                                 values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                               "V6" = "Hazard Ratio (95% CI)", "V7" = "P-value", "V8" = "FDR Corrected P-value",
                                               "V9" = "Hazard Ratio (95% CI)", "V10" = "P-value", "V11" = "FDR Corrected P-value"))
    flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Adjusted"), colwidths=c(1,1,1,1,1,3,3))
    flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,6))
  }
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  

  if(adj_only){
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the adjusted analyses", colwidths = 7)
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "*P-value < 0.2 after adjusting for multiple comparisons using the Benjamini-Hochberg procedure", colwidths = 7)
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "Hazard ratio could not be estimated for sitting without support since nearly all children had achieved this milestone before time of measurement", colwidths = 7)
    flextbl <- fontsize(flextbl, part = "all", size = 7)
    flextbl <- width(flextbl, 1:7, width=c(exp_col_size, out_col_size, .3, .55, .55, 1.2, .5))
    
  }else{
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the unadjusted analyses", colwidths = 11)
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "Hazard ratio could not be estimated for sitting without support since nearly all children had achieved this milestone before time of measurement", colwidths = 11)
    flextbl <- fontsize(flextbl, part = "all", size = 6)
    flextbl <- width(flextbl, 1:11, width=c(exp_col_size, out_col_size, .3, .55, .55, 1, .5, .7, 1, .5, .7))
  }
  
  flextbl
}

