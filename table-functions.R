library('flextable')
library('officer')

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
                      " " = character(), " " = character(), " " = character(), " " = character())
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Adjusted", " ", " ", " ", " "))
    tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                           "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", "FDR adjusted P-value"))
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
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped<-F
        }else {
          tbl <- rbind(tbl, list("", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","","",""))
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

growth_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj, adj_only=F){
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
    tbl <- data.table(matrix(nrow=0, ncol=10))
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
        if(j==1|skipped==T){
          tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
          skipped=F
        }else {
          tbl <- rbind(tbl, list(" ", out_var[j],  filtered_adj$N, round(filtered_adj$q1, 2), round(filtered_adj$q3, 2), 
                                 round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        }
      }
      if (i != length(exposure)) {
        tbl <- rbind(tbl, list("","","","","","","","","",""))
      }
    }
    flextbl<-flextable(tbl, col_keys=names(tbl))
    flextbl <- set_header_labels(flextbl,
                                 values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                               "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", 
                                               "V8" = "Coefficient (95% CI)", "V9" = "P-value", "V10" = "FDR Corrected P-value"))
    flextbl <- add_header_row(flextbl, values = c("","","","","", "Adjusted"), colwidths=c(1,1,1,1,1,5))
    flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,5))

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
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  if(adj_only){
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the adjusted analyses", colwidths = 10)
  }else{
    flextbl <- add_footer_row(flextbl, top=F, 
                              values = "N, 25th Percentile, and 75th Percentile are from the unadjusted analyses", colwidths = 15)
  }
  
  flextbl
}

