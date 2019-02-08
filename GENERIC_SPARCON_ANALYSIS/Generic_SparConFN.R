## generic SparCon FN

#normalize within well
#input: table of structure:
#wellID - ind_type - Measurement - comparison_var
norm.well <- function(dat)
{
  ###for ctrl means of EXP
  #get ctrl means
  mccd <- dat[, .SD]
  
  ctrl.means <- mccd[ind_type=="CTRL",
                     lapply(.SD, geometric.mean),
                     by=.(wellID),
                     .SDcols = "Measurement"]
  
  for (i in 1:dim(ctrl.means)[1])
  {
    iWell <- ctrl.means[i, wellID]
    print(paste("Normalized well", iWell))
    mccd[wellID==iWell, `:=`(Measurement = as.double(as.double(Measurement) / ctrl.means[wellID==iWell, Measurement]))]
  }
  return(mccd)
}



exp.plfn.ecdf <- function(gp, alpha = 0.5, size = 2.55, scale = 1)
{gp+
    stat_ecdf(geom = "line", lwd = 1) +
    stat_ecdf(geom="point", alpha=alpha, size = size, aes(fill = ind_type, color = ind_type), shape = 21) +
    theme_few()+
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 12*scale, color = "black"),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
  
}

exp.theme <- function(gp, scale = 1)
{gp+
    theme_few()+
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
  
}


exp.plfn.geomPOINT <- function(gp, alpha = 0.5, size = 2.55, jitterwidth = 0.2, scale = 1, meanwidth = 0.3,
                               sewidth = 0.1)
{gp+
    geom_point(alpha=alpha, size = size, aes(fill = ind_type), color = "black", shape = 21, 
               position = position_jitter(width=jitterwidth)) +
    #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", color = "black",width = meanwidth) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = sewidth)+
    theme_few()+
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12*scale, color = "black"),
      axis.ticks.x = element_blank(),
      strip.text = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale)
    )
}

#output a neat AD-test table by comparison_var
#input: db, var_list
ad_results <- function(db, var_label = "Measurement")
{
      ctrl <- (db[ind_type == "CTRL", Measurement ])
      mut <- (db[ind_type == "MUT", Measurement ])
      test_result <- ad.test(ctrl,mut)
      results_table <- data.frame(CTRL_grp_name = "CTRL",
                        MUT_grp_name = "MUT",
                        CTRL_n = length(ctrl),
                        MUT_n = length(mut),
                        Variable = var_label,
                        CTRL_mean = mean(ctrl),
                        CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                        MUT_mean = mean(mut),
                        MUT_sem = sd(mut)/sqrt(length(mut)),
                        Fold_change = mean(mut) / mean(ctrl),
                        t.AD = test_result$ad[2,2],
                        p_val = test_result$ad[2,3]
      )
  return(results_table)
}

fmt.lengthen <- function(nsm = 1){
  function(x) format(x,nsmall = nsm,scientific = FALSE)
}