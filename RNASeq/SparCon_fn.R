####SparCon Functions
#Author: Kirill Zaslavsky
#Note: some code inspired by the many snippets at StackOverflow. Thank you, SO members!

#functions to 
#1: normalize within-well for ephys and ICC data
#2: plotting functions setting the theme for most plots in the paper
#3: functions for calculating statistical t-test results
#4: miscelanneous functions

######
###1: Functions for within-well normalization
#Within-well Normalize ICC
#Will normalize to CTRL in EXP wells
#Will normalize to GFP in COLOR CONTROL wells
mcc.icc.norm.well <- function(dat)
{
  ###for ctrl means of EXP
  #get ctrl means
  mccd <- dat[, .SD]
  
  ctrl.means <- mccd[ind_type=="CTRL" & Setup == "EXP",
                     lapply(.SD, geometric.mean),
                     by=.(exp_ID, wellID),
                     .SDcols = c("Syn_Number","Dend_Length")]
  
  print(ctrl.means)
  for (i in 1:dim(ctrl.means)[1])
  {
    iWell <- ctrl.means[i, wellID]
    print(iWell)
    mccd[wellID==iWell, `:=`(Syn_Number = as.double(as.double(Syn_Number) / ctrl.means[wellID==iWell, Syn_Number]),
                         Dend_Length = Dend_Length / ctrl.means[wellID==iWell, Dend_Length])]
  }

  
  ###for GFP means of COLOR CONTROL
  ctrl.means2 <- mccd[Color=="GFP" & Setup == "COLOR_CONTROL",
                     lapply(.SD, geometric.mean),
                     by=.(exp_ID, wellID),
                     .SDcols = c("Syn_Number","Dend_Length")]
  print(ctrl.means2)
  wellIDs <- unique(mccd[Setup == "COLOR_CONTROL", wellID])
  print(wellIDs)
  
  for (i in 1:length(wellIDs))
  {
    iWell <- ctrl.means2[i, wellID]
    print(iWell)
    mccd[wellID==iWell, `:=`(Syn_Number = as.double(Syn_Number) / ctrl.means2[wellID==iWell, Syn_Number],
                         Dend_Length = Dend_Length / ctrl.means2[wellID==iWell, Dend_Length])]
  }
  
  return(mccd)
}

#Within-well Normalize ephys
#input: Unnormalized ephys dataset with sEPSC_Freq and sEPSC_Ampl values
#Will normalize to CTRL in EXP wells
#Will normalize to GFP in COLOR CONTROL wells
mcc.ephys.norm.well <- function(dat)
{
  mccd <- dat[, .SD]
  #get ctrl means
  
  mccd[, sEPSC_Ampl := abs(sEPSC_Ampl)] #can't get geometric mean of negative values
  print(mccd[1:10, sEPSC_Ampl])
  ctrl.means <- mccd[ind_type=="CTRL" & Setup == "EXP",
                     lapply(.SD, geometric.mean),
                     by=.(exp_ID, wellID),
                     .SDcols = c("sEPSC_Freq","sEPSC_Ampl")]
  
  #print(ctrl.means)
  for (i in 1:dim(ctrl.means)[1])
  {
    iWell <- ctrl.means[i, wellID]
    print(iWell)
    mccd[wellID==iWell, `:=`(sEPSC_Freq = sEPSC_Freq / ctrl.means[wellID==iWell, sEPSC_Freq],
                         sEPSC_Ampl = sEPSC_Ampl / ctrl.means[wellID==iWell, sEPSC_Ampl])]
  }
  
  ###for GFP means of COLOR CONTROL
  ctrl.means2 <- mccd[Color=="GFP" & Setup == "COLOR_CONTROL",
                      lapply(.SD, geometric.mean),
                      by=.(exp_ID, wellID),
                      .SDcols = c("sEPSC_Freq","sEPSC_Ampl")]
  print(ctrl.means2)
  wellIDs <- unique(mccd[Setup == "COLOR_CONTROL", wellID])
  print(wellIDs)
  
  for (i in 1:length(wellIDs))
  {
    
    iWell <- ctrl.means2[i, wellID]
    print(iWell)
    mccd[wellID==iWell, `:=`(sEPSC_Freq = sEPSC_Freq / ctrl.means2[wellID==iWell, sEPSC_Freq],
                         sEPSC_Ampl = sEPSC_Ampl / ctrl.means2[wellID==iWell, sEPSC_Ampl])]
  }
  
  
  return(mccd)
}


#Within-well normalize sholl
#Subtracts mean crossings of MUT cells from CTRL cells at a given radius within one well
mcc.sholl.norm.well <- function(mccd)
{
  #get ctrl means
  ctrl.means <- mccd[ind_type=="CTRL",
                     lapply(.SD, mean),
                     by=.(exp_ID, wellID, Radius),
                     .SDcols = c("Crossings")]
  
  # print(ctrl.means)
  for (i in 1:dim(ctrl.means)[1])
  {
    iWell <- ctrl.means[i, wellID]
    iRadius <- ctrl.means[i, Radius]
    #print(iWell)
    mccd[wellID==iWell & Radius == iRadius, 
         `:=`(deltaCrossings = Crossings - ctrl.means[wellID==iWell & Radius == iRadius, Crossings])]
  }

  return(mccd)
}



#plot fn for color controls
clr.ctrl.plfn <- function(gp)
{gp+
    geom_point(alpha=0.7, size = 3.5, aes(fill = Color), color = "black", shape = 21, position = position_jitter(width=0.2)) +
    stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
    theme_few()+
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
          axis.title.x = element_text(face = "bold", size = 16, margin = margin(10,0,10,0)),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.ticks.x = element_line(color = "black", size = 1),
          strip.text = element_blank(),
          panel.border = element_rect(fill=NA, color = "black", size = 2),
          axis.line.y = element_blank(),
          axis.line.x = element_blank()
#axis.line = element_blank()
    )
}



#plot ECDF for color controls
clr.ctrl.plfn.ecdf <- function(gp)
{gp+
    stat_ecdf(geom = "line", lwd= 2)+
    stat_ecdf(geom = "point", alpha=0.7, size = 3.5, aes(fill = Color), color = "black", shape = 21)+
   # geom_point(alpha=0.7, size = 3.5, aes(fill = Color), color = "black", shape = 21, position = position_jitter(width=0.2)) +
   # stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
   # stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
    theme_few()+
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
          axis.title.x = element_text(face = "bold", size = 16, margin = margin(10,0,10,0)),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.ticks.x = element_line(color = "black", size = 1),
          strip.text = element_blank(),
          panel.border = element_rect(fill=NA, color = "black", size = 1),
          axis.line.y = element_line(color = "black", size = 1),
          axis.line.x = element_line(color = "black", size = 1)
          #axis.line = element_blank()
    )
}

#


exp.plfn <- function(gp)
{gp+
    geom_point(alpha=0.5, size = 2.5, aes(fill = SHANK2_gtype), color = "black", shape = 21, 
               position = position_jitter(width=0.2)) +
    stat_summary(fun.y = mean, geom = "bar", fill = NA, color = "black", width = 0.3, size = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
    theme_few()+
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
          axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.ticks.x = element_blank(),
          strip.text = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size = 1))
  
}

exp.plfn.geomPOINT <- function(gp, alpha = 0.5, size = 2.55, jitterwidth = 0.2, scale = 1, meanwidth = 0.3,
                               sewidth = 0.1)
{gp+
    geom_point(alpha=alpha, size = size, aes(fill = SHANK2_gtype), color = "black", shape = 21, 
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

exp.plfn.geomPOINT.2 <- function(gp)
{gp+
    geom_point(alpha=0.4, size = 2.5, aes(fill = as.factor(exp_ID)), color = "black", shape = 21, 
               position = position_jitter(width=0.2)) +
    #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", color = "black",width = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
    theme_few()+
    theme(
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.ticks.x = element_blank(),
      strip.text = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14, face = "bold"),
      legend.title = element_text(colour="black", size = 16, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(hjust = 0.5))
  
}

exp.plfn.geomPOINT.3 <- function(gp)
{gp+
    geom_point(alpha=1, size = 1, aes(fill = SHANK2_gtype, color = SHANK2_gtype), shape = 21, 
               position = position_jitter(width=0.2)) +
    #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", color = "black",width = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
    theme_few()+
    theme(
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.ticks.x = element_blank(),
      strip.text = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14, face = "bold"),
      legend.title = element_text(colour="black", size = 16, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(hjust = 0.5))
  
}

exp.plfn.ecdf <- function(gp, alpha = 0.5, size = 2.55, scale = 1)
{gp+
    stat_ecdf(geom = "line", lwd = 1) +
    stat_ecdf(geom="point", alpha=alpha, size = size, aes(fill = SHANK2_gtype, color = SHANK2_gtype), shape = 21) +
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

#plotting function for sholl
plfn.sholl <- function(gp)
{gp+
    stat_summary(fun.data = mean_se, geom="errorbar", lwd = 1, alpha = 0.75)+
    stat_summary(fun.y = mean, geom = "line", lwd = 2, alpha = 0.75)+
    stat_summary(fun.y = mean, geom = "point", shape = 21, aes(fill = genotype), color = "black", size = 3, alpha = 0.75)+
    # stat_smooth(se=FALSE)+
    theme_classic() +
    theme(    
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.ticks.x = element_line(color="black", size = 1),
      strip.text = element_text(size = 16, face  = "bold"),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14, face = "bold"),
      legend.title = element_text(colour="black", size = 16, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1))
}





# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#output faceted ecdfs
plot_ecdf_facets <- function(db, var_list, var_names, pl_names, round = 0, margin=0.05)
{
  for (i in seq_along(var_list))
  {
    curr.plot <- ggplot(db, aes_string(x=var_list[[i]], 
                                       color = "SHANK2_gtype",
                                       group = "SHANK2_gtype")
    ) %>% exp.plfn.ecdf()+
      scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,1,2)])+
      scale_color_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,1,2)])+
      facet_grid(.~comparison_var)+
      labs(x =  var_names[[i]], y = "Cumulative Probability") +
      theme(legend.position = "none")+ 
      scale_x_continuous(
        limits = c(0,
                   (round(max(db[,.SD,.SDcols = var_list[[i]]])/2, round) * 2) + 
                     (round(max(db[,.SD,.SDcols = var_list[[i]]])/2, round) * 2) * margin),
        breaks = c(0,
                   round(max(db[,.SD,.SDcols = var_list[[i]]])/2,round),
                   round(max(db[,.SD,.SDcols = var_list[[i]]])/2,round) * 2),
        labels = fmt.round(0))
    plots.syn.cp[[i]] <- curr.plot
    
    ggsave(paste0(pl_names[[i]], ".pdf"), plots.syn.cp[[i]], device = "pdf",
           width = 9, height = 6)
    
    # 
  }
  return(plots.syn.cp)
}

#output dot plots by group
pl.grp <- function(db, comp_var, var_names, pl_names, round = -2, var.list)
{
  for (i in seq_along(var_list))
  {
    curr.plot <- ggplot(db[comparison_var == comp_var[i]],
                        aes_string(x="comparison_label", y = var_list[[i]], group = "SHANK2_gtype")
    ) %>% exp.plfn.geomPOINT(scale = 1.5)+
      labs(y = var_names[[i]])+
      scale_y_continuous(limits = c(0,
                                    round((max(db[,.SD,.SDcols = var_list[[i]]])+10)/2, round)*2),
                         breaks = c(0,
                                    round((max(db[,.SD,.SDcols = var_list[[i]]])+10)/2,round),
                                    round((max(db[,.SD,.SDcols = var_list[[i]]])+10)/2,round)*2))+
      #coord_trans(y="log2")+
      #coord_cartesian(ylim=c(0.1, 4))+
      scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,1,2)])+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 16),
            legend.position = "none")
    plots.pt[[i]] <- curr.plot
    
   # ggsave(paste0(pl_names[i], ".pdf"), curr.plot, device = "pdf",
    #       width = 3, height = 4)
  }
  return(plots.pt)
}



#plotting function for all batches to compare against one another
ecdf_plfn.allb <- function(db, var, lab, title, fname, scale = 1, margin = 0.3)
{
  vals <- abs(db[ind_type == "CTRL", .SD, .SDcols = var])
  min <- min(vals)
  max <- max(vals)
  print(max)
  print(max+margin)
  pl <- ggplot(db[ind_type == "CTRL"], aes_string(x=var, color = "exp_ID", group = "exp_ID"))
  
  out <- pl %>% exp.theme(scale = scale) + 
    stat_ecdf(geom="line", lwd=2, alpha = 0.6)+
    labs(y="Cumulative Probability", x = lab)+
    scale_color_manual("Batch", values = color2[c(3,1,2)])+
    scale_x_continuous(lim = c(0, max + margin))+
    theme(legend.position = "none")+
    ggtitle(title)
  
  ggsave(fname, out, device = pdf, width = 4, height = 4)
}

#plotting function for 1 batch for testing fitting
ecdf_plfn.1b <- function(db, var, lab, title, fname, scale=1)
{
  vals <- abs(db[ind_type == "CTRL", .SD, .SDcols = var])
  min <- min(vals)
  max <- max(vals)
  vals.dist <- fitdist(as.numeric(unlist(vals)), distr = "lnorm")
  #fit.curve = data.frame(a = rlnorm(10000, meanlog = vals.dist$estimate[1], sdlog=vals.dist$estimate[2]))
  pl <- ggplot(vals, aes_string(x=var))
  
  out <- pl %>% exp.theme(scale = scale) + 
#    stat_ecdf(data=fit.curve, aes(x=a), lwd=2, color="red", alpha = 0.75)+
    stat_function(geom = "line", lwd = 1, color = "red", 
                  fun = plnorm, args = list(meanlog = vals.dist$estimate[1], sdlog=vals.dist$estimate[2]))+
    stat_ecdf(geom = "point", shape = 19, color = "black", size = 3, alpha = 0.3)+
    labs(y="Cumulative Probability", x = lab)+
    scale_x_continuous(lim = c(0, max))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(fname, out, device = pdf, width = 4, height = 4)

}

#plotting function for 1 batch for testing fitting
dens_plfn.1b <- function(db, var, lab, title, fname, scale=1, bin_n = 20, margin = 0.1)
{
  vals <- abs(db[ind_type == "CTRL", .SD, .SDcols = var])
  min <- min(vals)
  max <- max(vals)
  vals.dist <- fitdist(as.numeric(unlist(vals)), distr = "lnorm")
  #fit.curve = data.frame(a = rlnorm(200000, meanlog = vals.dist$estimate[1], sdlog=vals.dist$estimate[2]))
  pl <- ggplot(vals, aes_string(x=var))
  
  out <- pl %>% exp.theme(scale = scale) + 
    #stat_density(data=fit.curve, aes(x=a), lwd=2, fill="red", alpha = 0.5)+
    stat_function(geom = "line", lwd = 1, color = "red", 
                  fun = dlnorm, args = list(meanlog = vals.dist$estimate[1], sdlog=vals.dist$estimate[2]))+
    geom_histogram(aes(y=..density..), bins = bin_n,
                   color = "black", fill = "black", size = 0.2, alpha = 0.2, hjust = 0)+
    labs(y="Density", x = lab)+
    scale_x_continuous(lim = c(0, max+margin))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(fname, out, device = pdf, width = 4, height = 4)
  
  return(out)
}

##########
#####3: Functions for outputting resutls of statistical tests
#output a neat AD-test table by comparison_var
#input: db, var_list
ad_results <- function(db, var_list)
{
  comp_groups <- as.character(unique(db[, comparison_var]))
  results_table <- data.frame(CTRL_grp_name = character(0),
                              MUT_grp_name = character(0),
                              CTRL_n = numeric(0),
                              MUT_n = numeric(0),
                              Variable = character(0),
                              CTRL_mean = numeric(0),
                              CTRL_sem = numeric(0),
                              MUT_mean = numeric(0),
                              MUT_sem = numeric(0),
                              Fold_Change = numeric(0),
                              t.AD = numeric(0),
                              p_val = numeric(0)
  )
  
  for (i in seq_along(comp_groups))
  {
    comp_group_var <- comp_groups[i]
    for (j in seq_along(var_list))
    {
      ctrl <- unlist(db[comparison_var == comp_group_var & ind_type == "CTRL", .SD, .SDcols = var_list[[j]] ])
      mut <- unlist(db[comparison_var == comp_group_var & ind_type == "MUT", .SD, .SDcols = var_list[[j]] ])
      test_result <- ad.test(ctrl,mut)
      row <- data.frame(CTRL_grp_name = (unique(db[comparison_var == comp_group_var & ind_type == "CTRL", (comparison_label)])),
                        MUT_grp_name = (unique(db[comparison_var == comp_group_var & ind_type == "MUT", (comparison_label)])),
                        CTRL_n = length(ctrl),
                        MUT_n = length(mut),
                        Variable = var_list[[j]],
                        CTRL_mean = mean(ctrl),
                        CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                        MUT_mean = mean(mut),
                        MUT_sem = sd(mut)/sqrt(length(mut)),
                        Fold_change = mean(mut) / mean(ctrl),
                        t.AD = test_result$ad[2,2],
                        p_val = test_result$ad[2,3]
      )
      results_table <- rbind(results_table, row)
    }
  }
  return(results_table)
}

#output a neat t-test table by comparison_var
#input: db, var_list
ttest_results <- function(db, var_list)
{
  comp_groups <- as.character(unique(db[, comparison_var]))
  results_table <- data.frame(CTRL_grp_name = character(0),
                              MUT_grp_name = character(0),
                              CTRL_n = numeric(0),
                              MUT_n = numeric(0),
                              Variable = character(0),
                              CTRL_mean = numeric(0),
                              CTRL_sem = numeric(0),
                              MUT_mean = numeric(0),
                              MUT_sem = numeric(0),
                              CTRL_estimate = numeric(0),
                              MUT_estimate = numeric(0),
                              CI_lower = numeric(0),
                              CI_upper = numeric(0),
                              df = numeric(0),
                              t = numeric(0),
                              p_val = numeric(0)
  )
  
  for (i in seq_along(comp_groups))
  {
    comp_group_var <- comp_groups[i]
    for (j in seq_along(var_list))
    {
      ctrl <- db[comparison_var == comp_group_var & ind_type == "CTRL", .SD, .SDcols = var_list[[j]] ]
      mut <- db[comparison_var == comp_group_var & ind_type == "MUT", .SD, .SDcols = var_list[[j]] ]
      test_result <- t.test(ctrl,mut)
      row <- data.frame(CTRL_grp_name = (unique(db[comparison_var == comp_group_var & ind_type == "CTRL", (comparison_label)])),
                        MUT_grp_name = (unique(db[comparison_var == comp_group_var & ind_type == "MUT", (comparison_label)])),
                        CTRL_n = length(ctrl),
                        MUT_n = length(mut),
                        Variable = var_list[[j]],
                        CTRL_mean = mean(ctrl),
                        CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                        MUT_mean = mean(mut),
                        MUT_sem = sd(mut)/sqrt(length(mut)),
                        CTRL_estimate = test_result[["estimate"]][[1]],
                        MUT_estimate = test_result[["estimate"]][[2]],
                        CI_lower = test_result[["conf.int"]][[1]],
                        CI_upper = test_result[["conf.int"]][[2]],
                        df = test_result[["parameter"]],
                        t = test_result[["statistic"]],
                        p_val = test_result[["p.value"]]
      )
      results_table <- rbind(results_table, row)
    }
  }
  return(results_table)
}

###########
###4: Miscellaneous functions
#Set color scheme
#magenta green blue
color2 <- brewer.pal(8, "Accent")[c(1,5,6)]
sfm2 <- scale_fill_manual(values = color2)


#Ensure at least one digit to the right of a decimal point
#Used in plotting with ggplot2 to make sure y-axis had the same margin
fmt.lengthen <- function(nsm = 1){
  function(x) format(x,nsmall = nsm,scientific = FALSE)
}

#format axis labels by rounding the number
fmt.round <- function(rnd = 1){
  function(x) round(x,digits = rnd)
}

#function returns substring of string from specific charcter towards the end of the string
#input: string, character
#output: string from beginning of character to end
char2end <- function(str, char){
  charLoc <- gregexpr(pattern=char, str)
  outStr <- substr(str,charLoc[[1]][1]+1, nchar(as.character(str))) 
  return(outStr)
}

#capitalize first letter of every word in a string
#useful when parsing
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#calculate power assuming normal dist
#input: mean of treatment group, mean of control group, standard deviation, signifiance level, sample size
power_calculator <- function(mu_t, mu_c, sigma, alpha=0.05, N){ 
  lowertail <- (abs(mu_t - mu_c)*sqrt(N))/(2*sigma) 
  uppertail <- -1*lowertail 
  beta <- pnorm(lowertail- qnorm(1-alpha/2), lower.tail=TRUE) + 1- pnorm(uppertail- qnorm(1-alpha/2), lower.tail=FALSE) 
  return(beta) 
} 


#function to name specific columns as factors in the datasets for 
#sparse seeding co-culture
#purpose: to facilitate data visualization
sparc_categorize <- function (dat, fac)
{
  for (i in seq_along(fac))
  {
    dat[[fac[[i]][[1]]]] = factor(dat[[fac[[i]][[1]]]], levels = fac[[i]][[2]], ordered = TRUE)
  }
  return(dat)
}


#function to extract specific data subset from db
#subset must be one of the values of setup, default is EXP, for experimental
sparc_extract <- function(dat, subset = "EXP")
{
  dat[Setup == subset]
}
