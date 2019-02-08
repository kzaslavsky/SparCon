### Analysis for the drug mini screen in SparCon
### Experimental setup: After reseeding for co-culture, drug was applied for weeks 5-9 with every feeding
### will produce plots in paper as well as output stats
### Authors: Fraser McReady, Kirill Zaslavsky

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

library(data.table)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(stringr)
library(grid)
library(gridExtra)


dir.create("MiniScreen_Output")

#plot for within-genotype normalized data (WT treatment to WT DMSO and Mutant treatment to Mutant DMSO)
WGen_norm <- fread("Within_genotype_normalized_data.csv")[,-1]
WGen_stats <- fread("Within_genotype_normalized_stats.csv")[,-1]

WGen_norm$Treatment <- factor(WGen_norm$Treatment, levels = c("DMSO", "IGF1", "RAPA", "DHPG", "TG003", "BDNF", "EFT508"), ordered = TRUE)
WGen_stats$Treatment <- factor(WGen_stats$Treatment, levels = c("DMSO", "IGF1", "RAPA", "DHPG", "TG003", "BDNF", "EFT508"), ordered = TRUE)
WGen_stats$DMSO_Mean <- factor(WGen_stats$DMSO_Mean, levels = unique(WGen_stats$DMSO_Mean))

plotmeans <- data.frame(ind_type = WGen_stats[Treatment == "DMSO", ind_type],
                        hline = as.numeric(as.character(WGen_stats[Treatment == "DMSO", DMSO_Mean])))

comp_labels <- c(
  'R841X-Correction'="R841X-C",
  'R841X-Mutant'="R841X"
)

WGen_plot <- ggplot(WGen_norm, aes(y=Normalized, x = Treatment)) + 
  geom_point(alpha = 0.5, size = 3, aes(fill = ind_type), 
             color = "black", shape = 21, position = position_jitter(width = 0.2)) +
  facet_grid(. ~ ind_type, labeller = labeller(ind_type = comp_labels)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", color = "black", width = 0.4, lwd = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=0.7, width = 0.1) +
  geom_text(data = WGen_stats, aes(x = Treatment, y = -0.5, label = FCR, size = 2, fontface = "bold"), color = "black") +
  geom_hline(aes(yintercept = hline), data = plotmeans, linetype = 2, size = 1, color = "black") +
  scale_fill_manual(values = color3[c(1,2)]) + 
  scale_color_manual(values = color3[c(1,2)]) +
  ylim(-0.5, 7) +
  theme(
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 16, color = "black"), 
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_line(color = "black", size = 1),
    plot.title = element_text(lineheight=.8, face="bold", size = 25),
    strip.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, colour = "black", size = 1)) + 
  ylab("Total Dendrite Length \n (Normalized to Genotype DMSO)")

pdf("MiniScreen_Output/WGen_plot.pdf", width = 9, height = 4)
WGen_plot
dev.off()



#plot for all data normalized to WT DMSO condition
WT_norm <- fread("WT_normalized_data.csv")[,-1]
Mut_stats <- fread("Mutant_normalized_stats.csv")[,-1]
WT_stats <- fread("WT_normalized_stats.csv")[,-1]

WT_norm$Treatment <- gsub("-", " ", WT_norm$Treatment)
Mut_stats$Treatment <- gsub("-", " " , Mut_stats$Treatment)
WT_stats$Treatment <- gsub("-", " " , WT_stats$Treatment)

WT_norm$Treatment <- str_wrap(WT_norm$Treatment, width = 6)
Mut_stats$Treatment <- str_wrap(Mut_stats$Treatment, width = 6)
WT_stats$Treatment <- str_wrap(WT_stats$Treatment, width = 6)

WT_norm$Treatment <- factor(WT_norm$Treatment, levels = unique(WT_norm$Treatment))
Mut_stats$Treatment <- factor(Mut_stats$Treatment, levels = unique(Mut_stats$Treatment))
WT_stats$Treatment <- factor(WT_stats$Treatment, levels = unique(WT_stats$Treatment))

Mut_subset <- WT_norm[Treatment == "DMSO\n(WT)" | Treatment == "DMSO",]
Mut_subset[ind_type == "R841X-Correction", ind_type := "R841X-C"]
Mut_subset[ind_type == "R841X-Mutant", ind_type := "R841X"]
Mut_subset$Treatment <- factor(Mut_subset$Treatment, levels = unique(Mut_subset$Treatment))
Mut_subset$ind_type <- factor(Mut_subset$ind_type, levels = unique(Mut_subset$ind_type), ordered = TRUE)




Mut_plot <- ggplot(Mut_subset, aes(y=Normalized, x = ind_type)) + 
  geom_point(alpha = 0.5, size = 3, aes(fill = ind_type), 
             color = "black", shape = 21, position = position_jitter(width = 0.2)) + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", color = "black", width = 0.4, lwd = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=0.7, width = 0.1) +
#  geom_hline(aes(yintercept = mean(Mut_subset[ind_type == "R841X-C",][[7]])), linetype = 2, size = 1.7, color = "black") +
#  geom_text(data = Mut_stats, aes(x = Treatment, y = -0.1, label = FCR, size = 2, fontface = "bold")) +
  geom_text(data = WT_stats[1:2], aes(x = c("R841X-C", "R841X"), y = -0.5, label = FCR, size = 2, fontface = "bold"), 
            color = "black") +
  scale_fill_manual(values = color3[c(1,2)]) + 
  scale_color_manual(values = color3[c(1,2)]) +
  ylim(-0.5, 7) +
  theme(
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 16, color = "black"), 
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_line(color = "black", size = 1),
    plot.title = element_text(lineheight=.8, face="bold", size = 25),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, colour = "black", size = 1)) + 
  ylab("Total Dendrite Length \n (Normalized to WT)")

ggsave("MiniScreen_Output/MUT_plot.pdf", plot = Mut_plot, device = "pdf", width = 2, height = 3.8)


#combine above dotplots
grid.arrange(WGen_plot, Mut_plot, nrow=2)
Combined_dotplots <- recordPlot()
dev.off()


#Combined probability and dotplots
full_data <- fread("Mutant_normalized_full.csv")[,-1]
treatments <- unique(full_data$Treatment)[-c(1:2)]


setwd("MiniScreen_Output")
for (i in 1:length(treatments)) {
  
  #subsetting drug treatment
  plotdf <- subset(full_data, Treatment == "DMSO-(WT)" | Treatment == "DMSO" | Treatment == treatments[i])
  
  #renaming treatment groups
  plotdf[which(plotdf$Treatment == "DMSO-(WT)"),2] <- "R841X-C"
  plotdf[which(plotdf$Treatment == "DMSO"),2] <- "R841X"
  plotdf[which(plotdf$Treatment == treatments[i] & plotdf$ind_type == "R841X-Correction"),2] <- paste("R841X-C + ", treatments[i], sep = "")
  plotdf[which(plotdf$Treatment == treatments[i] & plotdf$ind_type == "R841X-Mutant"),2] <- paste("R841X + ", treatments[i], sep = "")
  
  plotdf$Treatment <- factor(plotdf$Treatment, levels = c("R841X-C",paste("R841X-C + ", treatments[i], sep = ""),"R841X",
                                                          
                                                          
                                                          paste("R841X + ", treatments[i], sep = "") 
                                                          ))
  
  #cumulative probability plot
  plot1 <- ggplot(plotdf, aes(x=Normalized, group = Treatment, fill = Treatment, color = Treatment)) +
    stat_ecdf(geom = "line", lwd = 1) +
    stat_ecdf(geom="point", alpha=0.5, size = 2.7, aes(fill = Treatment, color = Treatment), shape = 21) +
    theme_few()+
    scale_fill_manual(values = c(color3[1], "violetred3", color3[2], "darkorchid")) + 
    scale_color_manual(values = c(color3[1],"violetred3", color3[2], "darkorchid")) +
    ylab("Cumulative Probability") +
    xlab("Normalized Total Dendrite Length") +
    ggtitle(treatments[i]) +
    theme(
      axis.title.y = element_text(face = "bold", size = 16, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      strip.text = element_text(size = 16, face = "bold"),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 12, face = "bold"),
      legend.title = element_text(colour="black", size = 16, face = "bold"),
      legend.position = "none",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1))
  
  ggsave(paste0(treatments[i], "_ecdf.pdf"), plot1, device = "pdf", width = 4, height =6)
  
  
  #Dotplot
  plot2 <- ggplot(plotdf, aes(y=Normalized, x = Treatment)) + 
    geom_point(alpha = 0.5, size = 2.5, aes(fill = Treatment), color = "black", shape = 21, position = position_jitter(width = 0.2)) + 
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", color = "black", width = 0.4, lwd = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=0.7, width = 0.1) +
    geom_hline(aes(yintercept = mean(plotdf[Treatment == "R841X-C",][[7]])), linetype = 2, size = 1, color = "black") +
    scale_fill_manual(values = c(color3[1], "violetred3", color3[2], "darkorchid")) + 
    scale_color_manual(values = c(color3[1],"violetred3", color3[2], "darkorchid")) +
    ylim(-0.5, 7) +
    theme(
      axis.title.y = element_text(size = 8, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold", size = 16, color = "black"), 
      legend.position = "none",
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 25),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      #panel.background = element_rect(fill = NA, colour = "black", size = 1)) + 
      panel.background = element_blank()) + 
    ylab("Normalized Total Dendrite Length")
  
  ggsave(paste0(treatments[i], "_dp.pdf"), plot2, device = "pdf", width = 2, height =3)

}



## Stats

results_table <- data.frame(CTRL_grp_name = character(0),
                            MUT_grp_name = character(0),
                            CTRL_n = numeric(0),
                            MUT_n = numeric(0),
                            CTRL_mean = numeric(0),
                            CTRL_sem = numeric(0),
                            MUT_mean = numeric(0),
                            MUT_sem = numeric(0),
                            Fold_Change = numeric(0),
                            t.AD = numeric(0),
                            p_val = numeric(0))



#first, compare everything to WT DMSO; n.b. use treatments to advantage
#ctrl dmso vs r841x dmso
ctrl <- as.numeric(unlist(full_data[ind_type == "R841X-Correction" & Treatment == "DMSO-(WT)", Normalized]))
mut <- as.numeric(unlist(full_data[ind_type == "R841X-Mutant" & Treatment == "DMSO", Normalized]))
test_result <- ad.test(ctrl,mut)
row <- data.frame(CTRL_grp_name = "R841X-C (DMSO)",
                  MUT_grp_name = "R841X (DMSO)",
                  CTRL_n = length(ctrl),
                  MUT_n = length(mut),
                  CTRL_mean = mean(ctrl),
                  CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                  MUT_mean = mean(mut),
                  MUT_sem = sd(mut)/sqrt(length(mut)),
                  Fold_change = mean(mut) / mean(ctrl),
                  t.AD = test_result$ad[2,2],
                  p_val = test_result$ad[2,3])
results_table <- rbind(results_table, row)

#CTRL vs CTRL treated
for (i in 1:length(treatments))
{
  mut <- as.numeric(unlist(full_data[ind_type == "R841X-Correction" & Treatment == treatments[i], Normalized]))
  test_result <- ad.test(ctrl,mut)
  row <- data.frame(CTRL_grp_name = "R841X-C (DMSO)",
                    MUT_grp_name = paste0("R841X-C + ", treatments[i]),
                    CTRL_n = length(ctrl),
                    MUT_n = length(mut),
                    CTRL_mean = mean(ctrl),
                    CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                    MUT_mean = mean(mut),
                    MUT_sem = sd(mut)/sqrt(length(mut)),
                    Fold_change = mean(mut) / mean(ctrl),
                    t.AD = test_result$ad[2,2],
                    p_val = test_result$ad[2,3])
  results_table <- rbind(results_table, row)
}

#CTRL DMSO vs CTRL treated
for (i in 1:length(treatments))
{
  mut <- as.numeric(unlist(full_data[ind_type == "R841X-Mutant" & Treatment == treatments[i], Normalized]))
  test_result <- ad.test(ctrl,mut)  
  row <- data.frame(CTRL_grp_name = "R841X-C (DMSO)",
                    MUT_grp_name = paste0("R841X + ", treatments[i]),
                    CTRL_n = length(ctrl),
                    MUT_n = length(mut),
                    CTRL_mean = mean(ctrl),
                    CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                    MUT_mean = mean(mut),
                    MUT_sem = sd(mut)/sqrt(length(mut)),
                    Fold_change = mean(mut) / mean(ctrl),
                    t.AD = test_result$ad[2,2],
                    p_val = test_result$ad[2,3])
  results_table <- rbind(results_table, row)
}


#second, compare treated WT vs MUT
for (i in 1:length(treatments))
{
  ctrl <- as.numeric(unlist(full_data[ind_type == "R841X-Correction" & Treatment == treatments[i], Normalized]))
  mut <- as.numeric(unlist(full_data[ind_type == "R841X-Mutant" & Treatment == treatments[i], Normalized]))
  test_result <- ad.test(ctrl,mut)  
  row <- data.frame(CTRL_grp_name = paste0("R841X-C + ", treatments[i]),
                    MUT_grp_name = paste0("R841X + ", treatments[i]),
                    CTRL_n = length(ctrl),
                    MUT_n = length(mut),
                    CTRL_mean = mean(ctrl),
                    CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                    MUT_mean = mean(mut),
                    MUT_sem = sd(mut)/sqrt(length(mut)),
                    Fold_change = mean(mut) / mean(ctrl),
                    t.AD = test_result$ad[2,2],
                    p_val = test_result$ad[2,3])
  results_table <- rbind(results_table, row)
}

#third, compare MUT DMSO to MUT treated
ctrl <- as.numeric(unlist(full_data[ind_type == "R841X-Mutant" & Treatment == "DMSO", Normalized]))
for (i in 1:length(treatments))
{
  mut <- as.numeric(unlist(full_data[ind_type == "R841X-Mutant" & Treatment == treatments[i], Normalized]))
  test_result <- ad.test(ctrl,mut)  
  row <- data.frame(CTRL_grp_name = "R841X (DMSO)",
                    MUT_grp_name = paste0("R841X + ", treatments[i]),
                    CTRL_n = length(ctrl),
                    MUT_n = length(mut),
                    CTRL_mean = mean(ctrl),
                    CTRL_sem = sd(ctrl)/sqrt(length(ctrl)),
                    MUT_mean = mean(mut),
                    MUT_sem = sd(mut)/sqrt(length(mut)),
                    Fold_change = mean(mut) / mean(ctrl),
                    t.AD = test_result$ad[2,2],
                    p_val = test_result$ad[2,3])
  results_table <- rbind(results_table, row)
}

write.table(results_table, "miniDS_stats.txt", sep = "\t", row.names = FALSE)


