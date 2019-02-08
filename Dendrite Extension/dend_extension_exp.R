#Dendrite Extension Experiment
#Author: Kirill Zaslavsky


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)




library(data.table)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(dplyr)
library(scales)
library(ggthemes)
library(psych)
library(purrr)
library(kSamples)
source("SparCon_fn.R") 

dir.create("Dendrite_Output")

dend_extend_baseline <- fread("Dend_extend_baseline.txt", header = TRUE)
dend_extend_normalized <- fread("Dend_extend_normalize.txt", header = TRUE)

dend_extend_baseline.2 <- gather(dend_extend_baseline, key = "Group", value = "Length", 1:2,na.rm=TRUE)
dend_extend_normalized.2 <- gather(dend_extend_normalized, key = "Group", value = "Length", 1:4, na.rm = TRUE)
dend_extend_normalized.2$Group <- factor(dend_extend_normalized.2$Group,
                                         levels = c("CTRL-NaCl", "CTRL-KCl", "R841X-NaCl", "R841X-KCl"),
                                         ordered = TRUE)

setwd("Dendrite_Output")
pl.baseline <- ggplot(dend_extend_baseline.2, aes(x=Group, y = Length, group = Group)) + 
  geom_point(alpha=0.75, size = 3, aes(fill = Group), shape = 21, 
                         position = position_jitter(width=0.2)) +
  #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", color = "black",width = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
  theme_classic()+
  theme(
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(0,10,0,10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.x = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    legend.background = element_rect(fill=NA, color = "black"),
    legend.text = element_text(colour="black", size = 14, face = "bold"),
    legend.title = element_text(colour="black", size = 16, face = "bold"),
    legend.position = "none",
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = 1.5, face = "bold", size = 12))+
  labs(y=expression(bold(paste("Dendrite Length (", mu, "m)", sep = ""))), x = "Group")+
  ggtitle("Starting Dendrite Length")+
  scale_y_continuous(limits = c(0,3000))+
  scale_fill_manual(breaks = c("CTRL", "R841X"), values = color3[c(4,2)])
ggsave("dend_extend_baseline.pdf", pl.baseline, device = "pdf", width = 2.75, height = 3.5)

color_stim <- c("#E5BAD2","#F0027F", "#B3C6B3", "#7FC97F")

pl.normalized <- ggplot(dend_extend_normalized.2, aes(x=Group, y = Length, group = Group)) + 
  geom_point(alpha=0.65, size = 3, aes(fill = Group), shape = 21, 
             position = position_jitter(width=0.2)) +
  #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", color = "black",width = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1)+
  theme_classic()+
  theme(
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(0,10,0,10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.x = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    legend.background = element_rect(fill=NA, color = "black"),
    legend.text = element_text(colour="black", size = 14, face = "bold"),
    legend.title = element_text(colour="black", size = 16, face = "bold"),
    legend.position = "none",
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12))+
  labs(y="Relative Dendrite Length", x = "Group")+
  ggtitle("Dendrite Length 500 min Post-Stim")+
  scale_y_continuous(limits = c(0,2.5))+
  #scale_fill_manual(breaks = c("CTRL-KCl", "R841X-KCl", "CTRL-NaCl", "R841X-NaCl"), values = color3[c(4,2, 3,1)])
  scale_fill_manual(breaks = c("CTRL-KCl", "R841X-KCl", "CTRL-NaCl", "R841X-NaCl"), values = color_stim)
ggsave("dend_extend_normalized.pdf", pl.normalized, device = "pdf", width = 4.5, height = 4.5)

