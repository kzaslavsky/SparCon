#Output normalized Sholl
#Author: Kirill Zaslavsky


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
#if not using R-studio, set working directory manually
#using setwd()

#algorithm to do normalized sholl
library(tidyr)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(dplyr)
library(WGCNA)
library(scales)
library(ggthemes)
source("SparCon_fn.R")
dir.create("Sholl_Plots")

#read unnorm data
sholldb <- fread("shollDB_unnorm.txt", header = TRUE)
#sholldb.norm <- fread("shollDB_norm.txt", header = TRUE)

#set factor orders to facilitate plotting
sholldb$genotype = with(sholldb, factor(genotype, 
                                          levels = c("WT", "HET", "NULL", "CORR"), ordered = TRUE))

sholldb$fam_desc_paper = with(sholldb, factor(fam_desc_paper, 
                                                levels = c("SHANK2 R841X", "SHANK2 DEL", "SHANK2 KO", "SHANK2 R841X CORR")))

sholldb$comparison_var = with(sholldb, factor(comparison_var, 
                                                levels = c("SHANK2 R841X", "SHANK2 DEL", "SHANK2 KO", "SHANK2 R841X CORR"),
                                                ordered = TRUE))

#normalize
sholldb.norm <- mcc.sholl.norm.well(sholldb) 
#write.table(sholldb.norm, "shollDB_norm.txt", sep = "\t", row.names = FALSE)

#plot
setwd("Sholl_Plots")


comp_labels <- c(
  'SHANK2 KO'="KO-/-",
  'SHANK2 DEL'="DEL+/-",
  'SHANK2 R841X'="R841X+/-",
  'SHANK2 R841X CORR'="R841X-C+/+"
)


#standardized figure containing everything
pl <- ggplot(sholldb.norm[], aes(x=Radius, y = deltaCrossings, group = genotype, color = genotype))
pl.out <- pl %>% plfn.sholl + 
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  scale_color_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  coord_fixed(ratio = 200) + 
  coord_cartesian(xlim=c(0,780), ylim=c(-0.5,3)) +
  scale_y_continuous(limits = c(-100,100), breaks = c(0, seq(0,5,1)))+
  labs(y="Change in Crossings", x=expression(paste("Distance from soma (", mu, "m)")))+
  theme(legend.position = "none")+
  facet_grid(.~comparison_var + Lawn, labeller = labeller(comparison_var=comp_labels))
ggsave("sholl-norm_comparison_var_ALL.pdf", plot = pl.out, device = "pdf", width = 8.5, height = 4)

#top part with R841X, DEL
pl <- ggplot(sholldb.norm[comparison_var == "SHANK2 R841X" | comparison_var == "SHANK2 DEL"], 
             aes(x=Radius, y = deltaCrossings, group = genotype, color = genotype))
pl.out <- pl %>% plfn.sholl + 
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  scale_color_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  coord_fixed(ratio = 200) + 
  coord_cartesian(xlim=c(0,780), ylim=c(-0.5,3)) +
  scale_y_continuous(limits = c(-100,100), breaks = c(0, seq(0,5,1)))+
  labs(y="Change in Crossings", x=expression(paste("Distance from soma (", mu, "m)")))+
  theme(legend.position = "none")+
  facet_grid(.~comparison_var+Lawn, labeller = labeller(comparison_var=comp_labels))
ggsave("sholl-norm_comparison_var_top-half.pdf", plot = pl.out, device = "pdf", width = 5.5, height = 4)

#bottom half
pl <- ggplot(sholldb.norm[comparison_var == "SHANK2 R841X CORR" | comparison_var == "SHANK2 KO"], 
             aes(x=Radius, y = deltaCrossings, group = genotype, color = genotype))
pl.out <- pl %>% plfn.sholl + 
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  scale_color_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2,3,1)])+
  coord_fixed(ratio = 200) + 
  coord_cartesian(xlim=c(0,780), ylim=c(-0.5,3)) +
  scale_y_continuous(limits = c(-100,100), breaks = c(0, seq(0,5,1)))+
  labs(y="Change in Crossings", x=expression(paste("Distance from soma (", mu, "m)")))+
  theme(legend.position = "none")+
  facet_grid(.~comparison_var+Lawn, labeller = labeller(comparison_var=comp_labels))
ggsave("sholl-norm_comparison_var_bottom-half.pdf", plot = pl.out, device = "pdf", width = 8.5, height = 4)


##statistics
sholldb.norm$Radius <- as.factor(sholldb$Radius)
sholl.subsetKO <- sholldb.norm[comparison_var == "SHANK2 KO"]
sholl.subset441 <- sholldb.norm[comparison_var == "SHANK2 R841X"]
sholl.subset217 <- sholldb.norm[comparison_var == "SHANK2 DEL"]
sholl.subset441CORRWT <- sholldb.norm[comparison_var == "SHANK2 R841X CORR" & Lawn == "CTRL1"]
sholl.subset441CORRMUT <- sholldb.norm[comparison_var == "SHANK2 R841X CORR" & Lawn == "R841X +/-"]

#two-way anova to compare distributions
anova(lm(Crossings ~ ind_type + Radius + ind_type * Radius, data = sholl.subsetKO))   #for SHANK2 KO
anova(lm(Crossings ~ ind_type + Radius + ind_type * Radius, data = sholl.subset217))  #for SHANK2 DEL
anova(lm(Crossings ~ ind_type + Radius + ind_type * Radius, data = sholl.subset441))  #for SHANK2 R841X
anova(lm(Crossings ~ ind_type + Radius + ind_type * Radius, data = sholl.subset441CORRWT))  #for SHANK2 R841X
anova(lm(Crossings ~ ind_type + Radius + ind_type * Radius, data = sholl.subset441CORRMUT))  #for SHANK2 R841X





