#Analysis of Sparse co-cultures for connectivity (SparCon)
#Author: Kirill Zaslavsky
#Purpose: Generate figures used in the paper (Zaslavsky and Zhang et al.)
#NOTE: to prevent overplotting, position_jitterdodge randomly spreads the points around 
#      in dot plots along the categorial (x-axis) to separate them. 
#      while the y-values remain the same, the x-axis position of the dots may change 
#      from one plot version to another.

### If all necessary packages are installed, this script will run without errors
### Uncomment the section immediately below to install the packages


#install & load necessary packages and functions
#uncomment if you do not have these installed
#source("https://bioconductor.org/biocLite.R")   #bioconductor to help install the packages below
#biocLite(c("data.table", "RColorBrewer", "ggplot2", "stringr", "dplyr", "scales", "ggthemes", "psych",
#           "purrr", "kSamples"))


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
                                                             #if not using R-studio, set working directory manually
                                                             #using setwd()


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
library(tidyr)
source("SparCon_fn.R")    #contains within-well normalization, plotting, and other miscellaneous functions                      

#create output directories
dir.create("Output")
dir.create("Output/Stats")
dir.create("Output/Plots")

#load data
SCdata.icc <- fread("SC_ICC_DB.txt", header=TRUE)
SCdata.ephys <- fread("SC_ephys_db.txt", header=TRUE)

SCdata.icc[, growth_condition := "AA"]

#write.table(SCdata.icc, file = "SC_ICC_DB_upd.txt", row.names = FALSE)
#write.table(SCdata.ephys, file = "SC_ephys_db_upd.txt", row.names = FALSE)

#####
#### Part 1: Data pre-processing
#### Set appropriate columns as keys, factors, normalize within-well
setkey(SCdata.icc, wellID)
setkey(SCdata.ephys, wellID)
SCdata.icc[, exp_ID := as.factor(exp_ID)]
SCdata.ephys[, exp_ID := as.factor(exp_ID)]


#Set Syn_Number to numeric; by default it is read in as an integer which causes problems with calculations
SCdata.icc[, Syn_Number := as.numeric(Syn_Number)]

#Set Amplitude to positive (simplifies calculations and plotting).
SCdata.ephys[, sEPSC_Ampl := abs(sEPSC_Ampl)]

#Normalize within well - mut to ctrl by geometric.mean
SCdata.icc.norm <- mcc.icc.norm.well(SCdata.icc)
SCdata.ephys.norm <- mcc.ephys.norm.well(SCdata.ephys)

#write.table(SCdata.icc.norm, "SC_ICC_DB_norm.txt", row.names = FALSE)
#write.table(SCdata.ephys.norm, "SC_ephys_DB_norm.txt", row.names = FALSE)

#####
# to assist with plotting and analysis, categorize to
# make factors out of specific variables and separate color controls from experimental controls

#prepare list of datasets
dat_list <- list(SCdata.icc, SCdata.icc.norm, SCdata.ephys, SCdata.ephys.norm)

#define factor ordering
levels_ind_desc <- c("CTRL1", "CTRL2", 
                     "CTRL3", "CTRL4",
                     "SHANK2 R841X",
                     "SHANK2 DEL","SHANK2 KO", "SHANK2 R841X-C"
)

levels_SHANK2_gtype <- c("WT", "CORR", "HET", "NULL" )

levels_comparison_var <- c("SHANK2 R841X","SHANK2 DEL","SHANK2 KO",
                           "SHANK2 R841X CORR"
                           
)



levels_comparison_label <- c("CTRL", "SHANK2 R841X-C",
                             "SHANK2 R841X",
                             "SHANK2 DEL","SHANK2 KO"
)

cb_list <- list(list("ind_desc", levels_ind_desc),
                list("SHANK2_gtype",levels_SHANK2_gtype),
                list("comparison_label",levels_comparison_label),
                list("comparison_var", levels_comparison_var))

#categorize all files and separate between EXP and color controls

dat_EXP <- dat_list %>% 
  map(sparc_extract, c("EXP", "MUTLAWN")) %>%
  map(sparc_categorize, fac = cb_list)

dat_CLRCTRL <- dat_list %>% 
  map(sparc_extract, "COLOR_CONTROL") %>%
  map(sparc_categorize, fac = cb_list)

dat_MUTLAWN <- dat_list %>%
  map(sparc_extract, "MUTLAWN") %>%
  map(sparc_categorize, fac = cb_list)

#if you don't want to separate
dat_ALL <- dat_list %>% 
  map(sparc_categorize, fac = cb_list)


#At this point: 
#dat_EXP[[1]] is unnormalized imaging dataset (Synapse Number & Dnedrite Length)
#dat_EXP[[2]] is normalized imaging dataset
#dat_EXP[[3]] is unnormalized electrophysiology dataset (sEPSC Frequency and Amplitude)
#dat_EXP[[4]] is normalized electrophysiology dataset

############# end of Part I ########################
## you can now run Distributions.R and Power_Sims.R
############# end of Part I ########################


#####
##### Part 2: Statistics  

setwd("Output/Stats")

#Anderson-Darling Tests for Ephys
var_list <- list("sEPSC_Freq", "sEPSC_Ampl")

unnorm.sepsc.ad.res <- ad_results(dat_EXP[[3]], var_list)
norm.sepsc.ad.res <- ad_results(dat_EXP[[4]], var_list)

write.table(unnorm.sepsc.ad.res, "ephys_adtest_unnorm.txt", row.names = FALSE, sep = "\t")
write.table(norm.sepsc.ad.res, "ephys_adtest_norm.txt", row.names = FALSE, sep = "\t")


#Anderson-Darling Tests for ICC
var_list <- list("Syn_Number", "Dend_Length")

unnorm.icc.ad.res <- ad_results(dat_EXP[[1]], var_list)
norm.icc.ad.res <- ad_results(dat_EXP[[2]], var_list)

write.table(unnorm.icc.ad.res, "icc_adtest_unnorm.txt", row.names = FALSE, sep = "\t")
write.table(norm.icc.ad.res, "icc_adtest_norm.txt", row.names = FALSE, sep = "\t")


#Synapse
#MUT          CTRL        FOLD CHANGE
#R841X 2.2    CTRL 1.15   FC 1.91
#DEL   1.79   CTRL 1.17   FC 1.55
#KO    1.86   CTRL 1.12   FC 1.67

#Dendrite length
#MUT          CTRL        FOLD CHANGE
#R841X 1.64   CTRL 1.09   FC 1.50
#DEL   1.20   CTRL 1.06   FC 1.13
#KO    1.60   CTRL 1.06   FC 1.51


#####
##### Part 3:  Generate plots
setwd("../Plots")
#####
## Well-to-well and Batch-to-Batch variation
## Synapse Number by wellID

## well-to-well variation by batch for the first 3 batches used for power analyses
## Part of Supp Fig 7
syn.wwvar.nonorm.all <- ggplot(dat_EXP[[1]][exp_ID %in% 1:3], aes(x=ind_desc, y = Syn_Number,
                                                            fill = SHANK2_gtype, 
                                                            group = SHANK2_gtype))
pl.dotplot.OUT <- syn.wwvar.nonorm.all %>% exp.plfn.geomPOINT.3()+
  labs(x="Line", y = "Synapse Number")+
  scale_y_continuous(limits = c(1,256), breaks = c(0.1, 2^(seq(0,9,1))))+
  coord_trans(y="log2")+
  facet_grid(.~exp_ID+SHANK2_gtype+wellID, switch = "x") +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    strip.text = element_text(size=8),
    legend.position = "none")+
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"),values = c("magenta","green","blue"))+
  ggtitle("Well-to-well variation")+
  theme(axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
ggsave("wwvar_syn_unnorm_all_cslip.pdf", plot = pl.dotplot.OUT, device = "pdf", width = 20, height = 8)

## Synapse Number across one line by batch, unnormalized
## Part of Supp Fig 7
syn.bbvar.nonorm.oneLine <- ggplot(dat_EXP[[1]][ind_desc == "CTRL1"], aes(x=ind_desc, y = Syn_Number,
                                                                                      group = SHANK2_gtype))
pl.dotplot.OUT <- syn.bbvar.nonorm.oneLine %>% exp.plfn.geomPOINT.2()+
  labs(x="Individual", y = "Synapse Number")+
  scale_y_continuous(limits = c(1,256), breaks = c(0.1, 2^(seq(0,9,1))))+
  coord_trans(y="log2")+
  facet_grid(.~exp_ID) +
  theme(legend.position = "none") + 
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,3,3)]) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", size = 2, fill = NA))
ggsave("bbvar_syn_unnorm_CTRL1.pdf", plot = pl.dotplot.OUT, width = 4, height = 4)

## batch-to-batch for one line, normalized
## Part of Supp Fig 7
## Synapse Number
syn.bbvar.norm.oneLine <- ggplot(dat_EXP[[2]][ind_desc == "CTRL1"], aes(x=ind_desc, y = Syn_Number,
                                                                                   
                                                                                   group = SHANK2_gtype))
pl.dotplot.OUT <- syn.bbvar.norm.oneLine %>% exp.plfn.geomPOINT.2()+
  labs(x="Individual", y = "Normalized Synapse Number")+
  scale_y_continuous(limits = c(0.2,4), breaks = c(0.1, 2^(seq(0,9,1))))+
  coord_trans(y="log2")+
  facet_grid(.~exp_ID) +
  theme(legend.position = "none") +
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,3,3)]) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", size = 2, fill = NA))
ggsave("bbvar_syn_norm_CTRL1.pdf", plot = pl.dotplot.OUT, width = 4, height = 4)


## For detailed distribution plots & power simulations, refer to Power_Sims.R - Supp Figs 8-11


#####
## Part 3b: Color Control plots
## Color Control plots
color.ctrl <- brewer.pal(8, "Accent")[c(1,6)]
sfm.ctrl <- scale_fill_manual(values = color.ctrl)
scm.ctrl <- scale_color_manual(values = color.ctrl)

##ICCdata
ICC.clrctrl.syn <- ggplot(dat_CLRCTRL[[2]], aes(x=Color, y=(Syn_Number), group = Color))
ICC.clrctrl.syn.ecdf <- ggplot(dat_CLRCTRL[[2]], aes(x=Syn_Number, group = Color, color = Color))
ICC.clrctrl.dend <- ggplot(dat_CLRCTRL[[2]], aes(x=Color, y=Dend_Length, group = Color))
ICC.clrctrl.dend.ecdf <- ggplot(dat_CLRCTRL[[2]], aes(x=Dend_Length, group = Color, color = Color))

#p1 - syn number
p1 <- ICC.clrctrl.syn %>% clr.ctrl.plfn()+
  facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0.03,8), breaks = c(0.1, 2^(seq(0,8,1))),labels = fmt.lengthen())+
  labs(y="Normalized Synapse Number")+
  coord_trans(y="log2")+
  sfm.ctrl

#p2 - dend length
p2 <- ICC.clrctrl.dend %>% clr.ctrl.plfn()+
  facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0,3.1),breaks = c(0.0, seq(0,3,1)),labels = fmt.lengthen())+
  labs(y="Normalized Dendrite Length")+
  sfm.ctrl

#p3 - summary syn number (GFP vs mKO2)
p3 <- ICC.clrctrl.syn.ecdf %>% clr.ctrl.plfn.ecdf()+
  #scale_y_continuous(limits = c(0.03,8), breaks = c(0.1, 2^(seq(0,8,1))), labels = fmt.lengthen())+
  labs(x="Normalized Synapse Number", y="Cumulative Probability")+
  theme(panel.border = element_blank(),
        axis.line = element_line(color="black", size = 1))+
  #coord_trans(y="log2")+
  sfm.ctrl+
  scm.ctrl

#p4 - summary dend legnth (GFP vs mKO2)
p4 <- ICC.clrctrl.dend.ecdf %>% clr.ctrl.plfn.ecdf()+
  #scale_y_continuous(limits = c(0,3),breaks = c(0.0, seq(0,3,1)), labels = fmt.lengthen())+
  labs(x = "Normalized Dendrite Length", y="Cumulative Probability")+
  theme(panel.border = element_blank(),
        axis.line = element_line(color="black", size = 1))+
  sfm.ctrl+
  scm.ctrl

p5 <- ICC.clrctrl.syn %>% clr.ctrl.plfn()+
  #facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0.01,8), breaks = c(0.1, 2^(seq(0,8,1))),labels = fmt.lengthen())+
  labs(y="Normalized Synapse Number")+
  coord_trans(y="log2")+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1))+
  sfm.ctrl

p6 <- ICC.clrctrl.dend %>% clr.ctrl.plfn()+
  scale_y_continuous(limits = c(0,3.1),breaks = c(0.0, seq(0,3,1)),labels = fmt.lengthen())+
  labs(y="Normalized Dendrite Length")+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1))+
  sfm.ctrl

#theme black red blue
pdf("color_ICC_controls_ecdf.pdf", height = 9, width = 12)
multiplot(p1,p2,p3,p4, cols =2)
dev.off() #color controls

pdf("clr_syn_grp.pdf", width =3, height =4)
p5
dev.off()

pdf("clr_dend_grp.pdf", width =3, height =4)
p6
dev.off()

## for ephys
ephys.clrctrl.sEPSC_Freq <- ggplot(dat_CLRCTRL[[4]], aes(x=Color, y=sEPSC_Freq, group = Color))
ephys.clrctrl.sEPSC_Freq.ecdf <- ggplot(dat_CLRCTRL[[4]], aes(x=sEPSC_Freq, group = Color, color = Color))
ephys.clrctrl.sEPSC_Ampl <- ggplot(dat_CLRCTRL[[4]], aes(x=Color, y=sEPSC_Ampl, group = Color))
ephys.clrctrl.sEPSC_Ampl.ecdf <- ggplot(dat_CLRCTRL[[4]], aes(x=sEPSC_Ampl, group = Color, color = Color))

#p1 - sEPSC_Freq number
p1 <- ephys.clrctrl.sEPSC_Freq %>% clr.ctrl.plfn()+
  facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0.03,8), breaks = c(0.1, 2^(seq(0,8,1))),labels = fmt.lengthen())+
  labs(y="Normalized sEPSC Frequency")+
  coord_trans(y="log2")+
  sfm.ctrl

#p2 - sEPSC_Ampl length
p2 <- ephys.clrctrl.sEPSC_Ampl %>% clr.ctrl.plfn()+
  facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0,3),breaks = c(0.0, seq(0,3,1)),labels = fmt.lengthen())+
  labs(y="Normalized sEPSC Amplitude")+
  sfm.ctrl

#p3 - summary sEPSC_Freq number (GFP vs mKO2)
p3 <- ephys.clrctrl.sEPSC_Freq.ecdf %>% clr.ctrl.plfn.ecdf()+
  #scale_y_continuous(limits = c(0.03,8), breaks = c(0.1, 2^(seq(0,8,1))), labels = fmt.lengthen())+
  labs(x="Normalized sEPSC Frequency", y="Cumulative Probability")+
  theme(panel.border = element_blank(),
        axis.line = element_line(color="black", size = 1))+
  #coord_trans(y="log2")+
  sfm.ctrl+
  scm.ctrl

#p4 - summary sEPSC_Ampl legnth (GFP vs mKO2)
p4 <- ephys.clrctrl.sEPSC_Ampl.ecdf %>% clr.ctrl.plfn.ecdf()+
  #scale_y_continuous(limits = c(0,3),breaks = c(0.0, seq(0,3,1)), labels = fmt.lengthen())+
  labs(x = "Normalized sEPSC Amplitude", y="Cumulative Probability")+
  theme(panel.border = element_blank(),
        axis.line = element_line(color="black", size = 1))+
  sfm.ctrl+
  scm.ctrl

p5 <- ephys.clrctrl.sEPSC_Freq %>% clr.ctrl.plfn()+
  #facet_wrap(~wellID, ncol =3)+
  scale_y_continuous(limits = c(0.01,8), breaks = c(0.1, 2^(seq(0,8,1))),labels = fmt.lengthen())+
  labs(y="Normalized sEPSC Frequency")+
  coord_trans(y="log2")+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1))+
  sfm.ctrl

p6 <- ephys.clrctrl.sEPSC_Ampl %>% clr.ctrl.plfn()+
  scale_y_continuous(limits = c(0,3),breaks = c(0.0, seq(0,3,1)),labels = fmt.lengthen())+
  labs(y="Normalized sEPSC Amplitude")+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1))+
  sfm.ctrl

pdf("color_ephys_controls_ecdf.pdf", height = 9, width = 12)
multiplot(p1,p2,p3,p4, cols =2)
dev.off() #color controls

pdf("clr_sEPSC_Freq_grp.pdf", width =3, height =4)
p5
dev.off()

pdf("clr_sEPSC_Ampl_grp.pdf", width =3, height =4)
p6
dev.off()


#####
## Part 3c: Experimental comparison plots

## Cumulative probability plots
## normalized data > Figs 3 and 4
comp_labels <- c(
  'SHANK2 KO'="KO-/-",
  'SHANK2 DEL'="DEL+/-",
  'SHANK2 R841X'="R841X+/-",
  'SHANK2 R841X CORR'="R841X-C+/+"
)

#Normalized Synapse Number and Dendrite Length
var_list <- list("Syn_Number", "Dend_Length")
var_names <- list("Normalized Synapse Number", 
                  "Normalized Dendrite Length")
scale_x_list <- list(scale_x_continuous(breaks = c(0,3,6,9), limits = c(0,12)),
                     scale_x_continuous(lim = c(0,6)))
trans_list <- list(coord_trans(), coord_trans())
pl_names <- c("ecdf_norm_syn", "ecdf_norm_dend")
plot_ecdf_facets(dat_EXP[[2]], var_list, var_names, pl_names, scale_x_list, trans_list, comp_labels = comp_labels, size = 3)

#Normalized sEPSC Frequency and Amplitude
var_list <- list("sEPSC_Freq", "sEPSC_Ampl")
var_names <- list("Normalized sEPSC Frequency", 
                  "Normalized sEPSC Amplitude")
scale_x_list <- list(scale_x_continuous(trans = "log10", breaks = c(0.1,  10, 1000), limits = c(0.01, 10000),
                                        labels = fmt.round(1)),
                     scale_x_continuous(limits = c(-0.2, 4.2)))
trans_list <- list(coord_trans(), coord_trans())
pl_names <- c("ecdf_norm_sepsc_freq",
              "ecdf_norm_sepsc_ampl")
plot_ecdf_facets(dat_EXP[[4]][growth_condition == "AA"], 
                 var_list, var_names, pl_names, scale_x_list, trans_list, comp_labels = comp_labels, margin = 0.05)

#Normalized sEPSC Freq and Ampl for BDNF experiment
pl_names <- c("ecdf_norm_sepsc_freq_bdnf",
              "ecdf_norm_sepsc_ampl_bdnf")
sfm <- scale_fill_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2)])
scm <- scale_color_manual("Genotype", breaks=c("WT","HET","NULL", "CORR"), values = color3[c(4,2)])
plot_ecdf_facets(dat_EXP[[4]][growth_condition != "AA"], 
                 var_list, var_names, pl_names, scale_x_list, comp_labels = comp_labels, sfm = sfm, scm = scm, margin = 0.05, width = 5, height = 6)


#####
#### paired group dot plots - > insets for Figs 3 and 4

#Synapse and Dendrite Length
comp_groups <- as.character(unique(dat_EXP[[2]][, comparison_var]))
lawn_groups <- as.character(unique(dat_EXP[[2]][, Lawn]))
var_list <- list("Syn_Number", "Dend_Length")
scale_y_list <- list(scale_y_continuous(limits = c(0.03,16), breaks = c(0.1, 2^(seq(0,8,1)))),
                     scale_y_continuous(limits = c(0,6), breaks = c(0, seq(0,10,2)), labels = fmt.lengthen()))
trans_list <- list(coord_trans(y="log2"), coord_trans())
color_groups <- list(c(4,2),
                     c(4,3),
                     c(4,2),
                     c(1,2),
                     c(1,2))

plot_dp(dat_EXP[[2]], var_list, lawn_groups, comp_groups, scale_y_list, trans_list)

#sEPSC Freq and Ampl
var_list <- list("sEPSC_Freq", "sEPSC_Ampl")
scale_y_list <- list(scale_y_continuous(limits = c(0.01,256), breaks = c(2^(seq(0,9,2))), labels = fmt.round(1)),
                     scale_y_continuous(limits = c(0.0,4.3), breaks = c(0, seq(0,4,1)),labels = fmt.lengthen()))

plot_dp(dat_EXP[[4]][growth_condition == "AA"], var_list, lawn_groups, comp_groups, scale_y_list, trans_list, width = 1.05)
plot_dp(dat_EXP[[4]][growth_condition != "AA"], var_list, lawn_groups, comp_groups, scale_y_list, trans_list) #for BDNF exp



#####
##unnormalized > Supp Fig 12
## 

#ICC - ecdf
var_list <- list("Syn_Number", "Dend_Length")
var_names <- list("Synapse Number", 
                  expression(bold(paste("Dendrite Length (", mu, "m)", sep = ""))))

#ecdf
scale_x_list <- list(scale_x_continuous(breaks = c(0,200,400), limits = c(0,400)),
                     scale_x_continuous(breaks = c(0,2000,4000), limits = c(0,4100)))
trans_list <- list(coord_trans(), coord_trans())
pl_names <- c("ecdf_unnorm_syn_cp",
              "ecdf_unnorm_dend_cp")
plot_ecdf_facets(dat_EXP[[1]], var_list, var_names, pl_names, scale_x_list, trans_list, comp_labels = comp_labels, size = 3)

#dotplot
scale_y_list <- list(scale_y_continuous(breaks = c(2^seq(1,9,2)), limits = c(1,530)),
                     scale_y_continuous(breaks = c(0,2000,4000), limits = c(0,4100)))
trans_list <- list(coord_trans(y="log2"), coord_trans())
plot_dp(dat_EXP[[1]], var_list, lawn_groups, comp_groups, scale_y_list, trans_list,file_mod = "unnorm",width = 1.05)


#Ephys
var_list <- list("sEPSC_Freq", "sEPSC_Ampl")
var_names <- list("sEPSC Frequency (Hz)", 
                  "sEPSC Amplitude (mV)")

#ecdf
scale_x_list <- list(scale_x_continuous(trans = "log2", breaks = c(2^seq(-5,7,6)), limits = c(0.001, 1024),
                                        labels = fmt.round(2)),
                     scale_x_continuous(breaks = c(0,20,43), limits = c(0,43)))
trans_list <- list(coord_trans(), coord_trans())
pl_names <- c("ecdf_unnorm_sepsc_freq_cp_test",
              "ecdf_unnorm_sepsc_ampl_cp_test")
plot_ecdf_facets(dat_EXP[[3]][growth_condition == "AA"], var_list, var_names, pl_names, scale_x_list, trans_list, comp_labels = comp_labels, size = 3)

#dotplot
scale_y_list <- list(scale_y_continuous(breaks = c(2^seq(0,9,3)), limits = c(0.001, 512)
                                             ),
                          scale_y_continuous(breaks = c(0,20,43), limits = c(0,43)))
trans_list <- list(coord_trans(y="log2"), coord_trans())
plot_dp(dat_EXP[[3]][growth_condition == "AA"], var_list, lawn_groups, comp_groups, scale_y_list, trans_list,file_mod = "unnorm",width = 1.05)



#####
# ZD7288 exp
setwd("../..")
ZD.db <- fread("ZD7288_exp_final.txt", header = TRUE)
ZD.db2 <- melt(ZD.db, id.vars = 1, variable.name = "NeuronID", value.name = "AP_number")
ZD.db2 <- separate(ZD.db2, "NeuronID", into = c("exp_ID", "genotype","ind_type", "condition", "neuronID"), sep="_", remove=TRUE)
ZD.db2[ind_type == "GFP", ind_type := "CTRL"]
ZD.db2[ind_type == "mKO2", ind_type := "MUT"]
ZD.db2 <- unite(ZD.db2, col = "Group", c("genotype", "condition"), sep = " ", remove = FALSE)

color_stim <- c( "#B3C6B3", "#7FC97F", "#FFE7D4","#FDAE61")

pl <- ggplot(ZD.db2, 
             aes(x=Current, y = AP_number, group = Group, color = Group)) 
pl.out <- pl %>% plfn.sholl + 
  scale_fill_manual("Group",  values = color_stim)+
  scale_color_manual("Group",  values = color_stim)+
  scale_x_continuous(limits = c(0,34), breaks = seq(-5,30,5))+
  scale_y_continuous(breaks = seq(0,10,2))+
  labs(y="AP Number", x="Current (pA)")+ facet_grid(ind_type~.)+
  
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Output/Plots/ZD7288.pdf", plot = pl.out, device = "pdf", width = 7, height = 6)

anova(lm(AP_number ~ Current + condition + Current * condition, data = ZD.db2[exp_ID == "Batch1",]))       #R841X-C
summary(aov(AP_number ~ Current + condition + Current * condition, data = ZD.db2[exp_ID == "Batch1",]))
anova(lm(AP_number ~ Current + condition + Current * condition, data = ZD.db2[exp_ID == "Batch2",]))       #R841X
summary(aov(AP_number ~ Current + condition + Current * condition, data = ZD.db2[exp_ID == "Batch2",]))


