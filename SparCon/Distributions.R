### Characterizing data distributions of SparCon Data
### Author: Kirill Zaslavsky


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
#if not using R-studio, set working directory manually
#using setwd()

library(fitdistrplus)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(kSamples)
library(ggthemes)
source("SparCon_fn.R")
dir.create("Data Distributions")
setwd("Data Distributions")

### Comparing fits based on AIC for unnormalized data
getAICtable <- function(db, db_title="test", AIC.df = AIC.df)
{
  fit.lnorm <- fitdist(db, distr = "lnorm")
  fit.norm <- fitdist(db, distr = "norm")
  fit.gamma <- fitdist(db, distr = "gamma")
  fit.weibull <- fitdist(db, distr = "weibull")
  
  dt.tobind <- (data.table("Data" = db_title, 
           "LogNormal_AIC" = fit.lnorm$aic, 
           "Normal_AIC" = fit.norm$aic,
           "Gamma_AIC" = fit.gamma$aic,
           "Weibull_AIC" = fit.weibull$aic))
  
  return(rbind(AIC.df, dt.tobind))
}


l_db <- list(dat_EXP[[1]][ind_type == "CTRL", Syn_Number],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 1, Syn_Number],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 2, Syn_Number],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 3, Syn_Number],
             dat_EXP[[1]][ind_type == "CTRL", Dend_Length],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 1, Dend_Length],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 2, Dend_Length],
             dat_EXP[[1]][ind_type == "CTRL" & exp_ID == 3, Dend_Length],
             dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA"  & (exp_ID == 4 | exp_ID ==5 ), sEPSC_Freq],
             abs(dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA"  & (exp_ID == 4 | exp_ID ==5 ), sEPSC_Ampl]),
             dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID==4, sEPSC_Freq],
             abs(dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID==4, sEPSC_Ampl]),
             dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID==5, sEPSC_Freq],
             abs(dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID == 5, sEPSC_Ampl])
             
)

l_dbname <- list("Unnormalized Synapse Number ALL BATCHES",
                 "Unnormalized Synapse Number BATCH 1",
                 "Unnormalized Synapse Number BATCH 2",
                 "Unnormalized Synapse Number BATCH 3",
                 "Unnormalized Dendrite Length ALL BATCHES",
                 "Unnormalized Dendrite Length BATCH 1",
                 "Unnormalized Dendrite Length BATCH 2",
                 "Unnormalized Dendrite Length BATCH 3",
                 "Unnormalized sEPSC Frequency ALL BATCHES",
                 "Unnormalized sEPSC Amplitude ALL BATCHES",
                 "Unnormalized sEPSC Frequency BATCH 4",
                 "Unnormalized sEPSC Amplitude BATCH 4",
                 "Unnormalized sEPSC Frequency BATCH 5",
                 "Unnormalized sEPSC Amplitude BATCH 5"
)



AIC.df <- data.frame( "Data" = character(),
                      "LogNormal_AIC" = double(),
                      "Normal_AIC" = double(),
                      "Gamma_AIC" = double(),
                      "Weibull_AIC" = double())

for (i in 1:length(l_db))
{
  AIC.df <- getAICtable(l_db[[i]], l_dbname[[i]], AIC.df)
}

#conclusion: LogNormal best for most, Gamma 2nd best

write.table(AIC.df, file = "Distribution_Fitting_AIC.txt", sep = "\t")





#Fitting Unnormalized data

#unnorm Synapse Number
dens_plfn.1b(dat_EXP[[1]][exp_ID == 1], var = "Syn_Number", lab = "Synapse Number", title ="Synapse Number \n in Batch 1 in CTRL",
             fname = "UnNorm_Synapse.pdf", scale = 0.75, bin_n = 25, margin = 0.5)


#unnorm Dendrite Length
dend_name <- expression(bold(paste("Dendrite Length (", mu, "m)", sep = "")))
dens_plfn.1b(dat_EXP[[1]][exp_ID == 1], var = "Dend_Length", lab = dend_name, 
                  title ="Dendrite Length \n in Batch 1 in CTRL",
                  fname = "UnNorm_DendriteLength.pdf", scale = 0.75, bin_n = 25,margin = 0.5)

#unnorm sEPSC freq
dens_plfn.1b(dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID == 5], 
             var = "sEPSC_Freq", lab = "sEPSC Frequency (Hz)", title ="sEPSC Frequency \n in Batch 5 in CTRL",
             fname = "UnNorm_sEPSC_Freq.pdf", scale = 0.75, bin_n = 25, margin = 0.5)


#unnorm Amplitude
dens_plfn.1b(dat_EXP[[3]][ind_type == "CTRL" & growth_condition == "AA" & exp_ID == 5], 
                  var = "sEPSC_Ampl", lab = "sEPSC Amplitude (pA)", title ="sEPSC Amplitude \n in Batch 5 in CTRL",
                  fname = "UnNorm_sEPSC_Ampl.pdf", scale = 0.75, bin_n = 25,margin = 2)



###for one batch
#compare fits

dat_list <- rep(list(dat_EXP[[1]][as.numeric(exp_ID) == 1], 
                     dat_EXP[[2]][as.numeric(exp_ID) == 1], 
                     dat_EXP[[3]][growth_condition == "AA" & exp_ID == 5], 
                     dat_EXP[[4]][growth_condition == "AA" & exp_ID == 5]), 
                2)
var_list <- rep(list("Syn_Number", "sEPSC_Freq", "Dend_Length", "sEPSC_Ampl"), each = 2)

x_list <- list("Synapse Number ", 
               "Normalized Synapse Number",
               "sEPSC Frequency (Hz)",
               "Normalized sEPSC Frequency",
               dend_name, 
               "Normalized Dendrite Length",
               "sEPSC Amplitude (pA)",
               "Normalized sEPSC Amplitude")

title_list <- list("Synapse Number \nin Batch 1 in CTRL", 
                   "Normalized Synapse Number\nin Batch 1 in CTRL",
                   "sEPSC Frequency\nin Batch 5 in CTRL",
                   "Normalized sEPSC Frequency\nin Batch 5 in CTRL",
                   "Dendrite Length\nin Batch 1 in CTRL", 
                   "Normalized Dendrite Length\nin Batch 1 in CTRL",
                   "sEPSC Amplitude\nin Batch 5 in CTRL",
                   "Normalized sEPSC Amplitude\nin Batch 5 in CTRL")


fnames <- as.list(paste0("1batch/", c("unnorm_syn_ecdf_1b.pdf",
                                      "norm_syn_ecdf_1b.pdf",
                                      "unnorm_sEPSC_Freq_1b.pdf",
                                      "norm_sEPSC_Freq_1b.pdf",
                                      "unnorm_Dend_ecdf_1b.pdf",
                                      "norm_Dend_ecdf_1b.pdf",
                                      "unnorm_sEPSC_Ampl_1b.pdf",
                                      "norm_sEPSC_Ampl_1b.pdf"
)))

list_exp <- list(dat_list, var_list, x_list, title_list, fnames)

dir.create("1batch")
pmap(list_exp, ecdf_plfn.1b, scale = 0.75)


###for all batches
#compare batches between one another
dat_list <- rep(list(dat_EXP[[1]][exp_ID %in% 1:3], 
                     dat_EXP[[2]][exp_ID %in% 1:3], 
                     dat_EXP[[3]][growth_condition == "AA" & exp_ID %in% 4:5],
                     dat_EXP[[4]][growth_condition == "AA" & exp_ID %in% 4:5]), 2)
var_list <- rep(list("Syn_Number", "sEPSC_Freq", "Dend_Length", "sEPSC_Ampl"), each = 2)

x_list <- list("Synapse Number ", 
               "Normalized Synapse Number",
               "sEPSC Frequency",
               "Normalized sEPSC Frequency",
               dend_name, 
               "Normalized Dendrite Length",
               "sEPSC Amplitude (pA)",
               "Normalized sEPSC Amplitude")
title_list <- list("Synapse Number \nacross all batches in CTRL", 
                   "Normalized Synapse Number\nacross all batches in CTRL",
                   "sEPSC Frequency\nacross all batches in CTRL",
                   "Normalized sEPSC Frequency\nacross all batches in CTRL",
                   "Dendrite Length\nacross all batches in CTRL", 
                   "Normalized Dendrite Length\nacross all batches in CTRL",
                   "sEPSC Amplitude\nacross all batches in CTRL",
                   "Normalized sEPSC Amplitude\nacross all batches in CTRL")


fnames <- as.list(paste0("allbatches/", c("unnorm_syn_ecdf_allb.pdf",
                                          "norm_syn_ecdf_allb.pdf",
                                          "unnorm_sEPSC_Freq_allb.pdf",
                                          "norm_sEPSC_Freq_allb.pdf",
                                          "unnorm_Dend_ecdf_allb.pdf",
                                          "norm_Dend_ecdf_allb.pdf",
                                          "unnorm_sEPSC_Ampl_allb.pdf",
                                          "norm_sEPSC_Ampl_allb.pdf"
)))

list_exp <- list(dat_list, var_list, x_list, title_list, fnames)

dir.create("allbatches")
pmap(list_exp, ecdf_plfn.allb, scale = 0.75, margin = 0.05)