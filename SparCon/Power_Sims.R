#Power analysis of SparCon Data
#Author: Kirill Zaslavsky


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
#if not using R-studio, set working directory manually
#using setwd()

library(MASS)
library(fitdistrplus)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(kSamples)
library(ggthemes)
source("SparCon_fn.R")
source("Power_Sim_Fn.R") #has functions for simulating by N and by tau (FC in the paper)
dir.create("Power") 



setwd("Power")
dir.create("tables_byN")
dir.create("plots_testComparison")
dir.create("plots_byTau")
dir.create("tables_byTau")


#####
#### Power Simulations
####power plots - very N
#set constants
possible.ns <- seq(from=5, to=180, by=5) # The sample sizes we'll be considering 
powersN.ad <- rep(NA, length(possible.ns)) # Empty object to collect simulation estimates 
powersN.ks <- rep(NA, length(possible.ns))
powersN.t <- rep(NA, length(possible.ns))
alpha <- 0.05 # Standard significance level 
sims <- 1000 # Number of simulations to conduct for each N 
sig.list <- list()
dat_list <- rep(list(dat_EXP[[1]][exp_ID %in% 1:3], dat_EXP[[2]][exp_ID %in% 1:3], 
                     dat_EXP[[3]][growth_condition == "AA" & (exp_ID == 4 | exp_ID ==5 )], 
                     dat_EXP[[4]][growth_condition == "AA" & (exp_ID == 4 | exp_ID ==5 )]), 2)
var_list <- rep(list("Syn_Number", "sEPSC_Freq", "Dend_Length", "sEPSC_Ampl"), each = 2)
tau_list <- rep(list(1.5, 2, 1.25, 1.25),each = 2)
powersN <- list()


list_exp <- list(dat_list, var_list, tau_list)

#do the simulation and prepare for output and plotting
powersN <- pmap(list_exp, simulate_power_byN, possible.ns = possible.ns, sims = sims, alpha = 0.05)
powersN.2 <- map(powersN, cbind, possible.ns = possible.ns)            #output table
powersN.3 <- map(powersN.2, gather, Test, Power, pow.ad:pow.t)         #table for plotting
powersN.3 <- map(powersN.3, as.data.table)

###########
########use if you need to input from files
#l_names <- as.list(paste0("tables_byN/", c("unnorm_syn.txt",
#                                               "norm_syn.txt",
#                                               "unnorm_sepsc_freq.txt",
#                                              "norm_sepsc_freq.txt",
#                                               "unnorm_dend_length.txt",
#                                               "norm_dend_length.txt",
#                                               "unnorm_sepsc_ampl.txt",
#                                               "norm_sepsc_ampl.txt")))
# pwr <- map(l_names, read.table, header=TRUE)
# powersN.2 <- pwr
# powersN.3 <- map(powersN.2, gather, Test, Power, pow.ad:pow.t)         #table for plotting
# powersN.3 <- map(powersN.3, as.data.table)




#######
file_names <- as.list(paste0("tables_byN/", c("unnorm_syn.txt",
                "norm_syn.txt",
                "unnorm_sepsc_freq.txt",
                "norm_sepsc_freq.txt",
                "unnorm_dend_length.txt",
                "norm_dend_length.txt",
                "unnorm_sepsc_ampl.txt",
                "norm_sepsc_ampl.txt")))

map2(powersN.2, file_names, write.table, row.names=FALSE, sep = "\t")    #output the tables

########
###output the plots
titles <- c("Unnormalized \nSynapse Number",
            "Normalized \nSynapse Number",
            "Unnormalized \nsEPSC Frequency",
            "Normalized \nsEPSC Frequency",
            "Unnormalized \nDendrite Length",
            "Normalized \nDendrite Length",
            "Unnormalized \nsEPSC Amplitude",
            "Normalized \nsEPSC Amplitude")  

comp_titles <- c("Anderson-Darling on \n Normalized vs. Unnormalized \n Synapse Number",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n sEPSC Frequency",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n Dendrite Length",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n sEPSC Amplitude")

fnames <- paste0("plots_testComparison/", c("syn.pdf",
            "sepsc_freq.pdf",
            "dend.pdf",
            "sepsc_ampl.pdf"))

#plot loop
for (i in 1:length(fnames)){
  p1 <- ggplot(powersN.3[[(i*2)-1]], aes(x=possible.ns, y = Power, color = Test))%>%
    exp.theme()+
    theme(legend.position = "none")+
    geom_line(alpha = 0.8, lwd = 2)+
    geom_line(aes(y = 0.8), color = "black", lwd = 2) +
    labs(y = "Power", x = "Sample Size") + 
    ggtitle(titles[(i*2)-1]) +
    scale_y_continuous(breaks = c(seq(0,1,0.2)))+
    scale_x_continuous(limits = c(0,185), breaks = c(seq(0,180,30)))+
    scale_color_manual("Statistical test", labels=c("Anderson-Darling", "Kolmogorov-Smirnov", "T-test"), values = color2[c(3,1,2)])
  
  p2 <- ggplot(powersN.3[[i*2]], aes(x=possible.ns, y = Power, color = Test)) %>%
    exp.theme()+
    theme(legend.position = "none")+
    geom_line(alpha = 0.8, lwd = 2)+
    geom_line(aes(y = 0.8), color = "black", lwd = 2) +
    labs(y = "Power", x = "Sample Size") + 
    ggtitle(titles[i*2]) +
    scale_y_continuous(breaks = c(seq(0,1,0.2)))+
    scale_x_continuous(limits = c(0,185), breaks = c(seq(0,180,30)))+
    scale_color_manual("Statistical test", labels=c("Anderson-Darling", "Kolmogorov-Smirnov", "T-test"), values = color2[c(3,1,2)])
    
  
  
  p3 <- ggplot(powersN.3[[i*2]][Test == "pow.ad"], aes(x=possible.ns, y = Power)) %>% exp.theme()+
    geom_line(alpha = 0.8, lwd = 2, color = color2[3])+
    #geom_point(shape = 21, alpha = 0.8, color = "black", fill = color2[3], size = 3, stroke = 1.2)+
    geom_line(data = powersN.3[[(i*2)-1]][Test=="pow.ad"], aes(y=Power), color = color2[1], lwd=2)+
    #geom_point(data = powersN.3[[(i*2)-1]][Test=="pow.ad"],
    #           shape = 21, alpha = 0.8, color = "black", fill = color2[1], size = 3, stroke = 1.2)+
    geom_line(aes(y = 0.8), color = "black", lwd = 2) +
    labs(y = "Power", x = "Sample Size") + 
    ggtitle(comp_titles[i]) +
    scale_y_continuous(breaks = c(seq(0,1,0.2)))+
    scale_x_continuous(limits = c(0,185), breaks = c(seq(0,180,30)))+
    scale_color_manual("Statistical test", labels=c("Normalized", "Unnormalized"), values = color2[c(3,1)]) 
  
  pdf(fnames[i], width = 5.25, height = 14.75)
  multiplot(p1,p2,p3, cols =1)
  dev.off()
}



######

#### Power by Fold Change (tau)
########use if you need to input from files
# l_names <- as.list(paste0("tables_byTau/", c("unnorm_syn.txt",
  #                                         "norm_syn.txt",
  #                                         "unnorm_sepsc_freq.txt",
  #                                         "norm_sepsc_freq.txt",
  #                                         "unnorm_dend_length.txt",
  #                                         "norm_dend_length.txt",
  #                                         "unnorm_sepsc_ampl.txt",
  #                                         "norm_sepsc_ampl.txt")))
# pwr <- map(l_names, read.table, header=TRUE)
# powersT.2 <- pwr
# powersT.3 <- map(powersT.2, gather, Test, Power, pow.ad:pow.t)         #table for plotting
# powersT.3 <- map(powersT.3, as.data.table)

#power plots vary by tau, keep n defined
#i.e. what is the smallest possible fold-change at a given N that you can be powered to detect
tau <- c(seq(1,3,0.1))                      #List of possible fold changes
# tau_down <- sort(1/tau_up)                                  #mirror image; not worth doing
# tau <- c(tau_down, tau_up)
powersT.ad <- rep(NA, length(tau)) # Empty object to collect simulation estimates 
powersT.ks <- rep(NA, length(tau))
powersT.t <- rep(NA, length(tau))
alpha <- 0.05 # Standard significance level 
sims <- 1000 # Number of simulations to conduct for each tau
sig.list <- list()
dat_list <- rep(list(dat_EXP[[1]], dat_EXP[[2]], 
                     dat_EXP[[3]][growth_condition == "AA"], dat_EXP[[4]][growth_condition == "AA"]), 2)
var_list <- rep(list("Syn_Number", "sEPSC_Freq", "Dend_Length", "sEPSC_Ampl"), each = 2)
possible.ns <- rep(list(20,35, 30, 20), each = 2)  
powersN <- list()

list_exp <- list(dat_list, var_list, possible.ns)


powersT <- pmap(list_exp, simulate_power_byTau, sims = sims, alpha = 0.05, tau = tau)


powersT.2 <- map(powersT, cbind, tau = tau)                            #output table
powersT.2 <- map(powersT.2, fixTau)
powersT.3 <- map(powersT.2, gather, Test, Power, pow.ad:pow.t)         #table for plotting
powersT.3 <- map(powersT.3, as.data.table)



l_names <- as.list(paste0("tables_byTau/", c("unnorm_syn.txt",
                                       "norm_syn.txt",
                                       "unnorm_sepsc_freq.txt",
                                       "norm_sepsc_freq.txt",
                                       "unnorm_dend_length.txt",
                                       "norm_dend_length.txt",
                                       "unnorm_sepsc_ampl.txt",
                                       "norm_sepsc_ampl.txt")))

map2(powersT.2, l_names, write.table, row.names=FALSE, sep = "\t")    #output the tables


titles <- c("Unnormalized \nSynapse Number",
            "Normalized \nSynapse Number",
            "Unnormalized \nsEPSC Frequency",
            "Normalized \nsEPSC Frequency",
            "Unnormalized \nDendrite Length",
            "Normalized \nDendrite Length",
            "Unnormalized \nsEPSC Amplitude",
            "Normalized \nsEPSC Amplitude")  

comp_titles <- c("Anderson-Darling on \n Normalized vs. Unnormalized \n Synapse Number",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n sEPSC Frequency",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n Dendrite Length",
                 "Anderson-Darling on \n Normalized vs. Unnormalized \n sEPSC Amplitude")

fnames <- paste0("plots_byTau/", c("syn_Tau.pdf",
                              "sepsc_freq_Tau.pdf",
                              "dend_Tau.pdf",
                              "sepsc_ampl_Tau.pdf"))


for (i in 1:length(fnames)){
p3 <- ggplot(powersT.3[[i*2]][Test == "pow.ad"], aes(x=tau, y = Power)) %>% exp.theme()+
  geom_line(alpha = 0.8, lwd = 2, color = color2[3])+
  #geom_point(shape = 21, alpha = 0.8, color = "black", fill = color2[3], size = 3, stroke = 1.2)+
  geom_line(data = powersT.3[[(i*2)-1]][Test=="pow.ad"], aes(y=Power), color = color2[1], lwd=2)+
  #geom_point(data = powersT.3[[(i*2)-1]][Test=="pow.ad"],
  #           shape = 21, alpha = 0.8, color = "black", fill = color2[1], size = 3, stroke = 1.2)+
  geom_line(aes(y = 0.8), color = "black", lwd = 2) +
  labs(y = "Power", x = "Fold Change") + 
  #ggtitle(comp_titles[i]) +
  scale_y_continuous(breaks = c(seq(0,1,0.2)))+
  scale_x_continuous(limits = c(1,3), breaks = c(seq(1,3,0.25)))+
  scale_color_manual("Statistical test", labels=c("Normalized", "Unnormalized"), values = color2[c(3,1)]) 

ggsave(fnames[i], p3, device = "pdf", width = 5.25, height = 4)
}

#####

