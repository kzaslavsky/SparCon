#### Generic & very basic SPARC Analysis Pipeline
#### You can use this to analyze your own data
#### Sample data is provided in 'sampledata.txt'

#### Input: a table of the following format:
#### wellID - Group - Measurement
#### Output:
#### 1. Within-well normalized dataset
#### 2. Statistical comparison of the data using Anderson-Darling test
#### 3. Cumulative distribution plot
#### 4. Dot plot

#install & load necessary packages and functions
#uncomment if you do not have these installed
#source("https://bioconductor.org/biocLite.R")   #bioconductor to help install the packages below
#biocLite(c("data.table", "RColorBrewer", "ggplot2", "stringr", "dplyr", "scales", "ggthemes", "psych",
#           "purrr", "kSamples"))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
                                                             #if not using R-studio, set working directory manually
                                                             #using setwd()

#load packages
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(psych)
library(purrr)
library(kSamples)
source("Generic_SparConFN.R")    #contains within-well normalization, plotting, and other miscellaneous functions                      

#create output directories
dir.create("Sample_Output")
dir.create("Sample_Output/Stats")
dir.create("Sample_Output/Plots")

#read file & pre-process data
db <- fread("sampleData_mutantDend.txt", header = TRUE)
setkey(db, wellID)
db[, Measurement := as.numeric(abs(Measurement))]

#normalize & output normalized data table
db.norm <- norm.well(db)
write.table(db.norm, file = "sampleData_normalized.txt", row.names = FALSE, sep = "\t")

#stats
results <- ad_results(db.norm, "Measurement")
write.table(results, "Sample_Output/Stats/stats.txt", row.names = FALSE, sep = "\t")

#plots
color_sc <- c(brewer.pal(11, "RdYlGn")[4],brewer.pal(8, "Accent")[1]) 

#ecdf
pl <- ggplot(db.norm, aes(x=Measurement, group = ind_type, fill = ind_type, color = ind_type))
pdf("Sample_Output/Plots/ecdf.pdf", width = 4, height = 5)
pl %>% exp.plfn.ecdf() + labs(y = "Cumulative Probability", x = "Normalized Dendrite Length")+
# coord_trans(y="log2")+                                                    #uncomment to specify to log-axis for y
  scale_fill_manual("Genotype", breaks = c("CTRL", "MUT"), values = color_sc)+
  scale_color_manual(breaks = c("CTRL", "MUT"), values = color_sc)+
  theme(legend.position = "none") +
 scale_x_continuous(limits = c(-0.2, 10.0), breaks = seq(0,10,2.5))                                 #uncomment to specify linear axis limits
# scale_x_continuous(trans = "log10", breaks = c(0.1,  1, 10, 100), limits = c(0.01, 1000))
#                                                                           #uncomment to specify log-axis limits

dev.off()

#dotplot
pl <- ggplot(db.norm, aes(y=Measurement, x = ind_type))
pdf("Sample_Output/Plots/dotplot.pdf", width = 2.5, height = 5)
pl %>% exp.plfn.geomPOINT() +
#scale_y_continuous(limits = c(0.01,256), breaks = c(0.1, 2^(seq(0,9,2))))+ #uncomment to specify log-axis limits
#coord_trans(y="log2")+                                                     #uncomment to specify log-axis
  scale_y_continuous(limits = c(0.0,5), breaks = c(0, seq(0,5,1)))+        #uncomment to specify linear axis limits
  scale_fill_manual("Genotype", breaks = c("CTRL", "MUT"), values = color_sc) + 
  labs (y = "Normalized Dendrite Length") +                                                 #change 'Measurement' to desired label
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
dev.off()