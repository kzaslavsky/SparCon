#Analysis of  Activity-Dependent Signaling and ASD modules for RNASeq
#requires deseq2_rnaseq_kz.R to be run if there are no deseq2 output files
#Authors: Kirill Zaslavsky, Marat Mufteev


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
library(gplots)
library(dplyr)
library(tidyr)
library(VennDiagram)
library("GO.db")
library("AnnotationDbi")
source("GOfnKZ.R")
source("SparCon_fn.R")
source("RNASeq_FN.R")



source("http://bioconductor.org/biocLite.R") 



dir.create("RNASeq_Output")

######

###### input processed data from own DESeq output
data.4wk.proc <- fread("Output_DESeq2/data.4v4.txt")
data.9wk.proc <- fread("Output_DESeq2/data.9v9.txt")
data.counts <- fread("Output_DESeq2/norm_counts_ALL.txt")
data.counts.4wk <- data.counts[,1:9]
data.counts.9wk <- data.counts[,c(1,10:17)]

#####
#set constants
threshP = 0.05
threshLFC = log2(1.25)

#add significance column for coloring
data.4wk.proc[, Sig := ifelse(padj < threshP & abs(log2FoldChange) > threshLFC, TRUE, FALSE) ]
data.9wk.proc[, Sig := ifelse(padj < threshP & abs(log2FoldChange) > threshLFC, TRUE, FALSE) ]


#####


#####

#fix colnames in counts for mut
colnames(data.4wk.counts)[2:5] <- c("5_mut_B1_R1_4wk", "5_mut_B1_R2_4wk", "5_mut_B2_R1_4wk", "5_mut_B2_R2_4wk")
colnames(data.9wk.counts)[2:5] <- c("5_mut_B1_R1_9wk", "5_mut_B1_R2_9wk", "5_mut_B2_R1_9wk", "5_mut_B2_R2_9wk")

#####

#####
## activity-dependent signaling
genes.ads <- fread("Pruunslid_Bic_Stim_7wk.txt")
geneBox.adsignal <- data.9wk.proc[NAME %in% genes.ads$NAME,]

geneBox.adsignal.out <- ggplot(geneBox.adsignal, aes(x=reorder(NAME, -log2FoldChange), y = -log2FoldChange,
                                                           fill=padj)) %>% exp.theme() + 
  geom_bar(stat = "identity", color = "black") +
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 12, angle = 90),
        plot.title = element_text(size = 30, color = "black"),
        legend.position = "right")+
  labs(x="Gene Name", y = "log2 Fold Change", fill = "Count") +
  ggtitle("Activity-Dependent Genes") +
  #scale_fill_viridis(option = "magma", direction = -1, trans = "log")
  scale_fill_gradientn(limits=c(1e-16, 100), trans = "log",
                       breaks = c(10^seq(-16,0,3)), colors = cols.25) + labs(fill = "P-value")
  #scale_y_continuous(limits = lims[[i]]) +
  

ggsave("RNASeq_Output/adsignal_genebox.pdf", plot = geneBox.adsignal.out, device = "pdf", width = 12, height = 6)



#####
### ASD-related gene sets - Parikshak and Willsey modules
#input

ASD.modules <- fread("ASD_Modules.txt") #parikshak + willsey
ASD.modules <- ASD.modules[Module_Label != "-",]

threshP <- 0.05
threshLFC <- log(1)

#hyperGtest
#9wk
universe.num <- dim(data.9wk.proc)[1]
sig.genes.all <- dim(data.9wk.proc[padj < threshP & abs(log2FoldChange) > threshLFC, ])[1]
sig.genes.up <- dim(data.9wk.proc[padj < threshP & log2FoldChange < -threshLFC, ])[1]
sig.genes.down <- dim(data.9wk.proc[padj < threshP & log2FoldChange > threshLFC, ])[1]
Module_Labels <- unique(ASD.modules[, Module_Label])
n.labels <- length(Module_Labels)
phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.all, lower.tail = FALSE)

ORA.modules.9wk.all <- data.table(Module = character(), Pval = numeric())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.9wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i], HGNC] &
                                   padj< 0.05 , ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.9wk.all <- rbind(ORA.modules.9wk.all, 
                           data.frame(Module = Module_Labels[i], 
                           Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.all, lower.tail = FALSE)))
  if (Module_Labels[i] == "M16")
  {print(paste(n.overlap, gs.size, universe.num-gs.size, sig.genes.all))}
}

ORA.modules.9wk.all[,Pval.adj.Bonf:=Pval * n.labels]

#up
ORA.modules.9wk.up <- data.table(Module = character(), Pval = numeric())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.9wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i], HGNC] &
                                   padj< 0.05 & log2FoldChange > threshLFC, ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.9wk.up <- rbind(ORA.modules.9wk.up, 
                           data.frame(Module = Module_Labels[i], 
                                      Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.up, lower.tail = FALSE)))
}

ORA.modules.9wk.up[,Pval.adj.Bonf:=Pval * n.labels]
ORA.modules.9wk.up[, categ := "9wk-Up"]

#down
ORA.modules.9wk.down <- data.table(Module = character(), Pval = numeric())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.9wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i] , HGNC] &
                                   padj< 0.05 & log2FoldChange < -threshLFC, ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.9wk.down <- rbind(ORA.modules.9wk.down, 
                              data.frame(Module = Module_Labels[i], 
                                         Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.down, lower.tail = FALSE)))
}

ORA.modules.9wk.down[,Pval.adj.Bonf:=Pval * n.labels]
ORA.modules.9wk.down[, categ := "9wk-Down"]

#####
#####
#4wk 
universe.num <- dim(data.4wk.proc)[1]
sig.genes.all <- dim(data.4wk.proc[padj < threshP & abs(log2FoldChange) > threshLFC, ])[1]
sig.genes.up <- dim(data.4wk.proc[padj < threshP & log2FoldChange < -threshLFC, ])[1]
sig.genes.down <- dim(data.4wk.proc[padj < threshP & log2FoldChange > threshLFC, ])[1]
Module_Labels <- unique(ASD.modules[, Module_Label])
n.labels <- length(Module_Labels)
phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.all, lower.tail = FALSE)

ORA.modules.4wk.all <- data.table(Module = character(), Pval = double())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.4wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i], HGNC] &
                                   padj< 0.05 , ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.4wk.all <- rbind(ORA.modules.4wk.all, 
                           data.frame(Module = Module_Labels[i], 
                           Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.all, lower.tail = FALSE)))

}

ORA.modules.4wk.all[,Pval.adj.Bonf:=Pval * n.labels]
ORA.modules.4wk.all[, week := "4wk"]

ORA.modules.4wk.up <- data.table(Module = character(), Pval = double())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.4wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i], HGNC] &
                                   padj< 0.05 & log2FoldChange > threshLFC, ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.4wk.up <- rbind(ORA.modules.4wk.up, 
                              data.frame(Module = Module_Labels[i], 
                                         Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.up, lower.tail = FALSE)))

}

ORA.modules.4wk.up[,Pval.adj.Bonf:=Pval * n.labels]
ORA.modules.4wk.up[, categ := "4wk-Up"]

#down
ORA.modules.4wk.down <- data.table(Module = character(), Pval = numeric())

for (i in 1:n.labels)
{
  n.overlap <- dim(data.4wk.proc[NAME %in% ASD.modules[Module_Label == Module_Labels[i] , HGNC] &
                                   padj< 0.05 & log2FoldChange < -threshLFC, ])[1]
  gs.size <- length(ASD.modules[Module_Label == Module_Labels[i], HGNC])
  ORA.modules.4wk.down <- rbind(ORA.modules.4wk.down, 
                                data.frame(Module = Module_Labels[i], 
                                           Pval = phyper(n.overlap, gs.size, universe.num - gs.size, sig.genes.down, lower.tail = FALSE)))
}

ORA.modules.4wk.down[,Pval.adj.Bonf:=Pval * n.labels]
ORA.modules.4wk.down[, categ := "4wk-Down"]

#####
##make graph for ORA modules

ORA.modules.directional.combined <- rbind(ORA.modules.4wk.down, ORA.modules.4wk.up, ORA.modules.9wk.down, ORA.modules.9wk.up)
ORA.modules.directional.combined$Module <- factor(ORA.modules.directional.combined$Module, 
                                                  levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7",
                                                             "M8", "M9", "M10", "M11", "M12", "M13",
                                                             "M14", "M15", "M16", "M17", "M18",
                                                             "Periods3-5_PFC-MSC",
                                                             "Periods4-6_PFC-MSC",
                                                             "Periods8-10_MD-CBC"))
ORA.modules.directional.combined$categ <- factor(ORA.modules.directional.combined$categ, levels = c("9wk-Up", "9wk-Down", "4wk-Up", "4wk-Down"))
ORA.modules.directional.combined[,BH.pval := p.adjust(ORA.modules.directional.combined[,Pval], method = "BH")]



#plot only for the ASD modules
ASD.modules.names <- c("M2","M3", "M13","M16","M17", "Periods3-5_PFC-MSC",
                       "Periods4-6_PFC-MSC",
                       "Periods8-10_MD-CBC")


ORA.mod.plot <- ggplot(ORA.modules.directional.combined[Module %in% ASD.modules.names], aes(y=categ, x=Module, fill = BH.pval)) + 
  geom_tile(color="white", linewidth=2, width=.9, height = .9) +
  theme_minimal() + 
scale_fill_gradientn( trans = "log10", breaks = 0.05,
                      colors = cols.pval) + labs(fill = "P-value")


ggsave("RNASeq_Output/ORA.plot_Direction_ASD_only.pdf", plot = ORA.mod.plot %>% exp.theme() + 
         theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 9), 
               plot.title=element_text(size = 16, hjust=0.5)) + ggtitle("Developmental Modules"),
       device=pdf, width = 8, height = 3)


write.table(ORA.modules.directional.combined, "RNASeq_Output/ORA_modules_directional.txt", row.names = FALSE)


#ASD modules
ASD.modules <- fread("ASD_Modules.txt") #parikshak + willsey
ASD.modules <- ASD.modules[Module_Label != "-",]

#top 5 and bottom genes for modules
modules.9wk <- c("M16", "M17", "Periods3-5_PFC-MSC", "Periods4-6_PFC-MSC")
modules.4wk <- c("M3")

module_extract(seq.data.9v9, modules.9wk, ASD.modules)
modules.9wk.list <- list(M16, M17, `Periods3-5_PFC-MSC`, `Periods4-6_PFC-MSC`)

module_extract(seq.data.4v4, modules.4wk, ASD.modules)
modules.4wk.list <- list(M3)


#Top 5 and bottom 5 genes for each volcano plot
genebox.5up5down(seq.data.9v9, modules.9wk.list,modules.9wk)
genebox.5up5down(seq.data.4v4, modules.4wk.list, modules.4wk)
