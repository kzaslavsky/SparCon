#Volcano plot for RNASeq for comparisons at 9weeks
#authors: Fraser McCready, Kirill Zaslavsky

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
library(purrr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(forcats)
library(viridis)
library(grid)
library(gridExtra)
library(gtable)
library(biomaRt)
library(tools)
library(GSA)
source("RNASeq_FN.R")

dir.create("Volcano_Plots")

# function for plotting volcano plots
exp_volcano <- function (gp,scale = 1, alfa = 0.5, size = 3 )
{
  gp + 
    geom_point(size = size, alpha = alfa, scale = scale)+
    theme_classic()+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 12*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank(),
      #      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "none",
      #axis.line.y = element_line(),
      #axis.line.x = element_line(),
      plot.title = element_text(lineheight=.8, face="bold", size = 12*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )+
    scale_x_continuous(limits = c(-4,4))+
    geom_vline(xintercept=log2(1.25))+
    geom_vline(xintercept=-log2(1.25))+
    geom_hline(yintercept=-log10(0.05))+
    scale_color_manual(values= c( "gray40", "red"))+
    labs(x="Log2 Fold Change", y = "-log10(P)")
}

#set constants
threshP = 0.05
threshLFC = log2(1.25)

#Read RNAseq data and add significance column for coloring
seq.data.9v9 <- fread("Output_DESeq2/data.9v9.txt")
seq.data.9v9[, Sig := ifelse(padj < threshP & abs(log2FoldChange) > threshLFC, TRUE, FALSE) ]

#
seq.data.4v4 <- fread("Output_DESeq2/data.4v4.txt")
seq.data.4v4[, Sig := ifelse(padj < threshP & abs(log2FoldChange) > threshLFC, TRUE, FALSE) ]


#Loading GO accessions and biomart
geneset_list <- fread("geneset_list.csv")
GO.id <- as.character(geneset_list[[2]])
GO.name <- make.names((geneset_list[[1]]))

#load genesets
genesets <- GSA.read.gmt('Human_GOBP_AllPathways_no_GO_iea_May_01_2018_symbol.gmt')
gsetnames <- as.data.table(genesets$geneset.names)
gsetnames.2 <- separate(gsetnames, "V1", into = c('Name', 'Database','db_ID'), sep = "%")


#getting Genesets from .gmt file
for (i in 1:length(GO.id)) {
  
  temp.id <- GO.id[i]
  temp.name <- GO.name[i]
  
  gene.data <- genesets$genesets[[which(gsetnames.2$db_ID == GO.id[[i]])]]
  
  temp.exp <- seq.data.9v9[NAME %in% gene.data ,]
  
  assign(temp.name, temp.exp)
  
}


#Manually adding FMRP data set
fmrp <- fread("FMRP.csv")
fmrp <- seq.data.9v9[NAME %in% fmrp$NAME,]

GO.list <- list(synapse.assembly, chemical.synaptic.transmission, glutamate.receptor.signaling, neuron.projection.development,
                regulation.of.synaptic.plasticity, translation, fmrp)

GO.name <- c(GO.name, "FMRP.targets")

plot_cols <- c(rep("lightseagreen", 3), "mediumpurple3", rep("deepskyblue", 3))

setwd("Volcano_Plots")
#plotting volcano plots
for (j in 1:length(GO.list)){
  
  temp.df <- GO.list[[j]]
  temp.col <- plot_cols[j]
  
  temp.plot <- ggplot(data = temp.df, 
                      aes(x=log2FoldChange, y = -log10(padj), color = Sig)) %>% exp_volcano(size = 2, alfa = 0.6)+
    scale_y_continuous(limits = c(-0.5, 20))
  
  temp.plot <- temp.plot + scale_color_manual(values= c( "gray40", temp.col)) #+ ggtitle(toTitleCase(gsub("\\.", " ", GO.name[j])))
  
  assign(paste("volc", GO.name[j], sep = "_"), temp.plot)
  ggsave(paste(GO.name[j], ".pdf", sep = ""), plot = temp.plot, device = "pdf", width = 5, height = 5, units = "cm")
}

#top 5 and bottom 5 genes for each volcano plot
genebox.5up5down(seq.data.9v9, GO.list, GO.name)






