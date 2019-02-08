#### Cell Marker and Synaptic Gene Analysis ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)

#load packages
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(viridis)
library(forcats)

dir.create("Heatmap")

#Import data
deseq2_wk4 <- fread("Output_DESeq2/data.4v4.txt")
deseq2_wk9 <- fread("Output_DESeq2/data.9v9.txt")
markers <- fread("Neuron_Markers.txt")
counts <- fread("Output_DESeq2/norm_counts_ALL.txt")

colnames(markers)[1] <- "id"
colnames(deseq2_wk4)[1] <- "id"
colnames(deseq2_wk9)[1] <- "id"

#Making heatmap
counts4 <- counts[,c(1,2,4,6,8,10,12,14,16)]
counts9 <- counts[,c(1,3,5,7,9,11,13,15,17)]

colnames(counts4) <- c("id", "X_B1_R1", "X_B1_R2","X_B2_R1","X_B2_R2","C_B1_R1", "C_B1_R2","C_B2_R1","C_B2_R2")
colnames(counts9) <- c("id", "X_B1_R1", "X_B1_R2","X_B2_R1","X_B2_R2","C_B1_R1", "C_B1_R2","C_B2_R1","C_B2_R2")

counts4 <- merge(counts4[id %in% markers[, id]], markers, by = "id" )
counts9 <- merge(counts9[id %in% markers[, id]], markers, by = "id" )

wk.list <- list(counts4, counts9)
wk.name <- c("4", "9")

for (i in 1:length(wk.list)) {
  
  temp <- wk.list[[i]]
  
  temp$Marker_type <- factor(temp$Marker_type, levels = c("Layer", "Neuron", "NPC", "Glia", "Pre-Syn", "Post-Syn"))
  temp$Subtype <- factor(temp$Subtype, levels = c("1", "2-3", "2-4", "5", "5-6", "6", "MT", "", "GABA", "Pan", "Glu"))

  temp <- temp[order(temp[,10], temp[,11]),]
  
  temp_mark <- temp[c(1:23),]
  temp_syn <- temp[c(24:48),]
  
  mark_id <- temp_mark[[1]]
  syn_id <- temp_syn[[1]]
  
  temp_mark$id <- factor(temp_mark$id, levels = mark_id)
  temp_syn$id <- factor(temp_syn$id, levels = syn_id)
  
  temp_mark <- temp_mark[,c(1:9)]
  temp_syn <- temp_syn[,c(1:9)]
  
  temp_mark <- melt(temp_mark)
  temp_syn <- melt(temp_syn)
  
  temp_mark$id <- fct_rev(factor(temp_mark$id, levels = mark_id))
  temp_syn$id <- fct_rev(factor(temp_syn$id, levels = syn_id))
  
  temp_mark$id <- (temp_mark$id)
  temp_syn$id <- (temp_syn$id)
  
  count_mark <- length(grep("X", temp_mark$variable))
  count_syn <- length(grep("X", temp_syn$variable))
  
  temp_mark$facet <- c(rep(0,count_mark), rep(1,count_mark))
  temp_syn$facet <- c(rep(0,count_syn), rep(1,count_syn))
  
  mut_mark <- subset(temp_mark, facet == 0)
  cor_mark <- subset(temp_mark, facet == 1)
  mut_syn <- subset(temp_syn, facet == 0)
  cor_syn <- subset(temp_syn, facet == 1)
  
  mut_mark <- mut_mark[,-4]
  cor_mark <- cor_mark[,-4]
  mut_syn <- mut_syn[,-4]
  cor_syn <- cor_syn[,-4]
  
  assign(paste("mut_mark", wk.name[i], sep = "_"), mut_mark)
  assign(paste("cor_mark", wk.name[i], sep = "_"), cor_mark)
  assign(paste("mut_syn", wk.name[i], sep = "_"), mut_syn)
  assign(paste("cor_syn", wk.name[i], sep = "_"), cor_syn)

}

setwd("Heatmap")
#plotting heatmap blocks
p1 <- ggplot(cor_mark_4, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    legend.position = "none"
  )

p2 <- ggplot(mut_mark_4, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    legend.position = "none"
  )


p3 <- ggplot(mut_mark_9, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=12, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    legend.position = "none"
  )

p4 <- ggplot(cor_mark_9, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    legend.position = "none"
  )

p5 <- ggplot(cor_syn_4, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    legend.position = "none"
  )

p6 <- ggplot(mut_syn_4, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text= element_blank(),
    axis.text.x = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    panel.margin.x=unit(0.5, "cm"),
    panel.margin.y=unit(0.5, "cm"),
    legend.position="none"
  )


p7 <- ggplot(mut_syn_9, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    axis.text.x = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    panel.margin.x=unit(0.5, "cm"),
    panel.margin.y=unit(0.5, "cm"),
    legend.position = "none"
  )

p8 <- ggplot(cor_syn_9, aes(x = variable, y = id, fill = log(value))) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis(option = "magma", na.value = "black") +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    axis.text.x = element_blank(),
    panel.border=element_blank(),
    plot.title=element_text(hjust=0),
    strip.text=element_text(hjust=0),
    panel.margin.x=unit(0.5, "cm"),
    panel.margin.y=unit(0.5, "cm"),
    legend.position = "none"
  )

mark.full <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow = 2)
pdf("markers.pdf", width = 5, height = 6)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow = 2)
dev.off()

#barplots for heatmap

#subsetting marker and synapse genes
wk4_markers <- merge(deseq2_wk4[id %in% markers[, id]], markers, by = "id" )
wk9_markers <- merge(deseq2_wk9[id %in% markers[, id]], markers, by = "id" )

#extract layer/cell type markers
wk4_celltype <- subset(wk4_markers, Marker_type == "Layer" | Marker_type == "Glia" | Marker_type == "NPC" | Marker_type == "Neuron")
wk9_celltype <- subset(wk9_markers, Marker_type == "Layer" | Marker_type == "Glia" | Marker_type == "NPC" | Marker_type == "Neuron")

#formatting data for plotting
wk.list <- list(wk4_celltype, wk9_celltype)
df.names <- c("wk4_celltype", "wk9_celltype")

for (i in 1:length(wk.list)) {
  
  temp <- wk.list[[i]]
  
  temp$Marker_type <- factor(temp$Marker_type, levels = c("Layer", "Neuron", "NPC", "Glia"))
  temp$Subtype <- factor(temp$Subtype, levels = c("1", "2-3", "2-4", "5", "5-6", "6", "MT", ""))
  temp <- temp[order(temp[,8], temp[,9]),]
  
  levs <- unique(temp$id)
  levs <- rev(levs)
  
  temp$id <- factor(temp$id, levels = levs)
  
  #setting code for significance
  temp$plot_pval <- temp$padj
  temp[which(temp$plot_pval > 0.05),10] <- 1
  temp[which(temp$plot_pval <= 0.05),10] <- 2
  temp$plot_pval <- factor(temp$plot_pval, levels = c(1,2))
  
  assign(df.names[i], temp)
  
}

#extracting synapse markers
wk4_synapse <- subset(wk4_markers, Marker_type == "Pre-Syn" | Marker_type == "Post-Syn")
wk9_synapse <- subset(wk9_markers, Marker_type == "Pre-Syn" | Marker_type == "Post-Syn")

#formatting data for plotting
wk.list <- list(wk4_synapse, wk9_synapse)
df.names <- c("wk4_synapse", "wk9_synapse")

for (i in 1:length(wk.list)) {
  
  temp <- wk.list[[i]]
  
  temp$Marker_type <- factor(temp$Marker_type, levels = c("Pre-Syn", "Post-Syn"))
  temp$Subtype <- factor(temp$Subtype, levels = c("GABA", "Pan", "Glu"))
  temp <- temp[order(temp[,8], temp[,9]),]
  
  levs <- unique(temp$id)
  levs <- rev(levs)
  
  temp$id <- factor(temp$id, levels = levs)
  
  #setting code for significance
  temp$plot_pval <- temp$padj
  temp[which(temp$plot_pval > 0.05),10] <- 1
  temp[which(temp$plot_pval <= 0.05),10] <- 2
  temp$plot_pval <- factor(temp$plot_pval, levels = c(1,2))
  
  assign(df.names[i], temp)
  
}

#subset significant genes
wk4.cell.sig <- subset(wk4_celltype, padj <= 0.05)
wk9.cell.sig <- subset(wk9_celltype, padj <= 0.05)
wk4.syn.sig <- subset(wk4_synapse, padj <= 0.05)
wk9.syn.sig <- subset(wk9_synapse, padj <= 0.05)

#organize for plotting
wk4.cell.sig$week <- rep(4)
wk9.cell.sig$week <- rep(9)
wk4.syn.sig$week <- rep(4)
wk9.syn.sig$week <- rep(9)

cell.sig <- rbind(wk4.cell.sig, wk9.cell.sig)
cell.sig$week <- factor(cell.sig$week, levels = c(4,9))

syn.sig <- rbind(wk4.syn.sig, wk9.syn.sig)
syn.sig$week <- factor(syn.sig$week, levels = c(4,9))


p_cell <- ggplot(cell.sig, aes(x=id, y = log2FoldChange, fill = week))  +
  geom_bar(stat='identity', position = position_dodge()) +
  coord_flip() +
  scale_fill_manual(values = c("#bababa", "#ce1256")) +
  ggtitle("Week 4 - Cell Markers") +
  ylim(-4,4) +
  theme(
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 1.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

pdf("sig_markers_layers.pdf", width = 4, height = 6)
p_cell
dev.off()

p_syn <- ggplot(syn.sig, aes(x=id, y = log2FoldChange, fill = week))  +
  geom_bar(stat='identity', position = position_dodge()) +
  coord_flip() +
  scale_fill_manual(values = c("#bababa", "#ce1256")) +
  ggtitle("Week 4 - Cell Markers") +
  ylim(-4,4) +
  theme(
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 1.5),
    plot.title = element_blank(),
    legend.position = "none"
  )


pdf("sig_markers_syn.pdf", width = 4, height = 6)
p_syn
dev.off()