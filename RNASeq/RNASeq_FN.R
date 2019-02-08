#### Functions for RNA_Seq_Analysis
###  Author: Kirill Zaslavsky

#function to output top 5 and down expressed genes of a given gene sets/lists
#input:
#exp.db = DEG table from DESEQ2 output
#mod.list = list of modules, genesets, etc. each member of a list is a subset of DEG table from DESEQ2
#mod.name = character gector of module names, corresponding in order to mod.list
#output:
#plots of 5 up and 5 down regulated genes

genebox.5up5down <- function(exp.db, mod.list, mod.names)
{
  
  base.counts <- exp.db[,c(1:2)]
  
  for (i in 1:length(mod.list)) {
    
    temp.df <- mod.list[[i]]
    temp.genes <- temp.df[[1]]
    merge.df <- merge(temp.df, base.counts[which(base.counts$NAME %in% temp.genes),])
    
    sub.df <- subset(merge.df, Sig == TRUE)
    sub.df <- sub.df[order(-log2FoldChange),]
    
    count.pos <- length(which(sub.df$log2FoldChange > 0))
    count.neg <- length(which(sub.df$log2FoldChange < 0))
    
    fc.max <- max(sub.df$log2FoldChange)
    fc.min <- min(sub.df$log2FoldChange)
    
    #setting maximum x value
    if(fc.max <= 3){
      fc.max <- 3
    } else {
      fc.max <- ceiling(fc.max)
      fc.min <- -fc.max
      
    }
    
    #setting minimum x value
    if(fc.min >= -3){
      fc.min <- -3
    } else {
      fc.min <- floor(fc.min)
      fc.max <- -fc.min
      
    }
    
    brk <- c(floor(-fc.max/2), 0, ceiling(fc.max/2))
    
    
    #plotting
    if(count.pos >= 5 & count.neg >= 5){
      
      plot.df <- sub.df[c(1:5, (nrow(sub.df)-4):nrow(sub.df)),]
      plot.df$NAME <- factor(plot.df$NAME, levels = plot.df$NAME)
      plot.df$NAME <- fct_rev(plot.df$NAME)
      
      temp.plot <- ggplot(plot.df, aes(x=NAME, y = log2FoldChange, fill = log(baseMean)))  +
        geom_bar(stat='identity', position = position_dodge(), color = "black", size = 0.4) +
        scale_y_continuous(breaks = brk)+
        coord_flip(ylim = c(fc.min, fc.max)) +
        scale_fill_viridis(option = "magma", limits = c(0,11)) +
        
        theme(
          axis.text.x = element_text(face = "bold", size = 10, color = "black"),
          axis.text.y = element_text(face = "bold", size = 10, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        )
      
      assign(paste("bar", mod.names[i], sep = "_"), temp.plot)
      ggsave(paste(mod.names[i], "5UP5DOWN.pdf", sep = ""), device = "pdf", width = 5, height = 5, units = "cm")
      
    } else if (count.pos + count.neg < 10) {
      
      plot.df <- sub.df
      plot.df$NAME <- factor(plot.df$NAME, levels = plot.df$NAME)
      plot.df$NAME <- fct_rev(plot.df$NAME)
      
      temp.plot <- ggplot(plot.df, aes(x=NAME, y = log2FoldChange, fill = log(baseMean)))  +
        geom_bar(stat='identity', position = position_dodge(), color = "black", size = 0.4) +
        scale_y_continuous(breaks = brk)+
        coord_flip(ylim = c(fc.min, fc.max)) +
        scale_fill_viridis(option = "magma", limits = c(0,11)) +
        
        theme(
          axis.text.x = element_text(face = "bold", size = 10, color = "black"),
          axis.text.y = element_text(face = "bold", size = 10, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) + scale_y_continuous(limits <- c(fc.min, fc.max), breaks = brk)
      
      assign(paste("bar", mod.names[i], sep = "_"), temp.plot)
      ggsave(paste(mod.names[i], "5UP5DOWN.pdf", sep = ""), device = "pdf", width = 5, height = 5, units = "cm")
      
    } else {
      
      t <- count.neg - 1
      
      plot.df <- sub.df[c(1:(10 - count.neg), (nrow(sub.df)-t):nrow(sub.df)),]
      plot.df <- plot.df[!duplicated(plot.df),]
      plot.df$NAME <- factor(plot.df$NAME, levels = plot.df$NAME)
      plot.df$NAME <- fct_rev(plot.df$NAME)
      
      temp.plot <- ggplot(plot.df, aes(x=NAME, y = log2FoldChange, fill = log(baseMean)))  +
        geom_bar(stat='identity', position = position_dodge(), color = "black", size = 0.4) +
        scale_y_continuous(breaks = brk)+
        coord_flip(ylim = c(fc.min, fc.max)) +
        scale_fill_viridis(option = "magma", limits = c(0,11)) +
        theme(
          axis.text.x = element_text(face = "bold", size = 10, color = "black"),
          axis.text.y = element_text(face = "bold", size = 10, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        )+ scale_y_continuous(limits <- c(fc.min, fc.max), breaks = brk)
      
      assign(paste("bar", mod.names[i], sep = "_"), temp.plot)
      ggsave(paste(mod.names[i], "5UP5DOWN.pdf", sep = ""), device = "pdf", width = 5, height = 5, units = "cm")
    }
    
  }
}


#function to output expression data.table for each module
#input: DEG dataset, module list
#output: DEG table for module
module_extract <- function(exp_db, modules.list, ASD.gene.list)
{
  Module.list <- list()
  
  for (i in 1:length(modules.list)) {
    
    temp.id <- modules.list[i]
    temp.name <- modules.list[i]
    
    gene.data <- ASD.gene.list[Module_Label %in% temp.id,HGNC]
    
    temp.exp <- exp_db[NAME %in% gene.data ,]
    
    print(temp.name)
    assign(temp.name, temp.exp, inherits = TRUE)
    
  }
}

### make heatmap from gene list
gg.go.heat <- function(geneQuery ="neuron maturation", geneList=NA, dat = data.counts.merged){
  if(is.na(geneList))
  {geneList <- geneQuery.grep(geneQuery)}
  
  data.counts.picked <- data.counts.merged[grep(geneList,NAME)] #pick gene/genes
  data.counts.picked2 <- gather(data.counts.picked, key = "SampleID", value = "EXPR", c(2:9, 11:18))
  data.counts.picked2 <- data.counts.picked2[, c(-2,-3)]
  data.counts.picked3 <- as.data.table(separate(data.counts.picked2, col = 2, into = c("Line","Genotype", "Batch", "Rep", "Timepoint"),  sep = "_"))
  data.counts.norm <- RNASeq.norm.ctrl.4wk(data.counts.picked3)
  
  ggplot(data.counts.norm, aes(x=Timepoint, y=NAME)) + geom_tile(aes(fill=EXPR.norm))+
    scale_fill_gradient2(low="blue", mid = "white", high = "red") +
    facet_wrap(~Genotype)
}


### make linegraph from gene list
gg.go.line <- function(geneQuery ="neuron maturation", geneList=NA, dat = data.counts.merged, title = "Relative change in RNA abundance"){
  if(is.na(geneList))
  {geneList <- geneQuery.grep(geneQuery)}
  
  data.counts.picked <- data.counts.merged[grep(geneList,NAME)] #pick gene/genes
  data.counts.picked2 <- gather(data.counts.picked, key = "SampleID", value = "EXPR", c(2:9, 11:18))
  data.counts.picked2 <- data.counts.picked2[, c(-2,-3)]
  data.counts.picked3 <- as.data.table(separate(data.counts.picked2, col = 2, into = c("Line","Genotype", "Batch", "Rep", "Timepoint"),  sep = "_"))
  data.counts.norm <- RNASeq.norm.ctrl.4wk(data.counts.picked3)
  
  ggplot(data.counts.norm, aes(x=Timepoint, y = EXPR.norm, color = Genotype, group = NAME)) +
    #geom_point(alpha = 0.5) +
    #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "line", lwd = 1, aes(group = Genotype)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1, width = 0.2, aes(group = Genotype)) +
    ggtitle(title) +
    labs(y="Mean Counts", x="Timepoint")+
    theme(    
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.ticks.x = element_line(color="black", size = 1),
      strip.text = element_text(size = 16, face  = "bold"),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14, face = "bold"),
      legend.title = element_text(colour="black", size = 16, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(color="black", size = 18, face = "bold", hjust = 0.5))
}


#create row-normalized dataset
RNASeq.norm.ctrl <- function(dat)
{
  ###for ctrl means of EXP
  #get ctrl means
  mccd <- dat[, .SD]
  
  ctrl.means <- mccd[Genotype=="cor" & Timepoint == "4wk",
                     lapply(.SD, mean),
                     by=.(NAME),
                     .SDcols = c("EXPR")]
  
  #print(ctrl.means)
  for (i in 1:dim(ctrl.means)[1])
  {
    iGeneName <- ctrl.means[i, NAME]
    #   print(iGeneName)
    mccd[NAME==iGeneName, `:=`(EXPR.norm = as.double(as.double(EXPR) / ctrl.means[NAME==iGeneName, EXPR]))]
  }
  
  return(mccd)
}

#plotting for geneBoxes
geneBox_plfn <- function(gp)
{gp+
    theme(axis.title = element_text(color="black", size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(fill=NA, size = 2, color = "black"),
          axis.line = element_blank(),
          plot.title = element_text(color="black", size = 20, hjust=0.5))+
    scale_y_continuous(breaks = -10:10, limits = c(-5,5))
}

#nice color scales
# hack together a colourbar
cols.25 <- c(colorRampPalette(c("#ce472e", "#f05336","#e29421","#e29421","#eec73a" 
))(75), 
colorRampPalette(c("#eec73a", "#c9e2f6", "#e7f0fa"
))(25)
)

# hack together a colourbar
cols.25.rvs <- c(colorRampPalette(c(  "#e7f0fa", "#c9e2f6", "#eec73a"
))(25),
colorRampPalette(c( "#eec73a","#e29421", "#e29421","#f05336","#ce472e"  ),
                 bias=2)(75) 

)

cols.50 <- c(colorRampPalette(c("#ce472e", "#f05336","#e29421","#e29421","#eec73a" 
))(50), 
colorRampPalette(c("#eec73a", "#c9e2f6", "#e7f0fa"
))(50)
)

cols.50.rvs <- c(colorRampPalette(c(  "#e7f0fa", "#c9e2f6", "#eec73a"
))(50),
colorRampPalette(c( "#eec73a","#e29421", "#e29421","#f05336","#ce472e"  ),
                 bias=2)(50) 
)


# hack together a colourbar
cols.10 <- c(colorRampPalette(c("#ce472e", "#f05336","#e29421","#e29421","#eec73a" 
))(90), 
colorRampPalette(c("#eec73a", "#c9e2f6", "#e7f0fa"
))(10)
)

# discrete to highlight p-values
cols.pval <- c(colorRampPalette(c("#ce472e", "#f05336","#e29421","#e29421","#eec73a" 
))(95), 
colorRampPalette(c("#e7f0fa"
))(5)
)

# hack together a colourbar
cols.10.rvs <- c(colorRampPalette(c(  "#e7f0fa", "#c9e2f6", "#eec73a"
))(10),
colorRampPalette(c( "#eec73a","#e29421", "#e29421","#f05336","#ce472e"  ),
                 bias=2)(90) 

)

cols.75 <- c(colorRampPalette(c("#ce472e", "#f05336","#e29421","#e29421","#eec73a" 
))(25), 
colorRampPalette(c("#eec73a", "#c9e2f6", "#e7f0fa"
))(75)
)

# hack together a colourbar
cols.75.rvs <- c(colorRampPalette(c(  "#e7f0fa", "#c9e2f6", "#eec73a"
))(75),
colorRampPalette(c( "#eec73a","#e29421", "#e29421","#f05336","#ce472e"  ),
                 bias=2)(25))

cols.1lvl <- colorRampPalette(c( "#eec73a","#e29421", "#e29421","#f05336","#ce472e"  ), bias =2) (100) 


cols.try <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(50),
                              colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(50))

cols.try2 <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#eec73a"))(50),
               colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(50))

lm_eqn <- function(df){
  m <- lm(df[,1] ~ df[,2], df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

