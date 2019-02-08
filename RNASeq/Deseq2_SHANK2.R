## DESEQ2 analysis and preparation of files for GSEA
## author: KZ
## Output: DEG tables

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
#if not using R-studio, set working directory manually
#using setwd()

#biocLite("tximport")
#initialize the main file
library(DESeq2)
source("SparCon_fn.R")
source("RNASeq_FN.R")

# counts tables - import
data.4wk.counts <- fread("simple_coding_counts_4wk.txt")
setkey(data.4wk.counts, Geneid)
data.9wk.counts <- fread("simple_coding_counts_9wk.txt")
setkey(data.9wk.counts, Geneid)

# counts tables - merge
counts.merged <- as.data.table(merge(data.4wk.counts, data.9wk.counts, by = "Geneid"))
counts.merged <- as.matrix(counts.merged[,-1])
rownames(counts.merged) <- data.9wk.counts$Geneid

# make coldata file
coldata <- data.frame(condition = rep(rep(c("MUT", "CORR"),each = 4),2), 
                      Timepoint = rep(c("4wk", "9wk"), each = 8))
row.names(coldata) <- colnames(counts.merged)
coldata$condition <- as.factor(coldata$condition)
coldata$Timepoint <- as.factor(coldata$Timepoint)

# check whether coldata row.names = counts.merged colnames
all(row.names(coldata) == colnames(counts.merged))

##===================##
# make grouping variable for contrasts
coldata$group <- factor(paste(coldata$condition, coldata$Timepoint, sep = "_"))

#=======================================
### make dds
dds <- DESeqDataSetFromMatrix(countData = counts.merged,
                              colData = coldata,
                              #design = ~ condition + Timepoint + condition:Timepoint)       #if want to use 2x2 design
                              design = ~group)     


#pre-filtering as per manual
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#perform deseq
dds <- DESeq(dds)

#specific contrasts with lfcShrinkage - pairwise

##uncomment to output contrast of CTRL or MUT at 4 and 9 wk i.e. if you want to look at developmental changes
# res.4_9CORR <- results(dds, contrast = c("group", "CORR_9wk", "CORR_4wk"), alpha = 0.05)
# res.4_9CORR.lfc <- lfcShrink(dds, contrast = c("group", "CORR_9wk", "CORR_4wk"), res = res.4_9CORR, method = "apeglm") #if you want to compare CTRL at 4 and 9 wk
# res.4_9MUT <- results(dds, contrast = c("group", "MUT_9wk", "MUT_4wk"), alpha = 0.05)
# res.4_9MUT.lfc <- lfcShrink(dds, contrast = c("group", "MUT_9wk", "MUT_4wk"), res = res.4_9MUT, method = "apeglm") #if you want to compare MUT at 4 and 9 wk

res.4_4 <- results(dds, contrast = c("group", "MUT_4wk", "CORR_4wk"), alpha = 0.05)
res.4_4.lfc <- lfcShrink(dds, contrast = c("group", "MUT_4wk", "CORR_4wk"), res = res.4_4, method = "apeglm")
res.9_9 <- results(dds, contrast = c("group", "MUT_9wk", "CORR_9wk"), alpha = 0.05)
res.9_9.lfc <- lfcShrink(dds, contrast = c("group", "MUT_9wk", "CORR_9wk"), res = res.9_9, method = "apeglm")


## if using 2x2 factorial design
# res.4_4 <- results(dds,  name = "condition_MUT_vs_CORR", alpha = 0.05   )
# res.9_9 <- results(dds, list(c("condition_MUT_vs_CORR", "conditionMUT.Timepoint9wk")), alpha = 0.05)

#####
#Output
dir.create("Output_DESeq2")

#output norm counts table
norm.counts <- counts(dds, normalized = TRUE) #this is from my data
norm.counts.dt <- as.data.table(norm.counts)
norm.counts.dt$NAME <- row.names(norm.counts)
setcolorder(norm.counts.dt, c("NAME",colnames(norm.counts)))
write.table(norm.counts.dt, file = "Output_DESeq2/norm_counts_ALL.txt", sep ="\t", col.names = TRUE, row.names = FALSE)

#output contrasts
# write.table(base::append(res.4_9CORR.lfc, list(NAME = row.names(res.4_9CORR.lfc)))[,c("NAME", colnames(res.4_9CORR.lfc))],
#             file = "Output_DESeq2/data.4_9devCORR.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
# write.table(base::append(res.4_9MUT.lfc, list(NAME = row.names(res.4_9MUT.lfc)))[,c("NAME", colnames(res.4_9MUT.lfc))],
#             file = "Output_DESeq2/data.4_9devMUT.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(base::append(res.4_4.lfc, list(NAME = row.names(res.4_4.lfc)))[,c("NAME", colnames(res.4_4.lfc))],
            file = "Output_DESeq2/data.4v4.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(base::append(res.9_9.lfc, list(NAME = row.names(res.9_9.lfc)))[,c("NAME", colnames(res.9_9.lfc))],
            file = "Output_DESeq2/data.9v9.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

#prep GSEA
dir.create("GSEAinput")
ranked.4v4 <- as.data.table(res.4_4.lfc)[, `:=`(rank = -log10(padj)*sign(log2FoldChange),
                                                NAME = row.names(res.4_4.lfc))]
ranked.9v9 <- as.data.table(res.9_9.lfc)[, `:=`(rank = -log10(padj)*sign(log2FoldChange),
                                                NAME = row.names(res.4_4.lfc))]



ranked.4v4.rnk <- ranked.4v4[order(ranked.4v4[, rank], decreasing = TRUE), .(NAME, rank)]
ranked.4v4.rnk <- ranked.4v4.rnk[!is.na(rank),]
write.table(ranked.4v4.rnk, "GSEAinput/4wk-MUToverCTRL.rnk", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

ranked.9v9.rnk <- ranked.9v9[order(ranked.9v9[, rank], decreasing = TRUE), .(NAME, rank)]
ranked.9v9.rnk <- ranked.9v9.rnk[!is.na(rank),]
write.table(ranked.9v9.rnk, "GSEAinput/9wk-MUToverCTRL.rnk", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


