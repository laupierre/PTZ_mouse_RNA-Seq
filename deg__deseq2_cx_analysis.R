## 230515_A00558_0215_BH2VJLDSX7

library (openxlsx)
library (DESeq2)
library (ggplot2)
library(ggrepel)


## See Github RNA-Seq_mouse/gene_annotation.R
#system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)


## normal STAR results from the RNA-Seq IIT pipeline
a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("Gene|bam", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))


counts <- annot <- a

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

row.names (counts) <- counts$Geneid
counts <- counts[ ,grep ("IIT", colnames (counts))]


colnames (counts) <- gsub ("IIT_RPZ_", "", colnames (counts))
colnames (counts) <- gsub ("_S.*", "", colnames (counts))
samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "condition", "area")
samples$sample <- colnames (counts)
samples$condition <- gsub ("_.*", "", colnames (counts))
samples$area <- gsub (".*_", "", colnames (counts))
row.names (samples) <- colnames (counts)


## CX area

samples.s <- samples[samples$area == "CX", ]
counts.s <- counts[ ,colnames (counts) %in% row.names (samples.s)]
idx <- match (samples.s$sample, colnames (counts.s))
samples.s <- samples.s[idx, ]
samples.s$condition [grepl ("Z", samples.s$condition)] <- "PTZ"
samples.s$condition [grepl ("C", samples.s$condition)] <- "CT"
stopifnot (samples.s$sample == colnames (counts.s))


## DESeq2 
dds <- DESeqDataSetFromMatrix(countData = round (counts.s), colData = samples.s, design = ~ condition)
                                 
# keep <- rowSums(counts(dds)) >= 30
keep <- rowSums(counts(dds) >= 30) >= dim (counts.s)[2]/2
dds <- dds[keep,]
dds


dds <- DESeq(dds)
resultsNames(dds)
# condition_PTZ_vs_CT

res <- results(dds, contrast=list("condition_PTZ_vs_CT"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- resa <- res[order (res$padj), ]

write.xlsx (res, "PTZ_CX_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)



## PCA plot

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed () + geom_label_repel (aes(label = name))

ggsave ("PCA plot PTZ CX experiment.pdf")










