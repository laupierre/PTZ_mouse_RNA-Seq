## 230515_A00558_0215_BH2VJLDSX7

library (openxlsx)
library (DESeq2)
library (ggplot2)


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


## HPC area

samples.s <- samples[samples$area == "HPC", ]
counts.s <- counts[ ,colnames (counts) %in% row.names (samples.s)]
idx <- match (samples.s$sample, colnames (counts.s))
samples.s <- samples.s[idx, ]
stopifnot (samples.s$sample == colnames (counts.s))













