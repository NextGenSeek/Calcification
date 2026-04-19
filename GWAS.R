# Use Biomart or similar to obtain gene set

geneset <- read.csv(file="/PATH/TO/file.csv", header = TRUE)

Cardio <- read.table("/PATH/TO/CARDIoGRAM_GWAS_RESULTS.txt", header = TRUE)
library("tidyr")
Cardio2 <- separate(data=Cardio, col = "chr_pos_.b36.", into = c("chr", "start"), sep = ":")
Cardio2$chr <- gsub("chr","",as.character(Cardio2$chr))
colnames(geneset) [3] <- "start"
colnames(geneset) [4] <- "end"
colnames(geneset) [5] <- "chr"
Cardio2$start2 <- as.numeric(Cardio2$start)
Cardio2$end <- Cardio2$start2 + 50000
Cardio2$start <- Cardio2$start2 - 50000

# Has to be numeric
Cardio2$start <- as.numeric(Cardio2$start)
Cardio2$end <- as.numeric(Cardio2$end)
hsa04022$start <- as.numeric(geneset$start)
hsa04022$end <- as.numeric(geneset$end)

# Start and end have to be the last two columns
Cardio2 <- Cardio2[,-c(4:6,8:14)] 
Cardio3 <- Cardio2[,c(1:2, 4, 3, 5)] 
genesetB <- geneset[,c(1,2,5,6,3,4)] 

library(data.table)
setDT(Cardio3, key = c("chr","start","end"))
setDT(genesetB, key = c("chr","start","end"))

GWAS <- foverlaps(x= genesetB, y= Cardio3, type = "within", nomatch = 0)

# We're only interested in SNPs that reached a significance level of 0.001
GWAS2 <- GWAS[GWAS$pvalue < 0.001,]


### optional:
# This gives multiple Rs for the same gene, select only one per gene
GWAS3 <- GWAS2[!duplicated(GWAS2$Gene.name),] 

# Match significantly differentially expressed genes with this list.
sigGenes <- subset(res, padj < 0.05)

idx1 <- match(GWAS3$Gene.name, sigGenes$hgnc_symbol)
GWAS3$`log2FC` <- sigGenes$log2FoldChange [idx1]
GWAS3$`padj` <- sigGenes$padj [idx1]

GWAS4 <- na.omit(GWAS3)
