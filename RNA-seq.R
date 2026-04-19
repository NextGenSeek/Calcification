 # Load libraries (install first if needed)
 
 #### Analysis to be performed on Salmon quantified files
 
 library("tximeta")
 library(GenomicAlignments)
 library("DESeq2")
 library("vsn")
 library("dplyr")
 library("ggplot2")
 library("pheatmap")
 library("RColorBrewer")
 library( "biomaRt")
 library("genefilter")
 library("AnnotationDbi")
 library("org.Hs.eg.db")
 library(GO.db)
 library(edgeR)
 library(pathview)
 library(gage)
 library(gageData)

 setwd(dir)
 list.files(dir)
 list.files(file.path(dir, "quants"))
 csvfile <- file.path(dir, "SamplesAthero.csv")
 coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
 coldata
 coldata <- coldata[1:8,] # Adjust for number of samples
 coldata$names <- coldata$Sample
 coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
 file.exists(coldata$files)

 se <- tximeta(coldata)
 dim(se)
 head(rownames(se))
 gse <- summarizeToGene(se)
 dim(gse)
 head(rownames(gse))

 assayNames(gse)
 head(assay(gse), 3)
 colSums(assay(gse))
 rowRanges(gse)
 seqinfo(rowRanges(gse))
 colData(gse)

# just check
 gse$Patient
 gse$Condition

# Check the millions of mapped reads round to number of decimal points
 round(colSums(assay(gse))/1e6, 1)

# Construct a DESeqDataSet object from the SummarizedExperiment object with an appropriate design for the analysis

 dds <- DESeqDataSet(gse, design = ~Patient + Condition)

########## Explore analysis and visualisation

### Pre-filtering the dataset: remove all the lines with only zeros
 nrow(dds)
 keep <- rowSums(counts(dds)) > 1
 dds <- dds[keep,]
 nrow(dds)

 ### Look at the variance
 lambda <- 10^seq(from = -1, to = 2, length = 1000)
 cts <- matrix(rpois(1000*100, lambda), ncol = 100)
 meanSdPlot(cts, ranks = FALSE)

 ### Estimate size factors
 dds <- estimateSizeFactors(dds)

 df <- bind_rows(
   as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
          mutate(transformation = "log2(x + 1)"),
   as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
   as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

 colnames(df)[1:2] <- c("x", "y")  

 ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
   coord_fixed() + facet_grid( . ~ transformation)  

 # Check sample distances

 sampleDists <- dist(t(assay(vsd)))
 sampleDists
 sampleDistMatrix <- as.matrix( sampleDists )
 rownames(sampleDistMatrix) <- paste( vsd$age, vsd$Condition, sep = " - " )
 colnames(sampleDistMatrix) <- NULL
 colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
 pheatmap(sampleDistMatrix,
          clustering_distance_rows = sampleDists,
          clustering_distance_cols = sampleDists,
          col = colors)

 ### PCA plot
 plotPCA(vsd, intgroup = c("Patient", "Condition"))

 ### Differential expression analysis
 dds <- DESeq(dds)

 ### Building the results table
 res <- results(dds)
 res
 mcols(res, use.names = TRUE)
 summary(res)

 ### Add gene names
 res <- results(dds)
 res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
 ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
 genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl)
 idx <- match( res$ensembl, genemap$ensembl_gene_id )
 res$entrez <- genemap$entrezgene[ idx ]
 res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
 head(res,4)

 # Example 20 genes with highest variance across samples (working on VST data)
 topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)

### Heatmap

 mat  <- assay(vsd)[ topVarGenes, ]
 mat  <- mat - rowMeans(mat)
 anno <- as.data.frame(colData(vsd)[, c("Condition","Patient")])
 pheatmap(mat, annotation_col = anno, labels_row = res$hgnc_symbol)
 Heatmap fold2change


 sigGenes <- subset(res, padj < 0.05)
 sigGenes <- head(order(sigGenes$log2FoldChange), decreasing = TRUE,25)
 mat <- assay(rld) [sigGenes,]

 # Heatmap
 mat <- mat - rowMeans(mat)
 anno <- as.data.frame(colData(rld)[, c("Condition","Patient")])
 pheatmap(mat, annotation_col = anno, labels_row = res$hgnc_symbol)

 # Annotating and exporting results

 columns(org.Hs.eg.db)

 ens.str <- substr(rownames(res), 1, 15)
 res$symbol <- mapIds(org.Hs.eg.db,
                      keys=ens.str,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
 res$entrez <- mapIds(org.Hs.eg.db,
                      keys=ens.str,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

 ### GO Pathway analysis 

 res.05 <- results(dds, alpha = 0.05)
 table(res.05$padj < 0.05)
 SigGenes <- head(res[order(res$padj), ],increasing = TRUE, 50)
 ListSigGenes <- SigGenes$entrez
 go <- goana(ListSigGenes, species="Hs")
 topGO<-topGO(go,n=50)
 topGO

 ### KEGG pathway analysis
 keg <- kegga(ListSigGenes, species="Hs")
 topKEGG<- topKEGG(keg, n=10, truncate=60)
 topKEGG

 ###  Test pathway

 foldchanges = res$log2FoldChange
 names(foldchanges) = res$entrez
 head(foldchanges)
 kg.hsa=kegg.gsets(species="hsa")
 kegg.sigmet.gs<-kg.hsa$kg.sets[kg.hsa$sigmet.idx]
 kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

 # Run enrichment analysis on all log fc
 fc.kegg.sigmet.p <- gage(foldchanges, gsets = kegg.sigmet.gs)
 fc.kegg.dise.p <- gage(foldchanges, gsets = kegg.dise.gs)

 # convert the kegg results to data frames
 fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
 fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

 # Check most significant
 head(fc.kegg.sigmet.p.up, 20) 
 head(fc.kegg.dise.p.up, 10) 

 # View the most significant pathway from the pathway analysis
 fc.kegg.sigmet.p.up[grepl("hsaNUMBER", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

 # Overlay the expression data onto this pathway
 pathview(gene.data=foldchanges, species="hsa", pathway.id="hsaNUMBER")
 # figure now in working directory 

 # View the pathways as diagrams
 fc.kegg.sigmet.p.up[grepl("hsaNUMBER", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]
 pathview(gene.data=foldchanges, species="hsa", pathway.id="hsaNUMBER", kegg.native=FALSE)
