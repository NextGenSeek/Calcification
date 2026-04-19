##### Single-cell RNA-seq analysis - QC

# Load libraries

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(stringr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(data.table)
library(gridExtra)
library(enrichplot)
library(SingleR)
library(scuttle)
library(scater)
library("GO.db")
library("AnnotationDbi")
library("org.Mm.eg.db")

# Create a Seurat object for each sample

files1 <- list.files(path = "/PATH/", pattern=".txt.gz")

for (file in files1)
{
seurat_data <-data.table::fread(file)
seurat_data <- data.frame(seurat_data, row.names = 1)
seurat_obj <- CreateSeuratObject(counts = seurat_data,
min.cells = 10,
min.features = 20,
project = file)
assign(file, seurat_obj)
}

# Filtering (mt/MT)

VlnPlot(matrix.txt.gz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
Sample <- subset(matrix.txt.gz, subset = nFeature_RNA > "NUMBER" & nFeature_RNA < "NUMBER" & percent.mt < 10)
Sample <- NormalizeData(Sample, normalization.method = "LogNormalize", scale.factor = "NUMBER")

# merge datasets

merged_seurat <- merge(x= Sample1, y = c(Sample2,
Sample3),
add.cell.ids = c("SampleID1",
"SampleID2",
"SampleID3"),
project = "name")

merged_seurat$sampleName <- rownames(merged_seurat@meta.data)

# Create metadata dataframe

metadata <- merged_seurat@meta.data

# Create sample column

merged_seurat$sample <- NA
merged_seurat$sample[which(str_detect(merged_seurat$sampleName, "^text"))] <- "SampleDescription"

# Create Weeks column

merged_seurat$Weeks <- NA
merged_seurat$Weeks[which(str_detect(merged_seurat$sampleName, "0"))] <- "0"
merged_seurat$Weeks[which(str_detect(merged_seurat$sampleName, "8"))] <- "8"
merged_seurat$Weeks[which(str_detect(merged_seurat$sampleName, "16"))] <- "16"
merged_seurat$Weeks[which(str_detect(merged_seurat$sampleName, "26"))] <- "26"

# Identification of highly variable features

merged_seurat <- FindVariableFeatures(object= merged_seurat, selection.method = "mvp", dispersion.cutoff = c(1.5,20), mean.cutoff = c(0.01,10))
merged_seurat <- ScaleData(object = merged_seurat)
merged_seurat <- RunPCA(object = merged_seurat)
merged_seurat<- FindNeighbors(object = merged_seurat, dims = 1:50)
merged_seurat<- FindClusters(object= merged_seurat, resolution = 0.4)
merged_seurat_<- RunUMAP(object = merged_seurat, dims = 1:20)

# Create plot
DimPlot(merged_seurat_Ldlr, reduction = "umap", group.by = "Weeks", cols = c("red","green","blue","yellow"))
DimPlot(merged_seurat_Ldlr, reduction = "umap",group.by = "Pos", cols = c("grey", "dark green"))

DimHeatmap(merged_seurat_Ldlr, dims =1:6, cells = 500, balanced = TRUE)

# Show genes that are differentially expressed due to condition

cluster0 <- FindConservedMarkers(merged_seurat, ident.1 = 0, assay = "RNA", grouping.var = "Weeks")
cluster1 <- FindConservedMarkers(merged_seurat, ident.1 = 1, assay = "RNA", grouping.var = "Weeks")

# If cluster consist of multiple cell types, split:

merged_seurat <- FindSubCluster(merged_seurat, cluster = "10",graph.name = "RNA_snn", subcluster.name = "cluster10", resolution = 0.2, algorithm = 1)
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:20)

DimPlot(merged_seurat, reduction = "umap", group.by = "cluster10", label =TRUE)
merged_seurat <- SetIdent(merged_seurat, value = merged_seurat@meta.data$cluster10) # Necessary to change names
merged_seurat <- RenameIdents(merged_seurat, "cluster number" = "cell type")

DimPlot(merged_seurat, reduction = "umap",label = TRUE, label.size = 7, cols = c("OC" = "orange2",
"Macrophage1" = "forestgreen",
"Macrophage2" = "cyan",
"Macrophage3" = "darkturquoise",
"SEM"="darkgoldenrod",
"FC" = "darkmagenta",
"Minor SMC"="darkolivegreen3",
"SMC" = "deepskyblue",
"EC" = "deeppink1",
"T-cell" = "deepskyblue3",
"Fibroblast 1" = "purple2",
"Fibroblast 2"= "gold4"))

########################################

# Gene ontology

ens.str <- row.names(OCvsSMC)
OCvsSMC$symbol <- mapIds(org.Mm.eg.db,
keys=ens.str,
column="ENSEMBL",
keytype="SYMBOL",
multiVals="first")
OCvsSMC$entrez <- mapIds(org.Mm.eg.db,
keys=ens.str,
column="ENTREZID",
keytype="SYMBOL",
multiVals="first")
GOExp <- OCvsSMC
GOExp$names <- row.names(GOExp)

df <- GOExp
original_gene_list <- df$padj
names(original_gene_list) <- df$symbol
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

m <- matrix(0, ncol = 0, nrow = 1254)
data <- data.frame(m)
data$ENTREZID <- df$entrez
data$ENSEMBL <- df$symbol
data$SYMBOL <- rownames(df)
data$padj <- df$padj
data$log2FoldChange <- df$log2FoldChange
data$baseMean <- df$baseMean
data$pvalue <- df$pvalue
data <- na.omit(data)

data2 <- sort(data$log2FoldChange, decreasing = TRUE)
names(data2) <- data$ENSEMBL

gse <- gseGO(geneList = data2,
ont = "ALL",
keyType = "ENSEMBL",
nPerm = 10000,
minGSSize = 3,
maxGSSize = 800,
pvalueCutoff = 0.05,
verbose = TRUE,
OrgDb = org.Mm.eg.db,
pAdjustMethod = "none")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

####################################

# Gene expression plots

Features <- c("gene1","gene2","gene3")
FeaturePlot(merged_seurat, features = Features, min.cutoff = "q10")

DoHeatmap(object = SMC_OC, slot = "data", features = top200)

####################################

# Create dot plot

DotPlotFeatures <- c("gene1","gene2","gene3")
DotPlot(merged_seurat, features = DotPlotFeatures) + RotatedAxis()

###########################################################

# Perform integration to correct for batch effects

# select integration features

seurat_list <- SplitObject(merged_seurat, split.by = "species")
seurat_list <- lapply(X=seurat_list, FUN = SCTransform)

seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)

reference_dataset <- which(names(seurat_list) == "Mouse")

anchors <- FindIntegrationAnchors(object.list = seurat_list,
normalization.method = "SCT",
anchor.features = seurat_features,
reference = reference_dataset)

# integrate data

seurat_integrated <- IntegrateData(anchorset = anchors)

# Scale, run PCA and visualize integrated data

seurat_integrated <- ScaleData(object = seurat_integrated)
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:20)

#############################################

### Overlap two different species

# Make sure gene names are similar for both species (see separate script)

Mouse_sce <- as.SingleCellExperiment(Mouse)
plotExpression(Mouse_sce, features ="gene")
plotPCA(Mouse_sce, colour_by = "ident")

Combined <- as.SingleCellExperiment(Human)

# run SingleR

results_OC <- SingleR(test = Combined, ref = Mouse_sce, labels = Mouse_sce$ident)

Human$SingleR.pruned.calls <- results_OC$pruned.labels
Human$SingleR.calls <- results_OC$labels

Human2 <- SetIdent(Human, value = Human@meta.data$SingleR.calls) # to transfer names

DimPlot(Human, reduction = 'umap', group.by = "SingleR.calls", cols = c("0C"= "red1","SMC" = "tomato1", "FC"="lightsalmon1","SEM"="brown4"))
DimPlot(Human, reduction = 'umap', group.by = "SingleR.calls", label = TRUE, cols = c("0C"= "red"))
