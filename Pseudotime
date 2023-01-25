# Pseudotime

library(SeuratWrappers)
library(monocle3)

Obj <- subset(merged_seurat, idents = c("SMC","OC","FC",'SEM'))
DimPlot(Obj, reduction = "umap")

# cds <- as.cell_data_set("dataset")
# cds <- as.cell_data_set(SMC_OC_FC_SEM)
cds <- as.cell_data_set(Obj, assay = "RNA")
cds

# to get cell metadata
colData(cds)
# to get gene metadata
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# Cluster cells using Seurat's UMAP

# assign partitions
recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

# assign the cluster info
list_cluster <- Obj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- Obj@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")
cluster.before.trajectory

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == "OC"]))
plot_cells(cds = cds,
           color_cells_by = "pseudotime",
           label_branch_points = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
colData(cds)
data.pseudo <- as.data.frame(colData(cds))

gplot(data.pseudo, aes(monocle3_pseudotime, ident , fill = ident)) +
  +          geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident , monocle3_pseudotime, median), fill = ident)) +
         geom_boxplot()

# Finding genes that change as a function of pseudotime
deg_OC <- graph_test(cds, neighbor_graph = "principal_graph", cores = 7)
view(deg_OC)
deg_OC_ids <- row.names(subset(deg_OC, q_value < 0.05))

deg_OC %>%
  arrange(q_value) %>%
  dplyr::filter(status == "OK") %>%
  head()

FeaturePlot(Obj, features = c("genes of interest"))

# visualizing pseudotime in seurat
Obj$pseudotime <- pseudotime(cds)
FeaturePlot(Obj, features = "pseudotime", label = TRUE)

# Finding genes that change as a function of pseudotime
#cds <- Exp_cds[,Matrix::colSums(exprs(cds)) != 0]
cds <-estimate_size_factors(cds)
ExpGenes <- c("genes of interest")
Exp_cds <- subset(cds, gene_short_name %in% ExpGenes)

plot_genes_in_pseudotime(Exp_cds,
                         color_cells_by = "ident")
