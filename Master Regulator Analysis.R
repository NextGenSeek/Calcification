library(ggplot2)
library(ggpubr)
library(viper)
library(pheatmap)
library(RColorBrewer)
library(MUDAN)
library(umap)
library(SingleR)
library(celldex)

### Master regulators analysis
# To infer master regulators (MRs), we applied metaVIPER on the samples of interest, following the pipeline designed for MRs inference in scRNA-seq
# (https://github.com/califano-lab/single-cell-pipeline).

# Briefly, we first built a regulatory network using ARACNe on all single cells. ARACNe was run separately on each set of regulators: 1458 transcription factors, 821 transcriptional cofactors, and 3193 signaling proteins in mouse. 

raw.mat <- read.table('/PATH/TO/matrix.txt',header = T,row.names = 1)
mur.mt <-  read.table('/PATH/TO/murine_mito_genes.csv', header = T, sep = ',', stringsAsFactors = FALSE) # http://asia.ensembl.org/biomart/martview/b2b79163486788b39f300cc958d90092
mur.mt <- mur.mt$mur.symb
mt.mat <- MTFilter(raw.mat,mur.mt)
QCPlots(raw.mat, mur.mt)

filt.mat <- QCTransform(mt.mat) if needed
QCPlots(filt.mat, mur.mt)

# run these if data was not normalized or transformed
cpm.mat <- CPMTransform(filt.mat)
rank.mat <- RankTransform(cpm.mat)
rank.mat = as.matrix(raw.mat)
ARACNeTable(rank.mat, "pbmc-cpm2.tsv")

############ Run script in terminal ############
bash ./run_aracne.sh
################################################

netw_file = "/PATH/TO/merged_net.tsv"
RegProcess(netw_file, rank.mat, out.dir ='/PATH/TO/ARACNA_out_2/', out.name = 'pbmc_r1-net-')
r1.net <- readRDS('/PATH/TO/ARACNA_out/pbmc_r1-net-pruned.rds')
r1.pAct <- viper(rank.mat, r1.net, method = 'mad')

# We then inferred protein activity based on the combined regulatory network and gene expression signature.

# Partitioning around medoids (PAM) clustering was applied to distance matrix built on protein activity with the number of clusters ranging from 2 to 10.

r1.viperDist <- as.dist(viperSimilarity(r1.pAct))
r1.clusts <- PamKRange(r1.viperDist, kmin = 2, kmax = 10)
r1.clustSil <- SilScoreEval(r1.clusts, r1.viperDist)
plot.dat <- data.frame('k' = 2:10, 'Silhouette.Scores' = r1.clustSil)
ggplot(plot.dat, aes(x = k, y = Silhouette.Scores)) + geom_point() + geom_line() +
  ggtitle('1.1 Clustering Silhouette Scores') + theme_bw()

# The optimal number of clusters was determined by silhouette score. In this study, we have identified 3 clusters in mouse data based on protein activity. 

# Next, we generated cluster-specific ARACNe networks using meta-cell profiles built from integrating gene expression profiles of cells with similar protein activity profiles within each cluster. Counts from 5 nearest cells were integrated to generate each meta-cell.

r1.clustMats <- MakeCMfA(rank.mat, r1.viperDist, clustering = r1.clusts$k4, out.dir = '/PATH/TO/', out.name = 'pbmc-r1-clusts')

# Protein activity was inferred again using cluster-specific networks as input to VIPER.


c1.net <- RegProcess(netw_file , r1.clustMats[[1]], 
                     out.dir = '/PATH/TO/', out.name = 'pbmc-r2-c1_')
c2.net <- RegProcess(netw_file , r1.clustMats[[2]], 
                     out.dir = '/PATH/TO/', out.name = 'pbmc-r2-c2_')
c3.net <- RegProcess(netw_file , r1.clustMats[[3]], 
                     out.dir = '/PATH/TO/', out.name = 'pbmc-r2-c3_')
# load in networks
c1.net <- readRDS('/PATH/TO//pbmc-r2-c1_pruned.rds')
c2.net <- readRDS('/PATH/TO//pbmc-r2-c2_pruned.rds')
c3.net <- readRDS('/PATH/TO//pbmc-r2-c3_pruned.rds')

# infer protein activity
r2.pAct <- viper(rank.mat, list('c1' = c1.net, 'c2' = c2.net, 'c3' = c3.net), method = 'none')

#T-test was conducted on transcription factor and cofactor activities between SMC and SEM populations with 100 bootstraps.

r2.cbcMRs <- CBCMRs(r2.pAct) # identify the most representative proteins
r2.pAct.cbc <- r2.pAct[r2.cbcMRs ,] # filter the protein activity matrix
r2.louvain <- LouvainClust(r2.pAct.cbc) # perform clustering analysis

r2.MRs <- BTTestMRs(r2.pAct, r2.louvain, bootstrap.num = 100)

# Top activated and repressed MRs and their target genes are visualized using Cytoscape v3.3.0.

Top.MRS = MR_UnWrap(r2.MRs, top = 10)

## Data visualisation and annotation

r2.cbcUMAP <- CustomUMAP(r2.pAct.cbc)
plot=ClusterScatter(r2.cbcUMAP, r2.louvain, 'Viper Clustering (Louvain)')
plot

ClusterHeatmap(r2.pAct[ MR_UnWrap(r2.MRs, top = 10) , ], clust = r2.louvain, plotTitle = 'Louvain Clustering: Differentially Activated Proteins')

markers <- c('MARKER1', 'MARKER2', 'MARKER3') # markers tissue of interest
MarkerGrid(r2.cbcUMAP, r2.louvain, r2.pAct, markers, 'PBMC Marker Activity')

ref.data <- MouseRNAseqData(cell.ont="nonna")
predictions <- SingleR(test=rank.mat, assay.type.test=1, 
    ref=ref.data, labels=ref.data$label.main)
table(predictions$labels)

plot=ClusterScatter(r2.cbcUMAP, predictions$label, 'Viper Clustering (Louvain)')
