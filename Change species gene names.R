# Mouse and Human have different gene names !

library(biomaRt)
library(tidyr)

files1 <- list.files(path = "/PATH/TO", pattern=".txt.gz")

matrix <- read.table("matrix.txt", header = TRUE)
matrix_list <- matrix$gene

human = useMart(host = "https://dec2021.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mouse = useMart(host = "https://dec2021.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

genesV2_matrix = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol",
                   values = matrix_list,
                   mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows = TRUE)

 # some don't exist, so match the genes that exist
idx1 <- match(matrix$gene, genesV2_matrix$HGNC.symbol)
matrix$Mgene <- genesV2_matrix$MGI.symbol [idx1]

# remove all the NAs 
matrix_b <- matrix[!is.na(matrix$Mgene), ]

# Replace gene column with Mgene column
matrix_b["gene"] <- matrix_b["Mgene"]
matrix_b <- matrix_b[order(matrix_b[,'gene']),]
matrix_b <- matrix_b[!duplicated(matrix_b$gene), ]
