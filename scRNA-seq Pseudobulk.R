#Differential gene expression analysis/GO terms/Volcano plots:
# Extract raw counts and metadata to create SingleCellExperiment object
library(Seurat)

data[["RNA"]] <- JoinLayers(data[["RNA"]])
counts <- data[["RNA"]]$counts

counts<-Seurat::GetAssayData(data, assay = "SCT", layer = "counts")
metadata <- data@meta.data
metadata <- metadata[match(colnames(counts), rownames(metadata)), , drop = FALSE]

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(data@active.ident)
metadata$samplefac<-factor(data@meta.data$sample)

# Create single cell experiment object
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts),colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample")]
head(groups)

# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(DESeq2)

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(sce$samplefac))

# Total number of samples
ns <- length(sids)

## Determine the number of cells per sample
table(sce$samplefac)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$samplefac))

## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$samplefac)

#manually making the input metadata
group_id <- c("wt_F","wt_F","wt_M","wt_M","piebald_F","piebald_F","piebald_F","piebald_M","piebald_M","cfpss_M","cfpss_M","cfpss_F","cfpss_F")
sample_id <- c("wt1","wt2","wt3","wt4","piebald1","piebald2","piebald3","piebald4","piebald5","cfpss1","cfpss2","cfpss3","cfpss4")

newdata <- data.frame(sample_id, group_id, n_cells)

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

library(Matrix.utils)
# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)),groupings = groups, fun = "sum")

pb[1:13, 1:13]
head(pb)

rownames(pb)
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb),pattern = "_", n = 2),`[`, 1)

library(magrittr)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb,factor(splitf)) %>%lapply(function(u)
    set_colnames(t(u),stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

#by this point pb has lost the sample_id info and only retained group_id info and I am not sure why
# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$samplefac)

# Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){pb[[x]] %>%colnames()}

library(purrr)
de_samples <- map(1:length(kids), get_sample_ids) %>%unlist()

# Get cluster IDs for each of the samples
samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x],
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

library(dplyr)
# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(cluster_id = de_cluster_ids, sample_id = de_samples)
gg_df <- left_join(gg_df, newdata[, c("sample_id", "group_id")])
metadata <- gg_df %>%dplyr::select(cluster_id, sample_id,group_id)
    


# Generate vector of cluster IDs
clusters <- unique(metadata$cluster_id)

#selecting a specific cluster to do differential analysis on
clusters[1] #"MSCs"
clusters[2] #"ENS"
clusters[3] #"Mesenchyme"
clusters[4] #"Epithelium"
clusters[5] #"Mesenchyme"
clusters[6] #"ENS Progenitors"

# Subset the metadata to only the neuronal cells
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[2]), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id

# Subset the counts to only the neuronal
counts2 <- pb[[clusters[2]]]
cluster_counts <- data.frame(counts2[, which(colnames(counts2) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))      

library(DESeq2)
dds <- DESeqDataSetFromMatrix(cluster_counts,colData = cluster_metadata,design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "group_id")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
library(pheatmap)
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])    

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

#source for this part of the code: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
resultsNames(dds) # lists the coefficients
res <- results(dds, name="group_id_piebald_M_vs_piebald_F")

#relevel
dds$group_id = relevel( dds$group_id, "piebald_F")

write.csv(res, "~/Desktop/results.csv")
          