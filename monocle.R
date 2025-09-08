#pseudotime analysis with monocle 3

library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(SeuratWrappers)
library(garnett)

data<-readRDS("~/Desktop/Ret Ednrb Manuscript 2022/scRNAseq/WT_CFPss_ss_only/E14.5_females_ENSonly_3genotypes.rds")
DefaultAssay(data) <- "SCT"

#can subset cells to look at, note seurat_clusters may be != idents
data2<- data[, data$seurat_clusters %in% c("8", "7", "6","0")]

cds <- as.cell_data_set(data)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, color_cells_by = "seurat_clusters")

#to plot genes along the pseudotime
genenames <- c("Ret","Ednrb","Fabp7","Vamp2","Pdgfra","Col6a1","Tubb3","Epcam")
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
plot_cells(cds, genes=genenames,label_cell_groups=FALSE, show_trajectory_graph=FALSE,color_cells_by = "seurat_clusters")

###pseudotime, first specify nodes then plot###
cds <- order_cells(cds)

#function to automatically define roots
get_earliest_principal_node <- function(cds, cell_phenotype="broad_cell_type", root_types){
  root_pr_nodes <- lapply(root_types, function(root_type) {
    cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    return(root_pr_nodes)
  })
  return(unlist(root_pr_nodes))
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, cell_phenotype="seurat_clusters", root_types=c("1","4")))
#plot
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

##3d pseudotime##
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)

cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds_3d, cell_phenotype="seurat_clusters", root_types=c("1","4")))

plot_cells_3d(cds_3d, color_cells_by="seurat_clusters")
#may need to increase colors, use color_palette = colors
colors<-colors <- brewer.pal(11, "Set3")

