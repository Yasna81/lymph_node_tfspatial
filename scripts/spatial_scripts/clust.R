sample[["TF"]] <-CreateAssayObject(counts = t(tf_data)) 
DefaultAssay(sample) <- "TF"
sample <- ScaleData(sample,assay = "TF" ,verbose = FALSE)
sample <- RunPCA(sample,assay = "TF" ,features = rownames(sample[["TF"]]), verbose = FALSE)
ElbowPlot(sample)
sample <- FindNeighbors(sample,assay = "TF" ,dims = 1:4)
sample <- FindClusters(sample,assay = "TF",resolution = 0.3)
SpatialDimPlot(sample, label = TRUE)
avg_tf <- AverageExpression(sample,assays = "TF",features = rownames(sample[["TF"]]))
# z score across clusters and make a heatmap 
avg_tf_mat <- avg_tf$TF
zscore_tf <- t(scale(t(avg_tf_mat)))
library(pheatmap)
pheatmap(zscore_tf,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         color = colorRampPalette(c("blue","white","red"))(50),
         main = "Z-scored TF activity across clusters")


#cell type enrichment
library(dplyr)
sample$TF_cluster <- sample$seurat_clusters
dominant_celltype <- sample@meta.data %>%
    group_by(TF_cluster) %>%
    count(predicted.id) %>%
    slice_max(n,n=1)

SpatialDimPlot(sample, group.by = "predicted.id",label = TRUE)
SpatialDimPlot(sample, group.by = "TF_cluster", label =  TRUE)

library(Seurat)
library(dplyr)
library(fgsea)
library(msigdbr)

# parameters
cluster_col <- "seurat_clusters" # change if you used a different cluster column
assay_to_use <- "RNA" # use "RNA" counts/normalized for expression-based FC
min_spots <- 5 # minimum cluster size to consider
n_top_terms <- 20 # for later reporting

# 1) average expression per cluster 
avg_expr <- AverageExpression(sample, assays = assay_to_use, slot = "data", group.by = cluster_col)
avg_mat <- avg_expr[[assay_to_use]] # genes x clusters

clusters <- colnames(avg_mat)

# 2) compute fold-change vector for each cluster: cluster vs other clusters 
avg_dense <- as.matrix(avg_mat)
fold_changes <- list()

for (cl in clusters) {
    
    cl_idx <- which(colnames(avg_dense) == cl)
    other_idx <- setdiff(seq_len(ncol(avg_dense)), cl_idx)
    
    # convert both reference and "others" to numeric
    cl_values <- as.numeric(avg_dense[, cl_idx])
    others_mat <- as.matrix(avg_dense[, other_idx, drop = FALSE])
    others_mean <- rowMeans(others_mat)
    
    fc_vec <- cl_values - others_mean
    
    # name the vector
    names(fc_vec) <- rownames(avg_mat)
    
    fold_changes[[cl]] <- fc_vec
}

# 3)  removing tiny clusters from interpretation
cluster_sizes <- table(sample[[cluster_col]])
keep_clusters <- colnames(avg_mat)  # ignore size filter completely

# 4) prepare gene sets 
msig <- msigdbr(species = "Homo sapiens", category = "H") # change species if needed
pathways <- split(msig$gene_symbol, msig$gs_name)

# 5) running fgsea for each  cluster using ranked vector (descending)
fgsea_results <- list()

for (cl in keep_clusters) {
    message("Running fgsea for cluster: ", cl)
    
    #  extract the pre-computed fold-change vector
    fc_vec <- fold_changes[[cl]]
    
    if (is.null(fc_vec) || length(fc_vec) == 0) {
        warning("No fold-change data for cluster ", cl)
        next
    }
    
   
    if (!is.numeric(fc_vec)) fc_vec <- as.numeric(fc_vec)
    
    # Crucial: keep the original gene names!
    gene_names <- names(fc_vec)
    if (is.null(gene_names) || length(gene_names) != length(fc_vec)) {
        warning("Gene names missing or mismatched for cluster ", cl, " – trying to recover from avg_mat")
        gene_names <- rownames(avg_mat)
    }
    
    # Remove any genes with NA fc
    keep <- !is.na(fc_vec)
    fc_vec <- fc_vec[keep]
    gene_names <- gene_names[keep]
    
    # Rank: highest fold-change = top
    ranked <- sort(fc_vec, decreasing = TRUE)
    names(ranked) <- gene_names
    
    # Run fgsea
    fg <- fgsea(pathways = pathways,
                stats    = ranked,
                minSize  = 15,
                maxSize  = 500,
                nperm    = 10000)  # 5000 is ok, 10000 is better
    
    fgsea_results[[cl]] <- fg %>% arrange(padj)
}


fgsea_results[["g6"]][1:n_top_terms, c("pathway", "padj", "NES", "leadingEdge")]




#
library(ggplot2)
library(dplyr)
library(tidyr)  

# Extract top 6 pathways per cluster 
top_n <- 6   

top_pathways <- lapply(names(fgsea_results), function(cl) {
    df <- fgsea_results[[cl]]
                        if (nrow(df) == 0) return(NULL)
                        df %>%
                            arrange(padj) %>%
                            head(top_n) %>%
                            mutate(cluster = cl,
                                   pathway = gsub("HALLMARK_", "", pathway),     
                                   pathway = gsub("_", " ", pathway),
                                   pathway = stringr::str_to_title(pathway)) %>%
                            select(cluster, pathway, padj, NES)
}) %>% bind_rows()


if (nrow(top_pathways) == 0) stop("No significant pathways found – try increasing nperm or lowering padj cutoff")

# Ordering clusters  (g0, g1, g2, ...) and pathways by first appearance
top_pathways$cluster <- factor(top_pathways$cluster,
                               levels = paste0("g", 0:20))  

top_pathways$pathway <- factor(top_pathways$pathway,
                               levels = rev(unique(top_pathways$pathway[order(top_pathways$cluster, top_pathways$padj)])))

#Create the dotplot
p <- ggplot(top_pathways, aes(x = cluster, y = pathway)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_color_gradient2(low = "blue", mid = "grey90", high = "red",
                          midpoint = 0, limits = c(-max(abs(top_pathways$NES)), max(abs(top_pathways$NES))),
                          name = "NES") +
    scale_size_continuous(name = expression("-log"[10]*"(adj. p)"),
                          range = c(2, 8),
                          breaks = c(1, 2, 3, 5),
                          labels = c("1", "2", "3", "≥5")) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_blank(),
          panel.grid.major.y = element_line(color = "grey90"),
          legend.position = "right") +
    labs(title = "Top Hallmark Pathways per TF Cluster",
         caption = paste("Top", top_n, "pathways (padj-ranked) • fgsea on avg logFC vs rest"))

print(p)


# to validate our assumption we want to check the markers : / between 2 b cell clusters/
library(dplyr)
library(Seurat)
library(ggplot2)
library(tibble)
clust_a <- "4"
clust_b <- "5"
Idents(sample) <- "TF_cluster"
markers_A_vs_B <- FindMarkers(sample,
                              ident.1 = clust_a,
                              ident.2 = clust_b,      
                              assay = "RNA",
                              slot = "data",
                              test.use = "wilcox",
                              min.pct = 0.1,
                              logfc.threshold = 0) %>%   
    rownames_to_column("gene") %>%
    arrange(p_val_adj) %>%
    mutate(cluster = paste0(clust_a, "_vs_", clust_b))

top_markers <- markers_A_vs_B %>%
    group_by(cluster) %>%
    arrange(p_val_adj, .by_group = TRUE) %>%
    slice_head(n = 10)