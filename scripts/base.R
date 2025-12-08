#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("Seurat")
library(Seurat)
#-----------------------------------------------
#if (!requireNamespace("remotes",quietly = TRUE))
   # install.packages("remotes")
library(remotes)
#remotes::install_github("satijalab/seurat-data")
#library(SeuratData)
#InstallData("pbmc3k")
#options(timeout = 600)
pbmc.data <- Read10X(data.dir = "~/practice/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(count = pbmc.data , project ="PBMC3K",min.cells = 3 ,min.features = 200)
#fixing dashes for mit genes :
rownames(pbmc) <- gsub("_","-",rownames(pbmc))
#QC :
#first mitochondria 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^MT-")
VlnPlot(pbmc,features = c ("nFeature_RNA","nCount_RNA","percent.mt"))
# a basic filtering :
pbmc <- subset(pbmc , subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#norm :
pbmc <- NormalizeData(pbmc)
#featureselection
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc)
# 15 pc sounds good for downstream analysis
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc , dim = 1:15 )
DimPlot(pbmc , reduction = "umap")


#----------------------------------------------------------------------------
#downstream_analysis 
#---------------------------------------------------------------------------
library(Seurat)
FeaturePlot(pbmc,features = c("CD3D","CD8A","MS4A1","MKG7","CD14","FCER1A","PPBP"), reduction = "umap",pt.size = 0.5)
#
FeaturePlot(pbmc,features = c("IL7R","CCR7","S100A4","MS4A1","CD79A","NKG7","GNLY"), reduction = "umap",pt.size = 0.5)
FeaturePlot(pbmc,features = c("GZMB","NKG7","MKI67","HLA-DRA","PDCD1","CD4"), reduction = "umap",pt.size = 0.5)
#annots 
corrected_annotations <- c(
    "Memory CD4 T",#0
    "Naive CD4 T",#1
    "CD14+ Monocytes",#2
    "B Cell",#3
    "NK Cell",#4
    "Activated Memory CD4 T",#5
    "CD8 T",#6
    "DC",#7
    "Platalets" #8
)
names(corrected_annotations) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, corrected_annotations)
DimPlot(pbmc, reduction =  "umap", label = TRUE , pt.size = 0.5,repel = TRUE) + ggplot2::ggtitle("PBMC anotation")

#Idents(pbmc) <- pbmc$seurat_clusters
library(fgsea)
library(msigdbr)
#comparing markers between 2 cell types :
de_genes <- FindMarkers(pbmc,
                        ident.1 = "Activated Memory CD4 T",   
                        ident.2 = "Naive CD4 T",
                        min.pct = 0.25,
                        logfc.threshold = 0.25,
                        test.use = "wilcox")
de_genes <- de_genes[!is.na(de_genes$avg_log2FC), ]
gene_ranks <- de_genes$avg_log2FC
names(gene_ranks) <- rownames(de_genes)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

fgsea_results <- fgsea(pathways = pathways,
                       stats = gene_ranks,
                       minSize = 15,      # Minimum genes in pathway
                       maxSize = 500)

fgsea_results <- fgsea_results[order(fgsea_results$NES, decreasing = TRUE), ]
#tops :
top_up <- head(fgsea_results[fgsea_results$NES > 0 & fgsea_results$padj < 0.05, ], 10)
top_down <- head(fgsea_results[fgsea_results$NES < 0 & fgsea_results$padj < 0.05, ], 10)
sig_pathways <- fgsea_results[padj < 0.05]
ggplot(sig_pathways, aes(x = NES, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Pathway Enrichment", x = "NES", y = NULL) +
    theme_minimal() +
    theme(legend.position = "right")

#-------------------------------------------------------------------------

#computing cell cell comunication :
#BiocManager::install("CellChat")
library(remotes)
#remotes::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
data.input <- GetAssayData(pbmc, assay = "RNA",slot = "data")
labels <- Idents(pbmc)
meta <- data.frame(group = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input,meta= meta ,group.by = "group")
#choose species :
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
#subset to relvent signaling :
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#overall network strength
netVisual_circle(cellchat@net$count,vertex.weight = TRUE, weight.scale = TRUE)
#top ligand receptor pairs
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(3,4))
####What signals do myeloid cells provide to the immune system?

myeloid_senders <- c("CD14+ Monocytes", "DC")

cellchat <- computeCommunProb(cellchat)

# Then visualize
netVisual_chord_gene(cellchat, 
                     sources.use = myeloid_senders,
                     targets.use = c("Naive CD4 T", "Memory CD4 T", "CD8 T", "B cells"))

#pathway specific network : which clusters communicate via this pathway 
netVisual_aggregate(cellchat,signaling = "CD99",layout = "circle")
netVisual_aggregate(cellchat,signaling = "MIF",layout = "circle")
netVisual_aggregate(cellchat,signaling = "GALECTIN",layout = "circle")
#which ligands and receptors contribute to a given pathway
netVisual_chord_gene(cellchat,signaling = "GALECTIN")
#the most influential ligand-receptors
netAnalysis_contribution(cellchat,signaling = "GALECTIN")
####----------------------------------------------------------------------------------
#trajectory_analysis :
#remotes::install_github("kstreet13/slingshot")
#remotes::install_github("statOmics/tradeSeq")
remotes::install_github('cole-trapnell-lab/monocle3')
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(matrixStats)
sce <- as.SingleCellExperiment(pbmc)
sce$clusters <- Idents(pbmc)
sce <- slingshot(
    sce,
    clusterLabels = "clusters",
    reducedDim = "UMAP"
)
summary(sce$slingPseudotime_1)
library(RColorBrewer)
colors <- brewer.pal(length(unique(sce$seurat_clusters)), "Set1")
plot(reducedDims(sce)$UMAP, col = sce$clusters,pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
##--------------------------------------------------
immune_genes <- c(
    # Naive markers
    "CCR7", "SELL", "TCF7", "LEF1", "IL7R",
    # Early activation
    "CD69", "CD25", "CD71", "CD38",
    # Housekeeping (controls)
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B"
)

# Filter 
immune_genes <- immune_genes[immune_genes %in% rownames(sce)]
sce_subset <- sce[immune_genes, ]

sce_full <- fitGAM(sce_subset, nknots = 4, parallel = TRUE)

# Find genes that significantly change along trajectory
pattern_results <- patternTest(sce_full)
significant_genes <- pattern_results[pattern_results$pvalue < 0.05, ]
significant_genes <- significant_genes[order(significant_genes$waldStat, decreasing = TRUE), ]

plotSmoothers(sce_full,assays(sce_subset)$counts, gene = "CD69", alpha = 0.6)
plotSmoothers(sce_full,assays(sce_subset)$counts, gene = "CCR7", alpha = 0.6)
plotSmoothers(sce_full,assays(sce_subset)$counts, gene = "CD4", alpha = 0.6)
#Pyscenic prerequisites 
expr <- as.matrix(GetAssayData(pbmc, assay ="RNA",layer = "data"))
rownames(expr) <- sub("\\..*","",rownames(expr))
expr_df <- data.frame(Gene = rownames(expr),expr,check.names = FALSE)
expr_hgnc <- expr_genes[expr_genes %in% tf_list]
write.csv(expr_df,"exp-clean.csv",row.names = FALSE,quote = FALSE)
tf_list <- readLines("~/pacakges-set/allTFs_hg38.txt")
expr_genes <- rownames(expr)
sum(tf_list %in% expr_genes)
head()