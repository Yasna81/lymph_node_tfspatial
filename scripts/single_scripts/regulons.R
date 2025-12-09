library(Seurat)
library(dplyr)
library(pheatmap)
auc <- read.csv("~/pacakges-set/auc_mtx.csv", row.names = 1 , check.names = FALSE)
dim(auc); summary(as.numeric(as.matrix(auc)))
# we have tfs as cols and cells as rows so we do need a transformation again!
expr <- t(auc)
pbmc[["SCENIC_AUC"]] <- CreateAssayObject(counts = as.matrix(expr))
auc_mat <- GetAssayData(pbmc, assay = "SCENIC_AUC",slot = "data")
top_regulon <- rownames(auc_mat)[order(apply(auc_mat,1,var),decreasing = TRUE)] [1:10]
tops <-c( "MAFF(+)",   "BRCA1(+)" ,  "NFE2L2(+)", "FOSB(+)" , "ZNF263(+)" , "ZNF646(+)"
          , "EBF1(+)"  , "ZNF384(+)", "ATF3(+)" ,  "IRF4(+)")
DefaultAssay(pbmc) <- "SCENIC_AUC"
FeaturePlot(pbmc,features = c("MAFF(+)" ,  "BRCA1(+)" , "NFE2L2(+)", "FOSB(+)", "ZNF263(+)"))
VlnPlot(pbmc,features = tops , pt.size = 0)
# wee keep markers with high biologically plausibility by excluding
# BRCA1 and ZNFs as house keeping or cell cycle regulators that are less informative.
FeaturePlot(pbmc, features = c("MAFF(+)",   "NFE2L2(+)", "FOSB(+)" 
                               , "EBF1(+)", "ATF3(+)" ,  "IRF4(+)"))
VlnPlot(pbmc,features = c("MAFF(+)",   "NFE2L2(+)", "FOSB(+)" 
                          , "EBF1(+)", "ATF3(+)" ,  "IRF4(+)") , pt.size = 0)
##
#pathway enrichment :
DefaultAssay(pbmc) <- "RNA"
library(Seurat)
library(dplyr)
de_becells <- FindMarkers(pbmc, ident.1="B Cell", min.ptc = 0.25,logfc.threshold = 0.25)
de_becells <- de_becells %>% arrange(desc(avg_log2FC))
ranked_genes <- de_becells$avg_log2FC
names(ranked_genes) <- rownames(de_becells)
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
library(msigdbr)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
msig_h_list <- hallmark_sets %>%
    dplyr::select(gs_name, gene_symbol)
gseaRes <- GSEA(geneList = ranked_genes,
                TERM2GENE = msig_h_list,
                pvalueCutoff = 0.05,
                )
dotplot(gseaRes, showCategory = 10) + ggtitle("B cell GSEA")
