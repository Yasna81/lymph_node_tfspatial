keep <- c("pbmc","auc","auc_mat")
rm(list = setdiff(ls(),keep))
#loading data :
library(Seurat)
counts <- Read10X(data.dir = "~/practice/spatial/data/filtered_feature_bc_matrix")
sample <- CreateSeuratObject(counts = counts)
#adding image
spatial_image <- Read10X_Image(image.dir = "~/practice/spatial/data/spatial")
spatial_image <- subset(spatial_image, cells = colnames(sample))

sample[["slice1"]] <- spatial_image
sample[["slice1"]]@assay <- "RNA"
sample[["percent.mt"]] <- PercentageFeatureSet(sample,pattern = "^MT-")
VlnPlot(sample, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
sample <- subset(sample,subset = nFeature_RNA >3000 & nFeature_RNA < 7500 & percent.mt < 2)
sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample)
sample<- ScaleData(sample)
sample <- RunPCA(sample)
ElbowPlot(sample)
#19 pc
sample <- RunUMAP(sample,dim = 1:19)
sample <- FindNeighbors(sample, dims = 1:19)
sample <- FindClusters(sample,resolution = 0.5)
SpatialDimPlot(sample,label = TRUE)
#single cell,  cell types in tissue
anchors <- FindTransferAnchors(reference = pbmc, query = sample , dims = 1:19,
                               reference.assay = "RNA",query.assay = "RNA")
prediction <- TransferData(anchorset = anchors,
                           refdata = pbmc@active.ident,
                           weight.reduction = sample[["pca"]],
                           dims = 1:19)
sample <- AddMetaData(sample,prediction)
SpatialDimPlot(sample,
               group.by = "predicted.id",
               label = TRUE)

#TF
AUC <- GetAssayData(pbmc[["SCENIC_AUC"]], slot = "data")
#
tf_var <- apply(AUC,1,var)
top_10 <- names(sort(tf_var, decreasing = TRUE))[1:10]
achors_1 <- FindBridgeTransferAnchors(reference = pbmc ,
    query = sample ,
    dims = 1:19 )

tf_matrix <- as.matrix(AUC[top_10,])

tf_matrix <- t(tf_matrix)
#colnames(tf_matrix) <- rownames(pbmc@assays)
tf_preds <- TransferData(
    anchorset = anchors ,
    refdata = tf_matrix,
    dims = 1:19)

tf_preds_mat <- GetAssayData(tf_preds,slot = "data")
tf_preds_t <- t(tf_preds_mat)
sample <- AddMetaData(sample , tf_preds_t)
library(ggplot2)
for (tf in colnames(tf_preds_t)){
    p <- SpatialFeaturePlot(sample , features = tf) +
        ggtitle(paste("spatial activity:",tf))
    print(p)
}
