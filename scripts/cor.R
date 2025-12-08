# checking our hypothesis on tf activity 
tf_names <- c("IRF4(+)","ATF3(+)","FOSB(+)","MAFF(+)","NFE2L2(+)")

tf_data <- sample@meta.data[,tf_names]

head(tf_data)
cor_sp <- cor(tf_data,method = "spearman")
cor_sp
corrplot::corrplot(cor_sp,method = "color" , addCoef.col = "black" ,tl.cex = 1.2)
SpatialFeaturePlot(sample,features = c("NFE2L2(+)","IRF4(+)"))

SpatialFeaturePlot(sample,features = c("ATF3(+)","FOSB(+)"))
cor_mat <- saveRDS(cor_sp, "cor_sp.rds")
