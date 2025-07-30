### read and merge
rm(list = ls());gc()
setwd("~/projects/immune_scRNA_TaoLab/analysis/")

mouse1 = CreateSeuratObject(counts = Read10X("~/projects/immune_scRNA_TaoLab/data/20240412rawdata/mtx_res_20240315_10x/4HPAA/filtered_feature_bc_matrix/"),
                           project = "4HPAA",min.cells = 3,min.features = 20)
mouse2 = CreateSeuratObject(counts = Read10X("~/projects/immune_scRNA_TaoLab/data/20240412rawdata/mtx_res_20240315_10x/HFD/filtered_feature_bc_matrix/"),
                            project = "HFD",min.cells = 3,min.features = 20)
mouse=merge(mouse1,mouse2)

mouse$treatment = mouse$orig.ident
mouse$condition = ifelse(mouse$treatment=="HFD","MOCK","4HPAA")

VlnPlot(mouse,group.by = "treatment",features = c("nFeature_RNA","nCount_RNA"))

### QC
mouse$percent_mt = PercentageFeatureSet(mouse,"^mt-")
VlnPlot(mouse,"percent_mt")
mouse[["percent_ribo"]]<-PercentageFeatureSet(mouse, "^Rp[sl]")
mouse[["percent_hb"]]<-PercentageFeatureSet(mouse, "^Hb[^(p)]-")
# ribo-gene exp
VlnPlot(mouse,group.by = "cell_type_refine",features = rownames(mouse)[grep("^Rp[sl]",rownames(mouse))],
        split.by = "condition",split.plot = T,stack = T,flip = T)+ggtitle("Ribo-gene expression")

# feature scatter
FeatureScatter(mouse,feature1 = "nCount_RNA",feature2 = "percent_mt",group.by = "treatment")+
  FeatureScatter(mouse,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "treatment")
# decide filter cutoff
hist(mouse@meta.data[["nFeature_RNA"]], breaks=1000)
abline(v=c(300, 4000), col="red")
mouse = subset(mouse,nFeature_RNA>=300 & nFeature_RNA<=4000)
VlnPlot(mouse, group.by = "treatment",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"),
        ncol = 3, pt.size = 0)

### clustering
DefaultAssay(mouse)="RNA"

mouse = NormalizeData(object = mouse)
mouse = FindVariableFeatures(mouse,selection.method = "vst",nfeatures = 2000,assay = "RNA")
# plot variable features with and without labels
# CombinePlots(plots = list(VariableFeaturePlot(mouse), 
#                           LabelPoints(VariableFeaturePlot(mouse),points = head(VariableFeatures(mouse),10),repel = TRUE)))
LabelPoints(VariableFeaturePlot(mouse),points = head(VariableFeatures(mouse),10),repel = F,xnudge = 0, ynudge = 0)
mouse = ScaleData(mouse,features = rownames(mouse))
mouse = RunPCA(mouse,features = VariableFeatures(mouse))
ElbowPlot(mouse,ndims = 40)
mouse <- FindNeighbors(mouse,dims = 1:35)
mouse <- FindClusters(mouse, resolution = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.8, 0.9, 1))
mouse = RunUMAP(mouse,assay = "PCA",dims = 1:35)
##https://www.jianshu.com/p/8ed5541b9da5   about assay using
## resolutions?   0.3
clustree(mouse@meta.data, prefix = "RNA_snn_res.")   # 0.3
DimPlot(mouse,label = T,group.by = "RNA_snn_res.0.3",split.by = "treatment")
DimPlot(mouse,group.by = c("RNA_snn_res.0.3","treatment"),label = T)



save(mouse1,mouse2,mouse,file = "afterClu.RData")













