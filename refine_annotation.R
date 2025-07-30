rm(list = ls());gc()
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
mouse = readRDS("afterAnnp_SeuratObj.RDS")
mouse$condition = ifelse(mouse$treatment=="HFD","MOCK","4HPAA")


Idents(mouse) = mouse$RNA_snn_res.1
mouse$cell_type = as.character(mouse$cell_type)

### Neutrophil?
FeaturePlot(mouse,label = T,order = T,features = c("Itgam","Csf3r","Hp","F13a1","Ly6g","Mpo","Mmp9"))
VlnPlot(mouse,stack = T,flip = T,features = c("Itgam","Csf3r","Hp","F13a1","Ly6g","Mpo","Mmp9"))  # 25?
### eosinophil
FeaturePlot(mouse,features = c("Siglecf","Il5ra","Ccr3","Epx"),label = T,order = T)
VlnPlot(mouse,stack = T,flip = T,features = c("Siglecf","Il5ra","Ccr3","Epx"))  # nothing
### Mono or Macro
VlnPlot(mouse,stack = T,flip = T,features = c("Cd14","S100a8","S100a9","Fcgr3"))  # 10, 17, 25?

VlnPlot(mouse,stack = T,flip = T,features = c("Cd14","S100a8","S100a9","Fcgr3",
                                              "Adgre1","Cd86","Cd80","Cd68","Cd163","Mrc1","Marco","C1qa","C1qb"))
DotPlot(mouse,features = c("Cd14","S100a8","S100a9","Fcgr3",
                           "Adgre1","Cd86","Cd80","Cd68","Cd163","Mrc1","Marco","C1qa","C1qb"))
DotPlot(subset(mouse,RNA_snn_res.1 %in% c(6,10,17)),features = c("Adgre1","Ly6c2","Ccr2"))
VlnPlot(subset(mouse,RNA_snn_res.1 %in% c(6,10,17)),stack = T,flip = T,features = c("Adgre1","Ly6c2","Ccr2"))  ## 17=mono
mouse$cell_type[mouse$RNA_snn_res.1 == 17]="Monocyte"

### subset ILC and Macrophage 
##### ILC
ILC = subset(mouse,cell_type == "ILC")
DefaultAssay(ILC)="RNA"

ILC = NormalizeData(object = ILC)
ILC = FindVariableFeatures(ILC,selection.method = "vst",nfeatures = 2000,assay = "RNA")
LabelPoints(VariableFeaturePlot(ILC),points = head(VariableFeatures(ILC),10),repel = F,xnudge = 0, ynudge = 0)
ILC = ScaleData(ILC,features = rownames(ILC))
ILC = RunPCA(ILC,features = VariableFeatures(ILC))
ElbowPlot(ILC,ndims = 40)
ILC <- FindNeighbors(ILC,dims = 1:35)
ILC <- FindClusters(ILC, resolution = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.8, 0.9, 1))
ILC = RunUMAP(ILC,assay = "PCA",dims = 1:35)
##https://www.jianshu.com/p/8ed5541b9da5   about assay using
## resolutions?   0.3/0.15
clustree(ILC@meta.data, prefix = "RNA_snn_res.")   # 0.3
Idents(ILC) = ILC$RNA_snn_res.0.5
FeaturePlot(ILC,features = c("nCount_RNA","percent_mt","percent_hb","nFeature_RNA"),label = T,order = T)

### Markers from https://www.nature.com/articles/s41467-022-35347-6
VlnPlot(ILC,features = c("Il7r","Il2ra",
                         "Ncr1","Tbx21","Ifng","Gzma","Gzmb","Id3",
                         "Il1rl1","Klrg1","Gata3","Il4","Il5","Il13",
                         "Rorc","Il22","Ccr6","Il10","Il17a"),
        group.by = "RNA_snn_res.0.5",stack = T,flip = T,fill.by = "ident")+
  DimPlot(ILC,group.by = "RNA_snn_res.0.5",label = T)       ####### 4=ILC1  1=ILC3 0,2,3=ILC2

VlnPlot(ILC,features = c("Il7r","Il2ra",
                         "Tbx21",
                         "Il1rl1","Klrg1","Gata3","Il4","Il5",
                         "Il17a"),
        group.by = "RNA_snn_res.0.5",stack = T,flip = T,fill.by = "ident")+
  DimPlot(ILC,group.by = "RNA_snn_res.0.5",label = T)
FeaturePlot(ILC,features = c("S100a4","Ly6a","Ccr2","Cxcr6","Il17a"),order = T)

## verify clu5 doublets with B cell?
VlnPlot(ILC,features = c("Cd19","Cd79a","Ms4a1","Cd79b"),stack = T)

ILC_noclu5 = subset(ILC,RNA_snn_res.0.5 != 5 )
DefaultAssay(ILC_noclu5)="RNA"

ILC_noclu5 = NormalizeData(object = ILC_noclu5)
ILC_noclu5 = FindVariableFeatures(ILC_noclu5,selection.method = "vst",nfeatures = 2000,assay = "RNA")
LabelPoints(VariableFeaturePlot(ILC_noclu5),points = head(VariableFeatures(ILC_noclu5),10),repel = F,xnudge = 0, ynudge = 0)
ILC_noclu5 = ScaleData(ILC_noclu5,features = rownames(ILC_noclu5))
ILC_noclu5 = RunPCA(ILC_noclu5,features = VariableFeatures(ILC_noclu5))
ElbowPlot(ILC_noclu5,ndims = 40)
ILC_noclu5 <- FindNeighbors(ILC_noclu5,dims = 1:35)
ILC_noclu5 <- FindClusters(ILC_noclu5, resolution = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.8, 0.9, 1))
ILC_noclu5 = RunUMAP(ILC_noclu5,assay = "PCA",dims = 1:35)

####### 4=ILC1  1=ILC3 0,2,3=ILC2
mouse$cell_type_refine = mouse$cell_type
mouse$cell_type_refine[colnames(subset(ILC,RNA_snn_res.0.5 == 4))] = "ILC-1"
mouse$cell_type_refine[colnames(subset(ILC,RNA_snn_res.0.5 == 1))] = "ILC-3"
mouse$cell_type_refine[colnames(subset(ILC,RNA_snn_res.0.5 %in% c(0,2,3)))] = "ILC-2"
mouse$cell_type_refine[colnames(subset(ILC,RNA_snn_res.0.5 == 5))] = "Doublet"


## markers from Lingjun
VlnPlot(ILC,features = c("Tbx21","Itga1","Gata3","E2f1","Rorc","Ccr6"),stack = T,flip = T)
DotPlot(ILC,features = c("Tbx21","Itga1","Gata3","E2f1","Rorc","Ccr6"))

ILC$cell_type_refine = mouse$cell_type_refine[colnames(ILC)]
DotPlot(ILC,features = c("Tbx21","Itga1","Gata3","E2f1","Rorc","Ccr6","Cd79a"),group.by = "RNA_snn_res.0.3")+ggtitle("ILC subtypes marker")+
  DimPlot(ILC,group.by = "RNA_snn_res.0.3",label = T)+DimPlot(ILC,group.by = "cell_type_refine",label = T)



##### Macro
Macro = subset(mouse,cell_type == "Macrophage")
DefaultAssay(Macro)="RNA"

Macro = NormalizeData(object = Macro)
Macro = FindVariableFeatures(Macro,selection.method = "vst",nfeatures = 2000,assay = "RNA")
LabelPoints(VariableFeaturePlot(Macro),points = head(VariableFeatures(Macro),10),repel = F,xnudge = 0, ynudge = 0)
Macro = ScaleData(Macro,features = rownames(Macro))
Macro = RunPCA(Macro,features = VariableFeatures(Macro))
ElbowPlot(Macro,ndims = 40)
Macro <- FindNeighbors(Macro,dims = 1:35)
Macro <- FindClusters(Macro, resolution = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.8, 0.9, 1))
Macro = RunUMAP(Macro,assay = "PCA",dims = 1:35)
##https://www.jianshu.com/p/8ed5541b9da5   about assay using
## resolutions?   0.3/0.15
clustree(Macro@meta.data, prefix = "RNA_snn_res.")   # 0.1
Idents(Macro) = Macro$RNA_snn_res.0.4
FeaturePlot(Macro,features = c("nCount_RNA","percent_mt","percent_hb","nFeature_RNA"),label = T,order = T)

### cluster 4 and 5 maybe doublets
FeaturePlot(Macro,features = c("Cd3g","Cd79a"),order = T)
VlnPlot(Macro,features = c("C1qa","Cd3g","Cd79a","nCount_RNA","nFeature_RNA","percent_mt"),stack = T,flip = T)

mouse$cell_type_refine[colnames(subset(Macro,RNA_snn_res.0.4 %in% 4:5))] = "Doublet"

DotPlot(Macro,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40",
                "Cd163","Mrc1","Cd79a","Cd3g"))+ggtitle("Macrophage subtypes marker")+
  DimPlot(Macro,label = T,group.by = "RNA_snn_res.0.4")+DimPlot(Macro,label = T,group.by = "cell_type_refine")


Macro = subset(mouse,cell_type_refine == "Macrophage")
DefaultAssay(Macro)="RNA"

Macro = NormalizeData(object = Macro)
Macro = FindVariableFeatures(Macro,selection.method = "vst",nfeatures = 2000,assay = "RNA")
LabelPoints(VariableFeaturePlot(Macro),points = head(VariableFeatures(Macro),10),repel = F,xnudge = 0, ynudge = 0)
Macro = ScaleData(Macro,features = rownames(Macro))
Macro = RunPCA(Macro,features = VariableFeatures(Macro))
ElbowPlot(Macro,ndims = 40)
Macro <- FindNeighbors(Macro,dims = 1:35)
Macro <- FindClusters(Macro, resolution = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5, 0.8, 0.9, 1))
Macro = RunUMAP(Macro,assay = "PCA",dims = 1:35)
##https://www.jianshu.com/p/8ed5541b9da5   about assay using
## resolutions?   0.3/0.15
clustree(Macro@meta.data, prefix = "RNA_snn_res.")   # 0.1
DimPlot(Macro,group.by = "RNA_snn_res.0.5",label = T)

Idents(Macro) = Macro$RNA_snn_res.0.4
# rerun clustering and umap
VlnPlot(Macro,features = c("Ninj1","Gas6","Mrc1","Cd163","Igf1"),stack = T,flip = T)
VlnPlot(Macro,features = c("Cd14","Cst3","C5ar1","Stab1","Trem2"),stack = T,flip = T)
VlnPlot(Macro,features = c("Ccr2","Cd52","H2-Eb1","H2-Ab1","S100a6"),stack = T,flip = T)

FindAllMarkers(Macro,only.pos = T) %>% group_by(cluster) %>% slice_head(n=10) %>% as.data.frame()
VlnPlot(Macro,features = c("Il1b","Cxcl10","Ccl4","Cd36",
                           "Il10","Cd163","Cd206a"),stack = T,flip = T)

FeaturePlot(mouse,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40","H2-Ab1"))
DotPlot(Macro,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40","H2-Ab1"))+DimPlot(Macro,label = T)

FeaturePlot(mouse,c("Cd163","Mrc1"),order = T)
DotPlot(Macro,c("Cd163","Mrc1"))+DimPlot(Macro,label = T)

### markers https://pmc.ncbi.nlm.nih.gov/articles/PMC11657537/#Sec24
VlnPlot(mouse,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40","H2-Ab1",
                "Cd163","Mrc1"),stack = T)


mouse$cell_type_refine[mouse$RNA_snn_res.1 == 10 & mouse$cell_type_refine != "Doublet"] = "Macro-1"
mouse$cell_type_refine[mouse$RNA_snn_res.1 == 6 & mouse$cell_type_refine != "Doublet"] = "Macro-2"

mouse$cell_type_refine[mouse$RNA_snn_res.1 == 14 & mouse$cell_type_refine == "Macrophage"] = "Dendritic"
mouse$cell_type_refine[mouse$RNA_snn_res.1 == 3 & mouse$cell_type_refine == "Macrophage"] = "Fibroblast"

## additionally
Idents(Macro) = Macro$RNA_snn_res.0.5
DotPlot(Macro,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40","Cd163","Mrc1","Cd79a","Cd3g"))
mouse$cell_type_refine[colnames(subset(Macro,RNA_snn_res.0.5 == 3))] = "Macro-1"
Macro$cell_type_refine[colnames(subset(Macro,RNA_snn_res.0.5 == 3))] = "Macro-1"

DotPlot(Macro,c("Cd86","Il1b","Cxcl10","Ccl5","Tnf","Cd40","Cd163","Mrc1","Cd79a","Cd3g"))+
DimPlot(Macro,label = T,group.by = c("RNA_snn_res.0.5","cell_type_refine"))


### Neutro and eosnophil
FeaturePlot(mouse,features = c("Siglecf","Il5ra","Ccr3","Epx"),order = T)
FeaturePlot(mouse,features = c("Itgam","Csf3r","Hp","F13a1","Ly6g","Mpo","Mmp9"),order = T)
mouse$cell_type_refine[mouse$cell_type == "Granulocyte"] = "Neutrophil"


DimPlot(mouse,group.by = "cell_type_refine",label = T,repel = T)
mouse$cell_type_refine = factor(mouse$cell_type_refine,levels = c("GC B","naive B","Plasma B","Dendritic","Neutrophil","Macro-1","Macro-2","Mast","Monocyte",
                                                                  "CD4+ effective" ,"CD4+ memory","CD4+ Treg","CD8+ cytotoxic","Proliferate T","DNT","Th17",
                                                                  "ILC-1","ILC-2","ILC-3","NK","NKT",
                                                                  "Endothelial","Fibroblast","Doublet","LowQ"))
ggarrange(DimPlot(mouse,group.by = "cell_type_refine",label = T,repel = T,cols = paletteer::paletteer_d("ggsci::default_igv"))+
            guides(color = guide_legend(ncol = 1,override.aes = list(size = 3)))+
            theme(legend.direction = "horizontal",legend.position = "right"),
          VlnPlot(mouse,group.by = "cell_type_refine",stack = T,flip = T,
                  features = c("Cd79a","Mki67","Ighm","Sdc1","Itgax","Hdc","C1qa","Tnf","Cd163","Mcpt8",
                               "Cd4","Icos","Ccr7","Foxp3","Cd8a","Gzma","Ikzf2", "Il17a",
                               "Tbx21","Gata3","Rorc","Nkg7","Il23r",
                               "Fabp4","Col1a2","nCount_RNA","percent_mt"),fill.by = "ident")+
            scale_fill_paletteer_d("ggsci::default_igv")+NoLegend()+ggtitle("Marker gene expression"),
          ncol = 2)


## Cell composition bar
as.data.frame(prop.table(table(mouse$cell_type_refine,mouse$condition),margin = 2)) %>% ggplot(aes(x=Var2,y=Freq*100,fill = Var1))+
  geom_col()+paletteer::scale_fill_paletteer_d("ggsci::default_igv")+theme_bw()+ggtitle("Cell proportion")+
  xlab("condition")+ylab("Cell proportion")+labs(fill = "cell_type")+guides(fill = guide_legend(ncol = 1,override.aes = list(size = 3)))+
  theme(legend.title = element_text(size=13),plot.title = element_text(face = "bold",hjust = 0.5,size = 16))
ggarrange(DimPlot(mouse,group.by = "cell_type_refine",label = T,repel = T,cols = paletteer::paletteer_d("ggsci::default_igv"))+
            NoLegend(),
          VlnPlot(mouse,group.by = "cell_type_refine",stack = T,flip = T,
                  features = c("Cd79a","Mki67","Ighm","Sdc1","Itgax","Hdc","C1qa","Tnf","Cd163","Mcpt8",
                               "Cd4","Icos","Ccr7","Foxp3","Cd8a","Gzma","Ikzf2", "Il17a",
                               "Tbx21","Gata3","Rorc","Nkg7","Il23r",
                               "Fabp4","Col1a2","nCount_RNA","percent_mt"),fill.by = "ident")+
            scale_fill_paletteer_d("ggsci::default_igv")+NoLegend()+ggtitle("Marker gene expression"),
          ncol = 2)
DimPlot(mouse,group.by = "cell_type_refine",label = T,repel = T,cols = paletteer::paletteer_d("ggsci::default_igv"),split.by = "condition")

save.image("detailly_anno.RData")
saveRDS(mouse,file = "refinedAnno_SeuratObj.RDS")

# mouse = readRDS("refinedAnno_SeuratObj.RDS")
mouse$cell_type_refine = plyr::mapvalues(mouse$cell_type_refine, 
                                         from = c("Plasma B", "CD4+ Treg","CD4+ effective","Proliferate T","CD8+ cytotoxic"), 
                                         to = c("Memory B", "Treg","CD4+ T","CD4+ T","CD8+ T"))
saveRDS(mouse,file = "refinedAnno_SeuratObj.RDS")




