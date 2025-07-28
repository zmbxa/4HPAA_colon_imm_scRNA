### annotation
rm(list = ls());gc()
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
load("afterClu.RData")

mouse_joined = JoinLayers(mouse)
FeaturePlot(mouse_joined,features = c("Ptprc","Sdc1"),label = T,order = T)

Idents(mouse_joined) = mouse_joined$RNA_snn_res.0.5

FindMarkers(mouse_joined,ident.1 = c(2,24,18)) %>% head() # Col1a2,Col1a1,Lama2,Dpt,Scn;  fibroblast?
### verify fibroblast markers
FeaturePlot(mouse_joined,features = c("Fbln1","Fbln2","Col5a1","Loxl1","Lum"),order = T)
VlnPlot(mouse_joined,features = c("Fbln1","Fbln2","Col5a1","Loxl1","Lum"),stack = T,flip = T)

mouse_merge = mouse
mouse = mouse_joined

### annotation
mouse$cell_class = as.character(mouse$RNA_snn_res.0.5)

mouse$cell_class[mouse$cell_class %in% c(2,24,18)] = "Fibroblast"


FeaturePlot(mouse,features = c("Cd3g","Cd3e","Cd19","Cd79a","Cst3","Lyz2"),order = T,label = T)

mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(3,5,7:9,14,16,23)] = "T"
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(0,4,14:16,21)] = "B"
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(1,11:12,19)] = "Mye"
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(10)] = "NK"
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(6)] = "ILC"

VlnPlot(mouse,c("Il22","Cd3d","Klrb1b","Clnk"),stack = T,flip = T)
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(13)] = "NKT"
FeaturePlot(mouse,features = c("Cma1","Cpa3","Car8"))
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(17)] = "Mye"

VlnPlot(mouse,features = c("Bank1","Ms4a1","Cd79a",
                           "Cd3d","Cd3e","Cd3g"),stack = T,flip = T)
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(16)] = "Doublet"

mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(22)] = "Endothelial"
mouse$cell_class[mouse$RNA_snn_res.0.5 %in% c(20)] = "LowQ"

DimPlot(mouse,group.by = "cell_class",label = T)+paletteer::scale_color_paletteer_d("ggthemes::Miller_Stone")
DimPlot(mouse,group.by = "cell_class",label = T,split.by = "treatment")+paletteer::scale_color_paletteer_d("ggthemes::Miller_Stone")

# save(mouse,file = "rough_anno.RData")

### sub-type anno
mouse$cell_type = mouse$cell_class

FeaturePlot(mouse,order = T,features = c("Cst3","Lyz2",
                                         "Cd14","S100a8","S100a9","Fcgr3",
                                         "Adgre1","Cd86","Cd80","Cd68","Cd163","Mrc1","Macro","C1qa","C1qb"))

### try singleR ImmGen
library(SingleR) 
sce_for_SingleR <- GetAssayData(mouse, slot="data")
clusters=mouse@meta.data$RNA_snn_res.0.5
mouseImmu <- ImmGenData()
pred.mouseImmu <- SingleR(test = sce_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.main,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")

mouseRNA <- MouseRNAseqData()
pred.mouseRNA <- SingleR(test = sce_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType=data.frame(ClusterID=levels(mouse@meta.data$RNA_snn_res.0.5),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
mouse$singleR = cellType$mouseRNA[match(mouse$RNA_snn_res.0.5,cellType$ClusterID)]
mouse$singleR_Imm = cellType$mouseImmu[match(mouse$RNA_snn_res.0.5,cellType$ClusterID)]

DimPlot(mouse,group.by = c("singleR","singleR_Imm"),label=T)

mouse$singleR[mouse$cell_class != "Immune"] = mouse$cell_class[mouse$cell_class != "Immune"]
DimPlot(mouse,group.by = c("singleR"),label=T,split.by = "treatment")

##############
DimPlot(mouse,group.by = c("RNA_snn_res.0.5","cell_type"),label = T)

### detaily annotation
all_markers = FindAllMarkers(mouse,only.pos = T) 
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(17)] = "Mast"
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(4,14,21)] = "Plasma B"

##### B cells
DotPlot(mouse,features = c("Cd79a","Cd19","Igha","Ighm","Mki67","Aicda"))
VlnPlot(mouse,features = c("Cd79a","Cd19","Igha","Ighm","Mki67","Aicda"),stack = T,flip = T)
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(0)] = "naive B"
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(15)] = "GC B"

### Myeloid lineage
VlnPlot(mouse,features = c("Zbtb46","Cd209a","Xcr1","Itgax"),stack = T,flip = T)
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(11,12)] = "Dendritic"

VlnPlot(mouse,features = c("C1qa","C1qb","Lyz2","Mrc1"),stack = T,flip = T)
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(1)] = "Macrophage"
## M1 M2?
VlnPlot(subset(mouse,cell_type == "Macrophage"),flip = T,stack= T,group.by = "RNA_snn_res.1",
        features = c("Cd36","Cd86","Cxcl10","Stat6","Cxcr2","Ccl17"))
DotPlot(subset(mouse,cell_type == "Macrophage"),features = c("Il1a","Tlr2","Ifit1","Isg15","Oasl5","Ddx60","Mafb","Mrc1","Ear2","Tlr4"))

DotPlot(mouse,features = c("Hdc"))
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(19)] = "Granulocyte"


### T 
VlnPlot(mouse,features = c("Cd4","Cd8a","Gzmb","Foxp3","Mki67"),stack = T,flip=T)
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(9)] = "CD4+ Treg" # Foxp3
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(23)] = "Proliferate T" # Mki67
FeaturePlot(subset(mouse,cell_class == "T"),label = T,features = c("Il17a","Ifng","Foxp3"),order = T)

VlnPlot(subset(mouse,cell_class == "T"),features = c("Il17a","Foxp3","Ifng","Mki67"),stack = T,flip = T)
mouse$cell_type[mouse$RNA_snn_res.0.5 %in% c(7)] = "Th17" # Il17a

Idents(mouse) = mouse$RNA_snn_res.0.8
FindMarkers(mouse,ident.1 = 9) %>% head(20)
VlnPlot(mouse,c("Cd4","Cd8a","Cd28","Icos"),flip = T,stack = T)
FeaturePlot(mouse,c("Cd4","Cd8a","Cd28","Icos"),label = T,order = T)
mouse$cell_type[mouse$RNA_snn_res.0.8 %in% c(9)] = "CD4+ effective"
mouse$cell_type[mouse$RNA_snn_res.0.8 == 5 & mouse$cell_class!= "ILC"] = "CD4+ effective"

FeaturePlot(subset(mouse,cell_type ==  "T"),features = c("Cd4","Cd8a","Ccr7","Gzmk","Lef1"),label = T,order = T)
mouse$cell_type[mouse$RNA_snn_res.0.8 %in% c(7)] = "CD4+ memory"

DotPlot(mouse,features = c("Cd8a","Gzmk"))
mouse$cell_type[mouse$RNA_snn_res.0.8 %in% c(15)] = "CD8+ cytotoxic"

DotPlot(mouse,features = c("Cd8a","Ikzf2","Cd7","Gzmk","Ccl5"))
FeaturePlot(mouse,c("Cd8a","Cd4","Ccl5"),label = T,order = T)
VlnPlot(mouse,c("Cd8a","Cd4","Nkg7","Cd3d","Mki67","Foxp3"),group.by = "cell_type",stack = T,flip = T)
mouse$cell_type[mouse$cell_type == "T"] = "DNT"

DimPlot(mouse,group.by = "cell_type",label = T,repel = T)+paletteer::scale_color_paletteer_d("trekcolors::lcars_series")
DimPlot(mouse,group.by = "cell_type",label = T,repel = T,split.by = "treatment")+
  paletteer::scale_color_paletteer_d("trekcolors::lcars_series")+
  guides(color = guide_legend(ncol = 1,override.aes = list(size = 3)))

DimPlot(mouse,group.by = "cell_type",label = T,repel = F,split.by = "treatment")+
  paletteer::scale_color_paletteer_d("trekcolors::lcars_series")+guides(color = guide_legend(nrow = 3,override.aes = list(size = 3)))+
  theme(legend.direction = "horizontal",legend.position = "bottom",legend.byrow = T)
DimPlot(mouse,group.by = "cell_type",label = T,repel = F,split.by = "treatment")+
  paletteer::scale_color_paletteer_d("trekcolors::lcars_series")+guides(color = guide_legend(ncol = 1,override.aes = list(size = 3)))+
  theme(legend.direction = "horizontal",legend.position = "right")

## Cell composition bar
as.data.frame(prop.table(table(mouse$cell_type,mouse$treatment),margin = 2)) %>% ggplot(aes(x=Var2,y=Freq*100,fill = Var1))+
  geom_col()+paletteer::scale_fill_paletteer_d("trekcolors::lcars_series")+theme_bw()+ggtitle("Cell proportion")+
  xlab("treatment")+ylab("Cell proportion")+labs(fill = "cell_type")+guides(fill = guide_legend(ncol = 1,override.aes = list(size = 3)))+
  theme(legend.title = element_text(size=13),plot.title = element_text(face = "bold",hjust = 0.5,size = 16))

#### 
ggarrange(DimPlot(mouse,group.by = "cell_class",label = T)+paletteer::scale_color_paletteer_d("ggthemes::Miller_Stone")+
            guides(color = guide_legend(ncol = 1,override.aes = list(size = 3))),
           DimPlot(mouse,group.by = c("cell_type"),label = T)+paletteer::scale_color_paletteer_d("trekcolors::lcars_series")+
            guides(color = guide_legend(ncol = 1,override.aes = list(size = 3))),
          nrow = 1)
save.image("afterAnno.RData")
saveRDS(mouse,"afterAnnp_SeuratObj.RDS")

mouse$cell_type = factor(mouse$cell_type,levels = c("GC B","naive B","Plasma B","Dendritic","Granulocyte","Macrophage","Mast",
                                  "CD4+ effective" ,"CD4+ memory","CD4+ Treg","CD8+ cytotoxic","Proliferate T","DNT","Th17",
                                  "ILC","NK","NKT",
                                  "Endothelial","Fibroblast","Doublet","LowQ"))

VlnPlot(mouse,group.by = "cell_type",stack = T,flip = T,
        features = c("Cd79a","Mki67","Ighm","Sdc1","Itgax","Hdc","C1qa","C1qb","Mcpt8",
                     "Cd4","Icos","Ccr7","Foxp3","Cd8a","Gzma","Ikzf2", "Il17a",
                     "Gata3","Nkg7","Il23r",
                     "Fabp4","Col1a2","nCount_RNA","percent_mt"))+
  scale_fill_paletteer_d("colorBlindness::SteppedSequential5Steps")+NoLegend()

write.csv(all_markers %>% group_by(cluster) %>% slice_head(n=30),file = "top30_res0.5.csv",row.names = F,quote = F)


# VlnPlot(subset(mouse,cell_class=="Mye"),stack = T,flip = T,
#         features = c("Ly6c1","Ccr2","Cd14","Cd68","Cx3cr1","Cd36","Cd163",
#                      "Itgax","Xcr1","Clec9a","Sirpa","Cd209a","Tlr7","Irf8",
#                      "Mrc1","Arg1"))


