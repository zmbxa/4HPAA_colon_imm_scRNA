### heatmap with and without ribo-genes
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
# load("afterAnno.RData")
mouse = readRDS("refinedAnno_SeuratObj.RDS")
mouse$condition = ifelse(mouse$treatment=="HFD","MOCK","4HPAA")

## DEG for all Macrophage
Macro_DEG = FindMarkers(subset(mouse,cell_type_refine %in% c("Macro-1","Macro-2")),ident.1 = "4HPAA",ident.2 = "MOCK")
Macro_DEG = data.frame(cell_type = "Macrophage",gene=rownames(Macro_DEG),Macro_DEG,change = ifelse(Macro_DEG$p_val_adj < 0.05,ifelse(Macro_DEG$avg_log2FC>0,"Up","Down"),"Stable"))
ggarrange(DoHeatmap(subset(mouse,cell_type_refine %in% c("Macro-1","Macro-2")),features = filter(Macro_DEG,change!="Stable")$gene,
                    group.by = "type_treat",size=4, group.bar.height = 0.006)+ggtitle("Macrophage DEG, split subtypes heatmap","with ribo-genes"),
          DoHeatmap(subset(mouse,cell_type_refine %in% c("Macro-1","Macro-2")),features = filter(Macro_DEG,change!="Stable")$gene,
                    group.by = "condition",size=4, group.bar.height = 0.006)+ggtitle("Macrophage DEG, condition","with ribo-genes"))

ggarrange(DoHeatmap(subset(mouse,cell_type_refine %in% c("Macro-1","Macro-2")),features = filter(Macro_DEG,change!="Stable" & !grepl("^Rp[ls]",gene))$gene,
                    group.by = "type_treat",size=4, group.bar.height = 0.006)+ggtitle("Macrophage DEG, split subtypes heatmap","exclude ribo-genes"),
          DoHeatmap(subset(mouse,cell_type_refine %in% c("Macro-1","Macro-2")),features = filter(Macro_DEG,change!="Stable" & !grepl("^Rp[ls]",gene))$gene,
                    group.by = "condition",size=4, group.bar.height = 0.006)+ggtitle("Macrophage DEG, condition","exclude ribo-genes"))

