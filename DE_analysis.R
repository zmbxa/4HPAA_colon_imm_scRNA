###
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
# load("afterAnno.RData")
mouse = readRDS("refinedAnno_SeuratObj.RDS")

as.data.frame(prop.table(table(mouse$cell_class,mouse$treatment),margin = 2)) %>% ggplot(aes(x=Var2,y=Freq*100,fill = Var1))+
  geom_col()+paletteer::scale_fill_paletteer_d("ggthemes::Miller_Stone")+theme_bw()+ggtitle("Cell proportion")+
  xlab("treatment")+ylab("Cell proportion")+labs(fill = "cell_class")+guides(fill = guide_legend(ncol = 1,override.aes = list(size = 3)))+
  theme(legend.title = element_text(size=13),plot.title = element_text(face = "bold",hjust = 0.5,size = 16))


Idents(mouse) = mouse$treatment
DEG_cellTypes_wilcoxon = Reduce(bind_rows,sapply(unique(mouse$cell_type_refine), function(ct){
  print(ct)
  out = FindMarkers(subset(mouse,cell_type_refine == ct),ident.1 = "4HPAA",ident.2 = "HFD")
  out = data.frame(cell_type = ct,gene=rownames(out),out,change = ifelse(out$p_val_adj < 0.05,ifelse(out$avg_log2FC>0,"Up","Down"),"Stable"))
  return(out)
},simplify = F))

DEG_cellClass_wilcoxon = Reduce(bind_rows,sapply(unique(mouse$cell_class), function(ct){
  print(ct)
  out = FindMarkers(subset(mouse,cell_class == ct),ident.1 = "4HPAA",ident.2 = "HFD")
  out = data.frame(cell_type = ct,gene=rownames(out),out,change = ifelse(out$p_val_adj < 0.05,ifelse(out$avg_log2FC>0,"Up","Down"),"Stable"))
  return(out)
},simplify = F))

# volcano
### cell types
load("~/pipelines/RNAseq/analysis/get_edgeR_DE.RData")
ct_to_plot = names(table(filter(DEG_cellTypes_wilcoxon,change != "Stable")$cell_type)[table(filter(DEG_cellTypes_wilcoxon,change != "Stable")$cell_type)>200])
ggarrange(plotlist = sapply(ct_to_plot, function(name){
  plot_volcano(filter(DEG_cellTypes_wilcoxon,cell_type == name),FDR_c = "p_val_adj",logFC_c = "avg_log2FC",label_c = "gene",label_num = 20,title = name)
},simplify = F),common.legend = T,legend = "right")
### cell class
ct_to_plot = unique(mouse$cell_class)[1:7]
ggarrange(plotlist = sapply(ct_to_plot, function(name){
  plot_volcano(filter(DEG_cellClass_wilcoxon,cell_type == name),FDR_c = "p_val_adj",logFC_c = "avg_log2FC",label_c = "gene",label_num = 20,title = name)
},simplify = F),common.legend = T,legend = "right",ncol = 4,nrow = 2)

write.csv(DEG_cellClass_wilcoxon,"DEG_cellClass.csv",quote = F,row.names = F)
write.csv(DEG_cellTypes_wilcoxon,"DEG_subtypes.csv",quote = F,row.names = F)
save.image("DEG.RData")

### heatmap
mouse$class_treat = paste(mouse$cell_class,mouse$treatment,sep = "_")
mouse$type_treat = paste(mouse$cell_type_refine,mouse$treatment,sep = "_")

pb_classTreat = AverageExpression(mouse,group.by = c("cell_class","treatment")) %>% as.data.frame
pb_typeTreat = AverageExpression(mouse,group.by = c("cell_type_refine","treatment")) %>% as.data.frame
colnames(pb_typeTreat) = (table(mouse$cell_type_refine,mouse$treatment) %>% as.data.frame() %>% 
                            arrange(factor(Var1,levels = levels(mouse$cell_type_refine))) %>% mutate(name=paste(Var1,Var2,sep="_")))$name
colnames(pb_classTreat) = gsub("RNA.","",colnames(pb_classTreat))

type_heatmap = sapply(unique(DEG_cellTypes_wilcoxon$cell_type), function(ct){
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == ct &change!="Stable")) > 0){
    print(ct)
    return(DoHeatmap(subset(mouse,cell_type_refine == ct),features = filter(DEG_cellTypes_wilcoxon,cell_type == ct &change!="Stable")$gene,
            group.by = "type_treat",size=4))
    }
},simplify = F)
names(type_heatmap) = unique(DEG_cellTypes_wilcoxon$cell_type)


### pseudobulk heatmap,cell subtype
pheatmap(pb_typeTreat[unique(filter(DEG_cellTypes_wilcoxon,change!="Stable" & !grepl("Doub|LowQ",cell_type) & abs(avg_log2FC) > 0.8)$gene),-grep("LowQ|Doublet",colnames(pb_typeTreat))],
         cluster_cols = F,scale = "row",annotation_col = data.frame(treat = colnames(pb_typeTreat),row.names = colnames(pb_typeTreat)) %>% mutate(treat = gsub(".*_","",treat)),
         main = "DEGene heatmap (cell type pseudo-bulk)",color = colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100))
pheatmap(pb_typeTreat[unique(filter(DEG_cellTypes_wilcoxon,change!="Stable" & !grepl("Doub|LowQ",cell_type) & abs(avg_log2FC) > 0.8)$gene),-grep("LowQ|Doublet",colnames(pb_typeTreat))] %>% t,
         cluster_rows = F,scale = "column",annotation_row = data.frame(treat = colnames(pb_typeTreat),row.names = colnames(pb_typeTreat)) %>% mutate(treat = gsub(".*_","",treat)),
         main = "DEGene heatmap (cell type pseudo-bulk)",color = colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100))

# pheatmap(pb_typeTreat[unique(filter(DEG_cellTypes_wilcoxon,change != "Stable" & abs(avg_log2FC)>0.8)$gene),rownames(data.frame(name = colnames(pb_typeTreat),row.names = colnames(pb_typeTreat)) %>% mutate(name = gsub(".*_","",name)) %>% arrange(name))],
#          cluster_cols = F,scale = "row",annotation_col = data.frame(name = colnames(pb_typeTreat),row.names = colnames(pb_typeTreat)) %>% mutate(name = gsub(".*_","",name)),
#          main = "DEGene heatmap (cell type pseudo-bulk)",color = colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100),gaps_col = 24)

pheatmap(pb_classTreat[unique(filter(DEG_cellClass_wilcoxon,change!="Stable" & !grepl("Doub|LowQ",cell_type) & abs(avg_log2FC) > 0.8)$gene),-grep("LowQ|Doublet",colnames(pb_classTreat))],
         cluster_cols = F,scale = "row",annotation_col = data.frame(treat = colnames(pb_classTreat),row.names = colnames(pb_classTreat)) %>% mutate(treat = gsub(".*_","",treat)),
         main = "DEGene heatmap (cell class pseudo-bulk)",color = colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100))
