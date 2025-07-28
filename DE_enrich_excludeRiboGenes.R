###  heatmap without Ribosomal genes (cell type&conditions)
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
# mouse=readRDS("refinedAnno_SeuratObj.RDS")
load("DE_enrich_excludeRiboGenes.RData")

## DEG list, exclude ribo- genes
DEG_cellTypes_wilcoxon = Reduce(bind_rows,sapply(unique(mouse$cell_type_refine), function(ct){
  print(ct)
  out = FindMarkers(subset(mouse,cell_type_refine == ct),ident.1 = "4HPAA",ident.2 = "HFD",features =rownames(mouse)[-grep("^Rp[sl]",rownames(mouse))])
  out = data.frame(cell_type = ct,gene=rownames(out),out,change = ifelse(out$p_val_adj < 0.05,ifelse(out$avg_log2FC>0,"Up","Down"),"Stable"))
  return(out)
},simplify = F))

DEG_cellClass_wilcoxon = Reduce(bind_rows,sapply(unique(mouse$cell_class), function(ct){
  print(ct)
  out = FindMarkers(subset(mouse,cell_class == ct),ident.1 = "4HPAA",ident.2 = "HFD",features = rownames(mouse)[-grep("^Rp[sl]",rownames(mouse))])
  out = data.frame(cell_type = ct,gene=rownames(out),out,change = ifelse(out$p_val_adj < 0.05,ifelse(out$avg_log2FC>0,"Up","Down"),"Stable"))
  return(out)
},simplify = F))


### heatmap
mouse$class_treat = paste(mouse$cell_class,mouse$treatment,sep = "_")
mouse$type_treat = paste(mouse$cell_type_refine,mouse$treatment,sep = "_")

pb_classTreat = AverageExpression(mouse,group.by = c("cell_class","treatment")) %>% as.data.frame
pb_typeTreat = AverageExpression(mouse,group.by = c("cell_type_refine","treatment")) %>% as.data.frame
colnames(pb_typeTreat) = (table(mouse$cell_type_refine,mouse$treatment) %>% as.data.frame() %>% 
                            arrange(factor(Var1,levels = levels(mouse$cell_type_refine))) %>% mutate(name=paste(Var1,Var2,sep="_")))$name
colnames(pb_classTreat) = gsub("RNA.","",colnames(pb_classTreat))

## store heatmap into list 
type_heatmap = sapply(unique(DEG_cellTypes_wilcoxon$cell_type), function(ct){
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == ct &change!="Stable")) > 0){
    print(ct)
    return(DoHeatmap(subset(mouse,cell_type_refine == ct),features = filter(DEG_cellTypes_wilcoxon,cell_type == ct &change!="Stable")$gene,
                     group.by = "type_treat",size=4))
  }
},simplify = F)
names(type_heatmap) = unique(DEG_cellTypes_wilcoxon$cell_type)

### pseudobulk heatmap
pb_typeTreat = AverageExpression(mouse,group.by = c("cell_type_refine","treatment")) %>% as.data.frame
colnames(pb_typeTreat) = (table(mouse$cell_type_refine,mouse$treatment) %>% as.data.frame() %>% 
                            arrange(factor(Var1,levels = levels(mouse$cell_type_refine))) %>% mutate(name=paste(Var1,Var2,sep="_")))$name
pheatmap(pb_typeTreat[unique(filter(DEG_cellTypes_wilcoxon,change!="Stable" & !grepl("Doub|LowQ",cell_type) & abs(avg_log2FC) > 0.8)$gene),-grep("LowQ|Doublet",colnames(pb_typeTreat))] %>% t,
         cluster_rows = F,scale = "column",annotation_row = data.frame(treat = colnames(pb_typeTreat),row.names = colnames(pb_typeTreat)) %>% mutate(treat = gsub(".*_","",treat)),
         main = "DEGene heatmap (cell type pseudo-bulk)",color = colorRampPalette(colors = c('#11427C','white','#C31E1F'))(100))
# ggsave(filename = "cellTypes_DEG_heatmap_pseudobulk_excludeRibo.pdf",height = 20,width = 200,limitsize = F)


## for subtypes
type_heatmap$`ILC-2` + ggtitle("ILC-2 diff-gene heatmap")+theme(plot.title = element_text(hjust = 0.5,face='bold',size=15))

## save result
write.csv(DEG_cellClass_wilcoxon,"DEG_cellClass_noRibo.csv",quote = F,row.names = F)
write.csv(DEG_cellTypes_wilcoxon,"DEG_subtypes_noRibo.csv",quote = F,row.names = F)

##########
#### enrichment GO + KEGG
gene_SYMBOL = read.csv("~/annotations/mm10_geneID2SYMBOL.csv")
library(clusterProfiler)
cpGO_class = cpGO_type = list()

for (type in unique(DEG_cellClass_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")) > 0){
    enGO = enrichGO(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "all",readable = T) 
    enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    cpGO_class[[paste0(type,"_Up")]] = enGO@result
  }
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")) > 0){
    enGO = enrichGO(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "all",readable = T) 
    enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    cpGO_class[[paste0(type,"_Down")]] = enGO@result
  }
}
for (type in unique(DEG_cellTypes_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")) > 0){
    enGO = try(enrichGO(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "all",readable = T))
    if (!inherits(enGO, "try-error") & !is.null(enGO)){
      enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
      cpGO_type[[paste0(type,"_Up")]] = enGO@result
    }
  }
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")) > 0){
    enGO = try(enrichGO(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "all",readable = T) )
    if (!inherits(enGO, "try-error") & !is.null(enGO)){
      enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
      cpGO_type[[paste0(type,"_Down")]] = enGO@result
    }
  }
}

cpGOsig_class = Reduce(bind_rows,mapply(cpGO_class,names(cpGO_class), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F))%>% separate(diffType,into = c("cellClass","change"),remove = F)

cpGOsig_subtype = Reduce(bind_rows,mapply(cpGO_type,names(cpGO_type), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F)) %>% separate(diffType,into = c("cellType","change"),remove = F,sep = "_")

ggplot(filter(cpGOsig_subtype,change=="Up",!grepl("ribosom",Description,ignore.case = T)) %>% group_by(ONTOLOGY,diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradient2(mid="orange",high="firebrick")+ggtitle("Up-gene GO enrichment, subtype-level")+facet_grid(ONTOLOGY~.,scales = "free_y",space = "free")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

ggplot(filter(cpGOsig_subtype,change=="Down") %>% group_by(ONTOLOGY,diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradientn(colours = paletteer::paletteer_d("ggsci::blue_material")[-1])+ggtitle("Down-gene GO enrichment, subtype-level")+facet_grid(ONTOLOGY~.,scales = "free_y",space = "free")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

write.csv(cpGOsig_class,"GO_class_signif_noRibo.csv",quote = F,row.names = F)
write.csv(cpGOsig_subtype,"GO_subtype_signif_noRibo.csv",quote = F,row.names = F)
openxlsx::write.xlsx(cpGO_class,file = "GO_class_all_noRibo.xlsx")
openxlsx::write.xlsx(cpGO_type,file = "GO_subtype_all_noRibo.xlsx")


##### KEGG
library(org.Mm.eg.db)
gene_SYMBOL$ENTREZ = AnnotationDbi::mapIds(org.Mm.eg.db,keys = gene_SYMBOL$gene_id_trim,keytype = "ENSEMBL",column = c("ENTREZID"))

cpKEGG_class = cpKEGG_type = list()
# for cell class
for (type in unique(DEG_cellClass_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")) > 0){
    enKEGG = enrichKEGG(na.omit(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")$gene)$ENTREZ),organism = "mouse") 
    enKEGG@result$FoldEnrichment = apply(enKEGG@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    enKEGG@result$geneSYMBOLs = apply(enKEGG@result,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneID'],"\\/")),gene_SYMBOL$ENTREZ),"gene_name"],collapse = ", ")})
    cpKEGG_class[[paste0(type,"_Up")]] = enKEGG@result
  }
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")) > 0){
    enKEGG = enrichKEGG(na.omit(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")$gene)$ENTREZ),organism = "mouse") 
    enKEGG@result$FoldEnrichment = apply(enKEGG@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    enKEGG@result$geneSYMBOLs = apply(enKEGG@result,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneID'],"\\/")),gene_SYMBOL$ENTREZ),"gene_name"],collapse = ", ")})
    cpKEGG_class[[paste0(type,"_Down")]] = enKEGG@result
  }
}
# for cell types
for (type in unique(DEG_cellTypes_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")) > 0){
    enKEGG = try(enrichKEGG(na.omit(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")$gene)$ENTREZ),organism = "mouse"))
    if (!inherits(enKEGG, "try-error") & !is.null(enKEGG)){
      enKEGG@result$FoldEnrichment = apply(enKEGG@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
      enKEGG@result$geneSYMBOLs = apply(enKEGG@result,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneID'],"\\/")),gene_SYMBOL$ENTREZ),"gene_name"],collapse = ", ")})
      cpKEGG_type[[paste0(type,"_Up")]] = enKEGG@result
    }
  }
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")) > 0){
    enKEGG = try(enrichKEGG(na.omit(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")$gene)$ENTREZ),organism = "mouse"))
    if (!inherits(enKEGG, "try-error") & !is.null(enKEGG)){
      enKEGG@result$FoldEnrichment = apply(enKEGG@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
      enKEGG@result$geneSYMBOLs = apply(enKEGG@result,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneID'],"\\/")),gene_SYMBOL$ENTREZ),"gene_name"],collapse = ", ")})
      cpKEGG_type[[paste0(type,"_Down")]] = enKEGG@result
    }
  }
}
## signif terms
cpKEGGsig_class = Reduce(bind_rows,mapply(cpKEGG_class,names(cpKEGG_class), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F)) %>% separate(diffType,into = c("cellClass","change"),remove = F) %>% mutate(TermName = gsub(" -.*","",Description))

cpKEGGsig_subtype = Reduce(bind_rows,mapply(cpKEGG_type,names(cpKEGG_type), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F)) %>% separate(diffType,into = c("cellType","change"),remove = F,sep = "_") %>% mutate(TermName = gsub(" -.*","",Description))

### saving KEGG tables
openxlsx::write.xlsx(cpKEGG_class,file = "KEGG_class_all_noRibo.xlsx")
openxlsx::write.xlsx(cpKEGG_type,file = "KEGG_subtype_all_noRibo.xlsx")

### GO BP/MF/CC for each types
ggplot(filter(cpGOsig_subtype,cellType == type) %>% group_by(ONTOLOGY,change) %>% slice_head(n=15),aes(y=factor(Description,levels = rev(unique(filter(cpGOsig_subtype,cellType == type)$Description))),x=ifelse(change=="Down",log(p.adjust),-log(p.adjust)),fill=change))+
  geom_col(width = 0.8,mapping = aes(alpha=FoldEnrichment))+  facet_grid(ONTOLOGY ~ .,space = "free",scales = "free_y",labeller = labeller(ONTOLOGY = label_wrap_gen(width = 10)))+
  scale_x_continuous(labels = function(x) abs(x))+  scale_fill_manual(values = c(Down="royalblue",Up="brown"))+geom_vline(xintercept = 0,lty="dashed",color="grey37")+
  xlab("-log(p.adj)")+ylab("TermName")+  ggtitle(paste0(type," GO enriched Terms"))+theme_bw()+theme(plot.title = element_text(face="bold",size=14))
barplot_GO_type = sapply(unique(cpGOsig_subtype$cellType), function(type){
  ggplot(filter(cpGOsig_subtype,cellType == type) %>% group_by(ONTOLOGY,change) %>% slice_head(n=15),aes(y=factor(Description,levels = rev(unique(filter(cpGOsig_subtype,cellType == type)$Description))),x=ifelse(change=="Down",log(p.adjust),-log(p.adjust)),fill=change))+
    geom_col(width = 0.8,mapping = aes(alpha=FoldEnrichment))+  facet_grid(ONTOLOGY ~ .,space = "free",scales = "free_y",labeller = labeller(ONTOLOGY = label_wrap_gen(width = 10)))+
    scale_x_continuous(labels = function(x) abs(x))+  scale_fill_manual(values = c(Down="royalblue",Up="brown"))+geom_vline(xintercept = 0,lty="dashed",color="grey37")+
    xlab("-log(p.adj)")+ylab("TermName")+  ggtitle(paste0(type," GO enriched Terms"))+theme_bw()+theme(plot.title = element_text(face="bold",size=14))
},simplify = F)

# combining Dotplot GO
ggplot(filter(cpGOsig_subtype,change=="Up") %>% group_by(ONTOLOGY,diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradient2(mid="orange",high="firebrick")+ggtitle("Up-gene GO enrichment, subtype-level")+
  facet_grid(ONTOLOGY ~ .,space = "free",scales = "free_y",labeller = labeller(ONTOLOGY = label_wrap_gen(width = 10)))+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

ggplot(filter(cpGOsig_subtype,change=="Down") %>% group_by(ONTOLOGY,diffType) %>% slice_head(n=3) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradientn(colours = paletteer::paletteer_d("ggsci::blue_material")[-1])+ggtitle("Down-gene GO enrichment, subtype-level")+
  facet_grid(ONTOLOGY ~ .,space = "free",scales = "free_y",labeller = labeller(ONTOLOGY = label_wrap_gen(width = 10)))+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))


### KEGG for each types
KEGG_pw = read.csv("KEGG_terms/KEGG_pathway_ko_uniq.txt",sep = "\t",colClasses = rep("character",10))
KEGG_class = KEGG_pw %>% dplyr::select(1:4) %>% unique() %>% filter(level2_pathway_name %in% gsub("mmu","",(c(cpKEGGsig_class$ID,cpKEGGsig_subtype$ID))))
KEGG_class$ID = paste0("mmu",KEGG_class$level2_pathway_name)

cpKEGGsig_subtype_toplot = merge(cpKEGGsig_subtype,KEGG_class,by="ID")
cpKEGGsig_class_toplot = merge(cpKEGGsig_class,KEGG_class,by="ID")

ggplot(filter(cpKEGGsig_subtype_toplot,cellType == type),aes(y=factor(TermName,levels = rev(unique(filter(cpKEGGsig_subtype_toplot,cellType == type)$TermName))),x=ifelse(change=="Down",log(p.adjust),-log(p.adjust)),fill=change))+
  geom_col(width = 0.8,mapping = aes(alpha=FoldEnrichment))+geom_vline(xintercept = 0,lty="dashed",color="grey37")+xlab("-log(p.adj)")+
  facet_grid(level1_pathway_id ~ .,space = "free",scales = "free_y",labeller = labeller(level1_pathway_id = label_wrap_gen(width = 10)))+
  scale_fill_manual(values = c(Down="royalblue",Up="brown"))+ylab("TermName")+scale_x_continuous(labels = function(x) abs(x))+
  ggtitle(paste0(type," KEGG enriched Terms"))+theme_bw()+theme(plot.title = element_text(face="bold",size=14))

barplot_KEGG_type = sapply(unique(cpKEGGsig_subtype$cellType), function(type){
  ggplot(filter(cpKEGGsig_subtype_toplot,cellType == type),aes(y=factor(TermName,levels = rev(unique(filter(cpKEGGsig_subtype_toplot,cellType == type)$TermName))),x=ifelse(change=="Down",log(p.adjust),-log(p.adjust)),fill=change))+
    geom_col(width = 0.8,mapping = aes(alpha=FoldEnrichment))+geom_vline(xintercept = 0,lty="dashed",color="grey37")+xlab("-log(p.adj)")+
    facet_grid(level1_pathway_id ~ .,space = "free",scales = "free_y",labeller = labeller(level1_pathway_id = label_wrap_gen(width = 10)))+
    scale_fill_manual(values = c(Down="royalblue",Up="brown"))+ylab("TermName")+scale_x_continuous(labels = function(x) abs(x))+
    ggtitle(paste0(type," KEGG enriched Terms"))+theme_bw()+theme(plot.title = element_text(face="bold",size=14))
},simplify = F)


write.csv(cpKEGGsig_class_toplot,"KEGG_class_signif_noRibo.csv",quote = F,row.names = F)
write.csv(cpKEGGsig_subtype_toplot,"KEGG_subtype_signif_noRibo.csv",quote = F,row.names = F)
# openxlsx::write.xlsx(cpKEGG_class,file = "KEGG_class_all_noRibo.xlsx")
# openxlsx::write.xlsx(cpKEGG_type,file = "KEGG_subtype_all_noRibo.xlsx")

### combining Dotplot of KEGG
ggplot(filter(cpKEGGsig_subtype_toplot,change=="Up") %>% group_by(diffType) %>% slice_head(n=10) %>% as.data.frame(),aes(x=diffType,y=TermName,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradient2(mid="orange",high="firebrick")+ggtitle("Up-gene KEGG enrichment, subtype-level")+
  facet_grid(level1_pathway_id ~ .,space = "free",scales = "free_y",labeller = labeller(level1_pathway_id = label_wrap_gen(width = 10)))+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

ggplot(filter(cpKEGGsig_subtype_toplot,change=="Down") %>% group_by(diffType) %>% slice_head(n=10) %>% as.data.frame(),aes(x=diffType,y=TermName,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradientn(colours = paletteer::paletteer_d("ggsci::blue_material")[-1])+ggtitle("Down-gene KEGG enrichment, subtype-level")+
  facet_grid(level1_pathway_id ~ .,space = "free",scales = "free_y",labeller = labeller(level1_pathway_id = label_wrap_gen(width = 10)))+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))



save.image("DE_enrich_excludeRiboGenes.RData")

#############
### test word cloud clustering for GO terms
# testGO_ILC2 = simplifyEnrichment::GO_similarity((cpGO_type$`ILC-2_Up` %>% filter(pvalue < 0.05))$ID)
# sG = simplifyEnrichment::simplifyGO(testGO_ILC2)





