rm(list = ls());gc()
setwd("~/projects/immune_scRNA_TaoLab/analysis/")
load("DEG.RData")

### GO enrich
gene_SYMBOL = read.csv("~/annotations/mm10_geneID2SYMBOL.csv")

## clusterProfiler
library(clusterProfiler)
cpGO_class = cpGO_type = list()

for (type in unique(DEG_cellClass_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")) > 0){
    enGO = enrichGO(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",readable = T) 
    enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    cpGO_class[[paste0(type,"_Up")]] = enGO@result
  }
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")) > 0){
    enGO = enrichGO(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",readable = T) 
    enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
    cpGO_class[[paste0(type,"_Down")]] = enGO@result
  }
}

for (type in unique(DEG_cellTypes_wilcoxon$cell_type)){
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  # Up
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")) > 0){
    enGO = try(enrichGO(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",readable = T))
    if (!inherits(enGO, "try-error") & !is.null(enGO)){
      enGO@result$FoldEnrichment = apply(enGO@result, 1, function(x){eval(parse(text = x["GeneRatio"]))/eval(parse(text = x["BgRatio"]))})
      cpGO_type[[paste0(type,"_Up")]] = enGO@result
    }
  }
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")) > 0){
    enGO = try(enrichGO(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")$gene,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",readable = T) )
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

ggplot(filter(cpGOsig_subtype,change=="Up") %>% group_by(diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradient2(mid="orange",high="firebrick")+ggtitle("Up-gene GO enrichment, subtype-level")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

ggplot(filter(cpGOsig_subtype,change=="Down") %>% group_by(diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=Description,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradientn(colours = paletteer::paletteer_d("ggsci::blue_material")[-1])+ggtitle("Down-gene GO enrichment, subtype-level")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

write.csv(cpGOsig_class,"GO_class_signif.csv",quote = F,row.names = F)
write.csv(cpGOsig_subtype,"GO_subtype_signif.csv",quote = F,row.names = F)
openxlsx::write.xlsx(cpGO_class,file = "GO_class_all.xlsx")
openxlsx::write.xlsx(cpGO_type,file = "GO_subtype_all.xlsx")

### KEGG
gene_SYMBOL$ENTREZ = mapIds(org.Mm.eg.db,keys = gene_SYMBOL$gene_id_trim,keytype = "ENSEMBL",column = c("ENTREZID"))

cpKEGG_class = cpKEGG_type = list()
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

cpKEGGsig_class = Reduce(bind_rows,mapply(cpKEGG_class,names(cpKEGG_class), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F)) %>% separate(diffType,into = c("cellClass","change"),remove = F) %>% mutate(TermName = gsub(" -.*","",Description))

cpKEGGsig_subtype = Reduce(bind_rows,mapply(cpKEGG_type,names(cpKEGG_type), FUN=function(N,name){
  if(nrow(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4))>0 )
    data.frame(filter(N,p.adjust<0.05 & nchar(gsub("\\/.*","",BgRatio)) < 4),diffType=name)
},SIMPLIFY = F)) %>% separate(diffType,into = c("cellType","change"),remove = F,sep = "_") %>% mutate(TermName = gsub(" -.*","",Description))

ggplot(filter(cpKEGGsig_subtype,change=="Up") %>% group_by(diffType) %>% slice_head(n=10) %>% as.data.frame(),aes(x=diffType,y=TermName,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradient2(mid="orange",high="firebrick")+ggtitle("Up-gene KEGG enrichment, subtype-level")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))

ggplot(filter(cpKEGGsig_subtype,change=="Down") %>% group_by(diffType) %>% slice_head(n=5) %>% as.data.frame(),aes(x=diffType,y=TermName,fill=-log(p.adjust),size=FoldEnrichment))+
  geom_point(color="grey",shape=21)+scale_fill_gradientn(colours = paletteer::paletteer_d("ggsci::blue_material")[-1])+ggtitle("Down-gene KEGG enrichment, subtype-level")+
  theme_bw()+theme(plot.title = element_text(face="bold"),axis.text.x = element_text(angle = 30,hjust = 0.8,vjust = 0.8))


### for each subtype
barplot_GO = sapply(unique(cpGOsig_subtype$cellType), function(type){
  ggplot(filter(cpGOsig_subtype,cellType == type) %>% group_by(change) %>% slice_head(n=40),aes(y=factor(Description,levels = rev(unique(filter(cpGOsig_subtype,cellType == type)$Description))),x=-log(p.adjust),fill=change))+geom_col(width = 0.8)+
    facet_grid(change ~ .,space = "free",scales = "free")+scale_fill_manual(values = c(Down="royalblue",Up="brown"))+ylab("TermName")+
    ggtitle(paste0(type," GO enriched Terms"))
},simplify = F)

barplot_KEGG = sapply(unique(cpKEGGsig_subtype$cellType), function(type){
  ggplot(filter(cpKEGGsig_subtype,cellType == type),aes(y=factor(TermName,levels = rev(unique(filter(cpKEGGsig_subtype,cellType == type)$TermName))),x=-log(p.adjust),fill=change))+
    geom_col(width = 0.8)+facet_grid(change ~ .,space = "free",scales = "free")+scale_fill_manual(values = c(Down="royalblue",Up="brown"))+ylab("TermName")+
    ggtitle(paste0(type," KEGG enriched Terms"))
},simplify = F)

write.csv(cpKEGGsig_class,"KEGG_class_signif.csv",quote = F,row.names = F)
write.csv(cpKEGGsig_subtype,"KEGG_subtype_signif.csv",quote = F,row.names = F)
openxlsx::write.xlsx(cpKEGG_class,file = "KEGG_class_all.xlsx")
openxlsx::write.xlsx(cpKEGG_type,file = "KEGG_subtype_all.xlsx")
