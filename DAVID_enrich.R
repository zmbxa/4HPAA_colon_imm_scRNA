## 1. DAVID
library(reticulate)
library(dplyr)

# reticulate::source_python("~/pipelines/miniscripts/DAVIDChart4Rreticulate.py")
source_python("~/pipelines/miniscripts/DAVIDChart4Rreticulate.py")

DAVID_cellClass = list()
for (type in unique(DEG_cellClass_wilcoxon$cell_type)) {
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  #Up
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")) > 0){
    DAVID = david_gene_enrichment(gene_list = paste(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Up")$gene)$gene_id_trim,collapse = ", "))
    DAVID$geneSYMBOLs = apply(DAVID,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneIds'],", ")),gene_SYMBOL$gene_id_trim),"gene_name"],collapse = ", ")})
    DAVID_cellClass[[paste0(type,"_Up")]] = DAVID
  }
  if(nrow(filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")) > 0){
    DAVID = david_gene_enrichment(gene_list = paste(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellClass_wilcoxon,cell_type == type & change=="Down")$gene)$gene_id_trim,collapse = ", "))
    DAVID$geneSYMBOLs = apply(DAVID,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneIds'],", ")),gene_SYMBOL$gene_id_trim),"gene_name"],collapse = ", ")})
    DAVID_cellClass[[paste0(type,"_Down")]] = DAVID
  }
}
DAVIDsig_cellClass = Reduce(bind_rows,sapply(names(DAVID_cellClass),simplify = F, function(N){
  if(nrow(filter(DAVID_cellClass[[N]],EASEBonferroni < 0.05)) > 0)  data.frame(filter(DAVID_cellClass[[N]],EASEBonferroni < 0.05 & popHits < 1000)%>%  distinct(termName, .keep_all = TRUE),cluster = N)
}))

DAVID_cellType = list()
for (type in unique(DEG_cellTypes_wilcoxon$cell_type)) {
  if(grepl("LowQ|Doublet",type)) next
  print(type)
  #Up
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")) > 0){
    DAVID = try(david_gene_enrichment(gene_list = paste(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Up")$gene)$gene_id_trim,collapse = ", ")))
    if (!inherits(DAVID, "try-error") & nrow(DAVID) != 0 ){
      DAVID$geneSYMBOLs = apply(DAVID,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneIds'],", ")),gene_SYMBOL$gene_id_trim),"gene_name"],collapse = ", ")})
      DAVID_cellType[[paste0(type,"_Up")]] = DAVID
    }
  }
  if(nrow(filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")) > 0){
    DAVID = try(david_gene_enrichment(gene_list = paste(filter(gene_SYMBOL,gene_name %in% filter(DEG_cellTypes_wilcoxon,cell_type == type & change=="Down")$gene)$gene_id_trim,collapse = ", ")))
    if (!inherits(DAVID, "try-error") & nrow(DAVID) != 0){
      DAVID$geneSYMBOLs = apply(DAVID,1,function(x){paste(gene_SYMBOL[match(unlist(strsplit(x['geneIds'],", ")),gene_SYMBOL$gene_id_trim),"gene_name"],collapse = ", ")})
      DAVID_cellType[[paste0(type,"_Down")]] = DAVID
    }
  }
}
DAVIDsig_subtype = Reduce(bind_rows,sapply(names(DAVID_cellType),simplify = F, function(N){
  if(nrow(filter(DAVID_cellType[[N]],EASEBonferroni < 0.05 & popHits < 1000)) > 0)  
    return(data.frame(filter(DAVID_cellType[[N]],EASEBonferroni < 0.05 & popHits < 1000)%>%  distinct(termName, .keep_all = TRUE),cluster = N))
}))

save(DAVID_cellClass,DAVIDsig_cellClass,DAVID_cellType,DAVIDsig_subtype,file = "DAVID_enrichment.RData")
