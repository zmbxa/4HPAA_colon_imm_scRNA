devtools::install_github("arc85/celltalker")  

library(celltalker)
library(Seurat)
library(dplyr)
library(magrittr)

## Run celltalker
mouse_lr = read.delim("/mnt/transposon1/zhangyanxiaoLab/niuyuxiao/projects/immune_scRNA_TaoLab/analysis/cell_communication/mouse_lr_pair.txt")
mouse_pairs = setNames(mouse_lr[,c(2,3,1)],colnames(ramilowski_pairs))
# mouse_pairs = ramilowski_pairs
# mouse_pairs$ligand = stringr::str_to_title(mouse_pairs$ligand)
# mouse_pairs$receptor = stringr::str_to_title(mouse_pairs$receptor)
# mouse_pairs$pair = paste0(mouse_pairs$ligand,"_",mouse_pairs$receptor)

immune_talk <- celltalk(input_object=subset(mouse,cell_type_refine!="LowQ" & cell_type_refine!="Doublet"),metadata_grouping="cell_type_refine",
                        number_cells_required = 50, ligand_receptor_pairs = mouse_pairs)    

top_stats <- immune_talk  %>%
  mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
  filter(fdr < 0.05) %>%
  group_by(cell_type1) %>%
  top_n(3, interact_ratio) %>%
  ungroup()

library(RColorBrewer)   
n_colors <- max(6, length(unique(mouse$cell_type_refine)))  
palette <- rep(brewer.pal(n = 8, name = "Set2"), length.out = n_colors) 
colors_use <- palette[1:n_colors]  # 
circos_plot(ligand_receptor_frame = top_stats ,   #
            cell_group_colors = colors_use,
            ligand_color = "blue",
            receptor_color = "red",
            cex_outer = 1.4,
            cex_inner = 0.4)
write.csv(top_stats,file="cell_communication/cellCommunication_sig.csv",quote = F,row.names = F)

save.image("cell_communication/celltalker.RData")
