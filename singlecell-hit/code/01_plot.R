library(Seurat)
library(BuenColors)
library(dplyr)

so <- readRDS("../../../brain_hsv1_annotated_sct_obj.rds")
head(so@meta.data)

p0 <- DimPlot(so, group.by = "location", shuffle = TRUE, pt.size = 0.2) +
  theme_void() + ggtitle("") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "purple3")) +
  theme(legend.position = "none") 
cowplot::ggsave2(p0, file = "../plots/what-tissue.png", dpi = 300, width = 3.6, height = 3.6)
  

pXY <- DimPlot(so, group.by = "newgrp", shuffle = TRUE, pt.size = 0.2) +
  theme_void() + ggtitle("") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "purple3")) +
  theme(legend.position = "none") 
cowplot::ggsave2(pXY, file = "../plots/newgrp.png", dpi = 300, width = 3.6, height = 3.6)


p1 <- FeaturePlot(so, features = c("hsv_reads"),  order = TRUE,
           pt.size = 0.1, max.cutoff = "q95") + FontSize(main = 0.0001) + 
  theme_void() + ggtitle("") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)])) +
  theme(legend.position = "none") 
cowplot::ggsave2(p1, file = "../plots/hsv-expression-umap.png", dpi = 300, width = 3.6, height = 3.6)

p3 <- FeaturePlot(so, features = c("RORB"),  order = TRUE,
                  pt.size = 0.1, max.cutoff = "q95") + FontSize(main = 0.0001) + 
  theme_void() + ggtitle("") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none") 
cowplot::ggsave2(p3, file = "../plots/RORB.png", dpi = 300, width = 3.6, height = 3.6)

custom_color_vec <- as.character(jdb_palette("corona"))[1:10]
names(custom_color_vec) <- c("Glutamatergic", "Oligo", "GABAergic", "OPC", "Micro-PVM", 
                             "Astro", "Endo", "L2/3 IT", "Sst","VLMC") 
  
p2 <- DimPlot(so, group.by = "predicted.celltypes", shuffle = FALSE, pt.size = 0.2) +
  theme_void() + ggtitle("") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_color_manual(values = custom_color_vec) +
  theme(legend.position = "none") 

cowplot::ggsave2(p2, file = "plots/celltype-umap.png", dpi = 300, width = 3.6, height = 3.6)


DimPlot(so, group.by = c("newgrp", "predicted.celltypes"))
FeaturePlot(so, features =  c("RORB", ), max.cutoff = "q99", order = TRUE)




##### ------
# Make viral expression heatmap
##### ------
DefaultAssay(so) <- "viral"
obj <- JoinLayers(so)
mat <- obj[["viral"]]$counts 
log_mat <- t(log2(mat[,colSums(mat) > 10] + 1))


# Annotate for a combined heatmap
library(ComplexHeatmap)
library(data.table)

# annotate
hsv1_mat_t <- data.matrix(t(log_mat))
hsv1_mat_t <- hsv1_mat_t[rowSums(hsv1_mat_t) > 0, ]
annotat_df <- fread("hsv1_geneannotation.csv")
vec <- annotat_df$Kinetic_Class; names(vec) <- annotat_df$Gene
anno_df <- data.frame(gene = rownames(hsv1_mat_t), anno = vec[as.character(rownames(hsv1_mat_t))])
adf_go <- anno_df[complete.cases(anno_df),] %>%
  arrange((anno))
dim(adf_go)
# resum and sort
data.frame(
  cell = rownames(log_mat), 
  library = substr(rownames(log_mat), 11, 11), 
  total = rowSums(log_mat)
) %>% arrange(library, desc(total)) %>% pull(cell) -> cell_order

ha <- HeatmapAnnotation(
  gene_annotation = adf_go$anno,
  col = list(
    gene_annotation = c("early" = "dodgerblue3",
                        "late" = "dodgerblue4",
                        "aimmedate_early" = "dodgerblue")
  ),
  gp = gpar(col = "black"),
  show_legend = FALSE, annotation_label = ""
)

hm <- Heatmap(t(hsv1_mat_t[adf_go$gene,cell_order]), cluster_rows = FALSE, cluster_columns = FALSE,
              col = jdb_palette("solar_rojos"), top_annotation =  ha, 
              column_names_gp = grid::gpar(fontsize = 4),
              column_split  = adf_go$anno,
              row_names_gp = grid::gpar(fontsize = 0), show_heatmap_legend = FALSE)
hm

pdf(file="plots/hsv-gscrna-plot.pdf", width = 4, height = 2)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm
dev.off()
