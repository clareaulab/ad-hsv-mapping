library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
library(ComplexHeatmap)

# import HSV and rna-seq data
rna1 <- Read10X_h5("../data/SRR19792154_filtered_feature_bc_matrix.h5")
rna2 <- Read10X_h5("../data/SRR19792155_filtered_feature_bc_matrix.h5"); colnames(rna2) <- gsub("-1", "-2", colnames(rna2))
rna3 <- Read10X_h5("../data/SRR19792156_filtered_feature_bc_matrix.h5"); colnames(rna3) <- gsub("-1", "-3", colnames(rna3))
rna <- cbind(rna1,rna2,rna3)

import_kite_counts <- function(srr, add){
  mtx <- fread(paste0("../data/",srr,"_HSV1_final_kboutput.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  rn <-  paste0(fread(paste0("../data/",srr,"_HSV1_final_kboutput.barcodes.txt"), header = FALSE)[[1]], "-", add)
  cn <- paste0(fread(paste0("../data/",srr,"_HSV1_final_kboutput.genes.txt"), header = FALSE)[[1]])
  matx <- sparseMatrix(i = c(mtx[[1]],length(rn)),  j = c(mtx[[2]],length(cn)), x = c(mtx[[3]],0))
  rownames(matx) <- rn
  colnames(matx) <- cn
  return(t(matx))
}

hsv_mat <- cbind(import_kite_counts("SRR19792154", "1"), 
                 import_kite_counts("SRR19792155", "2"),
                 import_kite_counts("SRR19792156", "3")
)

# Make plot of viral expression
# annotate
annotat_df <- fread("../data/hsv1_geneannotation.csv", header = TRUE, sep = ",")
vec <- annotat_df$Kinetic_Class; names(vec) <- annotat_df$Gene
hsv_mat_ss <- hsv_mat[,colSums(hsv_mat) >= 100]
anno_df <- data.frame(gene = rownames(hsv_mat_ss), anno = vec[as.character(rownames(hsv_mat_ss))])
adf_go <- anno_df[complete.cases(anno_df),] %>%
  arrange((anno))

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

hm <- Heatmap(t(log1p(hsv_mat_ss[adf_go$gene,])), cluster_rows = FALSE, cluster_columns = FALSE,
              col = jdb_palette("solar_rojos"), top_annotation =  ha, 
              column_names_gp = grid::gpar(fontsize = 4),
              column_split  = adf_go$anno,
              row_names_gp = grid::gpar(fontsize = 0), show_heatmap_legend = FALSE)
hm

pdf(file="hsv-gex-plot.pdf", width = 4.2, height = 1.5)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm
dev.off()


# 

fix_na <- function(vec){
  unname(ifelse(is.na(vec), 0, vec))
}
n_HSV_umis = colSums(hsv_mat)[colnames(rna)] %>% fix_na
n_HSV_genes = colSums(hsv_mat>0)[colnames(rna)] %>% fix_na

# Set up Seurat object
# meta annotations from here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA851960&o=acc_s%3Aa
meta_df <- data.frame(
  id = c(rep("SRR19792154", dim(rna1)[2]), rep("SRR19792155", dim(rna2)[2]), rep("SRR19792156", dim(rna3)[2])),
  region = c(rep("Occipital_cortex", dim(rna1)[2]), rep("Hippocampus", dim(rna2)[2]), rep("Frontal_cortex", dim(rna3)[2])),
  n_HSV_umis, n_HSV_genes, log_HSV = log1p(n_HSV_umis)
); rownames(meta_df) <- colnames(rna)

so <- CreateSeuratObject(counts = rna, meta.data = meta_df)

# Do Seurat things
summary(so@meta.data$n_HSV_umis >= 10)
so_filt <- subset(so, nFeature_RNA >= 500 & nCount_RNA >= 1000)
summary(so_filt@meta.data$n_HSV_umis >= 10)
so_filt <- NormalizeData(so_filt) %>% FindVariableFeatures() %>% ScaleData() 
so_filt <- so_filt %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30) 
FeaturePlot(so_filt, c("RORB", "log_HSV"))

FeaturePlot(so_filt, c("RORB", "log_HSV", "MAP2", "NSE", "TUBB3", "DCX", "RBFOX3"))
so_filt$hsv_high <- case_when(
  so_filt$n_HSV_umis > 50 ~"high", 
  so_filt$n_HSV_umis > 5 ~"mid", 
  TRUE ~"low", 
)

FeaturePlot(so_filt, c("RORB", "log_HSV", "RHOJ", "FOS", "JUN", "CUX2", "NEUROD2", "SATB2", "MET", "CYP26A1"),
            order = TRUE, min.cutoff = 0, max.cutoff = "q98")

library(annotables)
x_genes <- grch38 %>% filter(chr == "X") %>% pull(symbol) %>% unique()
so_filt[["percent.X"]] <- PercentageFeatureSet(so_filt, features = intersect(rownames(so_filt), x_genes))
so_filt@meta.data %>%
  group_by(seurat_clusters) %>%  summarize(mean(percent.X), mean(n_HSV_genes > 10)) %>% data.frame()

# Y chromosome genes
msy_genes <- c("ZFY", "TBL1Y", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "KDM5D", "EIF1AY") # no RPS4Y1 in flex

FeaturePlot(so_filt, msy_genes,
            order = TRUE, min.cutoff = 0, max.cutoff = "q98")


# higher res clustering
so_filt <- FindClusters(so_filt, resolution = 3)
DimPlot(so_filt, group.by = c("seurat_clusters"), label = TRUE)

so_filt@meta.data %>%
  group_by(seurat_clusters) %>% summarize(round(mean(n_HSV_umis),1), prop = sum(n_HSV_umis > 10)/n()) %>% data.frame()
FindMarkers(so_filt, ident.1 = "28",  group.by = "seurat_clusters")
FeaturePlot(so_filt, c( "log_HSV",  "MECOM", "HOXA9", "KRT18", "ANXA2","KRT8"))


so_filt[["percent.mt"]] <- PercentageFeatureSet(so_filt, pattern = "^MT-")

FeaturePlot(so_filt, c("percent.mt", "log_HSV", "XIST", "MEF2C", "CADM2", "STMN1", "TMSB10"), max.cutoff = "q99", order = TRUE)

FeaturePlot(so_filt, c( "log_HSV", "XIST", "CLDN5"), max.cutoff = "q99", order = TRUE, split.by = "region")
