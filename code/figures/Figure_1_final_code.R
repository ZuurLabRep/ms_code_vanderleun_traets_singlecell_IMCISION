##################
#### FIGURE 1 ####
##################

# load functions, create folder
setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_1"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_1")

library(dittoSeq)


## load data
# clean_mat_id: filtered count matrix
# filt_mc_id: info on filtered meta cell generation
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group),]


#### Figure 1B ####
mc_colors <- read.table("./data/mc_colors.txt",header=T,sep="\t")

# plot 2D projection
plot_coord_mc_colors(mc2d_name=mc2d_id,mat_name=clean_mat_id,mc_name=filt_mc_id,mc_colors=mc_colors)



#### Figure 1C ####
## projection of genes
genes_oi = c("CD8B","CD4","CD14","CD3D","CD19","NCAM1")

# prep data, simple scaling
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]
umis = as.matrix(mat@mat[,rownames(submat)])
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 

# plot 2D projection (in high resolution png)
for (gene_oi in genes_oi){
  plot_coord_gene_expression(mc2d_name=mc2d_id,mat_name=clean_mat_id,umis_n=umis_n,mc_name=filt_mc_id,gene_oi=gene_oi)
}



#### Figure 1D ####
## heatmap with the gene expression of the feature genes per cell
data = as.matrix(mat@mat)
met.data = mat@cell_metadata
met.data = met.data[!is.na(met.data$mc_group),]
data <- data[,rownames(met.data)]
logexp <- log2(data + 1)

bulkSCE <- importDittoBulk(x <- list(counts = data,logcounts = logexp), metadata = met.data)

gene_list_ordered <- read.table("./data/Gene_list_heatmap_marker_genes_mc.txt",header=F)
gene_list_ordered <- gene_list_ordered$V1
  
bulkSCE$cell_type_nr <- bulkSCE$cell_type
bulkSCE@colData[bulkSCE@colData$cell_type_nr == "CD4",]["cell_type_nr"] <- "2_CD4"
bulkSCE@colData[bulkSCE@colData$cell_type_nr == "CD8",]["cell_type_nr"] <- "1_CD8"
bulkSCE@colData[bulkSCE@colData$cell_type_nr == "NK",]["cell_type_nr"] <- "3_NK"
bulkSCE@colData[bulkSCE@colData$cell_type_nr == "B",]["cell_type_nr"] <- "4_B"
bulkSCE@colData[bulkSCE@colData$cell_type_nr == "M",]["cell_type_nr"] <- "5_M"
bulkSCE@colData[bulkSCE@colData$mc_group == "Non-classical T cell",]["cell_type_nr"] <- "1_CD8"

colors_cells <- read.table("./data/mc_colors.txt",header=T)
colors_cells <- colors_cells[order(colors_cells$name),]

in_data <- gene_list_ordered[!(gene_list_ordered %in% rownames(bulkSCE))]

pdf(file=paste0(fig_path,"/heatmap_feature_genes_v2.pdf"), width=15, height=10) 
dittoHeatmap(bulkSCE,genes = gene_list_ordered,scaled.to.max=TRUE,
             annot.colors = c(colors_cells$color,c("#7f7f7f","#bfbfbf","#606060","#353535","#e5e5e5")),
             cluster_rows=FALSE, show_colnames = FALSE, cluster_cols = FALSE, annot.by = c("mc_group","cell_type_nr"),
             order.by = c("cell_type_nr","mc_group"),heatmap.colors.max.scaled = colorRampPalette(c("white", "darkorange", "#3A2C09"))(30))
dev.off()

