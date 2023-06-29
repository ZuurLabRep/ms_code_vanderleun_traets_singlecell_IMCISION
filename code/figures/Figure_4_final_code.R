##################
#### FIGURE 4 ####
##################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_4"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_4")

library(monocle3)
library(VISION)
library(ComplexHeatmap)
library(circlize)


#### Figure 4A ####
# monocle analysis Tregs
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)

# limit dimension reduction to list of genes
gene_set <- read.table("./data/HN_full_dataset_clean_feats_filtered.txt",header=T)
gene_set <- gene_set$x

md <- mat@cell_metadata[names(mc@mc),]
md <- md[md$cell_type == "CD4",]
md <- md[! (md$mc_group == "Granulocyte") ,] 

md <- md[which(rownames(md) %in%  colnames(mat@mat)),]
mat_data <- as.matrix(mat@mat)
mat_in_dat <- mat_data[,which(colnames(mat_data) %in%  rownames(md))]
md$orig.ident.v2 <- substr(md$orig.ident,1,5)

# input monocle3
cds <- new_cell_data_set(as(mat_in_dat, "sparseMatrix"),cell_metadata=md)

# norm and preprocess
cds <- preprocess_cds(cds, num_dim = 25,method="PCA", norm_method = c("log"),use_genes=gene_set) 

# remove possible batch effects
cds <- align_cds(cds, alignment_group = "orig.ident.v2",preprocess_method = "PCA")

# reduce dim
cds <- reduce_dimension(cds, reduction_method="UMAP",umap.metric="euclidean",umap.min_dist=0.1,umap.n_neighbors=50)

# recluster
cds <- cluster_cells(cds, k = 12)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

cds <- learn_graph(cds)
#order_cells(cds) # choose naive CD4 population

pdf(file=paste0(fig_path,"/CD4_UMAP_lineage_mc_group_genelist.pdf"), width=6.2, height=5)
plot_cells(cds,
           color_cells_by = "mc_group",
           label_cell_groups=FALSE,
           label_leaves=FALSE,trajectory_graph_color = "black",
           label_branch_points=FALSE,trajectory_graph_segment_size = 1.2,
           graph_label_size=1.5,label_roots = FALSE,
           cell_size = 0.7)+ scale_color_manual(values = c("#b3d4f0", "#a14f00","#107538","#0042b0","#59b56f"))
dev.off()

#saveRDS(cds,"./data/cds_TregGIMAP_TregTNF_all.rds")


#### Figure 4B ####
# load subset Tregs, from choose_cells
cds_subset = choose_cells(cds) # choose Treg cluster
#saveRDS(cds_subset,"./data/cds_subset_TregGIMAP_TregTNF_v2.rds")
#cds_subset <- readRDS(file = "./data/cds_subset_TregGIMAP_TregTNF_v2.rds")

## heatmap CD4 trajectory genes
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)
pr_deg_ids <- subset(pr_graph_test_res, q_value < 0.05)

pr_deg_ids <- subset(pr_graph_test_res, morans_I > 0.07)
pr_deg_ids <- pr_deg_ids[order(pr_deg_ids$q_value),]
#write.csv(pr_deg_ids,"./data/CD4_DEG_alongTrajectory_TregGIMAP_TregTNF_v2.csv",quote=F,row.names = F)

treg_genes <- pr_deg_ids[1:100,]$gene_short_name
treg_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% treg_genes]

test_top100 <- plot_genes_in_pseudotime(treg_lineage_cds,
                                        color_cells_by="mc_group",
                                        min_expr=0.4)
test_top100_data <- test_top100$data
test_top100_data <- test_top100_data[,c("Cell","patient","mc_group","response","condition","pseudotime","feature_label","adjusted_expression")]
test_top100_data_genes <- reshape2::dcast(data = test_top100_data,formula = Cell~feature_label,fun.aggregate = sum,value.var = "adjusted_expression")
test_top100_data_sel <- unique(test_top100_data[,-c(7,8)])
test_top100_data_all <- merge(test_top100_data_genes,test_top100_data_sel)

#write.table(test_top100_data,"./data/CD4_pseudotime_data_combine_top100.txt",sep="\t",quote=F,row.names = F)

# bin data 
start_gene <- colnames(test_top100_data_all)[2]
binned_mean <- tapply(test_top100_data_all[[start_gene]], cut(test_top100_data_all$pseudotime, 30), mean,na.rm = FALSE)
binned_mean_all <- as.data.frame(binned_mean)
colnames(binned_mean_all)[1] <- start_gene
for (gene in colnames(test_top100_data_all)[3:101]){
  binned_mean <- tapply(test_top100_data_all[[gene]], cut(test_top100_data_all$pseudotime, 30), mean)
  binned_mean <- as.data.frame(binned_mean)
  colnames(binned_mean)[1] <- gene
  binned_mean_all <- cbind(binned_mean_all,binned_mean)
}

# barplot bins
binned_mean <- test_top100_data_all
binned_mean$bins <- cut(test_top100_data_all$pseudotime, 30)
start_bin = levels(as.factor(binned_mean$bins))[1]

# barplot bins RE
binned_mean <- test_top100_data_all
binned_mean$bins <- cut(test_top100_data_all$pseudotime, 30)
binned_mean[binned_mean$patient=="Pat27",][["response"]] <- "RE"
binned_mean[binned_mean$patient=="Pat36",][["response"]] <- "RE"
binned_mean = binned_mean[binned_mean$response == "RE",]
binned_mean$rep_cond <- paste0(binned_mean$response,"_",binned_mean$condition)

binned_mean$test_pat <- paste0(binned_mean$response,"_",binned_mean$patient)

example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$rep_cond)),rep(0,length(levels(as.factor(binned_mean$rep_cond))))))
freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == start_bin,]$rep_cond))
freq_bin_df <- plyr::join(example_df,freq_bin)
#freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0
freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$rep_cond))))
freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
freq_bin_df$bin <- rep(start_bin,length(levels(as.factor(binned_mean$rep_cond))))
freq_bin_df_save <- freq_bin_df

for (bin in levels(as.factor(binned_mean$bins))[2:30]){
  print(bin)
  example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$rep_cond)),rep(0,length(levels(as.factor(binned_mean$rep_cond))))))
  freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == bin,]$rep_cond))
  freq_bin_df <- plyr::join(example_df,freq_bin)
  tryCatch({
    freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$rep_cond))))
  freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
  freq_bin_df$bin <- rep(bin,length(levels(as.factor(binned_mean$rep_cond))))
  freq_bin_df_save <- rbind(freq_bin_df_save,freq_bin_df)
}

freq_bin_df_save_RE <- freq_bin_df_save

# barplot bins NR
binned_mean <- test_top100_data_all
binned_mean$bins <- cut(test_top100_data_all$pseudotime, 30)
binned_mean[binned_mean$patient=="Pat27",][["response"]] <- "RE"
binned_mean[binned_mean$patient=="Pat36",][["response"]] <- "RE"
binned_mean = binned_mean[binned_mean$response == "NR",]
binned_mean$rep_cond <- paste0(binned_mean$response,"_",binned_mean$condition)
binned_mean$test_pat <- paste0(binned_mean$response,"_",binned_mean$patient)

example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$rep_cond)),rep(0,length(levels(as.factor(binned_mean$rep_cond))))))
freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == start_bin,]$rep_cond))
freq_bin_df <- plyr::join(example_df,freq_bin)
#freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0
freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$rep_cond))))
freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
freq_bin_df$bin <- rep(start_bin,length(levels(as.factor(binned_mean$rep_cond))))
freq_bin_df_save <- freq_bin_df

for (bin in levels(as.factor(binned_mean$bins))[2:30]){
  print(bin)
  example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$rep_cond)),rep(0,length(levels(as.factor(binned_mean$rep_cond))))))
  freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == bin,]$rep_cond))
  freq_bin_df <- plyr::join(example_df,freq_bin)
  tryCatch({
    freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$rep_cond))))
  freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
  freq_bin_df$bin <- rep(bin,length(levels(as.factor(binned_mean$rep_cond))))
  freq_bin_df_save <- rbind(freq_bin_df_save,freq_bin_df)
}

freq_bin_df_save_NR <- freq_bin_df_save

# barplot bins mc groups
binned_mean <- test_top100_data_all
binned_mean$bins <- cut(test_top100_data_all$pseudotime, 30)

example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$mc_group)),rep(0,length(levels(as.factor(binned_mean$mc_group))))))
freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == start_bin,]$mc_group))
freq_bin_df <- plyr::join(example_df,freq_bin)
freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0
freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$mc_group))))
freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
freq_bin_df$bin <- rep(start_bin,length(levels(as.factor(binned_mean$mc_group))))
freq_bin_df_save <- freq_bin_df

for (bin in levels(as.factor(binned_mean$bins))[2:30]){
  print(bin)
  example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$mc_group)),rep(0,length(levels(as.factor(binned_mean$mc_group))))))
  freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == bin,]$mc_group))
  freq_bin_df <- plyr::join(example_df,freq_bin)
  tryCatch({
    freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  freq_bin_df$Total <- rep(sum(freq_bin$Freq),length(levels(as.factor(binned_mean$mc_group))))
  freq_bin_df$Fraction <- freq_bin_df$Freq/freq_bin_df$Total
  freq_bin_df$bin <- rep(bin,length(levels(as.factor(binned_mean$mc_group))))
  freq_bin_df_save <- rbind(freq_bin_df_save,freq_bin_df)
}

recast_freq_bin_df_save <- reshape2::dcast(data = freq_bin_df_save,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save[is.na(recast_freq_bin_df_save)] <- 0

# heatmap with ggplot
oi_genes <- pr_deg_ids[1:50,]$gene_short_name
temp_data <- binned_mean_all[,oi_genes]

# transform data into the correct dims for heatmap
recast_freq_bin_df_save <- reshape2::dcast(data = freq_bin_df_save,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save[is.na(recast_freq_bin_df_save)] <- 0

recast_freq_bin_df_save_RE <- reshape2::dcast(data = freq_bin_df_save_RE,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save_RE[is.na(recast_freq_bin_df_save_RE)] <- 0
rownames(recast_freq_bin_df_save_RE) <- recast_freq_bin_df_save_RE$Var1
recast_freq_bin_df_save_NR <- reshape2::dcast(data = freq_bin_df_save_NR,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save_NR[is.na(recast_freq_bin_df_save_NR)] <- 0
rownames(recast_freq_bin_df_save_NR) <- recast_freq_bin_df_save_NR$Var1
recast_freq_bin_df_save_RE <- rev(as.data.frame(recast_freq_bin_df_save_RE))
recast_freq_bin_df_save_NR <- rev(as.data.frame(recast_freq_bin_df_save_NR))
recast_freq_bin_df_save <- rev(as.data.frame(recast_freq_bin_df_save))

#write.table(recast_freq_bin_df_save,"./data/CD4_pseudotime_barplot_mcgroups.txt",sep="\t",quote=F,row.names =F)

ha = HeatmapAnnotation(RE = anno_barplot(t(recast_freq_bin_df_save_RE[,1:30]), gp = gpar(fill = c("#07522B","#71BB87")), 
                                         bar_width = 1, height = unit(1.3, "cm")),
                       NR = anno_barplot(t(recast_freq_bin_df_save_NR[,1:30]), gp = gpar(fill = c("#901A1D","#E18283","#07522B","#71BB87")), 
                                         bar_width = 1, height = unit(1.3, "cm")),
                       foo = anno_barplot(t(recast_freq_bin_df_save[,1:30]), gp = gpar(fill =c("#B3D5EF", "#A05225","#10753D","#1F4A9F","#56B66F")), 
                                          bar_width = 1, height = unit(1, "cm")))

ha2 = rowAnnotation(foo = anno_mark(at =c(1,2,3,5,6,8,9,10,18,37,45), labels = rownames(t(temp_data))[c(1,2,3,5,6,8,9,10,18,37,45)]))

temp_data_t <- t(temp_data)
temp_data_t <- rev(as.data.frame(temp_data_t))
x <- t(scale(t(temp_data_t)))

col_fun <- colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
pdf(file=paste0(fig_path,"/CD4_pseudotime_heatmap_barplot_combine_top50_both_v3_NA.pdf"),width=6.7, height=5.4)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,right_annotation = ha2,
        row_names_side = "left",row_names_gp = gpar(fontsize = 4),show_row_names = FALSE)
dev.off()

# with zeros
binned_mean_all[is.na(binned_mean_all)] <- 0
temp_data <- binned_mean_all[,oi_genes]
temp_data_t <- t(temp_data)
temp_data_t <- rev(as.data.frame(temp_data_t))
x <- t(scale(t(temp_data_t)))

col_fun <- colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
pdf(file=paste0(fig_path,"/CD4_pseudotime_heatmap_barplot_combine_top50_both_v3.pdf"),width=6.7, height=5.4)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,right_annotation = ha2,
        row_names_side = "left",row_names_gp = gpar(fontsize = 4),show_row_names = FALSE)
dev.off()


#### Figure 4C ####
## top 5 signatures ssGSEA
annotation_ptime <- read.table("./data/CD4_pseudotime_data_combine_top100.txt",sep="\t",header=T)

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

submat <- mat@cell_metadata       
md <- submat[rownames(submat) %in% unique(annotation_ptime$Cell),] 

# ssGSEA
count_mat <- as.matrix(mat@mat[,rownames(md)])
n.umi <- colSums(count_mat)
scaled_counts <- t(t(count_mat) / n.umi) * median(n.umi)
norm_counts <- sweep(count_mat,2,colSums(count_mat), "/") * 1000
md$Cell <- rownames(md)

# top 5, ssGSEA with visium, selection of top 5 GO terms
vis <- Vision(scaled_counts,
              signatures = c("./data/top5_plus_trajectory.txt"),
              meta = md)
options(mc.cores = 2)
vis <- analyze(vis)

sig_scores_all <- getSignatureScores(vis)
sig_scores_all <- as.data.frame(sig_scores_all)
sig_scores_all$cell_id <- rownames(sig_scores_all)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(sig_scores_all,md,by="cell_id")

sig_scores_all_pre <- sig_scores_all[sig_scores_all$condition == "pre",]
sig_scores_all_on <- sig_scores_all[sig_scores_all$condition == "post",]
sig_scores_all_pre <- sig_scores_all

collect <- annotation_ptime[,c("Cell","mc_group","condition","response","pseudotime")]
collect <- collect[!duplicated(collect),]
collect_ptime <-merge(sig_scores_all_pre,collect,by="Cell")

# bin gene sets on pseudotime
binned_mean <- tapply(collect_ptime$GOMF_CYTOKINE_RECEPTOR_ACTIVITY, cut(collect_ptime$pseudotime, 30), mean,na.rm = FALSE)
binned_mean_all <- as.data.frame(binned_mean)
colnames(binned_mean_all)[1] <- "GOMF_CYTOKINE_RECEPTOR_ACTIVITY"
for (gene in colnames(collect_ptime)[3:6]){
  binned_mean <- tapply(collect_ptime[[gene]], cut(collect_ptime$pseudotime, 30), mean)
  binned_mean <- as.data.frame(binned_mean)
  colnames(binned_mean)[1] <- gene
  binned_mean_all <- cbind(binned_mean_all,binned_mean)
}

# saved from figure 4B
recast_freq_bin_df_save = read.table("./data/CD4_pseudotime_barplot_mcgroups.txt",sep="\t",header=F)
colnames(recast_freq_bin_df_save) = recast_freq_bin_df_save[1,]
recast_freq_bin_df_save = recast_freq_bin_df_save[-1,]
recast_freq_bin_df_save <- recast_freq_bin_df_save %>% 
  mutate_if(is.character, as.numeric)

ha = HeatmapAnnotation(foo = anno_barplot(t(recast_freq_bin_df_save[,1:30]), gp = gpar(fill =c("#B3D5EF", "#A05225","#10753D","#1F4A9F","#56B66F")), 
                                          bar_width = 1, height = unit(1.3, "cm")))

temp_data_t <- t(binned_mean_all[rownames(binned_mean_all),])
temp_data_t <- rev(as.data.frame(temp_data_t))
x <- t(scale(t(temp_data_t)))

col_fun <- colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
pdf(file=paste0(fig_path,"/CD4_pseudotime_GOterm_enrichement_top5.pdf"),width=6.7, height=4.5)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,
        row_names_side = "right",row_names_gp = gpar(fontsize = 4),show_row_names = TRUE)
dev.off()



#### Figure 4D ####
cds_subset = choose_cells(cds) # choose Treg cluster (square)
#saveRDS(cds_subset,"./data/cds_subset_TregGIMAP_TregTNF_all.rds")
#cds_subset <- readRDS(file = "./data/cds_subset_TregGIMAP_TregTNF_all.rds")

pl_set = plot_cells(cds_subset)
umap_1 = pl_set$data$data_dim_1
umap_2 = pl_set$data$data_dim_2
umap_mc = pl_set$data$mc_group
umap_cell = rownames(pl_set$data)

## UMAP genes
mc = scdb_mc(filt_mc_id)
mat = cds_subset@assays@data@listData$counts
submat = cds_subset@colData
umis = as.matrix(mat[,rownames(submat)])
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 

for (gene_oi in c("TNFRSF4","IL2RA","TNFRSF1B")){
  plot_UMAP_expression(cds_subset,umis_n,gene_oi,"treg")
}



