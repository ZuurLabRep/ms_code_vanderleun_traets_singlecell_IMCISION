##################
#### FIGURE 5 ####
##################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_5"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_5")

library(monocle3)
library(VISION)
library(ComplexHeatmap)
library(circlize)
library(AUCell)


#### Figures 5A ####
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

# limit dimension reduction to list of genes
gene_set = read.table("./data/Gene_list_heatmap_marker_genes_mc.txt")
gene_set <- gene_set$V1
gene_set <- rownames(gene_set)

md <- mat@cell_metadata[names(mc@mc),]
md = md[md$cell_type == "CD8",]
md = md[md$mc_group %in% c("Dysf GZMB", "Dysf ZNF683", "Naive-like CD8" ,"Transitional", "Cytotoxic"),] 

md<- md[which(rownames(md) %in%  colnames(mat@mat)),]
mat_data <- as.matrix(mat@mat)
mat_in_dat = mat_data[,which(colnames(mat_data) %in%  rownames(md))]
md$orig.ident <- substr(md$orig.ident,1,5)

# input monocle3
cds <- new_cell_data_set(as(mat_in_dat, "sparseMatrix"),cell_metadata=md)

# norm and preprocess
cds <- preprocess_cds(cds, num_dim = 40,method="PCA", norm_method = c("log"),use_genes=gene_set)
plot_pc_variance_explained(cds)

# remove possible batch effects
cds <- align_cds(cds, alignment_group = "orig.ident",preprocess_method = "PCA")

# reduce dim
cds <- reduce_dimension(cds, reduction_method="UMAP",umap.metric="euclidean",umap.n_neighbors=30) 

# recluster
cds = cluster_cells(cds, k = 19)

cds <- learn_graph(cds)
order_cells(cds) # naive as starting point

pdf(file=paste0(fig_path,"/UMAP_lineage_mc_group_norm_genelist_CD8.pdf"), width=6.2, height=5)
plot_cells(cds,
           color_cells_by = "mc_group",
           label_cell_groups=FALSE,
           label_leaves=FALSE,trajectory_graph_color = "black",
           label_branch_points=FALSE,trajectory_graph_segment_size = 1.2,
           graph_label_size=1.5,label_roots = FALSE,
           cell_size = 0.7)+ scale_color_manual(values = c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8"))
dev.off()

#saveRDS(cds_subset,"./data/cds_Dysf_all.rds")



#### Figure 5B ####
## trajectory heatmap
cds_subset = choose_cells(cds) # choose dysf cluster
#saveRDS(cds_subset, "./data/cds_subset_TNF_ZNF.rds")
#cds_subset <- readRDS(file = "./data/cds_subset_TNF_ZNF.rds")

## heatmap CD8 trajectory genes
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)
pr_deg_ids <- subset(pr_graph_test_res, morans_I > 0.07)
pr_deg_ids <- pr_deg_ids[order(pr_deg_ids$q_value),]
#write.csv(pr_deg_ids,"./data/CD8_DEG_alongTrajectory_dysTNF_dysZNF.csv",quote=F,row.names = F)

dysf_genes <- pr_deg_ids[1:100,]$gene_short_name
dysf_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% dysf_genes]

test_top100 <- plot_genes_in_pseudotime(dysf_lineage_cds,
                                        color_cells_by="mc_group",
                                        min_expr=0.4)

test_top100_data <- test_top100$data
test_top100_data <- test_top100_data[,c("Cell","mc_group","response","condition","pseudotime","feature_label","adjusted_expression","patient")]
test_top100_data_genes <- reshape2::dcast(data = test_top100_data,formula = Cell~feature_label,fun.aggregate = sum,value.var = "adjusted_expression")
test_top100_data_sel <- unique(test_top100_data[,-c(6,7)])
test_top100_data_all <- merge(test_top100_data_genes,test_top100_data_sel)

#write.table(test_top100_data,"./data/CD8_pseudotime_data_combine_top100.txt",sep="\t",quote=F,row.names = F)

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
binned_mean[binned_mean$patient == "Pat27",]["response"] <- "RE"
binned_mean[binned_mean$patient == "Pat36",]["response"] <- "RE"
binned_mean = binned_mean[binned_mean$response == "RE",]
binned_mean$rep_cond <- paste0(binned_mean$response,"_",binned_mean$condition)

example_df <- as.data.frame(cbind(Var1=levels(as.factor(binned_mean$rep_cond)),rep(0,length(levels(as.factor(binned_mean$rep_cond))))))
freq_bin <- as.data.frame(table(binned_mean[binned_mean$bins == start_bin,]$rep_cond))
freq_bin_df <- plyr::join(example_df,freq_bin)
freq_bin_df[is.na(freq_bin_df$Freq),]$Freq <- 0
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
binned_mean[binned_mean$patient == "Pat27",]["response"] <- "RE"
binned_mean[binned_mean$patient == "Pat36",]["response"] <- "RE"
binned_mean = binned_mean[binned_mean$response == "NR",]
binned_mean$rep_cond <- paste0(binned_mean$response,"_",binned_mean$condition)

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

#write.table(recast_freq_bin_df_save,"./data/CD8_pseudotime_barplot_mcgroups.txt",sep="\t",quote=F,row.names =F)

# heatmap with ggplot
oi_genes <- pr_deg_ids[1:50,]$gene_short_name # only plot 50
temp_data <- binned_mean_all[,oi_genes]

recast_freq_bin_df_save_RE <- reshape2::dcast(data = freq_bin_df_save_RE,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save_RE[is.na(recast_freq_bin_df_save_RE)] <- 0
rownames(recast_freq_bin_df_save_RE) <- recast_freq_bin_df_save_RE$Var1
recast_freq_bin_df_save_NR <- reshape2::dcast(data = freq_bin_df_save_NR,formula = Var1~bin,fun.aggregate = sum,value.var = "Fraction")
recast_freq_bin_df_save_NR[is.na(recast_freq_bin_df_save_NR)] <- 0
rownames(recast_freq_bin_df_save_NR) <- recast_freq_bin_df_save_NR$Var1

library(ComplexHeatmap)
ha = HeatmapAnnotation(RE = anno_barplot(t(recast_freq_bin_df_save_RE[,2:31]), gp = gpar(fill = c("#07522B","#71BB87")), 
                                         bar_width = 1, height = unit(1.3, "cm")),
                       NR = anno_barplot(t(recast_freq_bin_df_save_NR[,2:31]), gp = gpar(fill = c("#901A1D","#E18283","#07522B","#71BB87")), 
                                         bar_width = 1, height = unit(1.3, "cm")),
                       mc = anno_barplot(t(recast_freq_bin_df_save[,2:31]), gp = gpar(fill = c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8")), 
                                          bar_width = 1, height = unit(1, "cm")))

ha2 = rowAnnotation(foo = anno_mark(at =c(1,2,4,6,7,21,32,34,45,42,38), labels = rownames(t(temp_data))[c(1,2,4,6,7,21,32,34,45,42,38)]))

temp_data_t <- t(temp_data)
x <- t(scale(t(temp_data_t)))

col_fun = colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
col_fun(seq(-3, 3))
pdf(file=paste0(fig_path,"/CD8_pseudotime_heatmap_barplot_combine_top50_both_fix_v3_NA.pdf"),width=6.7, height=5.4)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,right_annotation = ha2,
        row_names_side = "left",row_names_gp = gpar(fontsize = 4),show_row_names = FALSE)
dev.off()

binned_mean_all[is.na(binned_mean_all)] <- 0
temp_data <- binned_mean_all[,oi_genes]
temp_data_t <- t(temp_data)
x <- t(scale(t(temp_data_t)))
col_fun = colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
col_fun(seq(-3, 3))
pdf(file=paste0(fig_path,"/CD8_pseudotime_heatmap_barplot_combine_top50_both_fix_v3.pdf"),width=6.7, height=5.4)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,right_annotation = ha2,
        row_names_side = "left",row_names_gp = gpar(fontsize = 4),show_row_names = FALSE)
dev.off()


#### Figure 5C ####
#cds_subset = choose_cells(cds) # choose dysf cluster (square)
#saveRDS(cds_subset, "./data/cds_subset_Dysf_all.rds")
#cds_subset <- readRDS(file = "./data/cds_subset_Dysf_all.rds")

cells_df = plot_cells(cds_subset)
umap_1 = cells_df$data$data_dim_1
umap_2 = cells_df$data$data_dim_2
umap_mc = cells_df$data$mc_group
umap_cell = rownames(cells_df$data)

## UMAP genes
mc = scdb_mc(filt_mc_id)
mat = cds_subset@assays@data@listData$counts
submat = cds_subset@colData
umis = as.matrix(mat[,rownames(submat)])
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 

for (gene_oi in c("LAG3","HAVCR2","CTLA4")){
  plot_UMAP_expression(cds_subset,umis_n,gene_oi,cell_type="dysf")
}



#### Figure 5D ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

tcr_sig = read.table("./data/NIHMS1790595-supplement-TABLE_4_signature_NeoTCR.csv",sep=",",header=T)

mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]

geneSets <- list(neoTCRCD4=tcr_sig$NeoTCR4,neoTCRCD8=tcr_sig$NeoTCR8)
genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Dysf GZMB","Dysf ZNF683","Naive-like CD8","Transitional","Cytotoxic") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("pre","post"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data = cells_AUC@assays@data@listData$AUC
check_data = as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_all$mc_group = factor(sig_scores_all$mc_group ,levels=c("Dysf GZMB","Dysf ZNF683","Naive-like CD8","Transitional","Cytotoxic"))
pdf(file=paste0(fig_path,"/neoTCR_CD8_all_sig.pdf"),width=4, height=4)
ggplot(sig_scores_all) +  geom_jitter(aes(x=mc_group,y=neoTCRCD8),width=0.2,size=0.6,alpha=0.7) + 
  geom_violin(aes(x=mc_group,y=neoTCRCD8,fill=mc_group)) +
  geom_boxplot(aes(x=mc_group,y=neoTCRCD8,fill=mc_group),width=0.2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 40, hjust = 1),legend.position="none",
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("#4C005C","#C20088","#b3d4f0","#c8a2c8","#E5D3B3")) +
  ylab("neoTCR_CD8") + 
  ggtitle("neoTCR_CD8 signature all") +
  ylim(c(0.03,0.35))
dev.off()

wilcox.test(sig_scores_all[sig_scores_all$mc_group == "Naive-like CD8",]$neoTCRCD8,
            sig_scores_all[sig_scores_all$mc_group == "Transitional",]$neoTCRCD8)



#### Figure 5E ####
# subset dysfunctional cells
#cds_subset = choose_cells(cds) # choose dysf cluster
#saveRDS(cds_subset, "./data/cds_subset_TNF_ZNF.rds")
#ds_subset <- readRDS(file = "./data/cds_subset_TNF_ZNF.rds")

# init
#pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)

## generate small subset for pseudotime extraction
#dysf_genes <- c("CTLA4")
#dysf_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% dysf_genes]

# get pseudotime of each cell
#get_graph <- plot_genes_in_pseudotime(dysf_lineage_cds,
#                                      color_cells_by="condition_response",
#                                      min_expr=0.5,horizontal_jitter=0.01,vertical_jitter=0.01) 
#graph_data <- get_graph$data
#write.table(graph_data,"./data/Pseuodotime_dysf.txt",sep="\t",quote=F,row.names = F)
graph_data = read.table("./data/Pseuodotime_dysf.txt",sep="\t",header=T)

## replace clone column with newest one
mat = scdb_mat(clean_mat_id)
check_meta <- mat@cell_metadata
check_meta$Cell <- rownames(check_meta)
graph_data <- merge(graph_data,check_meta[,c("Cell","unique_clone_ID")],by="Cell")
graph_data$pat_clone_newest <- paste0(graph_data$patient,"_",graph_data$unique_clone_ID.y)

# top clones from the data, RE
clone_list <- c("Pat15_10156" ,"Pat15_10393", "Pat17_446" ,  "Pat21_10195", "Pat22_11952" ,"Pat22_3329" , "Pat22_9012" , "Pat22_9775" ,
                "Pat27_10374" ,"Pat27_10984" ,"Pat27_4204",  "Pat27_527" ,  "Pat27_5411" , "Pat27_7644" , "Pat27_9566") 

plot_list <- list()
for (name_clone in clone_list){
  graph_data$clone_check <- rep("0",nrow(graph_data))
  graph_data[graph_data$pat_clone_newest == name_clone,]["clone_check"] <- "1"
  graph_data$clone_check <- as.numeric(graph_data$clone_check)
  
  set.seed(1)
  graph_data[graph_data$clone_check =="1",]["clone_check"] <- runif(nrow(graph_data[graph_data$clone_check=="1",]), min=0.5, max=0.5)
  
  plot_list[[name_clone]] = ggplot(graph_data, aes(x=pseudotime, y=clone_check)) + geom_point(aes(color=condition_response),size=2) + 
    scale_color_manual(values=c("#075227", "#075227","#73ba86","#73ba86")) + ylim(c(0.3,0.7)) + 
    theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",
                       axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    ggtitle(name_clone) + ylab("") + xlab("") 
}

graph_data$random_nr <- runif(nrow(graph_data), min=0.1, max=0.9)

plot_b = ggplot(graph_data, aes(x=pseudotime, y=random_nr)) + geom_point(aes(color=mc_group),size=1) + 
  scale_color_manual(values=c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8")) + ylim(c(0.1,0.9)) + 
  theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cell_states") + ylab("") + xlab("") 

graph_data_sel = graph_data[graph_data$pat_clone_newest %in% clone_list[1:8],]

plot_d = ggplot(graph_data_sel, aes(x = pseudotime, colour = condition, fill = condition)) + 
  geom_density(lwd = 0.7, alpha = 0.25) + scale_color_manual(values=c("#73ba86","#0b522a")) +
  scale_fill_manual(values=c("#73ba86","#0b522a")) + ylim(0,0.3) + xlim(c(min(graph_data$pseudotime),max(graph_data$pseudotime))) +
  theme_bw() + theme(plot.margin=unit(c(+0.4,1,-0.1,1),"cm"),legend.position = "none", panel.grid.major = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),axis.ticks.y = element_blank(),
                     panel.grid.minor = element_blank())

library(gridExtra)
pdf(file=paste0(fig_path,"/CD8_pseudotime_top_clones_Responders_update_v4.pdf"), width=3.4, height=4.1)
grid.arrange(plot_d,plot_list$Pat15_10156,plot_list$Pat15_10393,plot_list$Pat17_446,plot_list$Pat21_10195,
             plot_list$Pat22_11952,plot_list$Pat22_3329,
             plot_list$Pat22_9012,plot_list$Pat22_9775,plot_list$Pat27_10374,
             plot_list$Pat27_10984,plot_list$Pat27_4204,
             plot_list$Pat27_527,plot_list$Pat27_5411,
             plot_list$Pat27_7644,plot_list$Pat27_9566,
             plot_b,nrow=17,heights=c(1.7,rep(1,length(clone_list)), 3))
dev.off()

# top clones from the data, NR
clone_list <- c("Pat10_10018", "Pat10_10021", "Pat10_3556" , "Pat10_6592"  ,"Pat10_8312",  "Pat26_10269" ,"Pat26_11930", "Pat26_6593" ,
                "Pat26_7425" , "Pat26_9764" , "Pat30_5234" , "Pat33_8273" , "Pat37_11989", "Pat37_8630" , "Pat37_8821" )

graph_data <- get_graph$data

# check with newest table
check_meta <- mat@cell_metadata
check_meta$Cell <- rownames(check_meta)
graph_data <- merge(graph_data,check_meta[,c("Cell","unique_clone_ID")],by="Cell")
graph_data$pat_clone_newest <- paste0(graph_data$patient,"_",graph_data$unique_clone_ID.y)

plot_list <- list()
for (name_clone in clone_list){
  graph_data$clone_check <- rep("0",nrow(graph_data))
  graph_data[graph_data$pat_clone_newest == name_clone,]["clone_check"] <- "1"
  graph_data$clone_check <- as.numeric(graph_data$clone_check)
  
  set.seed(1)
  graph_data[graph_data$clone_check =="1",]["clone_check"] <- runif(nrow(graph_data[graph_data$clone_check=="1",]), min=0.5, max=0.5)
  
  plot_list[[name_clone]] = ggplot(graph_data, aes(x=pseudotime, y=clone_check)) + geom_point(aes(color=condition_response),size=2) + 
    scale_color_manual(values=c("#e08282", "#e08282","#910606","#910606")) + ylim(c(0.3,0.7)) + 
    theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",axis.title.y=element_blank(),
                       axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(size = 9),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    ggtitle(name_clone) + ylab("") + xlab("") 
}

graph_data$random_nr <- runif(nrow(graph_data), min=0.1, max=0.9)

plot_b = ggplot(graph_data, aes(x=pseudotime, y=random_nr)) + geom_point(aes(color=mc_group),size=1) + 
  scale_color_manual(values=c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8","red")) + ylim(c(0.1,0.9)) + 
  theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cell_states") + ylab("") + xlab("") 

graph_data_sel = graph_data[graph_data$pat_clone_newest %in% clone_list,]

plot_d = ggplot(graph_data_sel, aes(x = pseudotime, colour = condition, fill = condition)) + 
  geom_density(lwd = 0.7, alpha = 0.25) + scale_color_manual(values=c("#e08282","#910606")) +
  scale_fill_manual(values=c("#e08282","#910606")) + ylim(0,0.5) + xlim(c(min(graph_data$pseudotime),max(graph_data$pseudotime))) +
  theme_bw() + theme(plot.margin=unit(c(+0.4,1,-0.1,1),"cm"),legend.position = "none", panel.grid.major = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),axis.ticks.y = element_blank(),
                     panel.grid.minor = element_blank())

library(gridExtra)
pdf(file=paste0(fig_path,"/CD8_pseudotime_top_clones_NonResponders_newest_v3.pdf"),width=3.9, height=4.1)
grid.arrange(plot_d,plot_list$Pat10_10018,plot_list$Pat10_10021,plot_list$Pat10_3556,plot_list$Pat10_6592,
             plot_list$Pat10_8312,plot_list$Pat26_10269,
             plot_list$Pat26_11930,plot_list$Pat26_6593,plot_list$Pat26_7425,
             plot_list$Pat26_9764,plot_list$Pat30_5234,
             plot_list$Pat33_8273,plot_list$Pat37_11989,
             plot_list$Pat37_8630,plot_list$Pat37_8821,
             plot_b,nrow=17,heights=c(1.7,rep(1,length(clone_list)), 3)) #plot_list$Pat26_9954, plot_list$Pat37_7245,
dev.off()



#### Signature generation, monocle3 output ####
DEG_treg <- read.csv("./all/CD4_DEG_alongTrajectory_TregGIMAP_TregTNF_v2.csv")
DEG_treg_top <- DEG_treg[DEG_treg$p_value < 0.01 & DEG_treg$morans_test_statistic >8,] 

DEG_CD8 <- read.csv("./all/CD8_DEG_alongTrajectory_dysTNF_dysZNF.csv")
DEG_CD8_top <- DEG_CD8[DEG_CD8$p_value < 0.01 & DEG_CD8$morans_test_statistic >8,]

# extract early - late cells, Tregs
treg_genes <- DEG_treg_top$X
cds_subset = readRDS(file = "./all/cds_subset_TregGIMAP_TregTNF_v2.rds")
cds_subset = order_cells(cds_subset)
treg_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% "LAG3",] # random gene to extract pseudotime
extract_pseudotime <- plot_genes_in_pseudotime(treg_lineage_cds,color_cells_by="mc_group", min_expr=0.4)$data
pseudotime_qs = quantile(extract_pseudotime$pseudotime,probs = c(0.3, 0.7))
Treg_cells_early = extract_pseudotime[extract_pseudotime$pseudotime < pseudotime_qs[[1]],]$Cell
Treg_cells_late = extract_pseudotime[extract_pseudotime$pseudotime > pseudotime_qs[[2]],]$Cell

# extract early - late cells, dysf
dysf_genes <- DEG_CD8_top$X
cds_subset = readRDS(file = "./all/cds_subset_TNF_ZNF.rds")
dysf_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% "LAG3",] # random gene to extract pseudotime
extract_pseudotime <- plot_genes_in_pseudotime(dysf_lineage_cds,color_cells_by="mc_group", min_expr=0.4)$data
pseudotime_qs = quantile(extract_pseudotime$pseudotime,probs = c(0.3, 0.7))
dysf_cells_early = extract_pseudotime[extract_pseudotime$pseudotime < pseudotime_qs[[1]],]$Cell
dysf_cells_late = extract_pseudotime[extract_pseudotime$pseudotime > pseudotime_qs[[2]],]$Cell

# prepare mat
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
genes_oi <- unique(DEG_treg_top$X,DEG_CD8_top$X)

cc_mat = mat@mat[intersect(genes_oi, rownames(mat@mat)),]
cc_fraction = cc_mat/colSums(mat@mat)

submat = mat@cell_metadata[names(mc@mc),]
submat = submat[submat$patient != "Pat04" &
                  submat$patient != "Pat12" &
                  submat$patient != "Pat38", ]
cc_fraction = log2(cc_fraction * 10000)
cc_fraction[!is.finite(cc_fraction)] = 0

# pseudotime Tregs
table_cells_r = as.data.frame(matrix(ncol = 0, nrow = length(treg_genes)))
rownames(table_cells_r) = rownames(cc_fraction)
table_cells_r["early_Treg"] <- rowMeans(cc_fraction[,colnames(cc_fraction) %in% Treg_cells_early],na.rm = TRUE)
table_cells_r["late_Treg"] <- rowMeans(cc_fraction[,colnames(cc_fraction) %in% Treg_cells_late],na.rm = TRUE)

# pseudotime dysf
table_cells_r["early_dysf"] <- rowMeans(cc_fraction[,colnames(cc_fraction) %in% dysf_cells_early],na.rm = TRUE)
table_cells_r["late_dysf"] <- rowMeans(cc_fraction[,colnames(cc_fraction) %in% dysf_cells_late],na.rm = TRUE)

pseudo_time_up_Treg = DEG_treg_top[DEG_treg_top$X %in% rownames(table_cells_r[table_cells_r$early_Treg < table_cells_r$late_Treg,]),]
pseudo_time_up_dysf = DEG_CD8_top[DEG_CD8_top$X %in% rownames(table_cells_r[table_cells_r$early_dysf < table_cells_r$late_dysf,]),]

overlap_DEG <- intersect(pseudo_time_up_Treg$X,pseudo_time_up_dysf$X)

DEG_treg_top$group <- rep("Treg",nrow(DEG_treg_top))
DEG_CD8_top$group <- rep("CD8",nrow(DEG_CD8_top))

write.table(DEG_treg_top,"./data/DEG_treg_top_directional_v3.txt",quote=F,row.names = F,sep="\t")
write.table(DEG_CD8_top,"./data/DEG_CD8_top_directional_v3.txt",quote=F,row.names = F,sep="\t")
write.table(overlap_DEG,"./data/overlap_DEG_directional_v3.txt",quote=F,row.names = F,sep="\t")






#### Figure 5F ####
## venn diagram
library(VennDiagram)

DEG_CD8_top = read.table("./data/DEG_CD8_top_directional_v3.txt",header=T)
DEG_treg_top = read.table("./data/DEG_treg_top_directional_v3.txt",header=T)

set1 <- DEG_CD8_top$X
set2 <- DEG_treg_top$X

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel1")

dev.off()
pdf(file=paste0(fig_path,"/Treg_Dysf_venn_diagram.pdf"), width=3, height=2)
draw.pairwise.venn(area1=length(set2), area2=length(set1),cross.area=length(intersect(set1,set2)),
                   category=c("Treg","dysf"),fill=c(myCol[1],myCol[2]))
dev.off()


#### Figure 5G ####
## correlation
DEG_CD8_top = read.table("./data/DEG_CD8_top_directional_v3.txt",header=T)
DEG_treg_top = read.table("./data/DEG_treg_top_directional_v3.txt",header=T)

geneSets <- list(geneSet_CD8=DEG_CD8_top$X,geneSet_treg=DEG_treg_top$X)

genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Dysf ZNF683","Dysf GZMB","Treg GIMAP","Treg TNFRSF9") &
              submat$response %in% c("RE","NR") &
              #!submat$patient %in% c("Pat04", "Pat12", "Pat38", NA) &
              submat$condition %in% c("post","pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

dim(exprMatrix)
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data = cells_AUC@assays@data@listData$AUC
check_data = as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_sel = sig_scores_all[sig_scores_all$mc_group %in% c("Dysf ZNF683","Dysf GZMB"),]

pdf(file=paste0(fig_path,"/Treg_vs_Dysf_signature_onDysf_v3.pdf"), width=3.5, height=2.2)
ggplot(sig_scores_sel, aes(x = geneSet_treg, y = geneSet_CD8))+
  geom_point(aes(colour=mc_group))+ #ylim(0.1, 0.55) + xlim(0.1, 0.55) +
  scale_colour_manual(values = c("#4B1C5B", "#C11989"))+ #,"#c8a2c8"
  theme_bw() + ggtitle("Treg sig vs Dysf sig") +  geom_smooth(method='lm',col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

summary(lm(geneSet_treg~geneSet_CD8,data=sig_scores_sel))

sig_scores_sel = sig_scores_all[sig_scores_all$mc_group %in% c("Treg GIMAP","Treg TNFRSF9"),]

pdf(file=paste0(fig_path,"/Treg_vs_Dysf_signature_onTregs.pdf"), width=3.5, height=2.2)
ggplot(sig_scores_sel, aes(x = geneSet_treg, y = geneSet_CD8))+
  geom_point(aes(colour=mc_group))+
  scale_colour_manual(values = c("#1F4A9F", "#56B66F"))+ #,"#c8a2c8"
  theme_bw() + ggtitle("Treg sig vs Dysf sig") +  geom_smooth(method='lm',col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

summary(lm(geneSet_treg~geneSet_CD8,data=sig_scores_sel))





#### Figure 5H ####
# shared signature plotted on dysf and treg subsets
# Treg
overlap_DEG = read.table("./Data/overlap_DEG_directional_v3.txt")
geneSets <- list(geneSet1=overlap_DEG$V1)

genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Treg GIMAP","Treg TNFRSF9") &
              submat$response %in% c("RE","NR") &
              !submat$patient %in% c("Pat04", "Pat12", "Pat38", NA) &
              submat$condition %in% c("post","pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

dim(exprMatrix)
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data = cells_AUC@assays@data@listData$AUC
check_data = as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_sel = sig_scores_all[sig_scores_all$mc_group %in% c("Treg GIMAP","Treg TNFRSF9"),]

sig_scores_sel$cond_resp = paste0(sig_scores_sel$condition,sig_scores_sel$response)
sig_scores_sel$cond_resp = factor(sig_scores_sel$cond_resp,levels=c("preNR","preRE","postNR","postRE"))
pdf(file=paste0(fig_path,"/Treg_overlap_sig_AUCell_dir.pdf"), width=3.3, height=2.8)
ggplot(sig_scores_sel, aes(y = geneSet1, x = response,fill=cond_resp))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#DC7985","#89CC9B","#9B1F23","#116936"))+
  theme_bw() + ylab("Shared sig") + ggtitle("overlap sig Treg") + ylim(c(0.02,0.55)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "RE",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "RE",][["geneSet1"]])

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "NR",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "NR",][["geneSet1"]])

# Dysf
overlap_DEG = read.table("./Data/overlap_DEG_directional_v3.txt")
geneSets <- list(geneSet1=overlap_DEG$V1)

genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Dysf ZNF683","Dysf GZMB") &
              submat$response %in% c("RE","NR") &
              !submat$patient %in% c("Pat04", "Pat12", "Pat38", NA) &
              submat$condition %in% c("post","pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

dim(exprMatrix)
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data = cells_AUC@assays@data@listData$AUC
check_data = as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_sel = sig_scores_all[sig_scores_all$mc_group %in% c("Dysf ZNF683","Dysf GZMB"),]

sig_scores_sel$cond_resp = paste0(sig_scores_sel$condition,sig_scores_sel$response)
sig_scores_sel$cond_resp = factor(sig_scores_sel$cond_resp,levels=c("preNR","preRE","postNR","postRE"))
pdf(file=paste0(fig_path,"/Dysf_overlap_sig_AUCell_dir.pdf"), width=3.3, height=2.8)
ggplot(sig_scores_sel, aes(y = geneSet1, x = response,fill=cond_resp))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#DC7985","#89CC9B","#9B1F23","#116936"))+
  theme_bw() + ylab("Shared sig") + ggtitle("overlap sig dysf") + ylim(c(0.02,0.55)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "RE",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "RE",][["geneSet1"]])

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "NR",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "NR",][["geneSet1"]])

