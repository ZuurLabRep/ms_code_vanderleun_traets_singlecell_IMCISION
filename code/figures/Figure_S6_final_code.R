###################
#### FIGURE S6 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S6"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S6")

library(monocle3)
library(VISION)
library(ComplexHeatmap)
library(circlize)
library(AUCell)


#### Figures S6A ####
# monocle analysis Tregs, pseudotime, generated in figure 5
cds = readRDS("./data/cds_Dysf_all.rds")
cds <- cluster_cells(cds, k = 19)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
cds <- learn_graph(cds)
order_cells(cds) # select naive

pdf(file=paste0(fig_path,"/CD8_UMAP_lineage_pseudotime_genelist_v3.pdf"), width=4.3, height=3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,trajectory_graph_color = "black",
           label_branch_points=FALSE,trajectory_graph_segment_size = 1.2,
           graph_label_size=1.5,label_roots = FALSE,
           cell_size = 0.7)
dev.off()


#### Figure S5B ####
annotation_ptime <- read.table("./data/Pseuodotime_dysf.txt",sep="\t",header=T)

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

submat <- mat@cell_metadata       
md = submat[rownames(submat) %in% unique(annotation_ptime$Cell),] 

count_mat <- as.matrix(mat@mat[,rownames(md)])

n.umi <- colSums(count_mat)
scaled_counts <- t(t(count_mat) / n.umi) * median(n.umi)
norm_counts <- sweep(count_mat,2,colSums(count_mat), "/") * 1000

md$Cell <- rownames(md)

vis <- Vision(scaled_counts,
              signatures = c("./data/CD8_dysf_trajectory_topGO.txt"),
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

binned_mean <- tapply(collect_ptime$GOMF_CHEMOKINE_ACTIVITY, cut(collect_ptime$pseudotime, 30), mean,na.rm = FALSE)
binned_mean_all <- as.data.frame(binned_mean)
colnames(binned_mean_all)[1] <- "GOMF_CHEMOKINE_ACTIVITY"
for (gene in colnames(collect_ptime)[3:6]){
  binned_mean <- tapply(collect_ptime[[gene]], cut(collect_ptime$pseudotime, 30), mean)
  binned_mean <- as.data.frame(binned_mean)
  colnames(binned_mean)[1] <- gene
  binned_mean_all <- cbind(binned_mean_all,binned_mean)
}

# saved from figure 5B
recast_freq_bin_df_save = read.table("./data/CD8_pseudotime_barplot_mcgroups.txt",sep="\t",header=F)
colnames(recast_freq_bin_df_save) = recast_freq_bin_df_save[1,]
recast_freq_bin_df_save = recast_freq_bin_df_save[-1,]
recast_freq_bin_df_save <- recast_freq_bin_df_save %>% 
  mutate_if(is.character, as.numeric)
ha = HeatmapAnnotation(foo = anno_barplot(t(recast_freq_bin_df_save[,2:31]), gp = gpar(fill =c("#E5D3B3", "#4C005C","#C20088","#ff7f7f","#c8a2c8")), 
                                          bar_width = 1, height = unit(1.3, "cm")))

temp_data_t <- t(binned_mean_all[rev(rownames(binned_mean_all)),])
temp_data_t <- rev(as.data.frame(temp_data_t))
x <- t(scale(t(temp_data_t)))
col_fun = colorRamp2(c(-2, 0, 2), c("#236EAF", "white", "#BC2B33"))
pdf(file=paste0(fig_path,"/CD8_pseudotime_GOterm_enrichement_top5.pdf"),width=6, height=3)
Heatmap(x, cluster_columns = FALSE, col = col_fun,top_annotation = ha,
        row_names_side = "right",row_names_gp = gpar(fontsize = 4),show_row_names = TRUE)
dev.off()



#### Figure S5C ####
## GSEA on trajectory dysfunctionals
# subset dysfunctional cells
# load from folder
cds_subset <- readRDS(file = "./data/cds_subset_TNF_ZNF.rds")

# init
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)

### rnk file on dysfunctional cells trajectory
res_df_sig <- as.data.frame(cbind(name=rownames(pr_graph_test_res),sign= pr_graph_test_res$morans_test_statistic))

res_df_sig <- res_df_sig[order(res_df_sig$sign),]
res_df_sig <- res_df_sig[res_df_sig$name!= "",]
res_df_sig <- res_df_sig[res_df_sig$sign!= "",]

res_df_sig <- res_df_sig[!is.na(res_df_sig[,1]),]
res_df_sig <- res_df_sig[!is.na(res_df_sig[,2]),]

## ready for GSEA Java tools
write.table(res_df_sig, file="./data/Genes_sig_trajectory_dysf_v3.rnk", sep="\t", row.names=F, col.names=F,  quote=F)

## runned on GSEA Java tools
GSEA_pos <- read.csv("./data/Genes_sig_trajectory_dysf_GOterms/gsea_report_for_na_pos_1655114525918.tsv",sep="\t")
#GSEA_neg <-read.csv("./HN_full_dataset_figs_final_last/Genes_sig_trajectory_tregs_GOterms/gsea_report_for_na_neg_1647267851795.tsv",sep="\t")

GSEA_pos_top <- GSEA_pos[order(GSEA_pos$NOM.p.val),][1:20,]
all_GSEA <- GSEA_pos_top
all_GSEA <- all_GSEA[order(all_GSEA$FDR.q.val,decreasing = T),]
all_GSEA$NAME <- factor(all_GSEA$NAME,levels= all_GSEA$NAME)

pdf(file=paste0(fig_path,"/Trajectory_Dysf_moransteststat_barplot.pdf"), width=7, height=6) #
p = ggplot(all_GSEA) + geom_point(aes(y=NAME,x=-log10(FDR.q.val),size=SIZE)) + geom_vline(xintercept=0) +
  scale_size(range=c(2,7)) + geom_hline(yintercept = c(1:nrow(GSEA_pos)), colour="grey") + xlim(c(-3,3)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p
dev.off()



#### Figure S6D ####
# cell cycle genes
mel_gset = scdb_gset(lateral_gset_mel_id)
cc_genes = unique(c(names(mel_gset@gene_set[mel_gset@gene_set == 1])))

IFN_umis = mat@mat[intersect(cc_genes, rownames(mat@mat)),] ## CHANGE
IFN_fraction <<- colSums(IFN_umis)/colSums(mat@mat)
submat = mat@cell_metadata[names(mc@mc),]
submat = submat[!is.na(submat$patient) & (submat$patient != "Pat27") &
                  submat$mc_group %in% c("Dysf GZMB", "Dysf ZNF683"),]
submat$CXCL8 = IFN_fraction[rownames(submat)]

th = log2(0.025)
percentage = submat %>%
  group_by(response, condition) %>%
  summarize(total = n(),above = sum(log2(CXCL8) > th), fraction = (above / total)*100)
print(percentage)

submat$response_cond = paste0(submat$response,"_",submat$condition)
pdf(file=paste0(fig_path,"/Cell_cycle_dysf_RE_NR_v2.pdf"),width=4, height=5)
hist = ggplot(submat, aes(x = log2(CXCL8)))+
  geom_density(aes(group = condition, color = response_cond), alpha = 0.2, adjust = 4)+
  facet_wrap(~response)+
  geom_vline(xintercept = th, linetype = "dashed", color = "lightgrey")+
  xlim(-12.5, -2.5) +
  scale_color_manual(values=c("#901A1D","#E18283","#07522B","#71BB87")) + theme_bw() + 
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(hist)
dev.off()



set.seed(123)
#### Figure S6E ####
tcr_sig = read.table("./data/NIHMS1790595-supplement-TABLE_4_signature_NeoTCR.csv",sep=",",header=T)

mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]

geneSets <- list(neoTCRCD4=tcr_sig$NeoTCR4,neoTCRCD8=tcr_sig$NeoTCR8)
genesets_col = geneSets

remove = c("Pat04", "Pat12", "Pat38")
mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Dysf GZMB","Dysf ZNF683") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("pre","post") &
              !is.na(submat$unique_clone_ID) &
              !(submat$patient %in% remove),] 
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

# RE
sig_scores_all_RE_pre = sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "pre",]
sig_scores_all_RE_post = sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "post",]
sig_scores_all_RE_freq_pre = sig_scores_all_RE_pre %>% group_by(patient,pat_clone_ID) %>% summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
check_pre = sig_scores_all_RE_freq_pre[sig_scores_all_RE_freq_pre$n >1,]
sig_scores_all_RE_freq_post = sig_scores_all_RE_post %>% group_by(patient,pat_clone_ID) %>% summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
check_post = sig_scores_all_RE_freq_post[sig_scores_all_RE_freq_post$n >1,]

overlap_tcr = intersect(check_pre$pat_clone_ID,check_post$pat_clone_ID)
overlap_tcr_pat = sig_scores_all[sig_scores_all$pat_clone_ID %in% overlap_tcr,]

to_check= as.data.frame(overlap_tcr)
to_check$avg_neoTCR_pre = 0
to_check$avg_neoTCR_post = 0
to_check$avg_neoTCR_all = 0
for (tcr_x in overlap_tcr){
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_all = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x,]$neoTCRCD8)
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_pre = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x & sig_scores_all$condition == "pre",]$neoTCRCD8)
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_post = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x & sig_scores_all$condition == "post",]$neoTCRCD8)
}

to_check$patient = do.call(rbind,strsplit(to_check$overlap_tcr,"_"))[,1]
keep_tcrs = c()
for (pat_x in unique(to_check$patient)){
  per_50_top = quantile(to_check[to_check$patient == pat_x,]$avg_neoTCR_all)[4]
  keep_tcrs = c(keep_tcrs,to_check[to_check$patient == pat_x & to_check$avg_neoTCR_all > per_50_top,]$overlap_tcr)
}

to_plot = to_check[to_check$overlap_tcr %in% keep_tcrs,]

#per_75_top = quantile(to_check$avg_neoTCR_all)
#to_plot = to_check[to_check$avg_neoTCR_all >per_75_top[3],]

to_plot_mlt = reshape2::melt(to_plot[,c(1:3)]) 
pdf(file=paste0(fig_path,"/NeoTCR_score_top25_perpatient_min2_RE_sm.pdf"), width=2,height=3)
ggplot(to_plot_mlt,aes(y=value,x=variable,fill=variable)) + geom_violin() +
  ggpubr::stat_compare_means(paired=T)  + theme_bw() + geom_point(aes(y=value,x=variable)) +
  geom_line(aes(y=value,x=variable,group=overlap_tcr)) + 
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("top 25% clones RE") + scale_fill_manual(values=c("#71bc88","#0d542c"))+ ylim(c(0.1,0.26))
dev.off()

# NR
sig_scores_all_RE_pre = sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "pre",]
sig_scores_all_RE_post = sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "post",]
sig_scores_all_RE_freq_pre = sig_scores_all_RE_pre %>% group_by(patient,pat_clone_ID) %>% summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
check_pre = sig_scores_all_RE_freq_pre[sig_scores_all_RE_freq_pre$n >1,]
sig_scores_all_RE_freq_post = sig_scores_all_RE_post %>% group_by(patient,pat_clone_ID) %>% summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
check_post = sig_scores_all_RE_freq_post[sig_scores_all_RE_freq_post$n >1,]

overlap_tcr = intersect(check_pre$pat_clone_ID,check_post$pat_clone_ID)

overlap_tcr_pat = sig_scores_all[sig_scores_all$pat_clone_ID %in% overlap_tcr,]

to_check= as.data.frame(overlap_tcr)
to_check$avg_neoTCR_pre = 0
to_check$avg_neoTCR_post = 0
to_check$avg_neoTCR_all = 0
for (tcr_x in overlap_tcr){
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_all = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x,]$neoTCRCD8)
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_pre = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x & sig_scores_all$condition == "pre",]$neoTCRCD8)
  to_check[to_check$overlap_tcr == tcr_x,]$avg_neoTCR_post = mean(sig_scores_all[sig_scores_all$pat_clone_ID == tcr_x & sig_scores_all$condition == "post",]$neoTCRCD8)
}

to_check$patient = do.call(rbind,strsplit(to_check$overlap_tcr,"_"))[,1]
keep_tcrs = c()
for (pat_x in unique(to_check$patient)){
  per_50_top = quantile(to_check[to_check$patient == pat_x,]$avg_neoTCR_all)[4]
  keep_tcrs = c(keep_tcrs,to_check[to_check$patient == pat_x & to_check$avg_neoTCR_all > per_50_top,]$overlap_tcr)
}

to_plot = to_check[to_check$overlap_tcr %in% keep_tcrs,]

#per_75_top = quantile(to_check$avg_neoTCR_all)
#to_plot = to_check[to_check$avg_neoTCR_all >per_75_top[3],]

to_plot_mlt = reshape2::melt(to_plot[,c(1:3)])
pdf(file=paste0(fig_path,"/NeoTCR_score_top25_perpatient_min2_NR_sm.pdf"), width=2,height=3)
ggplot(to_plot_mlt,aes(y=value,x=variable,fill=variable)) + geom_violin() +
  ggpubr::stat_compare_means(paired=T)  + theme_bw() + geom_point(aes(y=value,x=variable)) +
  geom_line(aes(y=value,x=variable,group=overlap_tcr)) + 
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("top 25% clones NR") + scale_fill_manual(values=c("#dc7986","#921d1f"))+ ylim(c(0.1,0.26))
dev.off()



#### Figure S6F ####
# top clones
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

# non responders
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
#remove = c("None")
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "NR" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Dysf GZMB","Dysf ZNF683"),]

# determine top clones pre on in CD8 dysf
md$pat_cond_clone <- paste0(md$patient,"_",md$condition,"_",md$unique_clone_ID)
md_clones <- as.data.frame(table(md$pat_cond_clone))
sep_md <- as.data.frame(do.call(rbind, strsplit(as.character(md_clones$Var1), '_')))
md_clones <- cbind(md_clones,sep_md)
md_clones <- md_clones[order(-md_clones$Freq),]

md_clones_pre <- md_clones[md_clones$V2 == "pre",]
md_clones_post <- md_clones[md_clones$V2 == "post",]

colnames(md_clones_pre) <- c("all_pre","freq_pre","pat_id","cond_pre","clone_id_pre")
colnames(md_clones_post) <- c("all_post","freq_post","pat_id","cond_post","clone_id_post")

md_clones_pre$pat_clone_id <- paste0(md_clones_pre$pat_id,"_",md_clones_pre$clone_id_pre)
md_clones_post$pat_clone_id <- paste0(md_clones_post$pat_id,"_",md_clones_post$clone_id_post)

md_clones_all <- merge(md_clones_pre,md_clones_post,by="pat_clone_id")
md_clones_top <- md_clones_all[md_clones_all$freq_pre > 2 & md_clones_all$freq_post > 2,]
md_clones_top_id <- md_clones_top$pat_clone_id

# load everything again, all cell types CD8
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "NR" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional","Non-classical T cell","Dysf GZMB","Dysf ZNF683","Naive-like CD8","Cytotoxic"),]

md$pat_clone_ID <- paste0(md$patient,"_",md$unique_clone_ID)
md_toppers <- md[md$pat_clone_ID %in% md_clones_top_id,]
md_toppers_pre <- md_toppers[md_toppers$condition == "pre",]
md_toppers_post <- md_toppers[md_toppers$condition == "post",]

sum_md_toppers_post <- table(md_toppers_post[,c("pat_clone_ID","mc_group")]) # some TCRs also in CD4 or even myeloid!
sum_md_toppers_post <- cbind(sum_md_toppers_post, total = rowSums(sum_md_toppers_post))
sum_md_toppers_post <- as.data.frame.matrix(sum_md_toppers_post)
sum_md_toppers_post_sw_NR <- sweep(sum_md_toppers_post[,-4], 1, sum_md_toppers_post[,4], "/")

sum_md_toppers_pre <- table(md_toppers_pre[,c("pat_clone_ID","mc_group")]) # some TCRs also in CD4 or even myeloid!
sum_md_toppers_pre <- cbind(sum_md_toppers_pre, total = rowSums(sum_md_toppers_pre))
sum_md_toppers_pre <- as.data.frame.matrix(sum_md_toppers_pre)
sum_md_toppers_pre_sw_NR <- sweep(sum_md_toppers_pre[,-6], 1, sum_md_toppers_pre[,6], "/")

sum_md_toppers_post_sw_NR$`Naive-like CD8` <- rep(0.0,nrow(sum_md_toppers_post_sw_NR))
sum_md_toppers_post_sw_NR$`Transitional` <- rep(0.0,nrow(sum_md_toppers_post_sw_NR))
sum_md_toppers_post_sw_NR <- sum_md_toppers_post_sw_NR[,c(colnames(sum_md_toppers_pre_sw_NR))]

sum_md_fraction_change_NR <- sum_md_toppers_post_sw_NR-sum_md_toppers_pre_sw_NR

sum_md_toppers_post$`Naive-like CD8` <- rep(0.0,nrow(sum_md_toppers_post))
sum_md_toppers_post$`Transitional` <- rep(0.0,nrow(sum_md_toppers_post))
sum_md_toppers_post <- sum_md_toppers_post[,c(colnames(sum_md_toppers_pre_sw_NR))]
sum_md_toppers_pre <- sum_md_toppers_pre[,c(colnames(sum_md_toppers_pre_sw_NR))]
sum_md_toppers_post_pre <- sum_md_toppers_post + sum_md_toppers_pre

sum_md_toppers_post_pre$clone_id <- rownames(sum_md_toppers_post_pre)
temp <- as.data.frame(rowSums(sum_md_toppers_post_pre[,1:5]))
temp$clone_id <- rownames(temp)
order_sort <- rownames(temp[order(temp$`rowSums(sum_md_toppers_post_pre[, 1:5])`,decreasing = T),])
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[match(order_sort,sum_md_toppers_post_pre$clone_id),]
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[1:nrow(sum_md_toppers_post_pre),]
sum_md_toppers_post_pre_melt <- melt(sum_md_toppers_post_pre)

sum_md_fraction_change_NR$clone_id <- rownames(sum_md_fraction_change_NR)
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[match(order_sort,sum_md_fraction_change_NR$clone_id),]
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[1:nrow(sum_md_fraction_change_NR),]
sum_md_fraction_change_NR_melt <- melt(sum_md_fraction_change_NR)

sum_md_fraction_change_NR_melt$clone_id <- factor(sum_md_fraction_change_NR_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable,levels=c( "Non-classical T cell","Dysf GZMB", "Dysf ZNF683", "Naive-like CD8", "Transitional"))
sum_md_fraction_change_NR_melt$variable <- factor(sum_md_fraction_change_NR_melt$variable,levels=c( "Non-classical T cell","Dysf GZMB", "Dysf ZNF683", "Naive-like CD8", "Transitional"))

pdf(file=paste0(fig_path,"/Change_clones_cell_states_NR_dotplot_order_fix_v2_top15.pdf"),width=3.9, height=5.2)
ggplot(sum_md_fraction_change_NR_melt,aes(x=variable,y=clone_id)) + 
  geom_point(color="black",shape=21,aes(fill=value, size=sum_md_toppers_post_pre_melt$value)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", limits = c(-1, 1),name = "lfc") +
  scale_size_continuous(range = c(0,8),name = "clone_size_sum",limits=c(0,82)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# responders
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
#remove = c("None")
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "RE" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Dysf GZMB","Dysf ZNF683"),]

# determine top clones pre on in CD8 dysf
md$pat_cond_clone <- paste0(md$patient,"_",md$condition,"_",md$unique_clone_ID)
md_clones <- as.data.frame(table(md$pat_cond_clone))
sep_md <- as.data.frame(do.call(rbind, strsplit(as.character(md_clones$Var1), '_')))
md_clones <- cbind(md_clones,sep_md)
md_clones <- md_clones[order(-md_clones$Freq),]

md_clones_pre <- md_clones[md_clones$V2 == "pre",]
md_clones_post <- md_clones[md_clones$V2 == "post",]

colnames(md_clones_pre) <- c("all_pre","freq_pre","pat_id","cond_pre","clone_id_pre")
colnames(md_clones_post) <- c("all_post","freq_post","pat_id","cond_post","clone_id_post")

md_clones_pre$pat_clone_id <- paste0(md_clones_pre$pat_id,"_",md_clones_pre$clone_id_pre)
md_clones_post$pat_clone_id <- paste0(md_clones_post$pat_id,"_",md_clones_post$clone_id_post)

md_clones_all <- merge(md_clones_pre,md_clones_post,by="pat_clone_id")
md_clones_top <- md_clones_all[md_clones_all$freq_pre > 3 &md_clones_all$freq_post > 3,]
md_clones_top_id <- md_clones_top$pat_clone_id

# load everything again, all cell types CD8
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
#remove = c("None")
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "RE" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional","Non-classical T cell","Dysf GZMB","Dysf ZNF683","Naive-like CD8","Cytotoxic"),]

md$pat_clone_ID <- paste0(md$patient,"_",md$unique_clone_ID)
md_toppers <- md[md$pat_clone_ID %in% md_clones_top_id,]
md_toppers_pre <- md_toppers[md_toppers$condition == "pre",]
md_toppers_post <- md_toppers[md_toppers$condition == "post",]

sum_md_toppers_post <- table(md_toppers_post[,c("pat_clone_ID","mc_group")]) # some TCRs also in CD4 or even myeloid!
sum_md_toppers_post <- cbind(sum_md_toppers_post, total = rowSums(sum_md_toppers_post))
sum_md_toppers_post <- as.data.frame.matrix(sum_md_toppers_post)
sum_md_toppers_post_sw_RE <- sweep(sum_md_toppers_post[,-6], 1, sum_md_toppers_post[,6], "/")

sum_md_toppers_pre <- table(md_toppers_pre[,c("pat_clone_ID","mc_group")]) # some TCRs also in CD4 or even myeloid!
sum_md_toppers_pre <- cbind(sum_md_toppers_pre, total = rowSums(sum_md_toppers_pre))
sum_md_toppers_pre <- as.data.frame.matrix(sum_md_toppers_pre)
sum_md_toppers_pre_sw_RE <- sweep(sum_md_toppers_pre[,-6], 1, sum_md_toppers_pre[,6], "/")

sum_md_fraction_change_RE <- sum_md_toppers_post_sw_RE-sum_md_toppers_pre_sw_RE

# strip plot pre to on, per R and NR
sum_md_toppers_post <- sum_md_toppers_post[,c(colnames(sum_md_toppers_post))]
sum_md_toppers_pre <- sum_md_toppers_pre[,c(colnames(sum_md_toppers_post))]
sum_md_toppers_post_pre <- sum_md_toppers_post + sum_md_toppers_pre

sum_md_toppers_post_pre$clone_id <- rownames(sum_md_toppers_post_pre)
temp <- as.data.frame(rowSums(sum_md_toppers_post_pre[,1:5]))
temp$clone_id <- rownames(temp)
order_sort <- rownames(temp[order(temp$`rowSums(sum_md_toppers_post_pre[, 1:5])`,decreasing = T),])
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[match(order_sort,sum_md_toppers_post_pre$clone_id),]
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[1:15,]
sum_md_toppers_post_pre_melt <- melt(sum_md_toppers_post_pre[,-6])

sum_md_fraction_change_RE$clone_id <- rownames(sum_md_fraction_change_RE)
sum_md_fraction_change_RE <- sum_md_fraction_change_RE[match(order_sort,sum_md_fraction_change_RE$clone_id),]
sum_md_fraction_change_RE <- sum_md_fraction_change_RE[1:15,]
sum_md_fraction_change_RE_melt <- melt(sum_md_fraction_change_RE)

sum_md_fraction_change_RE_melt$clone_id <- factor(sum_md_fraction_change_RE_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_fraction_change_RE_melt$variable <- factor(sum_md_fraction_change_RE_melt$variable, levels=c( "Non-classical T cell","Dysf GZMB", "Dysf ZNF683", "Naive-like CD8", "Transitional"))
sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable, levels=c( "Non-classical T cell","Dysf GZMB", "Dysf ZNF683", "Naive-like CD8", "Transitional"))

pdf(file=paste0(fig_path,"/Change_clones_cell_states_RE_dotplot_order_legend_fix_v2.pdf"), width=3.9, height=5.2)
ggplot(sum_md_fraction_change_RE_melt,aes(x=variable,y=clone_id)) + 
  geom_point(color="black",shape=21,aes(fill=value, size=sum_md_toppers_post_pre_melt$value)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", limits = c(-1, 1),name = "lfc") +
  scale_size_continuous(range = c(0,8),name = "clone_size_sum",limits=c(0,80)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
