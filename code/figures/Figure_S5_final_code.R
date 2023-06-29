###################
#### FIGURE S5 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S5"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S5")

library(monocle3)
library(ComplexHeatmap)
library(circlize)


#### Figure S5A ####
# monocle analysis Tregs, pseudotime, generated in figure 4
cds = readRDS("./data/cds_TregGIMAP_TregTNF_all.rds")
cds <- cluster_cells(cds, k = 12)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
cds <- learn_graph(cds)
order_cells(cds) # select naive

pdf(file=paste0(fig_path,"/CD4_UMAP_lineage_pseudotime_genelist_v3.pdf"), width=4.3, height=3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,trajectory_graph_color = "black",
           label_branch_points=FALSE,trajectory_graph_segment_size = 1.2,
           graph_label_size=1.5,label_roots = FALSE,
           cell_size = 0.7)
dev.off()



#### Figure S5B ####
# correlation subpopulations
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

md = md[!is.na(md$cell_type) &
          md$mc_group %in% c("Treg TNFRSF9","Treg GIMAP","Tfh LAG3","Tfh NR3C1","Naive-like CD4","Naive-like CD8","Dysf GZMB","Dysf ZNF683","Transitional","Cytotoxic"),]

cells_oi = rownames(md)
umis =  mat@mat[,cells_oi]
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 
mat_oi = as.data.frame(as.matrix(umis_n))

gene_set <- read.table("./data/HN_full_dataset_clean_feats_filtered.txt")$V1

mat_oi = as.data.frame(t(mat_oi))
mat_oi_sel = mat_oi[,colnames(mat_oi) %in% gene_set]
mat_oi_sel = merge(mat_oi_sel,md[,c(1,94)],by=0)

mat_oi_sel[is.na(mat_oi_sel)] = 0

avg_exp_mat <- mat_oi_sel %>% group_by(mc_group) %>%
  summarise_each(funs(mean))

to_plot = as.data.frame(avg_exp_mat)
rownames(to_plot) = to_plot[,1]
to_plot = to_plot[,-c(1:2)]
to_plot[is.na(to_plot)] = 0
to_plot = as.data.frame(t(to_plot))
to_plot = log(to_plot)+10

to_plot_cor = cor(to_plot[,c("Naive-like CD4","Tfh LAG3","Tfh NR3C1","Treg GIMAP","Treg TNFRSF9")],method = "spearman")

pdf(file=paste0(fig_path,"/Correlation_heatmap_CD4_spearman.pdf"), width=3.5, height=3)
Heatmap(to_plot_cor,column_title="correlations CD4", name="spearman",
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8))
dev.off()



#### Figure S5C ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
meta_dat = as.data.frame(mat@cell_metadata)

remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Treg TNFRSF9","Treg GIMAP"),]

md$pat_clone <- paste0(md$patient,"_",md$unique_clone_ID)

md_pre <- md[md$condition == "pre",]
md_post <- md[md$condition == "post",]

shared_clones <- intersect(md_post$pat_clone,md_pre$pat_clone)
data_pre <- md_pre %>% group_by(patient,pat_clone) %>% dplyr::summarize(n=n()) %>% mutate(freq_pre=n/sum(n)) 
data_post <- md_post %>% group_by(patient,pat_clone) %>% dplyr::summarize(n=n()) %>% mutate(freq_post=n/sum(n)) 
shared_pre <- data_pre[data_pre$pat_clone %in% shared_clones,]
shared_post <- data_post[data_post$pat_clone %in% shared_clones,]
shared_both <- cbind(shared_pre,shared_post)

md_pre_RE <- md_pre[md_pre$response == "RE",]
md_pre_NR <- md_pre[md_pre$response == "NR",]
clones_baseline <- as.data.frame(table(md_pre_RE$pat_clone))
clones_baseline <- clones_baseline[order(clones_baseline$Freq,decreasing = T),]
clones_baseline_RE <- clones_baseline[clones_baseline$Freq >2,][1:10,]
clones_baseline <- as.data.frame(table(md_pre_NR$pat_clone))
clones_baseline <- clones_baseline[order(clones_baseline$Freq,decreasing = T),]
clones_baseline_NR <- clones_baseline[clones_baseline$Freq >2,][1:10,]

expanded_both <- shared_both[shared_both$n...3 > 1 & shared_both$n...7 > 1,]
expanded_post <- shared_both[shared_both$n...3 > 0 & shared_both$n...7 > 1,]
expanded_both_all <- shared_both[shared_both$n...3 > 0 & shared_both$n...7 > 0,]

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
meta_dat = as.data.frame(mat@cell_metadata)

remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Treg TNFRSF9","Treg GIMAP"),]

md$pat_clone <- paste0(md$patient,"_",md$unique_clone_ID)

md_pre <- md[md$condition == "pre",]
md_post <- md[md$condition == "post",]

## frequency plot
data_pre <- md_pre %>% group_by(patient,pat_clone) %>% dplyr::summarize(n=n()) %>% mutate(freq_pre=n/sum(n)) 
data_post <- md_post %>% group_by(patient,pat_clone) %>% dplyr::summarize(n=n()) %>% mutate(freq_post=n/sum(n)) 

## check clones Tregs in data
all_clones_tregs <- merge(data_pre,data_post, by= "pat_clone",all=T)
md = mat@cell_metadata
md$pat_clone = paste0(md$patient,"_",md$unique_clone_ID)
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$pat_clone %in% all_clones_tregs$pat_clone &
          md$cell_type == "CD4",]

md_pre <- md[md$condition == "pre",]
md_post <- md[md$condition == "post",]

# both NR and RE
data_pre <- md_pre[md_pre$pat_clone %in% c(clones_baseline_RE$Var1,clones_baseline_NR$Var1),] %>% group_by(patient,pat_clone,mc_group) %>% dplyr::summarize(n=n()) %>% mutate(freq_pre=n/sum(n),total = sum(n)) 
data_pre_freq <- data_pre[,c("pat_clone","total")]
data_pre_freq <- data_pre_freq[!duplicated(data_pre_freq$pat_clone),]
data_pre_freq$y <- 1

pdf(file=paste0(fig_path,"/Barplot_top_clones_Treg_shared_pre_RE_NR.pdf"), width=6.5,height=3)
ggplot() + geom_bar(data=data_pre,aes(x=pat_clone,y=freq_pre,fill=mc_group),stat="identity") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Treg top pre") +
  geom_text(data=data_pre_freq, aes(x=pat_clone, y=y, label=as.factor(total)),vjust=-0.1) +
  scale_fill_manual(values=c("#107538","#0042b0","#59b56f"))
dev.off()




#### Figure S5D ####
## TCR clones Tregs
annotation_ptime <-read.table("./data/Treg_trajectory_all.txt",sep="\t",header=T) 
annotation_ptime  = annotation_ptime[,-c(7,8)]
annotation_ptime = annotation_ptime[!duplicated(annotation_ptime$Cell),]

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
mat_meta = mat@cell_metadata

mat_meta$Cell <- rownames(mat_meta)
mat_cell <- merge(annotation_ptime,mat_meta,by="Cell")

mat_cell$binned_time <- cut(mat_cell$pseudotime, 30)

mat_cell$pat_clone_ID_v2 <- paste0(mat_cell$unique_clone_ID,"_",mat_cell$patient)
mat_cell <- mat_cell[!is.na(mat_cell$unique_clone_ID),]

cds_subset_pre <- mat_cell[mat_cell$condition.x == "pre",]
cds_subset_post <- mat_cell[mat_cell$condition.x == "post",]
check_tregs_pre <- as.data.frame(table(cds_subset_pre$pat_clone_ID_v2))
check_tregs_pre <- check_tregs_pre[order(check_tregs_pre$Freq,decreasing = T),]
check_tregs_post <- as.data.frame(table(cds_subset_post$pat_clone_ID_v2))
check_tregs_post <- check_tregs_post[order(check_tregs_post$Freq,decreasing = T),]

shared <- merge(check_tregs_post,check_tregs_pre, by = "Var1")
shared_good <- shared[shared$Freq.x >= 2 & shared$Freq.y >= 2,]
shared_expand <- shared[shared$Freq.x > shared$Freq.y,]

mat_cell$Treg_TCR <- rep("none",nrow(mat_cell))
mat_cell[mat_cell$pat_clone_ID_v2 %in% shared$Var1,]["Treg_TCR"] <- "shared"
mat_cell[mat_cell$pat_clone_ID_v2 %in% shared_expand$Var1,]["Treg_TCR"] <- "shared_expanded"

for (x_n in c(1:nrow(shared_expand))){
  mat_cell$TCR_to_plot = rep(0,nrow(mat_cell))
  mat_cell[mat_cell$pat_clone_ID_v2 == shared_expand$Var1[x_n],]$TCR_to_plot = 1
}

cds_subset_pre <- mat_cell[mat_cell$condition.x == "pre" & mat_cell$response.x =="RE",]
cds_subset_post <- mat_cell[mat_cell$condition.x == "post" & mat_cell$response.x =="RE",]
check_tregs_pre <- as.data.frame(table(cds_subset_pre$pat_clone_ID_v2))
check_tregs_pre <- check_tregs_pre[order(check_tregs_pre$Freq,decreasing = T),]
check_tregs_post <- as.data.frame(table(cds_subset_post$pat_clone_ID_v2))
check_tregs_post <- check_tregs_post[order(check_tregs_post$Freq,decreasing = T),]

shared <- merge(check_tregs_post,check_tregs_pre, by = "Var1")
shared_good <- shared[shared$Freq.x >= 2 & shared$Freq.y >= 2,]
shared_expand <- shared[shared$Freq.x > shared$Freq.y,]
shared_all <- rbind(shared_good,shared_expand)

clone_list <- shared_all$Var1

plot_list <- list()
for (name_clone in clone_list){
  mat_cell$clone_check <- rep("0",nrow(mat_cell))
  mat_cell[mat_cell$pat_clone_ID_v2 == name_clone,]["clone_check"] <- "1"
  mat_cell$clone_check <- as.numeric(mat_cell$clone_check)
  
  set.seed(1)
  mat_cell[mat_cell$clone_check =="1",]["clone_check"] <- runif(nrow(mat_cell[mat_cell$clone_check=="1",]), min=0.5, max=0.5)
  
  plot_list[[name_clone]] = ggplot(mat_cell, aes(x=pseudotime, y=clone_check)) + geom_point(aes(color=condition_response),size=2) + 
    scale_color_manual(values=c("#B4D4F0", "#075227","#73ba86","#73ba86")) + ylim(c(0.3,0.7)) + 
    theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",
                       axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    ggtitle(name_clone) + ylab("") + xlab("") 
}

mat_cell$random_nr <- runif(nrow(mat_cell), min=0.1, max=0.9)

plot_b = ggplot(mat_cell, aes(x=pseudotime, y=random_nr)) + geom_point(aes(color=mc_group.x),size=1) + 
  scale_color_manual(values=c("#E5D3B3", "#A05124","#12763C","#224BA0","#58B46E")) + ylim(c(0.1,0.9)) + 
  theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cell_states") + ylab("") + xlab("") 

graph_data_sel = mat_cell[mat_cell$pat_clone_ID_v2 %in% clone_list,]

plot_d = ggplot(graph_data_sel, aes(x = pseudotime, colour = condition_response, fill = condition_response)) + 
  geom_density(lwd = 0.7, alpha = 0.25) + scale_color_manual(values=c("#0b522a","#73ba86")) +
  scale_fill_manual(values=c("#0b522a","#73ba86")) + ylim(0,0.5) + xlim(min(mat_cell$pseudotime),max(mat_cell$pseudotime)) +
  theme_bw() + theme(plot.margin=unit(c(+0.4,1,-0.1,1),"cm"),legend.position = "none", panel.grid.major = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),axis.ticks.y = element_blank(),
                     panel.grid.minor = element_blank())

library(gridExtra)
pdf(file=paste0(fig_path,"/CD4_pseudotime_top_clones_Responders_update.pdf"), width=3.5, height=3.4)
grid.arrange(plot_d,plot_list$`12355_Pat21`,plot_list$`2794_Pat15`,plot_list$`672_Pat21`,plot_list$`9325_Pat39`,
             plot_list$`13733_Pat21`, plot_list$`3426_Pat21`,plot_list$`6887_Pat15`,plot_list$`7505_Pat39`,plot_list$`8748_Pat15`,
             plot_b,nrow=11,heights=c(1.7,rep(1,length(clone_list)), 3))
dev.off()

cds_subset_pre <- mat_cell[mat_cell$condition.x == "pre" & mat_cell$response.x =="NR",]
cds_subset_post <- mat_cell[mat_cell$condition.x == "post" & mat_cell$response.x =="NR",]
check_tregs_pre <- as.data.frame(table(cds_subset_pre$pat_clone_ID_v2))
check_tregs_pre <- check_tregs_pre[order(check_tregs_pre$Freq,decreasing = T),]
check_tregs_post <- as.data.frame(table(cds_subset_post$pat_clone_ID_v2))
check_tregs_post <- check_tregs_post[order(check_tregs_post$Freq,decreasing = T),]

shared <- merge(check_tregs_post,check_tregs_pre, by = "Var1")
shared_good <- shared[shared$Freq.x >= 2 & shared$Freq.y >= 2,]
shared_expand <- shared[shared$Freq.x > shared$Freq.y,]

clone_list <- shared_expand$Var1

plot_list <- list()
for (name_clone in clone_list){
  mat_cell$clone_check <- rep("0",nrow(mat_cell))
  mat_cell[mat_cell$pat_clone_ID_v2 == name_clone,]["clone_check"] <- "1"
  mat_cell$clone_check <- as.numeric(mat_cell$clone_check)
  
  set.seed(1)
  mat_cell[mat_cell$clone_check =="1",]["clone_check"] <- runif(nrow(mat_cell[mat_cell$clone_check=="1",]), min=0.5, max=0.5)
  
  plot_list[[name_clone]] = ggplot(mat_cell, aes(x=pseudotime, y=clone_check)) + geom_point(aes(color=condition_response),size=2) + 
    scale_color_manual(values=c("#911B1E", "#075227","#E18183","#73ba86")) + ylim(c(0.3,0.7)) + 
    theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",
                       axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    ggtitle(name_clone) + ylab("") + xlab("") 
}

mat_cell$random_nr <- runif(nrow(mat_cell), min=0.1, max=0.9)

plot_b = ggplot(mat_cell, aes(x=pseudotime, y=random_nr)) + geom_point(aes(color=mc_group.x),size=1) + 
  scale_color_manual(values=c("#E5D3B3", "#A05124","#12763C","#224BA0","#58B46E")) + ylim(c(0.1,0.9)) + 
  theme_bw() + theme(plot.margin=unit(c(-0.4,1,-0.1,1),"cm"),legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cell_states") + ylab("") + xlab("") 


graph_data_sel <- mat_cell[mat_cell$pat_clone_ID_v2 %in% clone_list,]

plot_d = ggplot(graph_data_sel, aes(x = pseudotime, colour = condition_response, fill = condition_response)) + 
  geom_density(lwd = 0.7, alpha = 0.25) + scale_color_manual(values=c("#910606","#e08282")) + 
  scale_fill_manual(values=c("#910606","#e08282")) + ylim(0,0.4) + xlim(min(mat_cell$pseudotime),max(mat_cell$pseudotime)) +
  theme_bw() + theme(plot.margin=unit(c(+0.4,1,-0.1,1),"cm"),legend.position = "none", panel.grid.major = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),axis.ticks.y = element_blank(),
                     panel.grid.minor = element_blank())

library(gridExtra)
pdf(file=paste0(fig_path,"/CD4_pseudotime_top_clones_NonResponders_update.pdf"), width=3.5, height=2.9)
grid.arrange(plot_d, plot_list$`11205_Pat10`,plot_list$`1326_Pat30`,plot_list$`1456_Pat30`,plot_list$`1917_Pat10`,plot_list$`2420_Pat30`,
             plot_list$`6411_Pat34`,plot_list$`8492_Pat30`,
             plot_b,nrow=9,heights=c(1.7,rep(1,length(clone_list)), 3))
dev.off()



#### Figure S5E ####
## GSEA
# load subset Tregs, from choose_cells
cds_subset <- readRDS(file = "./data/cds_subset_TregGIMAP_TregTNF_v2.rds")
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)

### rnk file on dysfunctional cells trajectory
res_df_sig <- as.data.frame(cbind(name=rownames(pr_graph_test_res),sign= pr_graph_test_res$morans_test_statistic))

res_df_sig <- res_df_sig[order(res_df_sig$sign),]

res_df_sig = res_df_sig[res_df_sig$name!= "",]
res_df_sig = res_df_sig[res_df_sig$sign!= "",]

res_df_sig <- res_df_sig[!is.na(res_df_sig[,1]),]
res_df_sig<-res_df_sig[!is.na(res_df_sig[,2]),]

## ready for GSEA Java tools
write.table(res_df_sig, file="./data/Genes_sig_trajectory_treg_v3.rnk", sep="\t", row.names=F, col.names=F,  quote=F)
## run .rnk on GSEA java application 

## folder generated by GSEA java application, top 5 in heatmap based on NES
GSEA_pos <- read.csv("./data/Genes_sig_trajectory_tregs_GOterms/gsea_report_for_na_pos_1647267851795.tsv",sep="\t")

GSEA_pos_top <- GSEA_pos[order(GSEA_pos$NOM.p.val),][1:20,]
GSEA_pos_other = GSEA_pos[GSEA_pos$NAME %in% c("GOBP_PYRUVATE_METABOLIC_PROCESS"),]

all_GSEA <- rbind(GSEA_pos_other,GSEA_pos_top)
all_GSEA = all_GSEA[order(all_GSEA$FDR.q.val,decreasing = T),]
all_GSEA$NAME = factor(all_GSEA$NAME,levels= all_GSEA$NAME)

pdf(file=paste0(fig_path,"/Trajectory_Treg_moransteststat_barplot.pdf"), width=7, height=6) #
p = ggplot(all_GSEA) + geom_point(aes(y=NAME,x=-log10(FDR.q.val),size=SIZE)) + geom_vline(xintercept=0) +
  scale_size(range=c(2,7)) + geom_hline(yintercept = c(1:nrow(GSEA_pos)), colour="grey") + xlim(c(-3,3)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p
dev.off()


