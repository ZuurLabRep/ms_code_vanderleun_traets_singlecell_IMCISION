###################
#### FIGURE S7 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S7"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S7")

library(AUCell)
library(GSEABase)
library(dplyr)
library(ggpubr)


#### Figure S7A ####
#dysf
overlap_DEG = read.table("./data/overlap_DEG_directional_v3.txt") #v3

mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

df_1 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Dysf GZMB"),])
df_2 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Dysf ZNF683") ,])

df_1_df <- mat@cell_metadata[intersect(df_1, colnames(mat_ds)), ]
df_2_df <- mat@cell_metadata[intersect(df_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(df_1_df,df_2_df)
oi_genes <- overlap_DEG$V1

pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb <- pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Shared_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "Shared_signature"

pb_dysf <- pb

# tregs
mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

df_1 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Treg TNFRSF9"),])

df_2 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Treg GIMAP") ,])

df_1_df <- mat@cell_metadata[intersect(df_1, colnames(mat_ds)), ]
df_2_df <- mat@cell_metadata[intersect(df_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(df_1_df,df_2_df)
oi_genes <- overlap_DEG$V1

pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb <- pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Shared_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "Shared_signature"

pb_treg <- pb

to_plot <- merge(pb_treg,pb_dysf,by="gene")

pdf(file=paste0(fig_path,"/Volcano_DEG_byhand_treg_shared_signature_2way_REV.pdf"), width=2.8, height=2.8)
ggplot() + geom_point(data=to_plot[to_plot$sig.y == "Not Sig" & to_plot$sig.x == "Not Sig",],aes(x=fc.x,y=fc.y),color="grey",size=1.5,alpha=0.5,shape=16,stroke = 0)  +
  geom_point(data=to_plot[to_plot$sig.y == "pval<0.05" & to_plot$sig.x == "pval<0.05",],aes(x=fc.x,y=fc.y),color="black",size=1.5,alpha=0.5,shape=16,stroke = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + geom_point(data=to_plot[to_plot$sig.y == "Shared_signature",],aes(x=fc.x,y=fc.y),color="red",size=1.5,shape=16,stroke = 0) +
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("2-way DEG shared signature") +
  xlab("Treg TNFRSF9 vs GIMAP") +
  ylab("Dysf GZMB vs ZNF683") +
  ylim(c(-4.5,4.5)) + xlim(c(-4.5,4.5)) 
dev.off()



#### Figure S7B-C ####
# shared signature expression per patient
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

overlap_DEG = read.table("./data/overlap_DEG_directional_v3.txt",header=T)
geneSets <- list(geneSet1=overlap_DEG$x)

genesets_col = geneSets

remove = c("Pat04", "Pat12", "Pat38", NA)

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Dysf ZNF683","Dysf GZMB") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("post","pre") &
              !(submat$patient %in% remove),] 
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

sig_scores_sel = sig_scores_all 

## per patient plot
all_merged = sig_scores_all

mean_bm3_list = c()
for (pat_x in unique(all_merged$pat_condition)){
  all_merged_sel = all_merged[all_merged$pat_condition == pat_x,]
  mean_bm3 = mean(all_merged_sel$geneSet1)
  mean_bm3_list = c(mean_bm3_list,mean_bm3)
}
to_plot <- as.data.frame(cbind(mean_bm3_list,unique(all_merged$pat_condition)))
to_plot <- cbind(to_plot,data.frame(do.call(rbind, strsplit(as.vector(to_plot$V2), split = "_"))))
resp_info = mat@cell_metadata[,c("patient","response")]
colnames(to_plot)[3] <- "patient"
resp_info <- resp_info[!duplicated(resp_info),]
to_plot <- merge(to_plot,resp_info,by="patient")

to_plot$cond_resp = paste0(to_plot$response,to_plot$X2)
to_plot$cond_resp = factor(to_plot$cond_resp, levels=c("NRpre","NRpost","REpre","REpost"))

pdf(file=paste0(fig_path,"/Dysf_overlap_siganture_patient_lines.pdf"), width=2.5, height=2.8)
to_plot$X2 = factor(to_plot$X2, levels=c("pre","post"))
ggplot(to_plot, aes(x=X2, y=as.numeric(mean_bm3_list), fill = cond_resp)) +
  geom_boxplot(outlier.shape = NA) + ylim(c(0.10,0.4)) +
  geom_point(aes(group = cond_resp,fill=cond_resp), pch = 21, alpha = 0.5,
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
  geom_line(aes(group = patient), color = "lightgrey", alpha = 0.5)+
  scale_fill_manual(values=c("#dc7986","#921d1f","#71bc88","#0d542c")) +
  theme_bw() + theme(legend.position= "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(~response,scales = "free") + ylab("overlap_sig") + ggtitle("Treg_overlap_sig") 
dev.off()

sig_scores_sel_a = to_plot[to_plot$response == "NR",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)

sig_scores_sel_a = to_plot[to_plot$response == "RE",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)


# Treg
mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Treg GIMAP","Treg TNFRSF9") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("post","pre") &
              !(submat$patient %in% remove),] 
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

sig_scores_sel = sig_scores_all 

## per patient plot
all_merged = sig_scores_all

mean_bm3_list = c()
for (pat_x in unique(all_merged$pat_condition)){
  all_merged_sel = all_merged[all_merged$pat_condition == pat_x,]
  mean_bm3 = mean(all_merged_sel$geneSet1)
  mean_bm3_list = c(mean_bm3_list,mean_bm3)
}
to_plot <- as.data.frame(cbind(mean_bm3_list,unique(all_merged$pat_condition)))
to_plot <- cbind(to_plot,data.frame(do.call(rbind, strsplit(as.vector(to_plot$V2), split = "_"))))
resp_info = mat@cell_metadata[,c("patient","response")]
colnames(to_plot)[3] <- "patient"
resp_info <- resp_info[!duplicated(resp_info),]
to_plot <- merge(to_plot,resp_info,by="patient")

to_plot$cond_resp = paste0(to_plot$response,to_plot$X2)
to_plot$cond_resp = factor(to_plot$cond_resp, levels=c("NRpre","NRpost","REpre","REpost"))

pdf(file=paste0(fig_path,"/Treg_overlap_siganture_patient_lines.pdf"), width=2.5, height=2.8)
to_plot$X2 = factor(to_plot$X2, levels=c("pre","post"))
ggplot(to_plot, aes(x=X2, y=as.numeric(mean_bm3_list), fill = cond_resp)) +
  geom_boxplot(outlier.shape = NA) + ylim(c(0.10,0.4)) +
  geom_point(aes(group = cond_resp,fill=cond_resp), pch = 21, alpha = 0.5,
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
  geom_line(aes(group = patient), color = "lightgrey", alpha = 0.5)+
  scale_fill_manual(values=c("#dc7986","#921d1f","#71bc88","#0d542c")) +
  theme_bw() + theme(legend.position= "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(~response,scales = "free") + ylab("overlap_sig") + ggtitle("Treg_overlap_sig")
dev.off()

sig_scores_sel_a = to_plot[to_plot$response == "NR",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)

sig_scores_sel_a = to_plot[to_plot$response == "RE",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)


# Tfh
mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Tfh LAG3","Tfh NR3C1") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("post","pre") &
              !(submat$patient %in% remove),] 
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

sig_scores_sel = sig_scores_all 

## per patient plot
all_merged = sig_scores_all

mean_bm3_list = c()
for (pat_x in unique(all_merged$pat_condition)){
  all_merged_sel = all_merged[all_merged$pat_condition == pat_x,]
  mean_bm3 = mean(all_merged_sel$geneSet1)
  mean_bm3_list = c(mean_bm3_list,mean_bm3)
}
to_plot <- as.data.frame(cbind(mean_bm3_list,unique(all_merged$pat_condition)))
to_plot <- cbind(to_plot,data.frame(do.call(rbind, strsplit(as.vector(to_plot$V2), split = "_"))))
resp_info = mat@cell_metadata[,c("patient","response")]
colnames(to_plot)[3] <- "patient"
resp_info <- resp_info[!duplicated(resp_info),]
to_plot <- merge(to_plot,resp_info,by="patient")

to_plot$cond_resp = paste0(to_plot$response,to_plot$X2)
to_plot$cond_resp = factor(to_plot$cond_resp, levels=c("NRpre","NRpost","REpre","REpost"))

pdf(file=paste0(fig_path,"/Tfh_overlap_siganture_patient_lines.pdf"), width=2.5, height=2.8)
to_plot$X2 = factor(to_plot$X2, levels=c("pre","post"))
ggplot(to_plot, aes(x=X2, y=as.numeric(mean_bm3_list), fill = cond_resp)) +
  geom_boxplot(outlier.shape = NA) + ylim(c(0.10,0.4)) +
  geom_point(aes(group = cond_resp,fill=cond_resp), pch = 21, alpha = 0.5,
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
  geom_line(aes(group = patient), color = "lightgrey", alpha = 0.5)+
  scale_fill_manual(values=c("#dc7986","#921d1f","#71bc88","#0d542c")) +
  theme_bw() + theme(legend.position= "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(~response,scales = "free") + ylab("overlap_sig") + ggtitle("Treg_overlap_sig")
dev.off()

sig_scores_sel_a = to_plot[to_plot$response == "NR",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)

sig_scores_sel_a = to_plot[to_plot$response == "RE",]
wilcox.test(as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "pre",]$mean_bm3_list),
            as.numeric(sig_scores_sel_a[sig_scores_sel_a$X2 == "post",]$mean_bm3_list),paired = T)


# Tfh, violin, all cells
overlap_DEG = read.table("./data/overlap_DEG_directional_v3.txt")
geneSets <- list(geneSet1=overlap_DEG$V1)

genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Tfh NR3C1","Tfh LAG3") &
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

sig_scores_sel = sig_scores_all[sig_scores_all$mc_group %in% c("Tfh NR3C1","Tfh LAG3"),]

sig_scores_sel$cond_resp = paste0(sig_scores_sel$condition,sig_scores_sel$response)
sig_scores_sel$cond_resp = factor(sig_scores_sel$cond_resp,levels=c("preNR","preRE","postNR","postRE"))
pdf(file=paste0(fig_path,"/Tfh_overlap_sig_AUCell_dir.pdf"), width=3.3, height=2.8)
ggplot(sig_scores_sel, aes(y = geneSet1, x = response,fill=cond_resp))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#DC7985","#89CC9B","#9B1F23","#116936"))+
  theme_bw() + ylab("Shared sig") + ggtitle("shared sig Treg") + ylim(c(0.02,0.55)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "RE",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "RE",][["geneSet1"]])

wilcox.test(sig_scores_sel[sig_scores_sel$condition=="pre" & sig_scores_sel$response == "NR",][["geneSet1"]],
            sig_scores_sel[sig_scores_sel$condition=="post" & sig_scores_sel$response == "NR",][["geneSet1"]])



#### Figure S7D ####
## enrichement scores GO terms and sigantures related to TNF and IFN
geneSets <- getGmt("./data/signatures/Final_GeneSets_To_test_220613.txt")
genesets_col = geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$cell_type %in% c("CD8","CD4","NK","B","M") &
              submat$response %in% c("NR") &
              submat$condition %in% c("post","pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) <- geneNames

dim(exprMatrix)
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)

# concat cell states
md$mc_group_concat <- md$mc_group
md[md$mc_group %in% c("Dysf GZMB","Dysf ZNF683"),]["mc_group_concat"] <- "Dysf"
md[md$mc_group %in% c("Treg GIMAP","Treg TNFRSF9"),]["mc_group_concat"] <- "Treg"
md[md$mc_group %in% c("Tfh LAG3","Tfh NR3C1"),]["mc_group_concat"] <- "Tfh"
md[md$mc_group %in% c("NK FGFBP2","NK KLRC1"),]["mc_group_concat"] <- "NK"
md[md$mc_group %in% c("B cell CD69","B cell CD27" ),]["mc_group_concat"] <- "Bcells"
md[md$mc_group %in% c("cDC CLEC9A","cDC LAMP3","cDC CD1C" ),]["mc_group_concat"] <- "DC"
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_all_pre <- sig_scores_all

df_to_plot <- as.data.frame(levels(as.factor(sig_scores_all_pre$mc_group_concat)))
for (gs_x in colnames(sig_scores_all_pre)[1:17]){
  list_p <- c()
  for (mc_x in levels(as.factor(sig_scores_all_pre$mc_group_concat))){
    temp_df <- sig_scores_all_pre[sig_scores_all_pre$mc_group_concat == mc_x,]
    
    temp_df$condition <- factor(temp_df$condition,levels=c("post","pre"))
    model <- lm(temp_df[[gs_x]] ~ as.factor(temp_df$condition),data=temp_df)
    sub_m <- summary(model)
    
    slope_sub <- model$coefficients[[2]]
    R2 <- sub_m$r.squared
    
    pval <- wilcox.test(temp_df[temp_df$condition == "post",][[gs_x]],temp_df[temp_df$condition == "pre",][[gs_x]])
    
    if (is.na(pval$p.value)){
      pval_s = 1
    }
    if (pval$p.value < 1){
      pval_s = 1
    }
    if (pval$p.value < 0.1){
      pval_s = 0.1
    }
    if (pval$p.value < 0.05){
      pval_s = 0.05
    }
    if (pval$p.value < 0.001){
      pval_s = 0.001
    }
    list_p <- c(list_p,pval_s)
  }
  df_to_plot <- cbind(df_to_plot,list_p)
  colnames(df_to_plot)[ncol(df_to_plot)] <- gs_x
}


df_to_plot_m <- as.data.frame(levels(as.factor(sig_scores_all_pre$mc_group_concat)))
for (gs_x in colnames(sig_scores_all_pre)[1:17]){
  list_p <- c()
  for (mc_x in levels(as.factor(sig_scores_all_pre$mc_group_concat))){
    temp_df <- sig_scores_all_pre[sig_scores_all_pre$mc_group_concat == mc_x,]
    
    temp_df$condition <- factor(temp_df$condition,levels=c("post","pre"))
    model <- lm(temp_df[[gs_x]] ~ as.factor(temp_df$condition),data=temp_df)
    sub_m <- summary(model)
    
    slope_sub <- model$coefficients[[2]]
    R2 <- sub_m$r.squared
    
    list_p <- c(list_p,((-slope_sub)/abs(slope_sub))*R2)
  }
  df_to_plot_m <- cbind(df_to_plot_m,list_p)
  colnames(df_to_plot_m)[ncol(df_to_plot_m)] <- gs_x
}

df_to_plotmelted <- reshape2::melt(df_to_plot,value.name="levels(as.factor(sig_scores_all_pre$mc_group_concat))")
df_to_plot_m_melted <- reshape2::melt(df_to_plot_m)

colnames(df_to_plotmelted) <- c("mc_group_concat","pathway","pval")
colnames(df_to_plot_m_melted) <- c("mc_group_concat","pathway","dir")

to_plot <- merge(df_to_plotmelted,df_to_plot_m_melted)
to_plot$pval <- factor(to_plot$pval, levels= c("1","0.1","0.05","0.001"))
save_to_plot_r <- to_plot

#save_to_plot_r[save_to_plot_r$dir < -0.10,]["dir"] <- -0.10
save_to_plot_r[save_to_plot_r$dir > 0.10,]["dir"] <- 0.10
save_to_plot_r <- save_to_plot_r[!is.na(save_to_plot_r$pathway),]
save_to_plot_r <- save_to_plot_r[!(save_to_plot_r$pathway %in% c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON","SANA_TNF_SIGNALING_UP","TNF_MAARTEN","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN",
                                                             "GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS","GOBP_POSITIVE_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY")),]
save_to_plot_r$pathway = droplevels(save_to_plot_r$pathway)
save_to_plot_r$pathway <- factor(save_to_plot_r$pathway,levels=c("GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA","GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA",
                                                                "GOBP_RESPONSE_TO_INTERFERON_ALPHA","GOBP_RESPONSE_TO_INTERFERON_BETA",
                                                                "GOBP_POSITIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION","GOBP_RESPONSE_TO_INTERFERON_GAMMA",
                                                                "IFN_AYERS","IFN_ANNE","PID_TNF_PATHWAY","REACTOME_TNF_SIGNALING"))

pdf(file=paste0(fig_path,"/Go_terms_final_signatures_incl_NR.pdf"), width=9.5, height=3.9)
ggplot(save_to_plot_r,aes(x=mc_group_concat,y=pathway,color=dir,size=pval)) + geom_point(aes(fill=dir),color="grey",shape=21,stroke = 0.4) + 
  scale_fill_gradient2(midpoint=0, low="#6c2fad", mid="white",high="#eda600", space ="Lab" ,
                       limits=c(-0.10,0.10)) + scale_size_manual(values=c(3,3.5,4.5,6.5,8)) + ggtitle("logfc_postvspre_NR") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## RE
geneSets <- getGmt("./data/signatures/Final_GeneSets_To_test_220613.txt")
genesets_col <- geneSets

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$cell_type %in% c("CD8","CD4","NK", "M", "B") &
              submat$response %in% c("RE") &
              submat$condition %in% c("post","pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) <- geneNames

cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(genesets_col,cell_rank)

check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)

# concat cell states
md$mc_group_concat <- md$mc_group
md[md$mc_group %in% c("Dysf GZMB","Dysf ZNF683"),]["mc_group_concat"] <- "Dysf"
md[md$mc_group %in% c("Treg GIMAP","Treg TNFRSF9"),]["mc_group_concat"] <- "Treg"
md[md$mc_group %in% c("Tfh LAG3","Tfh NR3C1"),]["mc_group_concat"] <- "Tfh"
md[md$mc_group %in% c("NK FGFBP2","NK KLRC1"),]["mc_group_concat"] <- "NK"
md[md$mc_group %in% c("B cell CD69","B cell CD27" ),]["mc_group_concat"] <- "Bcells"
md[md$mc_group %in% c("cDC CLEC9A","cDC LAMP3","cDC CD1C" ),]["mc_group_concat"] <- "DC"
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_all_pre <- sig_scores_all

df_to_plot <- as.data.frame(levels(as.factor(sig_scores_all_pre$mc_group_concat)))
for (gs_x in colnames(sig_scores_all_pre)[1:16]){
  list_p <- c()
  for (mc_x in levels(as.factor(sig_scores_all_pre$mc_group_concat))){
    temp_df <- sig_scores_all_pre[sig_scores_all_pre$mc_group_concat == mc_x,]
    
    temp_df$condition <- factor(temp_df$condition,levels=c("post","pre"))
    model <- lm(temp_df[[gs_x]] ~ as.factor(temp_df$condition),data=temp_df)
    sub_m <- summary(model)
    
    slope_sub <- model$coefficients[[2]]
    R2 <- sub_m$r.squared
    
    pval <- wilcox.test(temp_df[temp_df$condition == "post",][[gs_x]],temp_df[temp_df$condition == "pre",][[gs_x]])
    
    if (is.na(pval$p.value)){
      pval_s = 1
    }
    if (pval$p.value < 1){
      pval_s = 1
    }
    if (pval$p.value < 0.1){
      pval_s = 0.1
    }
    if (pval$p.value < 0.05){
      pval_s = 0.05
    }
    if (pval$p.value < 0.001){
      pval_s = 0.001
    }
    list_p <- c(list_p,pval_s)
  }
  df_to_plot <- cbind(df_to_plot,list_p)
  colnames(df_to_plot)[ncol(df_to_plot)] <- gs_x
}


df_to_plot_m <- as.data.frame(levels(as.factor(sig_scores_all_pre$mc_group_concat)))
for (gs_x in colnames(sig_scores_all_pre)[1:16]){
  list_p <- c()
  for (mc_x in levels(as.factor(sig_scores_all_pre$mc_group_concat))){
    temp_df <- sig_scores_all_pre[sig_scores_all_pre$mc_group_concat == mc_x,]
    
    temp_df$condition <- factor(temp_df$condition,levels=c("post","pre"))
    model <- lm(temp_df[[gs_x]] ~ as.factor(temp_df$condition),data=temp_df)
    sub_m <- summary(model)
    
    slope_sub <- model$coefficients[[2]]
    R2 <- sub_m$r.squared
    
    list_p <- c(list_p,((-slope_sub)/abs(slope_sub))*R2)
  }
  df_to_plot_m <- cbind(df_to_plot_m,list_p)
  colnames(df_to_plot_m)[ncol(df_to_plot_m)] <- gs_x
}

df_to_plotmelted <- reshape2::melt(df_to_plot,value.name="levels(as.factor(sig_scores_all_pre$mc_group_concat))")
df_to_plot_m_melted <- reshape2::melt(df_to_plot_m)

colnames(df_to_plotmelted) <- c("mc_group_concat","pathway","pval")
colnames(df_to_plot_m_melted) <- c("mc_group_concat","pathway","dir")

to_plot <- merge(df_to_plotmelted,df_to_plot_m_melted)
to_plot$pval <- factor(to_plot$pval, levels= c("1","0.1","0.05","0.001"))
save_to_plot_r <- to_plot

save_to_plot_r[save_to_plot_r$dir < -0.10,]["dir"] <- -0.10
#save_to_plot_r[save_to_plot_r$dir > 0.10,]["dir"] <- 0.10
save_to_plot_r <- save_to_plot_r[!(save_to_plot_r$pathway %in% c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON","SANA_TNF_SIGNALING_UP","TNF_MAARTEN","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN",
                                                                 "GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS","GOBP_POSITIVE_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY")),]
save_to_plot_r$pathway = droplevels(save_to_plot_r$pathway)
save_to_plot_r$pathway <- factor(save_to_plot_r$pathway,levels=c("GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA","GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA",
                                                                 "GOBP_RESPONSE_TO_INTERFERON_ALPHA","GOBP_RESPONSE_TO_INTERFERON_BETA",
                                                                 "GOBP_POSITIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION","GOBP_RESPONSE_TO_INTERFERON_GAMMA",
                                                                 "IFN_AYERS","IFN_ANNE","PID_TNF_PATHWAY","REACTOME_TNF_SIGNALING"))

pdf(file=paste0(fig_path,"/Go_terms_final_signatures_incl_RE.pdf"), width=9.5, height=3.9)
ggplot(save_to_plot_r,aes(x=mc_group_concat,y=pathway,color=dir,size=pval)) + geom_point(aes(fill=dir),color="grey",shape=21,stroke = 0.4) + 
  scale_fill_gradient2(midpoint=0, low="#6c2fad", mid="white",high="#eda600", space ="Lab" ,
                       limits=c(-0.10,0.10)) + scale_size_manual(values=c(3,3.5,4.5,6.5,8)) + ggtitle("logfc_postvspre_RE") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



#### Figure S7E ####
# waterfall plots, shared siganture on subsets
overlap_DEG = read.table("./data/overlap_DEG_directional_v3.txt",header=T)
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                mat@cell_metadata$response %in% c("RE") &
                                                                mat@cell_metadata$condition %in% c("post") &
                                                                mat@cell_metadata$mc_group %in% c("Treg TNFRSF9","Treg GIMAP")], 
                 nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                    mat@cell_metadata$response  %in% c("RE") &
                                                    mat@cell_metadata$condition %in% c("pre") &
                                                    mat@cell_metadata$mc_group %in% c("Treg TNFRSF9","Treg GIMAP")], 
                 geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),]
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

#de_2[de_2$enr > 1.5,]$enr <- 1.5
de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Treg_overlap_sig_RE_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Treg_overlap_postvspre_RE")
dev.off()


## dysf
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                mat@cell_metadata$response %in% c("RE") &
                                                                mat@cell_metadata$condition %in% c("post") &
                                                                mat@cell_metadata$mc_group %in% c("Dysf GZMB","Dysf ZNF683")], 
                 nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                    mat@cell_metadata$response  %in% c("RE") &
                                                    mat@cell_metadata$condition %in% c("pre") &
                                                    mat@cell_metadata$mc_group %in% c("Dysf GZMB","Dysf ZNF683")], 
                 geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),] ## CHANGE
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

#de_2[de_2$enr > 1.5,]$enr <- 1.5
#de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Dysf_overlap_sig_RE_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Dysf_overlap_postvspre_RE")
dev.off()


## tfh
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                mat@cell_metadata$response %in% c("RE") &
                                                                mat@cell_metadata$condition %in% c("post") &
                                                                mat@cell_metadata$mc_group %in% c("Tfh LAG3","Tfh NR3C1")], 
                 nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                    mat@cell_metadata$response  %in% c("RE") &
                                                    mat@cell_metadata$condition %in% c("pre") &
                                                    mat@cell_metadata$mc_group %in% c("Tfh LAG3","Tfh NR3C1")], 
                 geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),] ## CHANGE
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

de_2[de_2$enr > 1.5,]$enr <- 1.5
de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Tfh_overlap_sig_RE_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Tfh_overlap_postvspre_RE")
dev.off()


# non-responders
# treg
overlap_DEG = read.table("./data/overlap_DEG_directional_v3.txt",header=T)
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                 mat@cell_metadata$response %in% c("NR") &
                                                                 mat@cell_metadata$condition %in% c("post") &
                                                                 mat@cell_metadata$mc_group %in% c("Treg TNFRSF9","Treg GIMAP")], 
                  nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                     mat@cell_metadata$response  %in% c("NR") &
                                                     mat@cell_metadata$condition %in% c("pre") &
                                                     mat@cell_metadata$mc_group %in% c("Treg TNFRSF9","Treg GIMAP")], 
                  geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),]
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

#de_2[de_2$enr > 1.5,]$enr <- 1.5
#de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Treg_overlap_sig_NR_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Treg_overlap_postvspre_NR")
dev.off()


## dysf
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                 mat@cell_metadata$response %in% c("NR") &
                                                                 mat@cell_metadata$condition %in% c("post") &
                                                                 mat@cell_metadata$mc_group %in% c("Dysf GZMB","Dysf ZNF683")], 
                  nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                     mat@cell_metadata$response  %in% c("NR") &
                                                     mat@cell_metadata$condition %in% c("pre") &
                                                     mat@cell_metadata$mc_group %in% c("Dysf GZMB","Dysf ZNF683")], 
                  geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),] ## CHANGE
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

#de_2[de_2$enr > 1.5,]$enr <- 1.5
de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Dysf_overlap_sig_NR_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Dysf_overlap_postvspre_NR")
dev.off()


## tfh
mat_ds = scdb_mat(paste0(mat_id, "_ds"))

conditions <- unique(mat@cell_metadata$patient)
remove <- c("Pat04", "Pat12", "Pat38", NA)
conditions <- conditions[!conditions %in% remove]

de_1 <- diff_expr(mc, mat_ds, nms1=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                                 mat@cell_metadata$response %in% c("NR") &
                                                                 mat@cell_metadata$condition %in% c("post") &
                                                                 mat@cell_metadata$mc_group %in% c("Tfh LAG3","Tfh NR3C1")], 
                  nms2=rownames(mat@cell_metadata)[mat@cell_metadata$patient %in% conditions &
                                                     mat@cell_metadata$response  %in% c("NR") &
                                                     mat@cell_metadata$condition %in% c("pre") &
                                                     mat@cell_metadata$mc_group %in% c("Tfh LAG3","Tfh NR3C1")], 
                  geo_mean = T, min_max_umi = 10)

de_2 <- de_1[unique(overlap_DEG$x),] ## CHANGE
de_2 <- de_2[order(de_2$enr),]
de_2 <- de_2[,c("gene","enr")]
de_2 <- de_2[complete.cases(de_2),]

de_2$dir <- rep("up",nrow(de_2))
de_2[de_2$enr < 0,]["dir"] <- "down"
de_2 <- de_2[order(de_2$enr, decreasing=F),]
de_2$gene <- factor(de_2$gene,levels=de_2$gene)

#de_2[de_2$enr > 1.5,]$enr <- 1.5
#de_2[de_2$enr < -1.5,]$enr <- -1.5

pdf(file=paste0(fig_path,"/waterfall_Tfh_overlap_sig_NR_postvspre.pdf"), width=2.5, height=11.5)
ggplot(de_2,aes(x=gene,y=enr,fill=dir)) + geom_bar(stat="identity") + ylim(values=c(-1.5,1.5))+
  coord_flip() + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values=c("#237BB7","#F15F41")) + ggtitle("diff_exp_Tfh_overlap_postvspre_NR")
dev.off()



#### Figure S7F ####
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
mat_df = as.matrix(mat@mat)
umis = mat_df[,names(mc@mc)]
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 
submat = mat@cell_metadata

# NR
dt_cells_save <- as.data.frame(cbind( c(rep("pre",6),rep("post",6)),
                                     c("Pat10","Pat26","Pat30","Pat33","Pat34","Pat37",
                                       "Pat10","Pat26","Pat30","Pat33","Pat34","Pat37")))
colnames(dt_cells_save) <- c("condition","patient")

# bubble plot
for (cells_x in list(c("Dysf GZMB","Dysf ZNF683"),c("Treg TNFRSF9","Treg GIMAP"),
                     c("Tfh NR3C1","Tfh LAG3"))){
  remove = c("Pat04", "Pat12", "Pat38") 
  md = submat
  md = md[!md$patient %in% remove &
            md$response %in% "NR" &
            md$mc_group %in% cells_x,]
  
  dt_cells_number = md %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  cells_x_name = gsub(" ","_",cells_x)
  cells_x_name = paste(cells_x_name,collapse = ":")
  colnames(dt_cells_number)[3] = paste0("count_",cells_x_name)
  
  for (gene_x in c("LAG3","CTLA4","PDCD1","HAVCR2","TIGIT")){
    # per gene
    gene_oi = gene_x
    genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
    genes_cells <- merge(md,genes_cells,by=0)
    colnames(genes_cells)[ncol(genes_cells)] = "gene_oi"
    
    dt_genes = genes_cells %>% 
      group_by(condition,patient) %>% 
      summarize(average=mean(gene_oi))
    
    colnames(dt_genes)[3] = paste0(gene_x,":",cells_x_name)
    dt_cells_number = cbind(dt_cells_number,dt_genes[paste0(gene_x,":",cells_x_name)])
  }
  
  dt_cells_save = cbind(dt_cells_save,dt_cells_number)
}

dysfGZMB_ZNF <- colMeans(log2(dt_cells_save[,c(6:10)][7:12,]/dt_cells_save[,c(6:10)][1:6,]))
tregTNF_GIMAP <- colMeans(log2(dt_cells_save[,c(14:18)][7:12,]/dt_cells_save[,c(14:18)][1:6,]))
tfhLAG_NR3C1 <- colMeans(log2(dt_cells_save[,c(22:26)][7:12,]/dt_cells_save[,c(22:26)][1:6,]))

pval_x = list()
for (x in 6:10){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][7:12],dt_cells_save[,c(x)][1:6],paired=T)$p.value)
}
dysf_pval <- unlist(pval_x)

pval_x = list()
for (x in 14:18){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][7:12],dt_cells_save[,c(x)][1:6],paired=T)$p.value)
}
treg_pval <- unlist(pval_x)

pval_x = list()
for (x in 22:26){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][7:12],dt_cells_save[,c(x)][1:6],paired=T)$p.value)
}
tfh_pval <- unlist(pval_x)

to_plot <- cbind(dysfGZMB_ZNF,tregTNF_GIMAP,tfhLAG_NR3C1)
rownames(to_plot) <- do.call(rbind,strsplit(rownames(to_plot),":"))[,1]
to_plot_p <- cbind(dysf_pval,treg_pval,tfh_pval)
rownames(to_plot_p) <- rownames(to_plot)

to_plot <- melt(to_plot)
to_plot_p <- melt(to_plot_p)
to_plot$pvalue <- to_plot_p$value

to_plot$pvalue_cat <- ""
for (x in 1:nrow(to_plot)){
  pval = to_plot$pvalue[x]
  if (is.na(pval)){
    pval_s = 1
  }
  if (pval < 1){
    pval_s = 1
  }
  if (pval < 0.1){
    pval_s = 0.1
  }
  if (pval < 0.05){
    pval_s = 0.05
  }
  if (pval < 0.001){
    pval_s = 0.001
  }
  to_plot$pvalue_cat[x] = pval_s
}

to_plot$pvalue_cat <- factor(to_plot$pvalue_cat,levels=c("1","0.1","0.05","0.001"))
pdf(file=paste0(fig_path,"/Lfc_IR_post_pre_NR.pdf"),width=4, height=2.5)
ggplot(to_plot,aes(x=Var1,y=Var2,color=value,size = pvalue_cat)) + geom_point(aes(fill=value),color="grey",shape=21,stroke = 0.4)+
  scale_fill_gradient2(midpoint=0, low="#6c2fad", mid="white",high="#eda600", space ="Lab",limits=c(-1,1)) +
  scale_size_manual(values=c(3.5,4.5,6.5,8)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("") + xlab("") + ggtitle("lfc IR - post vs pre - NR")
dev.off()


# RE
dt_cells_save <- as.data.frame(cbind( c(rep("pre",9),rep("post",9)),
                                     c("Pat15","Pat17","Pat21","Pat22","Pat27","Pat29","Pat31","Pat36","Pat39",
                                       "Pat15","Pat17","Pat21","Pat22","Pat27","Pat29","Pat31","Pat36","Pat39")))
colnames(dt_cells_save) <- c("condition","patient")

# bubble plot
library(dplyr)
for (cells_x in list(c("Dysf GZMB","Dysf ZNF683"),c("Treg TNFRSF9","Treg GIMAP"),
                     c("Tfh NR3C1","Tfh LAG3"))){
  remove = c("Pat04", "Pat12", "Pat38") 
  md = submat
  md = md[!md$patient %in% remove &
            md$response %in% "RE" &
            md$mc_group %in% cells_x,]
  
  dt_cells_number = md %>% group_by(condition,patient) %>% summarise(count=n())
  cells_x_name = gsub(" ","_",cells_x)
  cells_x_name = paste(cells_x_name,collapse = ":")
  colnames(dt_cells_number)[3] = paste0("count_",cells_x_name)
  
  for (gene_x in c("LAG3","CTLA4","PDCD1","HAVCR2","TIGIT")){
    # per gene
    gene_oi = gene_x
    genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
    genes_cells <- merge(md,genes_cells,by=0)
    colnames(genes_cells)[ncol(genes_cells)] = "gene_oi"
    
    dt_genes = genes_cells %>% 
      group_by(condition,patient) %>% 
      summarize(average=mean(gene_oi))
    
    colnames(dt_genes)[3] = paste0(gene_x,":",cells_x_name)
    dt_cells_number = cbind(dt_cells_number,dt_genes[paste0(gene_x,":",cells_x_name)])
  }
  
  dt_cells_save = cbind(dt_cells_save,dt_cells_number)
}

dysfGZMB_ZNF <- colMeans(log2(dt_cells_save[,c(6:10)][10:18,]/dt_cells_save[,c(6:10)][1:9,]))
tregTNF_GIMAP <- log2(dt_cells_save[,c(14:18)][10:18,]/dt_cells_save[,c(14:18)][1:9,])
tregTNF_GIMAP[tregTNF_GIMAP == "Inf"] <- NA
tregTNF_GIMAP <- colMeans(tregTNF_GIMAP,na.rm = T)
tfhLAG_NR3C1 <- log2(dt_cells_save[,c(22:26)][10:18,]/dt_cells_save[,c(22:26)][1:9,])
tfhLAG_NR3C1[tfhLAG_NR3C1 == "Inf"] <- NA
tfhLAG_NR3C1 <- colMeans(tfhLAG_NR3C1,na.rm = T)

pval_x = list()
for (x in 6:10){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][10:18],dt_cells_save[,c(x)][1:9],paired=T)$p.value)
}
dysf_pval <- unlist(pval_x)

pval_x = list()
for (x in 14:18){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][10:18],dt_cells_save[,c(x)][1:9],paired=T)$p.value)
}
treg_pval <- unlist(pval_x)

pval_x = list()
for (x in 22:26){
  pval_x = c(pval_x, wilcox.test(dt_cells_save[,c(x)][10:18],dt_cells_save[,c(x)][1:9],paired=T)$p.value)
}
tfh_pval <- unlist(pval_x)

to_plot <- cbind(dysfGZMB_ZNF,tregTNF_GIMAP,tfhLAG_NR3C1)
to_plot[is.infinite(to_plot)] <- 0
rownames(to_plot) <- do.call(rbind,strsplit(rownames(to_plot),":"))[,1]
to_plot_p <- cbind(dysf_pval,treg_pval,tfh_pval)
rownames(to_plot_p) <- rownames(to_plot)

to_plot <- melt(to_plot)
to_plot_p <- melt(to_plot_p)
to_plot$pvalue <- to_plot_p$value

to_plot$pvalue_cat <- ""
for (x in 1:nrow(to_plot)){
  pval = to_plot$pvalue[x]
  if (is.na(pval)){
    pval_s = 1
  }
  if (pval < 1){
    pval_s = 1
  }
  if (pval < 0.1){
    pval_s = 0.1
  }
  if (pval < 0.05){
    pval_s = 0.05
  }
  if (pval < 0.001){
    pval_s = 0.001
  }
  to_plot$pvalue_cat[x] = pval_s
}

to_plot$pvalue_cat <- factor(to_plot$pvalue_cat,levels=c("1","0.1","0.05","0.001"))
pdf(file=paste0(fig_path,"/Lfc_IR_post_pre_RE.pdf"),width=4, height=2.5)
ggplot(to_plot,aes(x=Var1,y=Var2,color=value,size = pvalue_cat)) + geom_point(aes(fill=value),color="grey",shape=21,stroke = 0.4)+
  scale_fill_gradient2(midpoint=0, low="#6c2fad", mid="white",high="#eda600", space ="Lab",limits=c(-1,1)) +
  scale_size_manual(values=c(3.5,4.5,6.5,8)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("") + xlab("") + ggtitle("lfc IR - post vs pre - RE")
dev.off()


#### Figure S7G ####
for (cells_x in c("Tregs","dysf","Tfh")){
  print(cells_x)
  if (cells_x == "Tregs"){
    cells_mc_x <- c("Treg GIMAP","Treg TNFRSF9")
  }
  if (cells_x == "dysf"){
    cells_mc_x <- c("Dysf GZMB","Dysf ZNF683")
  }
  if (cells_x == "Tfh"){
    cells_mc_x <- c("Tfh LAG3","Tfh NR3C1")
  }
  remove = c("Pat04", "Pat12", "Pat38") 
  md = mat@cell_metadata
  md = md[!md$patient %in% remove &
            md$response %in% "NR" &
            md$mc_group %in% cells_mc_x,]
  
  # per gene
  gene_oi <- c("LAG3")
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells <- merge(md,genes_cells,by=0)
  colnames(genes_cells)[ncol(genes_cells)] <- gene_oi
  
  dt_cells_nr <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  dt_genes <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarize(average=mean(LAG3))
  dt_genes$cell_nr <- dt_cells_nr$count
  
  pdf(file=paste0(fig_path,"/",cells_x,"_LAG3_NR_perpatient_fullmat.pdf"),width=3, height=4)
  p1=ggplot(dt_genes,aes(x=condition,y=average)) + geom_point(aes(x=condition,y=average,color=patient,size=cell_nr)) +
    geom_line(aes(x=condition,y=average,group=patient,color=patient)) +   theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle("LAG3") +
    ggpubr::stat_compare_means(method="wilcox",paired=T,label = "p.format") +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320),
               labels = c(10,20,40,80,160,320), range=c(1,3))
  print(p1)
  dev.off()
  
  gene_oi <- c("CTLA4")
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells <- merge(md,genes_cells,by=0)
  colnames(genes_cells)[ncol(genes_cells)] <- gene_oi
  
  dt_cells_nr <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  dt_genes <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarize(average=mean(CTLA4))
  dt_genes$cell_nr <- dt_cells_nr$count
  
  pdf(file=paste0(fig_path,"/",cells_x,"_CTLA4_NR_perpatient_fullmat.pdf"),width=3, height=4)
  p2=ggplot(dt_genes,aes(x=condition,y=average)) + geom_point(aes(x=condition,y=average,color=patient,size=cell_nr)) +
    geom_line(aes(x=condition,y=average,group=patient,color=patient)) +   theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    ggtitle("CTLA4")+
    ggpubr::stat_compare_means(method="wilcox",paired=T,label = "p.format")+
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320),
               labels = c(10,20,40,80,160,320), range=c(1,3))
  print(p2)
  dev.off()
  
  gene_oi <- c("PDCD1")
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells <- merge(md,genes_cells,by=0)
  colnames(genes_cells)[ncol(genes_cells)] <- gene_oi
  
  dt_cells_nr <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  dt_genes <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarize(average=mean(PDCD1))
  dt_genes$cell_nr <- dt_cells_nr$count
  
  pdf(file=paste0(fig_path,"/",cells_x,"_PDCD1_NR_perpatient_fullmat.pdf"),width=3, height=4)
  p3=ggplot(dt_genes,aes(x=condition,y=average)) + geom_point(aes(x=condition,y=average,color=patient,size=cell_nr)) +
    geom_line(aes(x=condition,y=average,group=patient,color=patient)) +   theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    ggtitle("PDCD1")+
    ggpubr::stat_compare_means(method="wilcox",paired=T,label = "p.format") +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320),
               labels = c(10,20,40,80,160,320), range=c(1,3))
  print(p3)
  dev.off()
  
  gene_oi <- c("HAVCR2")
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells <- merge(md,genes_cells,by=0)
  colnames(genes_cells)[ncol(genes_cells)] <- gene_oi
  
  dt_cells_nr <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  dt_genes <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarize(average=mean(HAVCR2))
  dt_genes$cell_nr <- dt_cells_nr$count
  
  pdf(file=paste0(fig_path,"/",cells_x,"_HAVCR2_NR_perpatient_fullmat.pdf"),width=3, height=4)
  p4=ggplot(dt_genes,aes(x=condition,y=average)) + geom_point(aes(x=condition,y=average,color=patient,size=cell_nr)) +
    geom_line(aes(x=condition,y=average,group=patient,color=patient)) +   theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    ggtitle("HAVCR2")+
    ggpubr::stat_compare_means(method="wilcox",paired=T,label = "p.format")+
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320),
               labels = c(10,20,40,80,160,320), range=c(1,3))
  print(p4)
  dev.off()
  
  gene_oi <- c("TIGIT")
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  genes_cells <- merge(md,genes_cells,by=0)
  colnames(genes_cells)[ncol(genes_cells)] <- gene_oi
  
  dt_cells_nr <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarise(count=n())
  dt_genes <- genes_cells %>% 
    group_by(condition,patient) %>% 
    summarize(average=mean(TIGIT))
  dt_genes$cell_nr <- dt_cells_nr$count
  
  pdf(file=paste0(fig_path,"/",cells_x,"_TIGIT_NR_perpatient_fullmat.pdf"),width=3, height=4)
  p5=ggplot(dt_genes,aes(x=condition,y=average)) + geom_point(aes(x=condition,y=average,color=patient,size=cell_nr)) +
    geom_line(aes(x=condition,y=average,group=patient,color=patient)) +   theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    ggtitle("TIGIT")+
    ggpubr::stat_compare_means(method="wilcox",paired=T,label = "p.format")+
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320),
               labels = c(10,20,40,80,160,320), range=c(1,3))
  print(p5)
  dev.off()
  
  
  pdf(file=paste0(fig_path,"/",cells_x,"_IR_all_NR_perpatient_fullmat.pdf"),width=8, height=3,onefile=F)
  figure <- ggpubr::ggarrange(p1,p2,p3,p4,p5,
                      ncol = 5, nrow = 1,
                      common.legend = TRUE, legend = "right")
  annotate_figure(figure, top = text_grob(paste0(cells_x," - NR"), 
                                          color = "black", size = 14))
  print(figure)
  dev.off()
  
}




