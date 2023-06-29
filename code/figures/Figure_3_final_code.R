##################
#### FIGURE 3 ####
##################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures/", "Figure_3"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_3")

library(AUCell)
library(GSEABase)
library(ComplexHeatmap)
library(Seurat)
library(Startrac)


#### Figure 3A ####
# fold change pre > post, RE and NR
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

# define pre and post in RE and NR, T cells
pre = md[md$condition == "pre" & 
           md$cell_type %in% c("CD4","CD8","NK") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12" &
           !is.na(md$mc_group),]

post = md[md$condition == "post" & 
            md$cell_type %in% c("CD4","CD8","NK") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12" &
            !is.na(md$mc_group),]

# calculate fractions, pre and post
pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

# calculate fold change between pre and post
states <- colnames(pre_final)
samples <- rownames(pre_final)
fold_change <- post_final[samples, states]/pre_final[samples, states]
fold_change <- as.matrix(fold_change)

# calculate paired wilcoxon between pre and post per cell state 
select_md <- md[,c("patient","response")]
pval_states = as.data.frame(states)
pval_states$pvalue_NR <- rep("NA",nrow(pval_states))
pval_states$pvalue_RE <- rep("NA",nrow(pval_states))
for (state in states){
  pre_final_NR <- pre_final[rownames(pre_final) %in% select_md[select_md$response=="NR",]$patient,]
  post_final_NR <- post_final[rownames(post_final) %in% select_md[select_md$response=="NR",]$patient,]
  pval_s <- wilcox.test(post_final_NR[,state],pre_final_NR[,state], paired = TRUE)$p.value

  pval_states[pval_states$states == state,]$pvalue_NR <- pval_s
  pre_final_RE <- pre_final[rownames(pre_final) %in% select_md[select_md$response=="RE",]$patient,]
  post_final_RE <- post_final[rownames(post_final) %in% select_md[select_md$response=="RE",]$patient,]
  pval_s <- wilcox.test(post_final_RE[,state],pre_final_RE[,state], paired = TRUE)$p.value
  pval_states[pval_states$states == state,]$pvalue_RE <- pval_s
}
save_pvals_t_cells <- pval_states

# log2(fold change)
lfc_df = log2(fold_change)
lfc_df <- replace(lfc_df, is.infinite(lfc_df),NA)
lfc_df[is.na(lfc_df)] <- 0

lfc_df <- as.data.frame(lfc_df)
select_md <- md[,c("patient","response")]
lfc_df$patient <- rownames(lfc_df)
select_md <- unique(select_md)
lfc_df_all <- merge(lfc_df,select_md,by="patient")

lfc_df_all <- lfc_df_all[order(lfc_df_all$response, decreasing = T),]
rownames(lfc_df_all) <- lfc_df_all$patient

cells_order <- c( "Cytotoxic"   , "Dysf GZMB"  ,  "Dysf ZNF683" ,"Naive-like CD4","Naive-like CD8" ,
                  "Tfh LAG3"   ,  "Tfh NR3C1"   ,"Transitional",   "Treg GIMAP"  , "Treg TNFRSF9","Non-classical T cell","NK KLRC1","NK FGFBP2")

lfc_df_all <- lfc_df_all[,c("patient",cells_order,"response")  ]

## median fold change, and wilcox test
lfc_df_all_RE <- lfc_df_all[lfc_df_all$response == "RE",]
lfc_df_all_NR <- lfc_df_all[lfc_df_all$response == "NR",]

median_fc_RE <- robustbase::colMedians(as.matrix(lfc_df_all_RE[,2:(ncol(lfc_df_all_RE)-1)]))
median_fc_NR <- robustbase::colMedians(as.matrix(lfc_df_all_NR[,2:(ncol(lfc_df_all_NR)-1)]))

median_fc <- rbind(median_fc_RE,median_fc_NR)
median_fc <- as.data.frame(t(median_fc))

median_fc <- median_fc[match(cells_order,rownames(median_fc)),]
save_t_cells_medians <- median_fc # save T cell analysis


# M/B, non-T cell
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

# define pre and post in RE and NR, non-T cells, M and B
pre = md[md$condition == "pre" & 
           md$cell_type %in% c("M","B") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12" &
           !is.na(md$mc_group),]

post = md[md$condition == "post" & 
            md$cell_type %in% c("M","B") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12" &
            !is.na(md$mc_group),]

# calculate fractions, pre and post
pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

# calculate difference between pre and post
states <- colnames(pre_final)
samples <- rownames(pre_final)
fold_change <- post_final[samples, states]/pre_final[samples, states]
fold_change <- as.matrix(fold_change)

# calculate paired wilcoxon between pre and post per cell state 
select_md <- md[,c("patient","response")]
pval_states = as.data.frame(states)
pval_states$pvalue_NR = rep("NA",nrow(pval_states))
pval_states$pvalue_RE = rep("NA",nrow(pval_states))
for (state in states){
  pre_final_NR <- pre_final[rownames(pre_final) %in% select_md[select_md$response=="NR",]$patient,]
  post_final_NR <- post_final[rownames(post_final) %in% select_md[select_md$response=="NR",]$patient,]
  pval_s <- wilcox.test(post_final_NR[,state],pre_final_NR[,state], paired = TRUE)$p.value
  pval_states[pval_states$states == state,]$pvalue_NR <- pval_s
  
  pre_final_RE <- pre_final[rownames(pre_final) %in% select_md[select_md$response=="RE",]$patient,]
  post_final_RE <- post_final[rownames(post_final) %in% select_md[select_md$response=="RE",]$patient,]
  pval_s <- wilcox.test(post_final_RE[,state],pre_final_RE[,state], paired = TRUE)$p.value
  pval_states[pval_states$states == state,]$pvalue_RE <- pval_s
}

save_pvals_non_t_cells = pval_states 

# log2(fold change)
lfc_df <- log2(fold_change)
lfc_df <- replace(lfc_df, is.infinite(lfc_df),NA)
lfc_df[is.na(lfc_df)] <- 0

lfc_df <- as.data.frame(lfc_df)
select_md <- md[,c("patient","response")]
lfc_df$patient <- rownames(lfc_df)
select_md <- unique(select_md)
lfc_df_all <- merge(lfc_df,select_md,by="patient")

lfc_df_all <- lfc_df_all[order(lfc_df_all$response, decreasing = T),]
rownames(lfc_df_all) <- lfc_df_all$patient

cells_order <- c("B cell CD27", "B cell CD69",   "Plasma cell" , "cDC CD1C"    , "cDC CLEC9A"  , "cDC LAMP3"   ,
                 "Granulocyte",  "Mono-macro"  , "pDC" )

lfc_df_all <- lfc_df_all[,c("patient"    , cells_order   , "response" )  ]

## median fold change
lfc_df_all_RE <- lfc_df_all[lfc_df_all$response == "RE",]
lfc_df_all_NR <- lfc_df_all[lfc_df_all$response == "NR",]

median_fc_RE <- colMeans(as.matrix(lfc_df_all_RE[,2:(ncol(lfc_df_all_RE)-1)]))
median_fc_NR <- colMeans(as.matrix(lfc_df_all_NR[,2:(ncol(lfc_df_all_NR)-1)]))

median_fc <- rbind(median_fc_RE,median_fc_NR)
median_fc <- as.data.frame(t(median_fc))
median_fc <- median_fc[match(cells_order,rownames(median_fc)),]

save_non_t_cells_medians = median_fc # save non T cell analysis

# combina all data, plot median log2(fc)
medians_cells = rbind(save_t_cells_medians,save_non_t_cells_medians)
pvals_cells = rbind(save_pvals_t_cells,save_pvals_non_t_cells)

medians_cells$states = rownames(medians_cells)
prep_table <- merge(medians_cells,pvals_cells)

## RE plots
prep_table$sig_p <- rep("p>0.5",nrow(prep_table))
prep_table[prep_table$pvalue_RE < 0.5,][["sig_p"]] <- "p<0.5"
prep_table[prep_table$pvalue_RE < 0.1,][["sig_p"]] <- "p<0.1"
prep_table[prep_table$pvalue_RE < 0.05,][["sig_p"]] <- "p<0.05"
prep_table[prep_table$pvalue_RE < 0.01,][["sig_p"]] <- "p<0.01"

prep_table$sig_p <- factor(prep_table$sig_p,levels=c("p>0.5","p<0.5","p<0.1","p<0.05","p<0.01"))

colors_mc <- read.table("./data/mc_colors.txt",header=T)
prep_table$V1 <- gsub(" ", "_", prep_table$states)

order_cells <- c( "Naive-like_CD4",  "Treg_GIMAP"   ,  "Treg_TNFRSF9" ,   "Tfh_NR3C1"  ,    "Tfh_LAG3"   ,    "Cytotoxic"    ,  "Dysf_GZMB"   ,   "Dysf_ZNF683" ,   
                  "Transitional"  ,   "Naive-like_CD8",   "Non-classical_T_cell"        , "NK_KLRC1"    ,    "NK_FGFBP2"      ,   "B_cell_CD69"   ,
                  "B_cell_CD27"  ,  "Plasma_cell"  ,  "pDC" ,  "cDC_LAMP3"   ,   "cDC_CD1C"      ,"cDC_CLEC9A"   ,  "Mono-macro") 

colors_mc <- colors_mc[match(order_cells, colors_mc$name),]

prep_table$V1 <- factor(prep_table$V1,levels=colors_mc$name)
prep_table <- prep_table[prep_table$V1 != "Granulocyte",] # not enough cells 

prep_table$V1 <- factor(prep_table$V1,levels=(colors_mc$name))
prep_table <- prep_table[complete.cases(prep_table),]
prep_table$Y <- rep(0,nrow(prep_table))


T_prep_table <- prep_table[prep_table$V1 %in% c("Cytotoxic","Dysf_GZMB","Dysf_ZNF683","Non-classical_T_cell","Naive-like_CD4",
                                               "Naive-like_CD8","NK_FGFBP2","NK_KLRC1","Tfh_LAG3","Tfh_NR3C1",
                                               "Transitional","Treg_GIMAP","Treg_TNFRSF9"),]
MB_prep_table <- prep_table[prep_table$V1 %in% c("B_cell_CD27","B_cell_CD69","cDC_CD1C","cDC_CLEC9A","cDC_LAMP3","Mono-macro","pDC","Plasma_cell"),]

T_prep_table <- T_prep_table[order(T_prep_table$median_fc_RE, decreasing = F),]
T_prep_table$V1 <- factor(T_prep_table$V1,levels=T_prep_table$V1)

colors_mc_sub <- colors_mc[colors_mc$name %in% c("Cytotoxic","Dysf_GZMB","Dysf_ZNF683","classical_T_cell","Naive-like_CD4",
                                                "Naive-like_CD8","NK_FGFBP2","NK_KLRC1","Tfh_LAG3","Tfh_NR3C1",
                                                "Transitional","Treg_GIMAP","Treg_TNFRSF9"),]
colors_mc_sub[colors_mc_sub$name %in% prep_table[prep_table$pvalue_RE > 0.1,]$V1,]$color <- "grey"
colors_mc_sub <- colors_mc_sub[match(T_prep_table$V1,colors_mc_sub$name),]
pdf(file=paste0(fig_path,"/RE_median_fold_change_CD48.pdf"), width=4.6, height=3.5)
p1 = ggplot(T_prep_table,aes(x=median_fc_RE,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= (colors_mc_sub$color)) + 
  xlab("Median lfc upon treatment ") +ylab("") + scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  geom_segment(aes(y = V1, yend = V1,x = Y, xend = median_fc_RE),size = 0.8) + xlim(c(-2.5,2.5))+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p1)
dev.off()

MB_prep_table = MB_prep_table[order(MB_prep_table$median_fc_RE, decreasing = F),]
MB_prep_table$V1 = factor(MB_prep_table$V1,levels=MB_prep_table$V1)

colors_mc_sub = colors_mc[colors_mc$name %in% c("B_cell_CD27","B_cell_CD69","cDC_CD1C","cDC_CLEC9A","cDC_LAMP3","Mono-macro","pDC","Plasma_cell"),]
colors_mc_sub[colors_mc_sub$name %in% prep_table[prep_table$pvalue_RE > 0.1,]$V1,]$color = "grey"
colors_mc_sub = colors_mc_sub[match(MB_prep_table$V1,colors_mc_sub$name),]
pdf(file=paste0(fig_path,"/RE_median_fold_change_pre_post_MB.pdf"), width=4.3, height=2.7)
p2 = ggplot(MB_prep_table,aes(x=median_fc_RE,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= (colors_mc_sub$color)) + 
  xlab("Median lfc upon treatment ") +ylab("") + scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  geom_segment(aes(y = V1, yend = V1,x = Y, xend = median_fc_RE),size = 0.8) + xlim(c(-2.5,2.5))+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

pdf(file=paste0(fig_path,"/RE_median_fold_change_pre_post_CD48MB_pvalue.pdf"), width=4, height=6,onefile = F)
ggpubr::ggarrange(p1,p2,ncol=1,nrow=2,heights=c(1.5,1),common.legend = T,legend="right")
dev.off()

# NR plots
medians_cells <- rbind(save_t_cells_medians,save_non_t_cells_medians)
pvals_cells <- rbind(save_pvals_t_cells,save_pvals_non_t_cells)

medians_cells$states <- rownames(medians_cells)
prep_table <- merge(medians_cells,pvals_cells)

prep_table$sig_p <- rep("p>0.5",nrow(prep_table))
prep_table[prep_table$pvalue_NR < 0.5,][["sig_p"]] <- "p<0.5"
prep_table[prep_table$pvalue_NR < 0.1,][["sig_p"]] <- "p<0.1"
prep_table[prep_table$pvalue_NR < 0.05,][["sig_p"]] <- "p<0.05"
#prep_table[prep_table$pvalue_NR < 0.01,][["sig_p"]] <- "p<0.01"

prep_table$sig_p <- factor(prep_table$sig_p,levels=c("p>0.5","p<0.5","p<0.1","p<0.05","p<0.01"))

colors_mc <- read.table("./data/mc_colors.txt",header=T)

prep_table$V1 <- gsub(" ", "_", prep_table$states)
colors_mc <- colors_mc[match(order_cells, colors_mc$name),]

prep_table$V1 <- factor(prep_table$V1,levels=colors_mc$name)
prep_table <- prep_table[prep_table$V1 != "Granulocyte",] # not enough cells 

prep_table$V1 <- factor(prep_table$V1,levels=(colors_mc$name))
prep_table <- prep_table[complete.cases(prep_table),]
prep_table$Y <- rep(0,nrow(prep_table))

T_prep_table <- prep_table[prep_table$V1 %in% c("Cytotoxic","Dysf_GZMB","Dysf_ZNF683","classical_T_cell","Naive-like_CD4",
                                                "Naive-like_CD8","NK_FGFBP2","NK_KLRC1","Tfh_LAG3","Tfh_NR3C1",
                                                "Transitional","Treg_GIMAP","Treg_TNFRSF9"),]
MB_prep_table <- prep_table[prep_table$V1 %in% c("B_cell_CD27","B_cell_CD69","cDC_CD1C","cDC_CLEC9A","cDC_LAMP3","Mono-macro","pDC","Plasma_cell"),]

T_prep_table <- T_prep_table[order(T_prep_table$median_fc_NR, decreasing = F),]
T_prep_table$V1 <- factor(T_prep_table$V1,levels=T_prep_table$V1)

colors_mc_sub <- colors_mc[colors_mc$name %in% c("Cytotoxic","Dysf_GZMB","Dysf_ZNF683","classical_T_cell","Naive-like_CD4",
                                                 "Naive-like_CD8","NK_FGFBP2","NK_KLRC1","Tfh_LAG3","Tfh_NR3C1",
                                                 "Transitional","Treg_GIMAP","Treg_TNFRSF9"),]
colors_mc_sub[colors_mc_sub$name %in% prep_table[prep_table$pvalue_NR > 0.1,]$V1,]$color <- "grey"
colors_mc_sub <- colors_mc_sub[match(T_prep_table$V1,colors_mc_sub$name),]
pdf(file=paste0(fig_path,"/NR_median_fold_change_CD48.pdf"), width=4.3, height=3.5)
p1 = ggplot(T_prep_table,aes(x=median_fc_NR,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= (colors_mc_sub$color)) + 
  xlab("Median lfc upon treatment ") +ylab("") + scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  geom_segment(aes(y = V1, yend = V1,x = Y, xend = median_fc_NR),size = 0.8) + xlim(c(-2.5,2.5))+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p1)
dev.off()

MB_prep_table = MB_prep_table[order(MB_prep_table$median_fc_NR, decreasing = F),]
MB_prep_table$V1 = factor(MB_prep_table$V1,levels=MB_prep_table$V1)

colors_mc_sub = colors_mc[colors_mc$name %in% c("B_cell_CD27","B_cell_CD69","cDC_CD1C","cDC_CLEC9A","cDC_LAMP3","Mono-macro","pDC","Plasma_cell"),]
colors_mc_sub[colors_mc_sub$name %in% prep_table[prep_table$pvalue_NR > 0.1,]$V1,]$color = "grey"
colors_mc_sub = colors_mc_sub[match(MB_prep_table$V1,colors_mc_sub$name),]
pdf(file=paste0(fig_path,"/NR_median_fold_change_pre_post_MB.pdf"), width=4.3, height=2.7)
p2 = ggplot(MB_prep_table,aes(x=median_fc_NR,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= (colors_mc_sub$color)) + 
  xlab("Median lfc upon treatment ") +ylab("") + scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  geom_segment(aes(y = V1, yend = V1,x = Y, xend = median_fc_NR),size = 0.8) + xlim(c(-2.5,2.5))+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

pdf(file=paste0(fig_path,"/NR_median_fold_change_pre_post_CD48MB_pvalue.pdf"), width=4, height=6,onefile=F)
ggpubr::ggarrange(p1,p2,ncol=1,nrow=2,heights=c(1.5,1),common.legend = T,legend = "right")
dev.off()



#### Figure 3B ####
# line charts, T cells
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

pre = md[md$condition == "pre" & 
           md$cell_type %in% c("CD4","CD8","NK") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12" &
           !is.na(md$mc_group),]

post = md[md$condition == "post" & 
            md$cell_type %in% c("CD4","CD8","NK") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12" &
            !is.na(md$mc_group),]

# calculate fractions
pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

pre_final$condition <- rep("pre",nrow(pre_final))
post_final$condition <- rep("post",nrow(post_final))
pre_final$patient <- rownames(pre_final)
post_final$patient <- rownames(post_final)
prepost_final <- rbind(pre_final,post_final)

pat_response <- unique(md[,c("patient","response")])
NR_pat <- pat_response[pat_response$response == "NR",]$patient
RE_pat <- pat_response[pat_response$response == "RE",]$patient
prepost_final_NR <- prepost_final[prepost_final$patient %in% NR_pat,]
prepost_final_RE <- prepost_final[prepost_final$patient %in% RE_pat,]

prepost_final_RE$condition <- factor(prepost_final_RE$condition,levels=c("pre","post"))
prepost_final_NR$condition <- factor(prepost_final_NR$condition,levels=c("pre","post"))

## counts
post_frac_c <- post_frac[,c("patient","mc_group","count")]
pre_frac_c <- pre_frac[,c("patient","mc_group","count")]
pre_frac_c$condition <- rep("pre",nrow(pre_frac_c))
post_frac_c$condition <- rep("post",nrow(post_frac_c))
prepost_frac_c <- rbind(pre_frac_c,post_frac_c)

# plot with ggplot boxplot
for (cell_x in c("Treg GIMAP","Treg TNFRSF9","Dysf ZNF683","Dysf GZMB","Transitional","Naive-like CD8")){
  print(paste0("cell: ",cell_x))
  ylim_ax = c(0,0.5)
  make_line_boxplot(cell_x,ylim_ax,prepost_frac_c,prepost_final_RE,prepost_final_NR)
}


# line charts, non T cells
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

pre = md[md$condition == "pre" & 
           md$cell_type %in% c("M","B") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12" &
           !is.na(md$mc_group),]

post = md[md$condition == "post" & 
            md$cell_type %in% c("M","B") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12" &
            !is.na(md$mc_group),]

# calculate fractions
pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

pre_final$condition <- rep("pre",nrow(pre_final))
post_final$condition <- rep("post",nrow(post_final))
pre_final$patient <- rownames(pre_final)
post_final$patient <- rownames(post_final)
prepost_final <- rbind(pre_final,post_final)

pat_response <- unique(md[,c("patient","response")])
NR_pat <- pat_response[pat_response$response == "NR",]$patient
RE_pat <- pat_response[pat_response$response == "RE",]$patient
prepost_final_NR <- prepost_final[prepost_final$patient %in% NR_pat,]
prepost_final_RE <- prepost_final[prepost_final$patient %in% RE_pat,]

prepost_final_RE$condition <- factor(prepost_final_RE$condition,levels=c("pre","post"))
prepost_final_NR$condition <- factor(prepost_final_NR$condition,levels=c("pre","post"))

## counts
post_frac_c <- post_frac[,c("patient","mc_group","count")]
pre_frac_c <- pre_frac[,c("patient","mc_group","count")]

pre_frac_c$condition <- rep("pre",nrow(pre_frac_c))
post_frac_c$condition <- rep("post",nrow(post_frac_c))

prepost_frac_c <- rbind(pre_frac_c,post_frac_c)

# plot with ggplot boxplot
for (cell_x in c("Mono-macro")){
  print(paste0("cell: ",cell_x))
  ylim_ax = c(0,0.5)
  make_line_boxplot(cell_x,ylim_ax,prepost_frac_c,prepost_final_RE,prepost_final_NR)
}



#### Figure 3C ####
## Generate NK signatures
NK_sigs = generate_NKcell_signatures()

# CD56 bright signature
geneSets <- list(geneSet1=NK_sigs$CD56_bright)

# pre select NK cells
mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat = mat@cell_metadata       
md = submat[submat$mc_group %in% c("NK FGFBP2","NK KLRC1") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("post","pre") &
              submat$patient != "Pat38" &
              submat$patient != "Pat04" &
              submat$patient != "Pat12",]
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) <- geneNames

# UCell score
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(geneSets,cell_rank)
check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_all$cond_resp <- paste0(sig_scores_all$condition,sig_scores_all$response)
sig_scores_all$cond_resp <- factor(sig_scores_all$cond_resp,levels=c("preNR","preRE","postNR","postRE"))
pdf(file=paste0(fig_path,"/NK_CD56bright_violingplot.pdf"), width=3.6, height=2.8)
ggplot(sig_scores_all, aes(y = geneSet1, x = response,fill=cond_resp))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#DC7985","#89CC9B","#9B1F23","#116936"))+
  theme_bw() +  ylab("NK cd56bright") + ggtitle("NK cd56bright") + ylim(c(0.0,0.7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

wilcox.test(sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "pre",][["geneSet1"]],
            sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "post",][["geneSet1"]])

wilcox.test(sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "pre",][["geneSet1"]],
            sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "post",][["geneSet1"]])


# activate/mature signature
geneSets <- list(geneSet1=NK_sigs$active_mature)

# pre select NK cells
mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("NK FGFBP2","NK KLRC1") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("post","pre") &
              submat$patient != "Pat38" &
              submat$patient != "Pat04" &
              submat$patient != "Pat12",]
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) <- geneNames

# UCell score
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(geneSets,cell_rank)
check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

sig_scores_all$cond_resp <- paste0(sig_scores_all$condition,sig_scores_all$response)
sig_scores_all$cond_resp <- factor(sig_scores_all$cond_resp,levels=c("preNR","preRE","postNR","postRE"))
pdf(file=paste0(fig_path,"/NK_active_mature_violingplot.pdf"), width=3.6, height=2.8)
ggplot(sig_scores_all, aes(y = geneSet1, x = response,fill=cond_resp))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#DC7985","#89CC9B","#9B1F23","#116936"))+
  theme_bw() +  ylab("NK active") + ggtitle("NK_active") + ylim(c(0.1,0.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

wilcox.test(sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "pre",][["geneSet1"]],
            sig_scores_all[sig_scores_all$response == "RE" & sig_scores_all$condition == "post",][["geneSet1"]])

wilcox.test(sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "pre",][["geneSet1"]],
            sig_scores_all[sig_scores_all$response == "NR" & sig_scores_all$condition == "post",][["geneSet1"]])



#### Figure 3D ####
data = mat@cell_metadata[names(mc@mc),]
data = data[data$clonal_exp != "no" & 
              !is.na(data$response) &
              data$pat_condition != "Pat04_post" &
              data$pat_condition != "Pat12_post" &
              data$patient != "Pat38",]

trans = data[data$mc_group %in% c("Dysf GZMB","Dysf ZNF683","Transitional",
                                  "Naive-like CD8", "Cytotoxic","Naive-like CD4","Treg GIMAP","Treg TNFRSF9", "Tfh LAG3","Tfh NR3C1" ),]
trans <- as.data.frame(trans)
trans$Cell_Name <- rownames(trans)
trans$clone.id <- trans[,colnames(trans) == "unique_clone_ID"]
trans$majorCluster <- trans[,colnames(trans) == "mc_group"]
trans$loc <- trans[,colnames(trans) == "condition"] 
out_trans <- Startrac.run(trans, proj="IMCISION",verbose=F)

# get data out of startrac object
out_cluster <- as.data.frame(out_trans@cluster.data)
out_cluster <- out_cluster[!is.na(out_cluster$expa),]
out_cluster_tmp <- out_cluster[out_cluster$aid != "IMCISION",]

# heatmap
plot = as.data.frame(out_trans@pIndex.sig.tran)
plot = plot[plot$aid == "IMCISION",]
plot = plot[,c("majorCluster", "index", "value")]

plot_2 = spread(plot, key = "index", value = "value")
rownames(plot_2) = plot_2$majorCluster

# annotation colors mc groups
mc_colors <- read.table("./data/mc_colors.txt",header=T,sep="\t")
mc_colors$names_correct <- gsub("_"," ",mc_colors$name)
an_colors = mc_colors[,2]
names(an_colors) = mc_colors[,4]
row_ha = rowAnnotation(cellstate = rownames(plot_2), col = list("cellstate" = an_colors))
col_ha = HeatmapAnnotation(cellstate = colnames(plot_2)[-1], col = list("cellstate" = an_colors))

cols = BlueAndRed()
pdf(file=paste0(fig_path,"/Startrac_shared_clones_heatmap_CD4CD8.pdf"), width=5, height=4.4)
Heatmap(as.matrix(plot_2[,2:11]), col = cols,
        cluster_columns = T, cluster_rows = T,name = " ",bottom_annotation = col_ha,
        right_annotation = row_ha,
        column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10), 
        column_dend_height = unit(20, "mm"), row_dend_width = unit(20, "mm"))
dev.off()



#### Figure 3E ####
# Cell cycle genes in Dysf GZMB, naive-like CD8 and transitionals
lat_gset = read.table("./data/signatures/mel_lateral_genesets_CC.txt",header=T)
cc_genes = lat_gset$x

cc_umis <- mat@mat[intersect(cc_genes, rownames(mat@mat)),]
cc_fraction <<- colSums(cc_umis)/colSums(mat@mat)
submat <- mat@cell_metadata[names(mc@mc),]
submat <- submat[!is.na(submat$patient) &
                   submat$mc_group %in% c("Dysf GZMB", "Transitional", "Naive-like CD8") &
                   submat$response == "RE",]
submat$cell_cycle <- log2(cc_fraction[rownames(submat)]+1)

pdf(file=paste0(fig_path,"/Cell_cycle_Anne_avg_RE.pdf"),width=4.2, height=3)
lognorm = ggplot(submat, aes(x=condition, y=cell_cycle,fill=mc_group)) +
  geom_point(position=position_jitter(h=0, w=0.1), fill = "lightgrey", pch = 21, alpha = 0.5, size = 0.5)+
  geom_violin()+
  theme_bw() +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values= c("#4B1D5B","#F37D7F","#C9A3C8"))+
  facet_wrap(~ mc_group)
print(lognorm)
dev.off()

re <- submat[submat$mc_group == "Dysf GZMB",]
print(wilcox.test(re$cell_cycle ~ re$condition))
re <- submat[submat$mc_group == "Transitional",]
print(wilcox.test(re$cell_cycle ~ re$condition)) 
re <- submat[submat$mc_group == "Naive-like CD8",]
print(wilcox.test(re$cell_cycle ~ re$condition)) 




#### Figure 3F ####
# Top clones transitionals
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

# responders
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "RE" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional"),]

# determine top clones pre on in transitionals
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
md_clones_top <- md_clones_all[md_clones_all$freq_pre > 2 &md_clones_all$freq_post > 2,]
md_clones_top_id <- md_clones_top$pat_clone_id # top clones transitionals in RE

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

# determine fractions pre and post
sum_md_toppers_post <- table(md_toppers_post[,c("pat_clone_ID","mc_group")]) 
sum_md_toppers_post <- cbind(sum_md_toppers_post, total = rowSums(sum_md_toppers_post))
sum_md_toppers_post <- as.data.frame.matrix(sum_md_toppers_post)
sum_md_toppers_post_sw_NR <- sweep(sum_md_toppers_post, 1, sum_md_toppers_post[,7], "/")

sum_md_toppers_pre <- table(md_toppers_pre[,c("pat_clone_ID","mc_group")]) 
sum_md_toppers_pre <- cbind(sum_md_toppers_pre, total = rowSums(sum_md_toppers_pre))
sum_md_toppers_pre <- as.data.frame.matrix(sum_md_toppers_pre)
sum_md_toppers_pre_sw_NR <- sweep(sum_md_toppers_pre, 1, sum_md_toppers_pre[,7], "/")

sum_md_toppers_post_sw_NR <- sum_md_toppers_post_sw_NR[,c(colnames(sum_md_toppers_pre_sw_NR))]
sum_md_fraction_change_NR <- sum_md_toppers_post_sw_NR-sum_md_toppers_pre_sw_NR

sum_md_toppers_post <- sum_md_toppers_post[,colnames(sum_md_toppers_pre_sw_NR)[-7]]
sum_md_toppers_pre <- sum_md_toppers_pre[,colnames(sum_md_toppers_pre_sw_NR)[-7]]
sum_md_toppers_post_pre <- sum_md_toppers_post + sum_md_toppers_pre

sum_md_toppers_post_pre$clone_id <- rownames(sum_md_toppers_post_pre)
temp <- as.data.frame(rowSums(sum_md_toppers_post_pre[,1:6]))
temp$clone_id <- rownames(temp)
order_sort <- rownames(temp[order(temp$`rowSums(sum_md_toppers_post_pre[, 1:6])`,decreasing = T),])
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[match(order_sort,sum_md_toppers_post_pre$clone_id),]
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[1:15,]
sum_md_toppers_post_pre_melt <- reshape2::melt(sum_md_toppers_post_pre)

sum_md_fraction_change_NR$clone_id <- rownames(sum_md_fraction_change_NR)
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[match(order_sort,sum_md_fraction_change_NR$clone_id),]
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[1:15,][,-7]
sum_md_fraction_change_NR_melt <- reshape2::melt(sum_md_fraction_change_NR)

sum_md_fraction_change_NR_melt$clone_id <- factor(sum_md_fraction_change_NR_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-7])
sum_md_fraction_change_NR_melt$variable <- factor(sum_md_fraction_change_NR_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-7])

pdf(file=paste0(fig_path,"/Change_clones_onlytransitional_RE.pdf"),width=3, height=5)
p=ggplot(sum_md_fraction_change_NR_melt[sum_md_fraction_change_NR_melt$variable == "Transitional",],aes(x=variable,y=clone_id)) + 
  geom_point(color="black",shape=21,aes(fill=value, size=sum_md_toppers_post_pre_melt[sum_md_toppers_post_pre_melt$variable == "Transitional",]$value)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", limits = c(-1, 1),name = "lfc") +
  scale_size_continuous(range = c(0,8),name = "clone_size_sum",limits=c(0,82)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

# non-responders
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "NR" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional"),]

# determine top clones pre on in transitionals
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
md_clones_top <- md_clones_all[md_clones_all$freq_pre > 2 &md_clones_all$freq_post > 2,]
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

sum_md_toppers_post <- table(md_toppers_post[,c("pat_clone_ID","mc_group")]) 
sum_md_toppers_post <- cbind(sum_md_toppers_post, total = rowSums(sum_md_toppers_post))
sum_md_toppers_post <- as.data.frame.matrix(sum_md_toppers_post)
sum_md_toppers_post_sw_NR <- sweep(sum_md_toppers_post, 1, sum_md_toppers_post[,6], "/")

sum_md_toppers_pre <- table(md_toppers_pre[,c("pat_clone_ID","mc_group")]) 
sum_md_toppers_pre <- cbind(sum_md_toppers_pre, total = rowSums(sum_md_toppers_pre))
sum_md_toppers_pre <- as.data.frame.matrix(sum_md_toppers_pre)
sum_md_toppers_pre_sw_NR <- sweep(sum_md_toppers_pre, 1, sum_md_toppers_pre[,5], "/")

# fill missing cell states
sum_md_toppers_pre_sw_NR$`Naive-like CD8` <- rep(0.0,nrow(sum_md_toppers_pre_sw_NR))
sum_md_toppers_pre_sw_NR$`Dysf ZNF683` <- rep(0.0,nrow(sum_md_toppers_pre_sw_NR))
sum_md_toppers_post_sw_NR$`Non-classical T cell` <- rep(0.0,nrow(sum_md_toppers_post_sw_NR))
sum_md_toppers_post_sw_NR <- sum_md_toppers_post_sw_NR[,c(colnames(sum_md_toppers_pre_sw_NR))]

sum_md_fraction_change_NR <- sum_md_toppers_post_sw_NR-sum_md_toppers_pre_sw_NR

sum_md_toppers_pre$`Naive-like CD8` <- rep(0.0,nrow(sum_md_toppers_pre))
sum_md_toppers_pre$`Dysf ZNF683` <- rep(0.0,nrow(sum_md_toppers_pre))
sum_md_toppers_post$`Non-classical T cell` <- rep(0.0,nrow(sum_md_toppers_post))
sum_md_toppers_post <- sum_md_toppers_post[,colnames(sum_md_toppers_post_sw_NR)[-5]]
sum_md_toppers_pre <- sum_md_toppers_pre[,colnames(sum_md_toppers_pre_sw_NR)[-5]]
sum_md_toppers_post_pre <- sum_md_toppers_post + sum_md_toppers_pre

sum_md_toppers_post_pre$clone_id <- rownames(sum_md_toppers_post_pre)
temp <- as.data.frame(rowSums(sum_md_toppers_post_pre[,1:6]))
temp$clone_id <- rownames(temp)
order_sort <- rownames(temp[order(temp$`rowSums(sum_md_toppers_post_pre[, 1:6])`,decreasing = T),])
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[match(order_sort,sum_md_toppers_post_pre$clone_id),]
sum_md_toppers_post_pre <- sum_md_toppers_post_pre[1:9,]
sum_md_toppers_post_pre_melt <- reshape2::melt(sum_md_toppers_post_pre)

sum_md_fraction_change_NR$clone_id <- rownames(sum_md_fraction_change_NR)
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[match(order_sort,sum_md_fraction_change_NR$clone_id),]
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[1:9,][,-5]
sum_md_fraction_change_NR_melt <- reshape2::melt(sum_md_fraction_change_NR)

sum_md_fraction_change_NR_melt$clone_id <- factor(sum_md_fraction_change_NR_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-5])
sum_md_fraction_change_NR_melt$variable <- factor(sum_md_fraction_change_NR_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-5])

pdf(file=paste0(fig_path,"/Change_clones_onlytransitional_NR.pdf"),width=3, height=3.5)
p=ggplot(sum_md_fraction_change_NR_melt[sum_md_fraction_change_NR_melt$variable == "Transitional",],aes(x=variable,y=clone_id)) + 
  geom_point(color="black",shape=21,aes(fill=value, size=sum_md_toppers_post_pre_melt[sum_md_toppers_post_pre_melt$variable == "Transitional",]$value)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", limits = c(-1, 1),name = "lfc") +
  scale_size_continuous(range = c(0,8),name = "clone_size_sum",limits=c(0,82)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()



#### Figure 3G ####
## new clones
df = mat@cell_metadata[names(mc@mc),]
# select df with all cells
df = df[df$response %in% c("RE", "NR") &
          #df$cell_type %in% c("CD8", "CD4") &
          df$patient != "Pat04" &
          df$patient != "Pat12" &
          df$patient != "Pat38",]

df = df[!is.na(df$response) &
          !is.na(df$mc_group) &
          !is.na(df$condition) &
          !is.na(df$unique_clone_ID),]

df[is.na(df$shared_TCR), "shared_TCR"] = "non-overlapping"
tot_cells = df[!is.na(df$mc_group),]

# count nr of cells in either pre or post per cell state
cellnrs = tot_cells %>% 
  group_by(patient, condition) %>% 
  summarise(count = n())
ds = cellnrs %>% 
  group_by(patient) %>% 
  summarise(min = min(count))
ds = ds %>% 
  textshape::column_to_rownames("patient")

for (pat in rownames(ds)){
  tot_cells[tot_cells$patient == pat, "sample"] = ds[pat,"min"]
}

# downsample separately per patient
downsamp = tot_cells %>%
  group_by(patient, condition) %>%
  mutate(random_id = sample(n())) %>%
  ungroup() %>%
  dplyr::filter(random_id <= sample) %>%
  dplyr::select(-random_id)

# select unique clones
md = downsamp[downsamp$shared_TCR == "non-overlapping" & 
                !is.na(downsamp$unique_clone_ID) & downsamp$response == "RE",]
clones = unique(md[md$dom_state == "Transitional", "unique_clone_ID"])
md = md[md$unique_clone_ID %in% clones$unique_clone_ID,]

x = md %>% group_by(response, patient, condition, unique_clone_ID) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
x = na.omit(x)

x$unique_clone_ID = as.character(x$unique_clone_ID)
x$unique_clone_ID <- reorder(x$unique_clone_ID, x$count)
x$unique_clone_ID <- factor(x$unique_clone_ID, levels=rev(levels(x$unique_clone_ID)))

x$condition = factor(x$condition,levels=c("post","pre"))

pdf(file=paste0(fig_path,"/New_clones_onlytransitional_RE.pdf"),width=3, height=4)
b = ggplot(x, aes(x = condition, y = count))+
  geom_bar(stat = "identity", position = "stack", aes(fill = unique_clone_ID), color = "lightgrey", size = 0.1)+
  #scale_fill_manual(values = viridis::viridis(n=810))+
  facet_wrap(patient ~ ., nrow = 10,strip.position="right")+
  theme_bw() +
  theme(strip.text = element_text(size=8),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_flip()
print(b)
dev.off()




#### Figure 3H ####
# load Adapative data
working_dir = "~/Documents/Zuur_lab/Single_cell_data/IMCISION_shared_Joleen/Adaptive/Data/"
all_templates_samples = load_adpative_data_bulkTCR(working_dir)

remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
md = mat@cell_metadata
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response == "RE" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional"),]

md$pat_clone <- paste0(md$patient,"_",md$TRB_nt)
md_pre <- md[md$condition == "pre",]
md_post <- md[md$condition == "post",]

shared_clones <- intersect(md_post$pat_clone,md_pre$pat_clone)

# calculate frequencies TCRs
data_pre <- md_pre %>% 
  group_by(patient,pat_clone) %>% 
  dplyr::summarize(n=n()) %>% 
  mutate(freq_pre=n/sum(n)) 
data_post <- md_post %>% 
  group_by(patient,pat_clone) %>% 
  dplyr::summarize(n=n()) %>% 
  mutate(freq_post=n/sum(n)) 

shared_pre <- data_pre[data_pre$pat_clone %in% shared_clones,]
shared_post <- data_post[data_post$pat_clone %in% shared_clones,]

shared_both <- cbind(shared_pre,shared_post)

nonshared_pre <- data_pre[!(data_pre$pat_clone %in% shared_clones),]
nonshared_post <- data_post[!(data_post$pat_clone %in% shared_clones),] # of interest

# post treatment
# prepare df to save data in 
to_plot_all <- c(paste0(unique(nonshared_post$patient),"_no"),
                paste0(unique(nonshared_post$patient),"_yes"))
to_plot_all <- to_plot_all[order(to_plot_all)]
to_plot_all = as.data.frame(to_plot_all)
colnames(to_plot_all) = "pat_b"

# downsample 100 times, calculate frequency TCR in blood, nonshared_post
for (r_s in 1:100){
  ds_templates_samples <- as.data.frame(matrix(ncol=ncol(all_templates_samples),nrow=0))
  colnames(ds_templates_samples) <- colnames(all_templates_samples)
  for (sample_x in unique(all_templates_samples$test)){
    subset_df <- all_templates_samples[all_templates_samples$test == sample_x,]
    subset_df_ds <- subset_df[rownames(subset_df) %in% sample(rownames(subset_df),39723),] # minimum TCRs sample
    ds_templates_samples <- rbind(ds_templates_samples,subset_df_ds)
  } 
  
  ds_templates_samples$pat_clone <- paste0(ds_templates_samples$Patient,"_",ds_templates_samples$cdr3_rearrangement)
  pre_all_templates <- ds_templates_samples[ds_templates_samples$condition == "pre",]
  pre_all_templates$pat_clone <- paste0(pre_all_templates$Patient,"_",pre_all_templates$cdr3_rearrangement)
  nonshared_post$blood_pre <- "no"
  nonshared_post[nonshared_post$pat_clone %in% pre_all_templates$pat_clone,]$blood_pre <- "yes"
  post_all_templates <- ds_templates_samples[ds_templates_samples$condition == "post",]
  post_all_templates$pat_clone <- paste0(post_all_templates$Patient,"_",post_all_templates$cdr3_rearrangement)
  nonshared_post$blood_post <- "no"
  nonshared_post[nonshared_post$pat_clone %in% post_all_templates$pat_clone,]$blood_post <- "yes"
  
  to_plot_pre <- nonshared_post %>% 
    group_by(patient,blood_pre) %>% 
    dplyr::summarize(n=n()) %>% 
    mutate(freq_pre=n/sum(n)) 
  to_plot_pre = to_plot_pre %>% 
    ungroup %>% 
    complete(patient, blood_pre, fill = list(freq = 0))
  to_plot_post <- nonshared_post %>% 
    group_by(patient,blood_post) %>% 
    dplyr::summarize(n=n()) %>% 
    mutate(freq_post=n/sum(n)) 
  to_plot_post = to_plot_post %>% 
    ungroup %>% 
    complete(patient, blood_post, fill = list(freq = 0))
  
  to_plot_pre$pat_b <- paste0(to_plot_pre$patient,"_",to_plot_pre$blood_pre)
  to_plot_post$pat_b <- paste0(to_plot_post$patient,"_",to_plot_post$blood_post)
  to_plot <- merge(to_plot_post,to_plot_pre,by="pat_b")
  to_plot <- to_plot[,c("pat_b","freq_post","freq_pre")]
  colnames(to_plot)[2:3] <- c(paste0("freq_post_",r_s),paste0("freq_pre_",r_s))
  
  to_plot_all <- merge(to_plot_all,to_plot,by="pat_b")
}

col_odd <- seq_len(ncol(to_plot_all)-4) %% 2  
col_odd <- c(2,2,2,2,col_odd)

to_plot_all$mean_freq_post <- rowMeans(to_plot_all[,which(col_odd == 0)])
to_plot_all$mean_freq_pre <- rowMeans(to_plot_all[,which(col_odd == 1)])
to_plot_all$mean_sd_post <- apply(to_plot_all[,which(col_odd == 0)], 1, sd, na.rm=TRUE)
to_plot_all$mean_sd_pre <- apply(to_plot_all[,which(col_odd == 1)], 1, sd, na.rm=TRUE)
to_plot_all_ready <- to_plot_all[,c("pat_b","mean_freq_pre","mean_freq_post","mean_sd_post","mean_sd_pre")]
to_plot_all_ready <- as.data.frame(cbind(c(to_plot_all_ready$mean_freq_pre,to_plot_all_ready$mean_freq_post),
                                        c(to_plot_all_ready$mean_sd_pre,to_plot_all_ready$mean_sd_post),
                                        c(rep("pre",nrow(to_plot_all_ready)),rep("post",nrow(to_plot_all_ready))),
                                        c(to_plot_all_ready$pat_b,to_plot_all_ready$pat_b)))

colnames(to_plot_all_ready) <- c("Freq","sd_freq","condition","pat_b")
to_plot_all_ready$patient <- do.call(rbind,strsplit(to_plot_all_ready$pat_b,"_"))[,1]
to_plot_all_ready$blood <- do.call(rbind,strsplit(to_plot_all_ready$pat_b,"_"))[,2]
to_plot_all_ready$Freq <- as.numeric(to_plot_all_ready$Freq)
to_plot_all_ready$sd_freq <- as.numeric(to_plot_all_ready$sd_freq)

#  boxplot, per patient, blood/non-blood pre post, nonshared post
to_plot_tmp = to_plot_all_ready
to_plot_tmp = to_plot_tmp[to_plot_tmp$blood == "yes",]
to_plot_tmp$condition = factor(to_plot_tmp$condition, levels=c("pre","post"))
pdf(file=paste0(fig_path,"/Freq_unique_TCRs_transitionals_inblood_boxplot_TCRfraction.pdf"), width=3.5,height=4)
ggplot(to_plot_tmp,aes(x=condition,y=as.numeric(Freq),fill=condition)) + geom_boxplot() +
  geom_line(aes(group=patient)) + geom_point() +
  ggpubr::stat_compare_means(paired=T)+ ggtitle("Freq unique TCRs post in blood") +
  scale_fill_manual(values=c("#89CC9B","#116936")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# pre-treatment    
# prepare df to save data in 
to_plot_all <- c(paste0(unique(nonshared_pre$patient),"_no"),
                 paste0(unique(nonshared_pre$patient),"_yes"))
to_plot_all <- to_plot_all[order(to_plot_all)]
to_plot_all <- as.data.frame(to_plot_all)
colnames(to_plot_all) <- "pat_b"

# downsample 100 times, calculate frequency TCR in blood, nonshared_pre
for (r_s in 1:100){
  ds_templates_samples <- as.data.frame(matrix(ncol=ncol(all_templates_samples),nrow=0))
  colnames(ds_templates_samples) <- colnames(all_templates_samples)
  for (sample_x in unique(all_templates_samples$test)){
    subset_df <- all_templates_samples[all_templates_samples$test == sample_x,]
    subset_df_ds <- subset_df[rownames(subset_df) %in% sample(rownames(subset_df),39723),]
    ds_templates_samples <- rbind(ds_templates_samples,subset_df_ds)
  } 
  
  ds_templates_samples$pat_clone <- paste0(ds_templates_samples$Patient,"_",ds_templates_samples$cdr3_rearrangement)
  pre_all_templates <- ds_templates_samples[ds_templates_samples$condition == "pre",]
  pre_all_templates$pat_clone <- paste0(pre_all_templates$Patient,"_",pre_all_templates$cdr3_rearrangement)
  nonshared_pre$blood_pre <- "no"
  nonshared_pre[nonshared_pre$pat_clone %in% pre_all_templates$pat_clone,]$blood_pre <- "yes"
  post_all_templates <- ds_templates_samples[ds_templates_samples$condition == "post",]
  post_all_templates$pat_clone <- paste0(post_all_templates$Patient,"_",post_all_templates$cdr3_rearrangement)
  nonshared_pre$blood_post <- "no"
  nonshared_pre[nonshared_pre$pat_clone %in% post_all_templates$pat_clone,]$blood_post <- "yes"
  
  to_plot_pre <- nonshared_pre %>% 
    group_by(patient,blood_pre) %>% 
    dplyr::summarize(n=n()) %>% 
    mutate(freq_pre=n/sum(n)) 
  to_plot_pre = to_plot_pre %>% 
    ungroup %>% 
    complete(patient, blood_pre, fill = list(freq = 0))
  to_plot_post <- nonshared_pre %>% 
    group_by(patient,blood_post) %>% 
    dplyr::summarize(n=n()) %>% 
    mutate(freq_post=n/sum(n)) 
  to_plot_post <- to_plot_post %>% 
    ungroup %>% 
    complete(patient, blood_post, fill = list(freq = 0))
  
  to_plot_pre$pat_b <- paste0(to_plot_pre$patient,"_",to_plot_pre$blood_pre)
  to_plot_post$pat_b <- paste0(to_plot_post$patient,"_",to_plot_post$blood_post)
  to_plot <- merge(to_plot_post,to_plot_pre,by="pat_b")
  to_plot <- to_plot[,c("pat_b","freq_post","freq_pre")]
  colnames(to_plot)[2:3] <- c(paste0("freq_post_",r_s),paste0("freq_pre_",r_s))
  
  to_plot_all <- merge(to_plot_all,to_plot,by="pat_b")
}

col_odd <- seq_len(ncol(to_plot_all)-4) %% 2  
col_odd <- c(2,2,2,2,col_odd)

to_plot_all$mean_freq_post <- rowMeans(to_plot_all[,which(col_odd == 0)], na.rm=TRUE)
to_plot_all$mean_freq_pre <- rowMeans(to_plot_all[,which(col_odd == 1)], na.rm=TRUE)
to_plot_all$mean_sd_post <- apply(to_plot_all[,which(col_odd == 0)], 1, sd, na.rm=TRUE)
to_plot_all$mean_sd_pre <- apply(to_plot_all[,which(col_odd == 1)], 1, sd, na.rm=TRUE)
to_plot_all_ready <- to_plot_all[,c("pat_b","mean_freq_pre","mean_freq_post","mean_sd_post","mean_sd_pre")]
to_plot_all_ready <- as.data.frame(cbind(c(to_plot_all_ready$mean_freq_pre,to_plot_all_ready$mean_freq_post),
                                        c(to_plot_all_ready$mean_sd_pre,to_plot_all_ready$mean_sd_post),
                                        c(rep("pre",nrow(to_plot_all_ready)),rep("post",nrow(to_plot_all_ready))),
                                        c(to_plot_all_ready$pat_b,to_plot_all_ready$pat_b)))

colnames(to_plot_all_ready) <- c("Freq","sd_freq","condition","pat_b")
to_plot_all_ready$patient <- do.call(rbind,strsplit(to_plot_all_ready$pat_b,"_"))[,1]
to_plot_all_ready$blood <- do.call(rbind,strsplit(to_plot_all_ready$pat_b,"_"))[,2]
to_plot_all_ready$Freq <- as.numeric(to_plot_all_ready$Freq)
to_plot_all_ready$sd_freq <- as.numeric(to_plot_all_ready$sd_freq)

#  boxplot, per patient, blood/non-blood pre post, nonshared pre
to_plot_tmp = to_plot_all_ready
to_plot_tmp = to_plot_tmp[to_plot_tmp$blood == "yes",]
to_plot_tmp$condition = factor(to_plot_tmp$condition, levels=c("pre","post"))
pdf(file=paste0(fig_path,"/Freq_unique_PRE_TCRs_transitionals_inblood_boxplot_TCRfraction.pdf"), width=3.5,height=4)
ggplot(to_plot_tmp,aes(x=condition,y=as.numeric(Freq),fill=condition)) + geom_boxplot() +
  geom_line(aes(group=patient)) + geom_point() +
  ggpubr::stat_compare_means(paired=T)+ ggtitle("Freq pre unique TCRs post in blood") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("#89CC9B","#116936")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()







