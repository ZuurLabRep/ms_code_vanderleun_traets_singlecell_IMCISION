###################
#### FIGURE S4 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S4"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S4")

library(ggrepel)

#### Figure S4A ####
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata    

md <- submat[!is.na(submat$mc_group),]
md_sel <- md[,c("patient","response")]
md_sel <- md_sel[!duplicated(md_sel),]

x <- md %>% group_by(pat_condition) %>% 
  summarise(count = n())
x <- as.data.frame(x)
x$patient <- do.call(rbind,strsplit(x$pat_condition,"_"))[,1]
x$condition <- do.call(rbind,strsplit(x$pat_condition,"_"))[,2]
x$response <- md_sel$response[match(x$patient,md_sel$patient)]
x$cond_res <- paste0(x$response,x$condition)
x$condition <- factor(x$condition,levels=c("pre","post"))

x <- x[x$patient != "Pat04",]
pdf(file=paste0(fig_path,"/Total_cell_count_both.pdf"), width=2.8, height=3)
ggplot(x,aes(x=condition,y=count)) + geom_boxplot(aes(x=condition)) +
  geom_line(aes(group=patient,color=response)) +
  geom_point(aes(color=response)) +
  ggtitle("Total cells")+
  theme_bw()+
  scale_color_manual(values=c("#dc7986","#71bc88"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggpubr::stat_compare_means() + ylim(c(0,3000))
dev.off()


#### Figure S4B ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
mc_groups <- unique(md$mc_group)

pre <- md[md$condition == "pre" & 
           md$cell_type %in% c("CD4","CD8","NK") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12",]
pre <- pre[pre$mc_group %in% mc_groups,]

post <- md[md$condition == "post" & 
            md$cell_type %in% c("CD4","CD8","NK") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12",]
post <- post[post$mc_group %in% mc_groups,]

pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

states <- colnames(pre_final)
samples <- rownames(pre_final)
fold_change <- post_final[samples, states]/pre_final[samples, states]
fold_change <- as.matrix(fold_change)

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
save_pvals_t_cells = pval_states

# log2(fold change)
lfc_calc <- log2(fold_change)
lfc_calc <- replace(lfc_calc, is.infinite(lfc_calc),NA)
lfc_calc[is.na(lfc_calc)] <- 0

lfc_calc <- as.data.frame(lfc_calc)
select_md <- md[,c("patient","response")]
lfc_calc$patient <- rownames(lfc_calc)
select_md <- unique(select_md)
lfc_calc_all <- merge(lfc_calc,select_md,by="patient")

lfc_calc_all <- lfc_calc_all[order(lfc_calc_all$response, decreasing = T),]
rownames(lfc_calc_all) <- lfc_calc_all$patient
lfc_calc_all <- lfc_calc_all[,c("patient"    ,  "Cytotoxic"   , "Dysf GZMB"  ,  "Dysf ZNF683"   ,"Naive-like CD4","Naive-like CD8" ,
                        "Tfh LAG3"   ,  "Tfh NR3C1"   ,"Transitional",   "Treg GIMAP"  , "Treg TNFRSF9","Non-classical T cell" ,"NK KLRC1","NK FGFBP2", "response")  ]

library(colorspace)
cols1 <- diverging_hcl(palette="Berlin",n=255)

library(ComplexHeatmap)
to_plot = lfc_calc_all[,2:14]
col_ann =rowAnnotation(response=lfc_calc_all$response, col=list(response=c("RE"= "#339933","NR"="#DB4437")))

# set breaks
to_plot[to_plot > 3] = 3
to_plot[to_plot < -3] = -3

pdf(file=paste0(fig_path,"/Fold_change_heatmap_CD4CD8NK_v3_fix.pdf"), width=4.5, height=5.5)
out_h1 = Heatmap(as.matrix(to_plot),column_title="CD8/CD4/NK fraction fold change", name="zscore", left_annotation = col_ann,
                 row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
                 cluster_rows=F, cluster_columns=F,col = cols1)
out_h1
dev.off()

## median fold change, and wilcox test
lfc_calc_all_RE <- lfc_calc_all[lfc_calc_all$response == "RE",]
lfc_calc_all_NR <- lfc_calc_all[lfc_calc_all$response == "NR",]
median_fc_RE <- robustbase::colMedians(as.matrix(lfc_calc_all_RE[,2:(ncol(lfc_calc_all_RE)-1)]))
median_fc_NR <- robustbase::colMedians(as.matrix(lfc_calc_all_NR[,2:(ncol(lfc_calc_all_NR)-1)]))
median_fc <- rbind(median_fc_RE,median_fc_NR)

library(RColorBrewer)
cols2 <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

median_fc <- as.data.frame(t(median_fc))

cells_order <- c(  "Cytotoxic"   , "Dysf GZMB"  ,  "Dysf ZNF683"   ,"Naive-like CD4","Naive-like CD8" ,
                  "Tfh LAG3"   ,  "Tfh NR3C1"   ,"Transitional",   "Treg GIMAP"  , "Treg TNFRSF9","Non-classical T cell" ,"NK KLRC1","NK FGFBP2")

median_fc <- median_fc[match(cells_order,rownames(median_fc)),]

pdf(file=paste0(fig_path,"/Median_old_change_heatmap_CD4CD8NK_v3.pdf"), width=4.7, height=2.8)
out_h2 = Heatmap(t(median_fc),column_title="CD8/CD4/NK REvsNR", name="median",
                 row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
                 cluster_rows=F, cluster_columns=F,col = cols2)
out_h2
dev.off()

pdf(file=paste0(fig_path,"/Foldchange_Median_heatmap_CD4CD8NK_v3.pdf"), width=6.5, height=6)
ht_list = out_h1 %v% out_h2
draw(ht_list, column_km = 1)
dev.off()


## M, B non-T cell
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md <- mat@cell_metadata[names(mc@mc),]
mc_groups <- unique(md$mc_group)

pre <- md[md$condition == "pre" & 
           md$cell_type %in% c("M","B") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12",]
pre <- pre[pre$mc_group %in% mc_groups,]

post <- md[md$condition == "post" & 
            md$cell_type %in% c("M","B") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12",]
post <- post[post$mc_group %in% mc_groups,]

pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

post_clean <- post_frac[,c("patient", "mc_group", "fraction")]
post_final <- spread(post_clean, key = mc_group, value = fraction)
post_final <- post_final %>% tibble::column_to_rownames("patient")
post_final[is.na(post_final)] <- 0

states <- colnames(pre_final)
samples <- rownames(pre_final)
fold_change <- post_final[samples, states]/pre_final[samples, states]
fold_change <- as.matrix(fold_change)

select_md <- md[,c("patient","response")]
pval_states <- as.data.frame(states)
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

save_pvals_non_t_cells <- pval_states

# log2(fold change)
lfc_calc <- log2(fold_change)
lfc_calc <- replace(lfc_calc, is.infinite(lfc_calc),NA)
lfc_calc[is.na(lfc_calc)] <- 0

lfc_calc <- as.data.frame(lfc_calc)
select_md <- md[,c("patient","response")]
lfc_calc$patient <- rownames(lfc_calc)
select_md <- unique(select_md)
lfc_calc_all <- merge(lfc_calc,select_md,by="patient")

lfc_calc_all <- lfc_calc_all[order(lfc_calc_all$response, decreasing = T),]

rownames(lfc_calc_all) <- lfc_calc_all$patient

lfc_calc_all <- lfc_calc_all[,c("patient"    ,  "B cell CD27", "B cell CD69",   "Plasma cell" , "cDC CD1C"    , "cDC CLEC9A"  , "cDC LAMP3"   ,
                        "Granulocyte",  "Mono-macro"  , "pDC"      , "response" )  ]
library(colorspace)
cols1 <- diverging_hcl(palette="Berlin",n=255)

library(ComplexHeatmap)
to_plot = lfc_calc_all[,2:10]
col_ann =rowAnnotation(response=lfc_calc_all$response, col=list(response=c("RE"= "#339933","NR"="#DB4437")))

# set breaks
to_plot[to_plot > 3] = 3
to_plot[to_plot < -3] = -3

pdf(file=paste0(fig_path,"/Fold_change_heatmap_MB_v3_fix.pdf"), width=4.5, height=5.5)
out_h1 = Heatmap(as.matrix(to_plot),column_title="M/B fraction fold change", name="zscore", left_annotation = col_ann,
                 row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
                 cluster_rows=F, cluster_columns=F,col = cols1)
out_h1
dev.off()

## median fold change, and wilcox test
test_all_RE <- lfc_calc_all[lfc_calc_all$response == "RE",]
test_all_NR <- lfc_calc_all[lfc_calc_all$response == "NR",]
median_fc_RE <- colMeans(as.matrix(test_all_RE[,2:(ncol(test_all_RE)-1)]))
median_fc_NR <- colMeans(as.matrix(test_all_NR[,2:(ncol(test_all_NR)-1)]))
median_fc <- rbind(median_fc_RE,median_fc_NR)
median_fc <- as.data.frame(t(median_fc))

cols<- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
cells_order <- c("B cell CD27", "B cell CD69",   "Plasma cell" , "cDC CD1C"    , "cDC CLEC9A"  , "cDC LAMP3"   ,
                 "Granulocyte",  "Mono-macro"  , "pDC" )
median_fc <- median_fc[match(cells_order,rownames(median_fc)),]

pdf(file=paste0(fig_path,"/Median_old_change_heatmap_MB_v3.pdf"), width=4.7, height=2.8)
out_h2 = Heatmap(t(median_fc),column_title="M/B REvsNR", name="median",
                 row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
                 cluster_rows=F, cluster_columns=F,col = cols2)
out_h2
dev.off()

pdf(file=paste0(fig_path,"/Foldchange_Median_heatmap_MB_v3.pdf"), width=5.5, height=6)
ht_list = out_h1 %v% out_h2
draw(ht_list, column_km = 1)
dev.off()



#### Figures S4C, S4E, S4G ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
mc_groups = unique(md$mc_group)

pre <- md[md$condition == "pre" & 
           md$cell_type %in% c("CD4","CD8","NK") & 
           !is.na(md$patient) &
           md$patient != "Pat38" &
           md$patient != "Pat04" &
           md$patient != "Pat12",]
pre <- pre[pre$mc_group %in% mc_groups,]

post <- md[md$condition == "post" & 
            md$cell_type %in% c("CD4","CD8","NK") & 
            !is.na(md$patient) &
            md$patient != "Pat38" &
            md$patient != "Pat04" &
            md$patient != "Pat12",]
post <- post[post$mc_group %in% mc_groups,]

pre_frac <- pre %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

pre_clean <- pre_frac[,c("patient", "mc_group", "fraction")]
pre_final <- spread(pre_clean, key = mc_group, value = fraction)
pre_final <- pre_final %>% tibble::column_to_rownames("patient")
pre_final[is.na(pre_final)] <- 0

post_frac <- post %>% 
  group_by(patient, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

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
post_frac_c <- post_frac[,c("patient","mc_group","n")]
pre_frac_c <- pre_frac[,c("patient","mc_group","n")]

pre_frac_c$condition <- rep("pre",nrow(pre_frac_c))
post_frac_c$condition <- rep("post",nrow(post_frac_c))

prepost_frac_c <- rbind(pre_frac_c,post_frac_c)

for (cell_name in c("Treg GIMAP","Treg TNFRSF9","Dysf GZMB","Dysf ZNF683","Transitional","Naive-like CD8",
                    "Tfh LAG3","Tfh NR3C1","Naive-like CD4","NK FGFBP2", "NK KLRC1")){
  cell_x = cell_name
  print(cell_name)
  cell_x_counts <- prepost_frac_c[prepost_frac_c$mc_group == cell_x,]
  cell_x_counts$pat_cond <- paste0(cell_x_counts$patient,cell_x_counts$condition)
  prepost_final_RE$pat_cond <- paste0(prepost_final_RE$patient,prepost_final_RE$condition)
  prepost_final_RE_plot <- merge(prepost_final_RE,cell_x_counts,by= "pat_cond")
  pdf(file=paste0(fig_path,"/",cell_x,"_post_pre_RE.pdf"), width=2.5, height=2.3) #
  p = ggplot(prepost_final_RE_plot, aes(x=condition.x, y=prepost_final_RE_plot[[cell_x]], fill = condition.x)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group = condition.x, size = n, fill=condition.x), pch = 21, alpha = 0.5,
               position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
    geom_line(aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    ggtitle(paste0(cell_x, " ","RE")) +
    scale_fill_manual(values=c("#89CC9B","#116936")) +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320,640),
               labels = c(10,20,40,80,160,320,640),range=c(1,4))+
    ylab("fraction")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  print(wilcox.test(prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "pre",][[cell_x]],
              prepost_final_RE_plot[prepost_final_RE_plot$condition.y == "post",][[cell_x]]),paired=T)
  
  prepost_final_NR$pat_cond <- paste0(prepost_final_NR$patient,prepost_final_NR$condition)
  prepost_final_NR_plot <- merge(prepost_final_NR,cell_x_counts,by= "pat_cond")
  pdf(file=paste0(fig_path,"/",cell_x,"_post_pre_NR.pdf"), width=2.5, height=2.3) #
  p = ggplot(prepost_final_NR_plot, aes(x=condition.x, y=prepost_final_NR_plot[[cell_x]], fill = condition.x)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group = condition.x, size = n, fill=condition.x), pch = 21, alpha = 0.5,
               position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0))+ 
    geom_line(aes(group = patient.x), color = "lightgrey", alpha = 0.5)+
    ggtitle(paste0(cell_x, " ","NR")) +
    scale_fill_manual(values=c("#DC7985","#9B1F23")) +
    scale_size(name   = "# cells",
               breaks = c(10,20,40,80,160,320,640),
               labels = c(10,20,40,80,160,320,640),range=c(1,4))+
    ylab("fraction")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  print(wilcox.test(prepost_final_NR_plot[prepost_final_NR_plot$condition.y == "pre",][[cell_x]],
              prepost_final_NR_plot[prepost_final_NR_plot$condition.y == "post",][[cell_x]]),paired=T)
}


#### Figure S4D ####
## pre ratio versus on ratio in Tregs
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md <- mat@cell_metadata[names(mc@mc),]
md <- md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$patient != "Pat38" &
          md$patient != "Pat04" &
          md$patient != "Pat12" &
          md$cell_type %in% c("CD4"),]

data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

x <- data %>% 
  group_by(patient, response,mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

y = x[x$mc_group %in% c("Treg GIMAP", "Treg TNFRSF9"), ]

plot_treg_pre = pivot_wider(data = y, 
                            id_cols = c("patient", "response"), 
                            names_from = mc_group, 
                            values_from = c("n", "total", "fraction"))

plot_treg_pre$ratio_pre = log2(plot_treg_pre$`fraction_Treg TNFRSF9`/plot_treg_pre$`fraction_Treg GIMAP`)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "post" &
          md$patient != "Pat38" &
          md$patient != "Pat04" &
          md$patient != "Pat12" &
          md$cell_type %in% c("CD4"),]

data = md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

x <- data %>% 
  group_by(patient, response,mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

y = x[x$mc_group %in% c("Treg GIMAP", "Treg TNFRSF9"), ]

plot_treg_post = pivot_wider(data = y, 
                             id_cols = c("patient", "response"), 
                             names_from = mc_group, 
                             values_from = c("n", "total", "fraction"))

plot_treg_post$ratio_post = log2(plot_treg_post$`fraction_Treg TNFRSF9`/plot_treg_post$`fraction_Treg GIMAP`)

to_plot = plyr::join(plot_treg_post,plot_treg_pre, by="patient")

pdf(file=paste0(fig_path,"/ratio_pre_post_Treg.pdf"), width=3.5, height=2.5)
ggplot(to_plot, aes(x=ratio_pre, y=ratio_post,color=response...2)) +
  geom_point() + 
  ggtitle("Treg_ratio_TNFRSF9/GIMAP_PrevsPost") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#DC7985","#8BCB9B")) + 
  ylim(c(-4,6.5)) + xlim(c(-4,6.5)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0)
dev.off()

# Dysf
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$patient != "Pat38" &
          md$patient != "Pat04" &
          md$patient != "Pat12" &
          md$cell_type %in% c("CD8"),]

data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

x <- data %>% 
  group_by(patient, response,mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

y <- x[x$mc_group %in% c("Dysf GZMB", "Dysf ZNF683"), ]

plot_treg_pre <- pivot_wider(data = y, 
                            id_cols = c("patient", "response"), 
                            names_from = mc_group, 
                            values_from = c("n", "total", "fraction"))

plot_treg_pre$ratio_pre <- log2(plot_treg_pre$`fraction_Dysf GZMB`/plot_treg_pre$`fraction_Dysf ZNF683`)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "post" &
          md$patient != "Pat38" &
          md$patient != "Pat04" &
          md$patient != "Pat12" &
          md$cell_type %in% c("CD8"),]

data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

x <- data %>% 
  group_by(patient, response,mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n))

y = x[x$mc_group %in% c("Dysf GZMB", "Dysf ZNF683"), ]

plot_treg_post <- pivot_wider(data = y, 
                             id_cols = c("patient", "response"), 
                             names_from = mc_group, 
                             values_from = c("n", "total", "fraction"))

plot_treg_post$ratio_post <- log2(plot_treg_post$`fraction_Dysf GZMB`/plot_treg_post$`fraction_Dysf ZNF683`)
to_plot <- plyr::join(plot_treg_post,plot_treg_pre, by="patient")

pdf(file=paste0(fig_path,"/ratio_pre_post_Dysf.pdf"), width=3.5, height=2.5)
ggplot(to_plot, aes(x=ratio_pre, y=ratio_post,color=response...2)) +
  geom_point() + 
  ggtitle("Dysf_ratio_GZMB/ZNF683_PrevsPost") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#DC7985","#8BCB9B")) + 
  ylim(c(-4,6.5)) + xlim(c(-4,6.5)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0)
dev.off()



#### Figure S4F S4H ####
# NK signature volcano NR
NK_sigs = generate_NKcell_signatures()

mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

nk_1 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK KLRC1","NK FGFBP2") &
                                     mat@cell_metadata$condition == "post" &
                                     mat@cell_metadata$response == "NR",])

nk_2 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK KLRC1","NK FGFBP2") &
                                     mat@cell_metadata$condition == "pre" &
                                     mat@cell_metadata$response == "NR",])

nk_1_df <- mat@cell_metadata[intersect(nk_1, colnames(mat_ds)), ]
nk_2_df <- mat@cell_metadata[intersect(nk_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(nk_1_df,nk_2_df)

oi_genes = NK_sigs$active_mature

pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb <- pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Activation_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "Activation_signature"

oi_genes <- NK_sigs$CD56_bright
levels(pb$sig) <- c(levels(pb$sig), "CD56bright_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "CD56bright_signature"
pdf(file=paste0(fig_path,"/Volcano_DEG_byhand_NK_NR_postvspre_active.pdf"), width=3, height=3.6) #
volc = ggplot(pb, aes(fc, -log10(pval))) +
  geom_point(aes(col=as.factor(sig)),stroke=0.01) +
  scale_color_manual(values=c("black", "#9e9e9e","red","orange")) #508abf
volc+geom_text_repel(data=head(pb, 50), aes(label=pb$gene[1:50]),size=2,max.overlaps = getOption("ggrepel.max.overlaps", default = 12))+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5,linetype="dotted") +
  geom_vline(xintercept = -1.5,linetype="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype="dotted") +   
  labs(title  = paste0("NK cells, post vs pre, NR"),
       x = paste0("log2 fold change") ,
       y = paste0("-log10(p value)")) +
  geom_point(data=pb[pb$sig == "CD56bright_signature",],aes(x=fc, y= -log10(pval)),color="orange",stroke=0.01) +
  geom_point(data=pb[pb$sig == "Activation_signature",],aes(x=fc, y= -log10(pval)),color="red",stroke=0.01)
dev.off()

# RE
mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

nk_1 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK KLRC1","NK FGFBP2") &
                                     mat@cell_metadata$condition == "post" &
                                     mat@cell_metadata$response == "RE",])

nk_2 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK KLRC1","NK FGFBP2") &
                                     mat@cell_metadata$condition == "pre" &
                                     mat@cell_metadata$response == "RE",])

nk_1_df <- mat@cell_metadata[intersect(nk_1, colnames(mat_ds)), ]
nk_2_df <- mat@cell_metadata[intersect(nk_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(nk_1_df,nk_2_df)

oi_genes <- NK_sigs$active_mature

pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb <- pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Activation_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "Activation_signature"

oi_genes <- NK_sigs$CD56_bright
levels(pb$sig) <- c(levels(pb$sig), "CD56bright_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "CD56bright_signature"
pdf(file=paste0(fig_path,"/Volcano_DEG_byhand_NK_RE_postvspre_active.pdf"), width=3, height=3.6) #
volc = ggplot(pb, aes(fc, -log10(pval))) +
  geom_point(aes(col=as.factor(sig)),stroke=0.01) +
  scale_color_manual(values=c("black", "#9e9e9e","red","orange")) #508abf
volc+geom_text_repel(data=head(pb, 50), aes(label=pb$gene[1:50]),size=2,max.overlaps = getOption("ggrepel.max.overlaps", default = 12))+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5,linetype="dotted") +
  geom_vline(xintercept = -1.5,linetype="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype="dotted") +   
  labs(title  = paste0("NK cells, post vs pre, RE"),
       x = paste0("log2 fold change") ,
       y = paste0("-log10(p value)")) +
  geom_point(data=pb[pb$sig == "CD56bright_signature",],aes(x=fc, y= -log10(pval)),color="orange",stroke=0.01) +
  geom_point(data=pb[pb$sig == "Activation_signature",],aes(x=fc, y= -log10(pval)),color="red",stroke=0.01)
dev.off()

## NK signature volcano NK vs NK
mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

nk_1 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK KLRC1"),])

nk_2 <- rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("NK FGFBP2") ,])

nk_1_df <- mat@cell_metadata[intersect(nk_1, colnames(mat_ds)), ]
nk_2_df <- mat@cell_metadata[intersect(nk_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(nk_1_df,nk_2_df)

oi_genes <- NK_sigs$active_mature
pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb = pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Activation_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "Activation_signature"

oi_genes <- NK_sigs$CD56_bright
levels(pb$sig) <- c(levels(pb$sig), "CD56bright_signature")
pb[pb$gene %in% oi_genes,]["sig"] <- "CD56bright_signature"
pdf(file=paste0(fig_path,"/Volcano_DEG_byhand_NK_KLRC1_vs_FGFBP2_w.pdf"), width=3, height=3.6) #
volc = ggplot(pb, aes(fc, -log10(pval))) +
  geom_point(aes(col=as.factor(sig)),stroke=0.01) +
  scale_color_manual(values=c("black", "#9e9e9e","red","orange")) #508abf
volc+geom_text_repel(data=head(pb, 50), aes(label=pb$gene[1:50]),size=2,max.overlaps = getOption("ggrepel.max.overlaps", default = 12))+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5,linetype="dotted") +
  geom_vline(xintercept = -1.5,linetype="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype="dotted") +   
  labs(title  = paste0("NK cells, KLRC1 vs FGFBP2"),
       x = paste0("log2 fold change") ,
       y = paste0("-log10(p value)")) +
  geom_point(data=pb[pb$sig == "CD56bright_signature",],aes(x=fc, y= -log10(pval)),color="orange",stroke=0.01) +
  geom_point(data=pb[pb$sig == "Activation_signature",],aes(x=fc, y= -log10(pval)),color="red",stroke=0.01)
dev.off()



#### Figure S4I, S4J ####
# Transitional fractions TCR clones cell states, barplots
# Transitional clone fraction changes, strip plot
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
md_clones_top_id <- md_clones_top$pat_clone_id # top clones

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
sum_md_toppers_post_pre_melt <- melt(sum_md_toppers_post_pre)

sum_md_fraction_change_NR$clone_id <- rownames(sum_md_fraction_change_NR)
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[match(order_sort,sum_md_fraction_change_NR$clone_id),]
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[1:15,][,-7]
sum_md_fraction_change_NR_melt <- melt(sum_md_fraction_change_NR)

sum_md_fraction_change_NR_melt$clone_id <- factor(sum_md_fraction_change_NR_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-7])
sum_md_fraction_change_NR_melt$variable <- factor(sum_md_fraction_change_NR_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-7])


## barplots
to_plot_clones <- sum_md_fraction_change_NR_melt[sum_md_fraction_change_NR_melt$variable == "Transitional",]$clone_id

sum_md_toppers_pre$sum <- rowSums(sum_md_toppers_pre)
sum_md_toppers_post$sum <- rowSums(sum_md_toppers_post)

sum_md_toppers_pre_sel <- sum_md_toppers_pre[rownames(sum_md_toppers_pre) %in% to_plot_clones,]
sum_md_toppers_post_sel <- sum_md_toppers_post[rownames(sum_md_toppers_post) %in% to_plot_clones,]

sum_md_toppers_pre_sel$id <- rownames(sum_md_toppers_pre_sel)
sum_md_toppers_post_sel$id <- rownames(sum_md_toppers_post_sel)
m_pre <- reshape2::melt(sum_md_toppers_pre_sel[,-(ncol(sum_md_toppers_pre_sel)-1)],)
m_post <- reshape2::melt(sum_md_toppers_post_sel[,-(ncol(sum_md_toppers_post_sel)-1)],)

pdf(file=paste0(fig_path,"/Fractions_change_clones_onlytransitional_RE_post.pdf"),width=5.5, height=2.2)
ggplot(m_post, aes(x = id, y = value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#e5d3b3","#4b1d5b","#c01a88","#f37d7f","#fbcccb","#c9a3c8"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file=paste0(fig_path,"/Fractions_change_clones_onlytransitional_RE_pre.pdf"),width=5.5, height=2.2)
ggplot(m_pre, aes(x = id, y = value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#e5d3b3","#4b1d5b","#c01a88","#f37d7f","#fbcccb","#c9a3c8"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## line chart
sum_md_toppers_pre_sel$condition <- rep("pre",nrow(sum_md_toppers_pre_sel))
sum_md_toppers_post_sel$condition <- rep("post",nrow(sum_md_toppers_post_sel))
merged_plot <- rbind(sum_md_toppers_pre_sel,sum_md_toppers_post_sel)

merged_plot$condition <- factor(merged_plot$condition, levels=c("pre","post"))

merged_plot <- merged_plot[merged_plot$id %in% to_plot_clones,]

pdf(file=paste0(fig_path,"/Absolute_change_clones_onlytransitional_RE_post.pdf"),width=2.5, height=3)
ggplot(merged_plot, aes(x=condition, y=Transitional, fill = condition,group = condition)) +
  geom_point(pch = 21, alpha = 0.5) +
  geom_line(data=merged_plot,aes(group = id), color = "grey", alpha = 0.5)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()  

# non-responders
remove = c("Pat04", "Pat12", "Pat38") # because of low cell nrs in post sample, POST
#remove = c("None")
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          !md$patient %in% remove &
          md$response %in% "NR" &
          !is.na(md$unique_clone_ID) &
          md$mc_group %in% c("Transitional"),]

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
sum_md_toppers_post_pre_melt <- melt(sum_md_toppers_post_pre)

sum_md_fraction_change_NR$clone_id <- rownames(sum_md_fraction_change_NR)
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[match(order_sort,sum_md_fraction_change_NR$clone_id),]
sum_md_fraction_change_NR <- sum_md_fraction_change_NR[1:9,][,-5]
sum_md_fraction_change_NR_melt <- melt(sum_md_fraction_change_NR)

sum_md_fraction_change_NR_melt$clone_id <- factor(sum_md_fraction_change_NR_melt$clone_id, levels=rev(order_sort))
sum_md_toppers_post_pre_melt$clone_id <- factor(sum_md_toppers_post_pre_melt$clone_id, levels=rev(order_sort))

sum_md_toppers_post_pre_melt$variable <- factor(sum_md_toppers_post_pre_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-5])
sum_md_fraction_change_NR_melt$variable <- factor(sum_md_fraction_change_NR_melt$variable,levels=colnames(sum_md_toppers_pre_sw_NR)[-5])

## barplots
to_plot_clones= sum_md_fraction_change_NR_melt[sum_md_fraction_change_NR_melt$variable == "Transitional",]$clone_id

sum_md_toppers_pre$sum <- rowSums(sum_md_toppers_pre)
sum_md_toppers_post$sum <- rowSums(sum_md_toppers_post)

sum_md_toppers_pre_sel = sum_md_toppers_pre[rownames(sum_md_toppers_pre) %in% to_plot_clones,]
sum_md_toppers_post_sel = sum_md_toppers_post[rownames(sum_md_toppers_post) %in% to_plot_clones,]

sum_md_toppers_pre_sel$id = rownames(sum_md_toppers_pre_sel)
sum_md_toppers_post_sel$id = rownames(sum_md_toppers_post_sel)
m_pre = reshape2::melt(sum_md_toppers_pre_sel[,-(ncol(sum_md_toppers_pre_sel)-1)],)
m_post = reshape2::melt(sum_md_toppers_post_sel[,-(ncol(sum_md_toppers_post_sel)-1)],)

pdf(file=paste0(fig_path,"/Fractions_change_clones_onlytransitional_NR_post.pdf"),width=4.3, height=2.3)
ggplot(m_post, aes(x = id, y = value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#e5d3b3","#4b1d5b","#fbcccb","#c9a3c8","#f37d7f","#c01a88"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file=paste0(fig_path,"/Fractions_change_clones_onlytransitional_NR_pre.pdf"),width=4.3, height=2.3)
ggplot(m_pre, aes(x = id, y = value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#e5d3b3","#4b1d5b","#fbcccb","#c9a3c8","#f37d7f","#c01a88"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


## lineplot
sum_md_toppers_pre_sel$condition = rep("pre",nrow(sum_md_toppers_pre_sel))
sum_md_toppers_post_sel$condition = rep("post",nrow(sum_md_toppers_post_sel))
merged_plot = rbind(sum_md_toppers_pre_sel,sum_md_toppers_post_sel)

merged_plot$condition = factor(merged_plot$condition, levels=c("pre","post"))

pdf(file=paste0(fig_path,"/Absolute_change_clones_onlytransitional_NR_post.pdf"),width=2.5, height=3)
ggplot(merged_plot, aes(x=condition, y=Transitional, fill = condition,group = condition)) +
  geom_point(pch = 21, alpha = 0.5) +
  geom_line(data=merged_plot,aes(group = id), color = "grey", alpha = 0.5)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()  




