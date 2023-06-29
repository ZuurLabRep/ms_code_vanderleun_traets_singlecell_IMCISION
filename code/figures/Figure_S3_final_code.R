###################
#### FIGURE S3 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S3"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S3")

library(circlize)
library(ComplexHeatmap)
library(AUCell)


## first run 1.kurten_metacell_generation.R

#### Figure S3A ####
## plot original annotations on top
mat = scdb_mat(clean_mat_id_ferris)
mc = scdb_mc(filt_mc_id_ferris)
supmc_tab = read.table("./data/CD45_tumor_supmc_Kurten.txt",sep="\t")

mc2d <- scdb_mc2d(mc2d_id_ferris)
xco <- as.data.frame(mc2d@sc_x)
xco$cellnms <- rownames(xco)
yco <- as.data.frame(mc2d@sc_y)
yco$cellnms <- rownames(yco)
coo <- merge(xco, yco, by="cellnms")
coo <- na.omit(coo)
colnames(coo)[2:3] <- c("sc_x","sc_y")

cells <- rownames(mat@cell_metadata)

cells <- cells[cells %in% names(mc@mc)]
cell_codes <- intersect(cells, coo$cellnms)

coo_df <- as.data.frame(coo)
rownames(coo_df) <- coo_df$cellnms
coo_df <- coo_df[,-1]
coo_df <- coo_df[intersect(rownames(coo_df), cell_codes),]
coo_df[cell_codes, "orig.ident"] <- mat@cell_metadata[cell_codes, "orig.ident"]
coo_df[cell_codes, "mc_group"] <- mat@cell_metadata[cell_codes, "mc_group"]

supmc_tab = supmc_tab[order(supmc_tab$name),]
supmc_tab_f <- supmc_tab
supmc_tab_f[!(supmc_tab_f$name %in% c("cDC_CD1C","cDC_CLEC9A","cDC_LAMP3","Treg_GIMAP","Treg_TNF","pDC")),]$color <- "grey"
pdf(file=paste0(fig_path,"/2dprojection_hpvneg_mc_group_v2.pdf"), width=5.5, height=4)
p = ggplot()+
  geom_point(data = coo_df, aes(x=sc_x, y = sc_y, fill = mc_group),colour = "black",pch=21, size=1.4,stroke = 0.08)+
  guides(colour = guide_legend(override.aes = list(size=5)))+ scale_fill_manual(values=c(supmc_tab_f$color,"grey")) +
  theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
print(p)
dev.off()


#### Figure S3B ####
mc = scdb_mc(filt_mc_id)
submat = md
umis = mat@mat
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 
for (gene_oi in c("CLEC9A","FOXP3","CD1C","TNFRSF9","IL3RA","LAMP3")){
  genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])
  
  genes_cells$gene_exp <- rowSums(genes_cells)
  genes_cells$lognorm = log2(genes_cells$gene_exp + 1)
  
  genes_cells$cell_id <- rownames(genes_cells)
  coo_df$cell_id <- rownames(coo_df)
  gene_scores_all <- plyr::join(genes_cells,coo_df,by="cell_id")
  
  png(file=paste0(fig_path,"/2dprojection_gene_",gene_oi,".png"),width=1000, height=1000)
  p = ggplot()+geom_point(data = gene_scores_all[gene_scores_all$lognorm == 0,], aes(x=sc_x, y = sc_y),fill = "lightgrey",pch=21, size=6,stroke = 0,color='transparent') +
    geom_point(data = gene_scores_all[gene_scores_all$lognorm !=0,], aes(x=sc_x, y = sc_y, fill = lognorm),color='transparent',pch=21, size=6,stroke = 0) +
    scale_fill_gradient2(low="#ffe7ca", mid="#f48622", high="#3f2800",midpoint = max(gene_scores_all$lognorm)/2, 
                         limits = c(0, max(gene_scores_all$lognorm)))+ ggtitle(gene_oi) +
    theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()
}


#### Figure S3C ####
sig_tregs = read.table("./data/active_treg_vs_less_activate_treg_sig_DEG_signature.txt",header=T)
treg_active = sig_tregs[sig_tregs$fc > 1,]$gene
geneSets <- list(treg_active=treg_active)

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Treg_GIMAP","Treg_TNF"),]

exprMatrix <- as.matrix(mat@mat)
exprMatrix <- exprMatrix[,colnames(exprMatrix) %in% rownames(md)]
rownames(exprMatrix) <- geneNames

dim(exprMatrix)
cell_rank <- AUCell_buildRankings(exprMatrix, nCores=10, plotStats = F)
cells_AUC <- AUCell_calcAUC(geneSets,cell_rank)

check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
coo_df$cell_id <- rownames(coo_df)
sig_scores_all <- plyr::join(check_data,coo_df,by="cell_id")

pdf(file=paste0(fig_path,"/2dprojection_subset_Treg_pot_treg_active_violin.pdf"),width=3.7, height=3.5)
ggplot(sig_scores_all,aes(y=treg_active,x=mc_group,fill=mc_group))+ geom_jitter(width=0.1)+geom_violin()  +
  scale_fill_manual(values=c("#0042b0","#59b56f")) + ggtitle("treg_active") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggpubr::stat_compare_means(label = "p.format")
dev.off()


#### Figure S3D ####
## Tregs per patient
non_tcell = c("cDC_CD1C","cDC_LAMP3","cDC_CLEC9A","pDC","cell13","Undet","cell12")

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!(md$mc_group %in% non_tcell),]
data = md

# calculate fractions
x <- data %>% 
  group_by(orig.ident, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

# plot boxplots
fractions_x <- as.data.frame(x)
fractions_x_sub <- fractions_x[fractions_x$mc_group %in% c("Treg_GIMAP","Treg_TNF"),]
fractions_x_sub$orig.ident <- do.call(rbind, strsplit(fractions_x_sub$orig.ident,"_"))[,2]

treg_frac <- fractions_x_sub[fractions_x_sub$mc_group == "Treg_TNF",]$fraction/fractions_x_sub[fractions_x_sub$mc_group == "Treg_GIMAP",]$fraction
to_plot_ferris <- as.data.frame(treg_frac)

# combine with IMCISION data
to_plot_imcision = read.table("./data/Treg_ratio_table_IMCISION.txt",sep="\t")
to_plot_imcision = to_plot_imcision[,c("ratio","patient")]
to_plot_imcision = to_plot_imcision[to_check$patient != "Pat04",]
colnames(to_plot_imcision) = c("treg_frac","patient")
to_plot_imcision$condition = "IMCISION"
to_plot_ferris$condition = "Kürten"
to_plot_ferris$patient = "none"
to_plot_both = rbind(to_plot_imcision,to_plot_ferris)
pdf(file=paste0(fig_path,"/Boxplot_Treg_ratios.pdf"),width=2.2, height=3.2)
ggplot(to_plot_both,aes(x=condition,y=treg_frac,fill=condition)) + geom_boxplot() + geom_point() +
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("ratio TregTNFRSF9/TregGIMAP") +
  scale_fill_manual(values=c("lightgrey","dimgrey"))
dev.off()


#### Figure S3E ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!(md$mc_group %in% c("pDC","cDC_LAMP3","cDC_CLEC9A","cDC_CD1C","cell13","Undet","cell12")) &
          !is.na( md$mc_group),]
data = md[,c("orig.ident", "amp_batch_id","mc_group")]

# caclulate fractions
x <- data %>% 
  group_by(orig.ident, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

## myeloids, non-T cells
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]

md = md[md$mc_group %in% c("pDC","cDC_LAMP3","cDC_CLEC9A","cDC_CD1C","cell13","Undet","cell12") &
          !is.na( md$mc_group),]

data = md[,c("orig.ident", "amp_batch_id","mc_group")]

# caclulate fractions
x_2 <- data %>% 
  group_by(orig.ident, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

x_all <- rbind(x,x_2)
to_plot <- dcast(data = x_all,formula = orig.ident~mc_group,value.var = "fraction")
to_plot[is.na(to_plot)] <- 0
to_plot <- to_plot[,-1]

col_fun = colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(file=paste0("/DATA/j.traets/Zuur_lab/IMCISION_scRNA_Anne/Revisions/Figures_External_data_analysis_Ferris/cDC_correlations_Ferris_Tregs_all.pdf"), width=4.6, height=4)
Heatmap(cor(to_plot[,c("pDC","Treg_GIMAP","Treg_TNF","cDC_LAMP3","cDC_CLEC9A","cDC_CD1C")],method="spearman"),col=col_fun)
dev.off()

# IMCISION
# correlations tregs vs dcs
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("CD4","CD8","NK") ,]

data = md[,c("patient", "response", "condition", "cell_type","mc_group")]

# caclulate fractions
x <- data %>% 
  group_by(patient,response, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

x <- x[x$mc_group %in% c("Treg GIMAP","Treg TNFRSF9"),]
x_tregs <- x
x_tregs = dcast(x_tregs,patient~mc_group,value.var = "fraction")

## myeloids, non-T cells
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("M","B"),]

# remove cells with TCR
md <- md[is.na(md$unique_clone_ID),]
data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

# caclulate fractions
x <- data %>% 
  group_by(patient,response, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

x_DCs <- x[x$mc_group %in% c("cDC LAMP3","cDC CD1C","cDC CLEC9A","pDC"),]
x_DCs <- dcast(x_DCs,patient~mc_group,value.var = "fraction")
to_plot <- merge(x_DCs,x_tregs,by="patient")

library(circlize)
col_fun = colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(file=paste0("./Figures_final_rev/correlations_tregs_dcs_all.pdf"), width=4.6, height=4)
Heatmap(cor(to_plot[,-1],method="spearman"),col = col_fun)
dev.off()



#### Figure S3F ####
mat_ds = scdb_mat("HN_full_dataset_ds")
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

treg_1 = rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Treg TNFRSF9") &
                                       mat@cell_metadata$condition == "pre",])
treg_2 = rownames(mat@cell_metadata[ mat@cell_metadata$mc_group %in% c("Treg GIMAP") &
                                       mat@cell_metadata$condition == "pre",])

treg_1_df <- mat@cell_metadata[intersect(treg_1, colnames(mat_ds)), ]
treg_2_df <- mat@cell_metadata[intersect(treg_2, colnames(mat_ds)), ]

pb <- DE_wilcoxon(treg_1_df,treg_2_df) # differentially expressed genes between group 1 and 2, based on wilcoxon

treg_genes = read.table("./data/signatures/Treg_activation_genes.txt",header=T)$x

# plot volcano
pb$sig <- as.factor(mutate(pb, sig=ifelse(pb$pval<0.05, "pval<0.05", "Not Sig"))[,7])
pb = pb[order(pb$pval),]
levels(pb$sig) <- c(levels(pb$sig), "Activation_signature")
pb[pb$gene %in% treg_genes,]["sig"] <- "Activation_signature"
pdf(file=paste0(fig_path,"/Volcano_DEG_byhand_Tregs_pre_wilcox.pdf"), width=3.8, height=4.5) #
volc = ggplot(pb, aes(fc, -log10(pval))) +
  geom_point(aes(col=as.factor(sig)),stroke=0.01) +
  scale_color_manual(values=c("black", "#9e9e9e","orange")) #508abf
volc+geom_text_repel(data=head(pb, 50), aes(label=pb$gene[1:50]),size=2,max.overlaps = getOption("ggrepel.max.overlaps", default = 12))+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5,linetype="dotted") +
  geom_vline(xintercept = -1.5,linetype="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype="dotted") +   
  labs(title  = paste0("Tregs_GIMAPvsTNFRSF9"),
       x = paste0("log2 fold change") ,
       y = paste0("-log10(p value)")) +
  geom_point(data=pb[pb$sig == "Activation_signature",][1:11,],aes(x=fc, y= -log10(pval)),color="orange")
dev.off()



#### Figure S3G-J ####
# gene expression 
gene_oi = c("TNFRSF9") 

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
mc = scdb_mc(filt_mc_id)
umis = mat@mat
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 

genes_cells <- as.data.frame(umis_n[which(rownames(umis_n) %in% gene_oi),])

genes_cells$gene_exp <- rowSums(genes_cells)
genes_cells$lognorm = log2(genes_cells$gene_exp + 1)

genes_cells$cell_id = rownames(genes_cells)
md = mat@cell_metadata[names(mc@mc),]
md$cell_id = rownames(md)
genes_cells_md = merge(genes_cells,md,by="cell_id")
genes_cells_md$mc_group_pat = paste0(genes_cells_md$mc_group,"_",genes_cells_md$patient,"_",genes_cells_md$condition)

# signature
treg_genes = read.table("./data/signatures/Treg_activation_genes.txt",header=T)$x

mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]

geneSets <- list(geneSet1=treg_genes)

mat = scdb_mat(clean_mat_id)
geneNames = rownames(mat@mat)
submat <- mat@cell_metadata       
md = submat[submat$mc_group %in% c("Treg TNFRSF9","Treg GIMAP") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("pre"),] 
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(geneSets,cell_rank)

check_data <- cells_AUC@assays@data@listData$AUC
check_data <- as.data.frame(t(check_data))

# plot per cell state
check_data$cell_id <- rownames(check_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(check_data,md,by="cell_id")

# TNFRSF9 expression
to_add <- genes_cells_md[,c("cell_id","lognorm")]
sig_scores_all_exp <- merge(sig_scores_all,to_add,by="cell_id")

sig_scores_all_exp$TNFRSF9_pos <- "no_expression"
sig_scores_all_exp[sig_scores_all_exp$lognorm > 0,]$TNFRSF9_pos <- "pos_TNFRSF9"
sig_scores_all_exp_sub <- sig_scores_all_exp[sig_scores_all_exp$mc_group == "Treg TNFRSF9",]
sig_scores_all_exp_sub$TNFRSF9_pos <- factor(sig_scores_all_exp_sub$TNFRSF9_pos,levels=c("pos_TNFRSF9","no_expression"))
pdf(file=paste0(fig_path,"/TNFRSF9_expression_histogram_TNF.pdf"),width=4.5, height=2.7)
ggplot(sig_scores_all_exp_sub,aes(x=geneSet1,fill=TNFRSF9_pos)) + geom_histogram( alpha=1,position="stack") +
  scale_fill_manual(values=c("red","grey")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("expression_TNFSRF9_in_TregTNF") + xlab("Activation_sig")
dev.off()

sig_scores_all_exp$TNFRSF9_pos <- "no_expression"
sig_scores_all_exp[sig_scores_all_exp$lognorm > 0,]$TNFRSF9_pos <- "pos_TNFRSF9"
sig_scores_all_exp_sub <- sig_scores_all_exp[sig_scores_all_exp$mc_group == "Treg GIMAP",]
sig_scores_all_exp_sub$TNFRSF9_pos <- factor(sig_scores_all_exp_sub$TNFRSF9_pos,levels=c("pos_TNFRSF9","no_expression"))
pdf(file=paste0(fig_path,"/TNFRSF9_expression_histogram_GIMAP.pdf"),width=4.5, height=2.7)
ggplot(sig_scores_all_exp_sub,aes(x=geneSet1,fill=TNFRSF9_pos)) + geom_histogram( alpha=1,position="stack") +
  scale_fill_manual(values=c("red","grey")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("expression_TNFSRF9_in_TregGIMAP") + xlab("Activation_sig")
dev.off()

# single gene expression vs signature
pdf(file=paste0(fig_path,"/TNFRSF9_expression_vs_activation_signature_jitter_TNF.pdf"),width=5, height=3.5)
ggplot(sig_scores_all_exp[sig_scores_all_exp$mc_group == "Treg TNFRSF9",],aes(x=lognorm,y=geneSet1,color=mc_group)) + geom_point(alpha=0.5,position = position_jitter(w = 0.05, h = 0)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Activation sig") + xlab("exp_TNFRSF9") +
  ggtitle("Treg TNFRSF9 baseline") +
  scale_color_manual(values=c("#59b56f"))
dev.off()

pdf(file=paste0(fig_path,"/TNFRSF9_expression_vs_activation_signature_jitter_GIMAP.pdf"),width=4.8, height=3.5)
ggplot(sig_scores_all_exp[sig_scores_all_exp$mc_group == "Treg GIMAP",],aes(x=lognorm,y=geneSet1,color=mc_group)) + geom_point(alpha=0.5,position = position_jitter(w = 0.05, h = 0)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Activation sig") + xlab("exp_TNFRSF9") +
  ggtitle("Treg GIMAP baseline") +
  scale_color_manual(values=c("#0042b0"))
dev.off()


# correlation plot 1, IHC
IHC_all = read.csv("./data/41BB_FoxP3_AnalysisSettings_Results_newannotations_JT.csv",header=T)
IHC_all$Patient_ID <- paste0("Pat",IHC_all$Patient_ID)
IHC_all$Patient_ID[1] <- "Pat04"

pdf(file=paste0(fig_path,"/ratio_single_cell_data_vs_ratio_Treg41BBpos_vs_Treg41BBneg_mm2.pdf"),width=2.7, height=2.7)
ggplot(IHC_all,aes(x=ratio_single_cell_data,y=ratio_Treg41BBpos_vs_Treg41BBneg_mm2)) + geom_point() + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("ratio_Treg41BBpos_vs_Treg41BBneg_mm2") + xlab("ratio_single_cell_data") +
  ggtitle("sc vs IHC excluded") +
  stat_smooth(method = "lm", col = "darkblue",se = FALSE)
dev.off()
cor.test(IHC_all$ratio_Treg41BBpos_vs_Treg41BBneg_mm2,IHC_all$ratio_single_cell_data,method="pearson")
summary(lm(data = IHC_all,formula=ratio_Treg41BBpos_vs_Treg41BBneg_mm2~ratio_single_cell_data))$r.squared

# correlation plot 2, expression
sig_scores_all_exp$TNFRSF9_pos <- "no_expression"
sig_scores_all_exp[sig_scores_all_exp$lognorm > 0,]$TNFRSF9_pos <- "pos_TNFRSF9"
freq_exp_patient <- sig_scores_all_exp %>% 
  group_by(patient,TNFRSF9_pos) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n), total=sum(n)) 

freq_exp_patient$Patient_ID = freq_exp_patient$patient
IHC_all_plot = merge(IHC_all,freq_exp_patient,by="Patient_ID")

IHC_all_plot_sub <- IHC_all_plot[IHC_all_plot$TNFRSF9_pos == "pos_TNFRSF9",]
IHC_all_plot_sub <- IHC_all_plot_sub[IHC_all_plot_sub$Overweeg_exclusie != "Yes",]
IHC_all_plot_sub$n_neg <- IHC_all_plot_sub$total-IHC_all_plot_sub$n
IHC_all_plot_sub$ratio_freq_pos <- IHC_all_plot_sub$n /IHC_all_plot_sub$n_neg
pdf(file=paste0(fig_path,"/Ratio_pos_TNFRSF9_expression_vs_ratio_IHC_scatter_excl.pdf"),width=2.7, height=2.7)
ggplot(IHC_all_plot_sub,aes(x=ratio_freq_pos,y=ratio_Treg41BBpos_vs_Treg41BBneg_mm2)) + geom_point() + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("ratio_Treg41BBpos_vs_Treg41BBneg_mm2") + xlab("ratio_TNFRSF9_pos_vs_neg_sc") +
  ggtitle("sc vs IHC excluded") +
  stat_smooth(method = "lm", col = "darkblue",se = FALSE)
dev.off()

cor.test(IHC_all_plot_sub$ratio_freq_pos,IHC_all_plot_sub$ratio_Treg41BBpos_vs_Treg41BBneg_mm2,method="pearson")
summary(lm(data = IHC_all_plot_sub,formula=ratio_freq_pos~ratio_Treg41BBpos_vs_Treg41BBneg_mm2))
