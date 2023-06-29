##################
#### FIGURE 2 ####
##################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_2"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_2")

library(AUCell)
library(GSEABase)
library(tidyverse)


## load data
# clean_mat_id: filtered count matrix
# filt_mc_id: info on filtered meta cell generation
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group),]



#### Figure 2A ####
# colors mc groups
mc_colors <- read.table("./Data/mc_colors.txt",header=T,sep="\t")

## two plots, 2D projection cells colored by cell state, baseline RE and NR
plot_coord_mc_colors_baseline_resp(mc2d_name=mc2d_id,mat_name=clean_mat_id,mc_name=filt_mc_id,mc_colors=mc_colors)



#### Figure 2B ####
# RE vs NR coefficient, per cell type, T cell and non-T cell
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[  md$condition == "pre" &
            md$cell_type %in% c("CD4","CD8","NK"),]
md = md[md$patient != "Pat04",] # exclude nivo patient

md$pat_response <- paste0(md$patient,md$response)
data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

fractions_x <- as.data.frame(x)
prep_table <- as.data.frame(cbind(unique(fractions_x$mc_group),rep("none",length(unique(fractions_x$mc_group))),rep("none",length(unique(fractions_x$mc_group)))))

counter = 1
for (cell_t in unique(prep_table$V1)){
  fractions_x_sub <- fractions_x[fractions_x$mc_group %in% c(cell_t),]
  fractions_x_sub$response <- factor(fractions_x_sub$response,levels=c("RE","NR"))
  model <- lm(fractions_x_sub$fraction ~ as.factor(fractions_x_sub$response),data=fractions_x_sub) # RE vs NR on baseline
  sub_m <- summary(model)
  
  slope_sub <- model$coefficients[[2]]
  R2 <- sub_m$r.squared
  test_p <- wilcox.test(fractions_x_sub[fractions_x_sub$response == "RE",][["fraction"]],
                        fractions_x_sub[fractions_x_sub$response == "NR",][["fraction"]])
  prep_table$V3[counter] <- test_p$p.value
  prep_table$V2[counter] <- ((-slope_sub)/abs(slope_sub))*R2 # direction based on linear model
  counter = counter+1
}

prep_table$V2 <- as.numeric(prep_table$V2)
prep_table$V3 <- as.numeric(prep_table$V3)
prep_table$Y <- rep(0,nrow(prep_table))
prep_table_save <- prep_table

## Predicitve coeff other compartments (B and M)
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("B","M"),] #
md = md[md$patient != "Pat04",]

# remove labeled cells with TCR
md <- md[is.na(md$unique_clone_ID),]
data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

fractions_x <- as.data.frame(x)
prep_table <- as.data.frame(cbind(unique(fractions_x$mc_group),rep("none",length(unique(fractions_x$mc_group))),rep("none",length(unique(fractions_x$mc_group)))))

counter = 1
for (cell_t in unique(prep_table$V1)){
  fractions_x_sub <- fractions_x[fractions_x$mc_group %in% c(cell_t),]
  fractions_x_sub$response <- factor(fractions_x_sub$response,levels=c("RE","NR"))
  model <- lm(fractions_x_sub$fraction ~ as.factor(fractions_x_sub$response),data=fractions_x_sub)
  sub_m <- summary(model)
  
  slope_sub <- model$coefficients[[2]]
  R2 <- sub_m$r.squared
  test_p <- wilcox.test(fractions_x_sub[fractions_x_sub$response == "RE",][["fraction"]],
                        fractions_x_sub[fractions_x_sub$response == "NR",][["fraction"]])
  prep_table$V3[counter] <- test_p$p.value
  prep_table$V2[counter] <- ((-slope_sub)/abs(slope_sub))*R2
  counter = counter+1
}

prep_table$V2 <- as.numeric(prep_table$V2)
prep_table$V3 <- as.numeric(prep_table$V3)
prep_table$Y <- rep(0,nrow(prep_table))

# combine T and non-T cells
prep_table_save2 <- rbind(prep_table_save,prep_table)
prep_table <- prep_table_save2 # overwrite

# add pvalue cutoffs and color cell states
prep_table$sig_p <- rep("p>0.5",nrow(prep_table))
prep_table[prep_table$V3 < 0.3,][["sig_p"]] <- "p<0.5"
prep_table[prep_table$V3 < 0.1,][["sig_p"]] <- "p<0.1"
prep_table[prep_table$V3 < 0.05,][["sig_p"]] <- "p<0.05"
#prep_table[prep_table$V3 < 0.01,][["sig_p"]] <- "p<0.01"

prep_table$sig_p <- factor(prep_table$sig_p,levels=c("p>0.5","p<0.5","p<0.1","p<0.05","p<0.01"))
colors_mc <- read.table("./data/mc_colors.txt",header=T)
prep_table$V1 <- gsub(" ", "_", prep_table$V1)

prep_table <- prep_table[order(prep_table$V2, prep_table$V3),]
prep_table <- prep_table[prep_table$V1 != "Granulocyte",] # population too small

colors_mc <- colors_mc[match(prep_table$V1,colors_mc$name),]
prep_table$V1 <- factor(prep_table$V1,levels=prep_table$V1)
prep_table <- prep_table[complete.cases(prep_table),]

colors_mc[colors_mc$name %in% prep_table[prep_table$V3 > 0.1,]$V1,]$color <- "grey"

prep_table_s <- prep_table[prep_table$V1 %in% c("pDC","cDC_LAMP3","cDC_CD1C","cDC_CLEC9A",
                                               "Plasma_cell","B_cell_CD69","B_cell_CD27","Mono-macro"),]
colors_mc_s = colors_mc[colors_mc$name %in% prep_table_s$V1,]
pdf(file=paste0(fig_path,"/Predictive_index_baseline_MB_ver.pdf"), width=3.8, height=2.1)
p1 = ggplot(prep_table_s,aes(x=V2,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= colors_mc_s$color) + 
  xlab("Predictive coefficient") +ylab("") +scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  xlim(c(-0.3,0.3))+ geom_segment(aes(y = V1, yend = V1,x = Y, xend = V2),size = 0.8) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p1)
dev.off()

prep_table_s = prep_table[!(prep_table$V1 %in% c("pDC","cDC_LAMP3","cDC_CD1C","cDC_CLEC9A",
                                                 "Plasma_cell","B_cell_CD69","B_cell_CD27","Mono-macro")),]
colors_mc_s = colors_mc[colors_mc$name %in% prep_table_s$V1,]
pdf(file=paste0(fig_path,"/Predictive_index_baseline_CD48_ver.pdf"), width=4.6, height=3.4)
p2 = ggplot(prep_table_s,aes(x=V2,y=V1,color=V1)) + geom_point(aes(size=sig_p)) +
  geom_vline(xintercept=0) + scale_color_manual(values= colors_mc_s$color) + 
  xlab("Predictive coefficient") +ylab("") +scale_size_manual(values=c(2,3.5,5.5,6.5,8)) + 
  xlim(c(-0.3,0.3))+ geom_segment(aes(y = V1, yend = V1,x = Y, xend = V2),size = 0.8) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

pdf(file=paste0(fig_path,"/Predictive_index_baseline_CD48_MB_pvalue.pdf"), width=4, height=6,onefile=FALSE)
ggpubr::ggarrange(p2,p1,ncol=1,nrow=2,heights=c(1.5,1),common.legend = TRUE,legend="right")
dev.off()




#### Figure 2C ####
## boxplots per cells oi
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("CD4","CD8","NK") ,]
md = md[md$patient != "Pat04",]

data <- md[,c("patient", "response", "condition", "cell_type","mc_group")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

# plot boxplots
fractions_x <- as.data.frame(x)
fractions_x_sub <- fractions_x[fractions_x$mc_group %in% c("Treg TNFRSF9","Treg GIMAP"),]
pdf(file=paste0(fig_path,"/Pretreatment_Treg_fractions_by_response_CD48NK.pdf"), width=3, height=2.9)
p = ggplot(fractions_x_sub, aes(x=mc_group, y=fraction)) +
  geom_boxplot(aes(fill = response), outlier.shape = NA)  +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0), 
             aes(group = response, size = count,fill=response), pch = 21)+
  ggtitle("Tregs baseline") + ylim(0,0.4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(name   = "# cells",
             breaks = c(10,20,40,80,160,320,640),
             labels = c(10,20,40,80,160,320,640), range=c(1,3)) +
  scale_fill_manual(values=c("#DC7985","#8BCB9B"))
print(p)
dev.off()

# color pat 27
fractions_x_sub$response_color = fractions_x_sub$response
fractions_x_sub[fractions_x_sub$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/Pretreatment_Treg_fractions_by_response_CD48NK_pat27.pdf"), width=3, height=2.9)
p = ggplot(fractions_x_sub, aes(x=mc_group, y=fraction)) +
  geom_boxplot(aes(fill = response), outlier.shape = NA)  +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0), 
             aes(group = response, size = count,fill=response_color), pch = 21)+
  ggtitle("Tregs baseline") + ylim(0,0.4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(name   = "# cells",
             breaks = c(10,20,40,80,160,320,640),
             labels = c(10,20,40,80,160,320,640), range=c(1,3)) +
  scale_fill_manual(values=c("#DC7985","#8BCB9B","yellow"))
print(p)
dev.off()

fractions_x_sub_sel <- fractions_x_sub[fractions_x_sub$mc_group == "Treg TNFRSF9",]
wilcox.test(fractions_x_sub_sel[fractions_x_sub_sel$response=="RE",][["fraction"]],fractions_x_sub_sel[fractions_x_sub_sel$response=="NR",][["fraction"]])

fractions_x_sub_sel <- fractions_x_sub[fractions_x_sub$mc_group == "Treg GIMAP",]
wilcox.test(fractions_x_sub_sel[fractions_x_sub_sel$response=="RE",][["fraction"]],fractions_x_sub_sel[fractions_x_sub_sel$response=="NR",][["fraction"]])


## myeloids/B cells, non-T cells
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("M","B"),]
md = md[md$patient != "Pat04",]

# remove cells with TCR
md <- md[is.na(md$unique_clone_ID),]

data <- md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

# plot boxplots
fractions_x <- as.data.frame(x)
fractions_x_sub <- fractions_x[fractions_x$mc_group %in% c("cDC CD1C","pDC"),] #,"cDC LAMP3"
pdf(file=paste0(fig_path,"/Pretreatment_cDCs_fractions_by_response_MB.pdf"), width=3, height=2.9)
p = ggplot(fractions_x_sub, aes(x=mc_group, y=fraction)) +
  geom_boxplot(aes(fill = response), outlier.shape = NA)  +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0), 
             aes(group = response, size = count,fill=response), pch = 21)+
  ggtitle("cDC CD1C baseline") + ylim(0,0.6) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(name   = "# cells",
             breaks = c(10,20,40,80,160,320,640),
             labels = c(10,20,40,80,160,320,640), range=c(1,3)) +
  scale_fill_manual(values=c("#DC7985","#8BCB9B"))
print(p)
dev.off()

# color pat 27
fractions_x_sub$response_color = fractions_x_sub$response
fractions_x_sub[fractions_x_sub$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/Pretreatment_cDCs_fractions_by_response_MB_pat27.pdf"), width=3, height=2.9)
p = ggplot(fractions_x_sub, aes(x=mc_group, y=fraction)) +
  geom_boxplot(aes(fill = response), outlier.shape = NA)  +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0), 
             aes(group = response, size = count,fill=response_color), pch = 21)+
  ggtitle("cDC CD1C baseline") + ylim(0,0.6) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(name   = "# cells",
             breaks = c(10,20,40,80,160,320,640),
             labels = c(10,20,40,80,160,320,640), range=c(1,3)) +
  scale_fill_manual(values=c("#DC7985","#8BCB9B","yellow"))
print(p)
dev.off()

fractions_x_sub_sel <- fractions_x_sub[fractions_x_sub$mc_group == "pDC",]
wilcox.test(fractions_x_sub_sel[fractions_x_sub_sel$response=="RE",][["fraction"]],fractions_x_sub_sel[fractions_x_sub_sel$response=="NR",][["fraction"]])

fractions_x_sub_sel <- fractions_x_sub[fractions_x_sub$mc_group == "cDC CD1C",]
wilcox.test(fractions_x_sub_sel[fractions_x_sub_sel$response=="RE",][["fraction"]],fractions_x_sub_sel[fractions_x_sub_sel$response=="NR",][["fraction"]])



#### Figure 2D ####
# ratio Tregs
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("CD4"),]
md = md[md$patient != "Pat04",]

data = md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

y = x[x$mc_group %in% c("Treg GIMAP", "Treg TNFRSF9"), ]

plot_treg = pivot_wider(data = y, 
                        id_cols = c("patient", "response"), 
                        names_from = mc_group, 
                        values_from = c("count", "total", "fraction"))

plot_treg$ratio = log2(plot_treg$`fraction_Treg TNFRSF9`/plot_treg$`fraction_Treg GIMAP`)

pdf(file=paste0(fig_path,"/Treg_ratio_baseline_boxplot.pdf"), width=2.1, height=2.8) #
p=ggplot(plot_treg, aes(x = response, y = ratio))+
  geom_boxplot(aes(fill = response), outlier.shape = NA)+
  geom_point(position=position_jitter(width = 0.1, height = 0), 
             aes(group = response,fill=response),pch = 21, size = 2)+
  scale_fill_manual(values=c("#DC7985","#8BCB9B")) + ylim(c(-1,5))+
  ylab("log2(Treg_TNFRSF9/Treg_GIMAP)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

# color pat 27
plot_treg$response_color = plot_treg$response
plot_treg[plot_treg$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/Treg_ratio_baseline_boxplot_pat27.pdf"), width=2.3, height=2.8) #
p=ggplot(plot_treg, aes(x = response, y = ratio))+
  geom_boxplot(aes(fill = response), outlier.shape = NA)+
  geom_point(position=position_jitter(width = 0.1, height = 0), 
             aes(group = response,fill=response_color),pch = 21, size = 2)+
  scale_fill_manual(values=c("#DC7985","#8BCB9B","yellow")) + ylim(c(-1,5))+
  ylab("log2(Treg_TNFRSF9/Treg_GIMAP)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

wilcox.test(plot_treg[plot_treg$response == "RE",]$ratio,plot_treg[plot_treg$response == "NR",]$ratio)



## Figure 2D
## correlation Tregs 
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)
md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$cell_type) &
          md$condition == "pre" &
          md$cell_type %in% c("CD4","CD8","NK") ,]
md = md[md$patient != "Pat04",]

data = md[,c("patient", "response", "condition", "cell_type","mc_group","mc_group_subdiv")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, mc_group) %>% 
  summarize(count=n()) %>% 
  mutate(fraction=count/sum(count), total=sum(count)) 

t_both <- cbind(as.data.frame(x[x$mc_group=="Treg GIMAP",]),
               as.data.frame(x[x$mc_group=="Treg TNFRSF9",]))

colnames(t_both)[c(5,11)] <- c("fraction_Treg_GIMAP","fraction_Treg_TNFRSF9")
t_both <- t_both[,c("response","fraction_Treg_GIMAP","fraction_Treg_TNFRSF9","patient")]

pdf(file=paste0(fig_path,"/ratio_correlations_tregs_plot.pdf"), width=3.4, height=2.8)
p = ggplot(t_both, aes(x=fraction_Treg_GIMAP, y=fraction_Treg_TNFRSF9,color=response)) +
  geom_point() + 
  ggtitle("Baseline_Treg_correlation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#DC7985","#8BCB9B"))
print(p)
dev.off()

# color pat 27
t_both$response_color = t_both$response
t_both[t_both$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/ratio_correlations_tregs_plot_pat27.pdf"), width=3.8, height=2.8)
p = ggplot(t_both, aes(x=fraction_Treg_GIMAP, y=fraction_Treg_TNFRSF9,color=response_color)) +
  geom_point() + 
  ggtitle("Baseline_Treg_correlation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#DC7985","#8BCB9B","yellow"))
print(p)
dev.off()



#### Figure 2E ####
# UCell score
treg_genes <- read.table("./data/signatures/Treg_activation_genes.txt",header=T) # Treg activation genes

mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]

geneSets <- list(geneSet1=treg_genes$x)
geneNames <- rownames(mat@mat)

md = submat[submat$mc_group %in% c("Treg TNFRSF9","Treg GIMAP") &
              submat$response %in% c("RE","NR") &
              submat$condition %in% c("pre"),] 
md = md[md$patient != "Pat04",]
exprMatrix <- as.matrix(mat@mat[,rownames(md)])
rownames(exprMatrix) = geneNames

cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(geneSets,cell_rank)

plot_data <- cells_AUC@assays@data@listData$AUC
plot_data <- as.data.frame(t(plot_data))

# plot per cell state
plot_data$cell_id <- rownames(plot_data)
md$cell_id <- rownames(md)
sig_scores_all <- plyr::join(plot_data,md,by="cell_id")

pdf(file=paste0(fig_path,"/Treg_activation_signature_aucell.pdf"), width=1.9, height=3.2)
p = ggplot(sig_scores_all, aes(y = geneSet1, x = mc_group,fill=mc_group))+
  geom_violin() + geom_boxplot(width=.1,position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0042b0","#59b56f","#9B1F23","116936"))+
  theme_bw() + ylab("Treg_activation_genes") + ggtitle("Treg_activation_genes") + ylim(c(0,0.8)) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
print(p)
dev.off()

sig_scores_sel_a = sig_scores_all
wilcox.test(sig_scores_sel_a[sig_scores_sel_a$mc_group == "Treg GIMAP",]$geneSet1,
            sig_scores_sel_a[sig_scores_sel_a$mc_group == "Treg TNFRSF9",]$geneSet1)






