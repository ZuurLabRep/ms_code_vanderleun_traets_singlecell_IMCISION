###################
#### FIGURE S2 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S2"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S2")

library(AUCell)
library(GSEABase)
library(ggrepel)
library(tidyverse)
library(Startrac)
library(circlize)
library(ComplexHeatmap)


#### Figure S2A ####
## plot CTLA4, PD1
mc_colors <- read.table("./data/mc_colors.txt",header=T,sep="\t")

# scale data
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)
submat = mat@cell_metadata[names(mc@mc),]
umis = as.matrix(mat@mat[,rownames(submat)])
umis_n = sweep(umis,2,colSums(umis), "/") * 10000 

gene_oi = c("CTLA4")
plot_violin_gene_expression(mat_name=clean_mat_id,umis_n=umis_n,mc_name=filt_mc_id,gene_oi=gene_oi,mc_colors=mc_colors)

gene_oi = c("PDCD1")
plot_violin_gene_expression(mat_name=clean_mat_id,umis_n=umis_n,mc_name=filt_mc_id,gene_oi=gene_oi,mc_colors=mc_colors)



#### Figure S2B ####
facs = read.csv("./data/20220112 FACS stats all samples R input.csv")
check = as.data.frame(cbind(mat@cell_metadata$patient,mat@cell_metadata$response))
check = check[!duplicated(check),]
facs[facs$patient %in% c(27,36),]$response = "RE" # fix old annotation
facs = facs[facs$patient != 4,]

facs = facs[facs$condition == "pre",]
pdf(file=paste0(fig_path,"/FACS_PD1_baseline_response_rev.pdf"), width=2, height=2.6) #
p=ggplot(facs,aes(x=response,y=PD1_expr,fill=response)) +  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = response, fill=response), pch = 21, alpha = 0.7,
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0)) +
  geom_point(data=facs[facs$patient == "27",],aes(x=response,y=PD1_expr,fill=response,group = response),color="yellow")+
  theme_bw() + scale_fill_manual(values = c("#dc7985","#8bcb9b")) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,2300)
print(p)
dev.off()

wilcox.test(facs[facs$response =="RE",]$PD1_expr,facs[facs$response =="NR",]$PD1_expr)



#### Figure S2C ####
## barplot compartments baseline
celltype_proportions(clean_mat_id,filt_mc_id,"pre","RE","CD8",c("CD8", "CD4", "NK", "B", "M"))
celltype_proportions(clean_mat_id,filt_mc_id,"pre","NR","CD8",c("CD8", "CD4", "NK", "B", "M"))


#### Figure S2D ####
## boxplots compartments baseline, NR vs RE
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

md <- mat@cell_metadata[names(mc@mc),]
md <- md[!is.na(md$mc_group) &
          md$condition == "pre",]
md <- md[md$patient != "Pat04",]
  
# remove M,B labeled cells with TCR
md_t <- md[md$cell_type %in% c("CD4","CD8","NK"),] 
md_m <- md[md$cell_type %in% c("M","B"),] 
md_m <- md_m[is.na(md_m$unique_clone_ID),]
md <- rbind(md_t,md_m)

data <- md[,c("patient", "response", "condition", "cell_type","mc_group")]

# calculate fractions
x <- data %>% 
  group_by(patient, response, cell_type) %>% 
  summarize(n=n()) %>% 
  mutate(fraction=n/sum(n), total=sum(n)) 

fractions_x <- as.data.frame(x)
fractions_x$response_color = fractions_x$response
fractions_x[fractions_x$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/Patients_boxplot_response_pre.pdf"), width=4, height=3)
p=ggplot(fractions_x, aes(x=cell_type, y=fraction)) +
  geom_boxplot(aes(fill = response), outlier.shape = NA)  +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0), 
             aes(group = response,fill=response_color), pch = 21)+
  ggtitle("cell_type_baseline") + ylim(0,0.8) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("#DC7985","#8BCB9B","yellow"))
print(p)
dev.off()

for (cell_x in c("NK","CD4","CD8","M","B")){
  print(cell_x)
  print(wilcox.test(fractions_x[fractions_x$cell_type == cell_x & fractions_x$response == "RE",]$fraction,
              fractions_x[fractions_x$cell_type == cell_x & fractions_x$response == "NR",]$fraction))
}


#### Figure S2E ####
# boxpot
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          md$cell_type %in% c("CD8") & 
          md$condition == "pre" &
          md$patient != "Pat04",]

## clonality startrac
data_pat = md[,c("patient","response")]
data_pat = data_pat[!duplicated(data_pat$patient),]

# Run Startrac on data
data <- as.data.frame(md)
data$Cell_Name <- rownames(data)
data$clone.id <- data[,colnames(data) == "unique_clone_ID"]
data$majorCluster <- data[,colnames(data) == "patient"]
data$loc <- data[,colnames(data) == "condition"] 
out <- Startrac.run(data, proj="IMCISION",verbose=F)

# get data out of startrac object
out_cluster <- as.data.frame(out@cluster.data)
out_cluster <- out_cluster[!is.na(out_cluster$expa),]
out_cluster_tmp <- out_cluster[out_cluster$aid != "IMCISION",]
out_cluster_tmp = out_cluster_tmp[order(out_cluster_tmp$expa,decreasing = T),]
out_cluster_tmp$majorCluster = factor(out_cluster_tmp$majorCluster, levels=out_cluster_tmp$majorCluster)

# combine with response
colnames(out_cluster_tmp)[2] <- "patient"
out_cluster_tmp <- merge(out_cluster_tmp,data_pat,by="patient")

pdf(file=paste0(fig_path,"/Startrac_expansion_index_boxplots_baseline_CD8.pdf"), width=2.2, height=3)
startrac = ggplot(out_cluster_tmp,aes(x=response,y=expa,fill=response)) + geom_boxplot(outlier.alpha = 0)+
  geom_jitter(width = 0.05)+ ylab("expansion index") + theme_bw() + 
  scale_fill_manual(values=c("#E18283","#71BB87"))+
  ylim(c(0.03,0.18)) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
startrac
dev.off()

out_cluster_tmp$response_color = out_cluster_tmp$response
out_cluster_tmp[out_cluster_tmp$patient == "Pat27",]$response_color = "RE_pat27"
pdf(file=paste0(fig_path,"/Startrac_expansion_index_boxplots_baseline_CD8_pat27.pdf"), width=2.2, height=3)
startrac = ggplot(out_cluster_tmp,aes(x=response,y=expa)) + geom_boxplot(outlier.alpha = 0,aes(fill=response))+
  geom_jitter(aes(color=response_color),width = 0.05)+ ylab("expansion index") + theme_bw() + 
  scale_color_manual(values=c("black","black","yellow"))+
  scale_fill_manual(values=c("#E18283","#71BB87"))+
  ylim(c(0.03,0.18)) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
startrac
dev.off()


# barplots
mc = scdb_mc(filt_mc_id)
mat = scdb_mat(clean_mat_id)

md = mat@cell_metadata[names(mc@mc),]
md = md[!is.na(md$mc_group) &
          md$cell_type %in% c("CD8") & #"CD4",
          md$condition == "pre" &
          md$patient != "Pat04" &
          !is.na(md$unique_clone_ID),]

md_freq <- md %>% group_by(patient,pat_clone_ID) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n), total=sum(n))
md_freq <- md_freq[order(md_freq$freq,decreasing = T),]

md_freq$top_clones <- "none"
for (x in unique(md_freq$patient)){
  md_freq[md_freq$patient == x,][1,]$top_clones <- "top_1"
  md_freq[md_freq$patient == x,][2:5,]$top_clones <- "top_2_5"
  md_freq[md_freq$patient == x,][6:15,]$top_clones <- "top_6_15"
}

md_freq <- merge(md_freq,data_pat,by="patient")

md_freq$top_clones = factor(md_freq$top_clones,levels=c("top_1","top_2_5","top_6_15","none"))
md_freq$patient = factor(md_freq$patient, levels=c("Pat30","Pat38","Pat10","Pat26","Pat34","Pat37","Pat33",
                                                   "Pat12","Pat36","Pat21","Pat15","Pat31","Pat27","Pat22",
                                                   "Pat39","Pat17","Pat29")) # order figure S2C
pdf(file=paste0(fig_path,"/Top_most_n_TCRs_patient_baseline_CD8_order.pdf"), width=7.5, height=3.4)
g= ggplot(md_freq,aes(x=patient,y=n,fill=top_clones)) + 
  geom_bar(stat="identity") +
  theme_bw() + ggtitle("Most abundant TCRs, abs (CD8)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~response,scales = "free_x",space = "free_x",) +
  scale_fill_manual(values=c("lightblue","#c99dc1","#9dc9a0","grey"))
g
dev.off()

pdf(file=paste0(fig_path,"/Top_most_n_TCRs_patient_baseline_CD8_order_NR.pdf"), width=4.5, height=3.4)
g= ggplot(md_freq[md_freq$response == "NR",],aes(x=patient,y=n,fill=top_clones)) + 
  geom_bar(stat="identity") +
  theme_bw() + ggtitle("Most abundant TCRs, abs (CD8)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("lightblue","#c99dc1","#9dc9a0","grey"))
g
dev.off()

pdf(file=paste0(fig_path,"/Top_most_n_TCRs_patient_baseline_CD8_order_RE.pdf"), width=5.2, height=3.4)
g= ggplot(md_freq[md_freq$response == "RE",],aes(x=patient,y=n,fill=top_clones)) + 
  geom_bar(stat="identity") +
  theme_bw() + ggtitle("Most abundant TCRs, abs (CD8)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("lightblue","#c99dc1","#9dc9a0","grey"))
g
dev.off()

#### Figure S2F ####
# multiplex data
multi_imci = read.table("./data/Multiplex_IMCISION_stainings_Joris_paper.txt",sep="\t",header=T)

mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

annotations <- mat@cell_metadata[,c("patient","response")]
annotations <- annotations[!is.na(annotations$patient),]
annotations <- annotations[!duplicated(annotations$patient),]

# CD3+FOXP3+ / CD3+ - fix columns
multi_imci$patient <- multi_imci$ID
multi_imci$patient[1:68] <- paste0("Pat0",multi_imci$patient[1:68])
multi_imci$patient[69:252] <- paste0("Pat",multi_imci$patient[69:252])
multi_imci$patient_timepoint <- paste0(multi_imci$patient,"_",multi_imci$Time.point)
multi_imci_cast <- dcast(data = multi_imci,formula = patient_timepoint ~ Phenotype,fun.aggregate = sum,value.var = "Intratumoral.density..cells.mm2.")

multi_imci_cast$patient <- do.call(rbind,strsplit(multi_imci_cast$patient_timepoint, "_"))[,1]
multi_imci_cast$timepoint <- do.call(rbind,strsplit(multi_imci_cast$patient_timepoint, "_"))[,2]

# match patients, and only pre treatment
multiplex_im_anno <- merge(multi_imci_cast,annotations,by="patient")
to_plot <- multiplex_im_anno
to_plot <- to_plot[to_plot$timepoint == "Pre-treatment",]

rownames(to_plot) <- to_plot$patient
to_plot_scaled <- t(scale(to_plot[,c(4:6)])) 
col_ann <- HeatmapAnnotation(response = to_plot$response,
                           col=list( response=c("RE"= "#89CC9B","NR"="#DC7985")))

to_plot_nonscaled <- to_plot[,c(4:6)]
col_fun <- colorRamp2(c(0, 600), c("white", "red"))
to_plot_single <- as.data.frame(to_plot_nonscaled[,2])
colnames(to_plot_single) <- colnames(to_plot_nonscaled)[2]
pdf(file=paste0(fig_path,"/Multiplex_Baseline_selected_heatmap_split_CD8_CD3_abs.pdf"), width=6, height=2)
Heatmap(t(to_plot_single),column_title="Multiplex_baseline", name="abs", top_annotation = col_ann,
        col = col_fun)
dev.off()


