###################
#### FIGURE S1 ####
###################

setwd("/DATA/project")
source("./code/Figures_general_functions.R")
tbk_reload()

dir.create(file.path("./figures", "Figure_S1"), showWarnings = FALSE)
fig_path = paste0(.scfigs_base,"Figure_S1")

library(Startrac)

#### Figure S1A ####

#### Figure S1B-C ####
# in pre-processing part


#### Figure S1D ####
mat = scdb_mat(clean_mat_id)
mc = scdb_mc(filt_mc_id)

TCR_coverage(clean_mat_id, width = 2400, heigth = 600)


#### Figure S1E ####
data = mat@cell_metadata[names(mc@mc),]
data = data[data$clonal_exp != "no" & !is.na(data$response),]
data = data[data$mc_group %in% c("Dysf GZMB","Dysf ZNF683","Tfh LAG3","Tfh NR3C1","Transitional",
                                 "Naive-like CD8", "Cytotoxic","Naive-like CD4","Treg GIMAP","Treg TNFRSF9"),]

# Run Startrac on data
data <- as.data.frame(data)
data$Cell_Name <- rownames(data)
data$clone.id <- data[,colnames(data) == "unique_clone_ID"]
data$majorCluster <- data[,colnames(data) == "mc_group"]
data$loc <- data[,colnames(data) == "condition"] 
data_sub <- dplyr::select(data,"Cell_Name","clone.id","patient","majorCluster","loc")
out <- Startrac.run(data, proj="IMCISION",verbose=F)

# get data out of startrac object
out_cluster <- as.data.frame(out@cluster.data)
out_cluster <- out_cluster[!is.na(out_cluster$expa),]
out_cluster_tmp <- out_cluster[out_cluster$aid != "IMCISION",]

order_cells <- levels(as.factor(out_cluster_tmp$majorCluster))

# plot
mc_colors <- read.table("./data/mc_colors.txt",header=T,sep="\t")
mc_colors$name <- gsub("_"," ",mc_colors$name)
mc_colors <- mc_colors[match(order_cells, mc_colors$name),]

pdf(file=paste0(fig_path,"/Startrac_expansion_index_boxplots.pdf"), width=5, height=3.5)
startrac = ggplot(out_cluster_tmp, aes(x = majorCluster, y = expa))+
  geom_point(position=position_jitter(width = 0.1, height = 0), 
             aes(group = aid), fill = "lightgrey", 
             alpha = 0.8)+
  scale_fill_manual(values = mc_colors$color)+
  stat_boxplot(aes(fill = majorCluster), outlier.shape = NA)+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))+
  ylab("expansion index") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
startrac
dev.off()



